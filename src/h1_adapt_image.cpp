#include <hermes2d.h>
#include "shared/img2d.h"
#include "pix_by_pix_integrator.h"
#include "h1_adapt_image.h"

#define H2D_FEI_VALUE  0
#define H2D_FEI_DX     1
#define H2D_FEI_DY     2

H1AdaptImage::H1AdaptImage(const PixByPixInt* pix_by_pix_int, H1Space* space)
  : H1Adapt(space), image(pix_by_pix_int->get_image()), pix_by_pix_int(pix_by_pix_int) {}

H1AdaptImage::~H1AdaptImage() {}

scalar H1AdaptImage::eval_error(biform_val_t bi_fn, biform_ord_t bi_ord,
                    MeshFunction *sln1, MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2,
                    RefMap *rv1,        RefMap *rv2,        RefMap *rrv1,        RefMap *rrv2) {
  assert_msg(sln1 == sln2, "The (coarse) solution 1 is not equal to the (coarse) solution 2. Only one component is supported.");
  assert_msg(rsln1 == rsln2, "The reference solution 1 is not equal to the refence solution 2. Only one component is supported.");

  int sln_order = sln1->get_fn_order();
  Element* e = sln1->get_active_element();
  scalar error_squared = pix_by_pix_int->integrate_element(e, sln_order, sln_order, eval_element_error, sln1, true);
  return error_squared;
}

scalar H1AdaptImage::eval_norm(biform_val_t bi_fn, biform_ord_t bi_ord,
                               MeshFunction *rsln1, MeshFunction *rsln2, RefMap *rrv1, RefMap *rrv2) {
  assert_msg(rsln1 == rsln2, "The reference solution 1 is not equal to the refence solution 2. Only one component is supported.");

  Element* e = rsln1->get_active_element();
  scalar norm_squared = pix_by_pix_int->integrate_element(e, image->get_interpolation_h_order(), image->get_interpolation_v_order(), eval_element_norm, NULL, true); /** \todo Norm is a function of the image: an integral image can be used instead. */
  return norm_squared;
}

double H1AdaptImage::eval_element_error(const Vector2D& ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr) {
  Solution* sln = (Solution*)pars_ptr;

  //obtain values from image
  double img_value, img_dx, img_dy;
  img_value = image->get_sample(phys.x, phys.y, img_dx, img_dy);

  //obtain values from element
  double sln_value = sln->get_ref_value(element, ref.x, ref.y, 0, H2D_FEI_VALUE);
  double sln_dx = sln->get_ref_value(element, ref.x, ref.y, 0, H2D_FEI_DX) * diff_rescale.x; //remapping to physical domain
  double sln_dy = sln->get_ref_value(element, ref.x, ref.y, 0, H2D_FEI_DY) * diff_rescale.y;

  //calculate
  return (sqr(img_value - sln_value) + sqr(img_dx - sln_dx) + sqr(img_dy - sln_dy));
}

double H1AdaptImage::eval_element_norm(const Vector2D& ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr) {
  //obtain values from image
  double img_value, img_dx, img_dy;
  img_value = image->get_sample(phys.x, phys.y, img_dx, img_dy);

  return (sqr(img_value) + sqr(img_dx) + sqr(img_dy));
}

void H1AdaptImage::add_neighbours_nonadapt(Mesh* mesh, const int element_id, const int start_vertex_id, const int end_vertex_id, const Range<double>& range_x, const Range<double>& range_y) {
  int vertex_id = start_vertex_id;
  int prev_element_id = element_id;
  do {
    Element* e;
    for_all_active_elements(e, mesh) { //TODO: not efficient, K-D tree would be more fancy
      if (e->id != element_id && e->id != prev_element_id) { //skip the processed element        
        //find vertex that matches current vertex
        const int nvert = e->nvert;
        int inx_vert = 0;
        while(inx_vert < nvert && e->vn[inx_vert]->id != vertex_id)
          inx_vert++;
        if (inx_vert < nvert) { //find which vertex is the next one, if any
          int inx_vert_edge = -1;
          int inx_vert_next = (inx_vert+1) % nvert, inx_vert_prev = (inx_vert-1+nvert) % nvert;
          if (range_x.is_in_closed(e->vn[inx_vert_next]->x) && range_y.is_in_closed(e->vn[inx_vert_next]->y))
            inx_vert_edge = inx_vert_next;
          else if (range_x.is_in_closed(e->vn[inx_vert_prev]->x) && range_y.is_in_closed(e->vn[inx_vert_prev]->y))
            inx_vert_edge = inx_vert_prev;
          
          //if a vertex of an edge is found, add to a map and restart the search
          if (inx_vert_edge >= 0) {
            no_adapt_elem[e->id] = 1;
            prev_element_id = e->id;
            vertex_id = e->vn[inx_vert_edge]->id;
            break; //restart the search for elements
          }
        }
      }
    }
  } while (vertex_id != end_vertex_id && vertex_id != start_vertex_id); //stop if end vertex was found or no other vertex was found (element is on an edge)
}

void H1AdaptImage::add_neighbours_1px_to_adapt(const Mesh* mesh, const Element* elem_input) {
#define EDGE_CNT 4
  //assumes that all coordinates are integers despite that they are stored to doubles

  //edge pair indices
  int oriented_edge_pairs[EDGE_CNT][2] = { {0, 1}, {1, 2}, {3, 2}, {0, 3} };

  //build list of edges
  EdgeRange elem_input_edges[EDGE_CNT];
  bool elem_input_edge_unit[EDGE_CNT];
  for(int i = 0; i < EDGE_CNT; i++) {
    int inx = oriented_edge_pairs[i][0], inx_next = oriented_edge_pairs[i][1];
    elem_input_edges[i] = EdgeRange(Range<double>(elem_input->vn[inx]->x, elem_input->vn[inx_next]->x), Range<double>(elem_input->vn[inx]->y, elem_input->vn[inx_next]->y));
    double edge_len = std::max(elem_input_edges[i].x.upper() - elem_input_edges[i].x.lower(), elem_input_edges[i].y.upper() - elem_input_edges[i].y.lower());
    if (std::abs(edge_len - 1) < 1E-13)
      elem_input_edge_unit[i] = true;
    else
      elem_input_edge_unit[i] = false;
  }
 
  //inspect active elements
  Element* elem_active;
  for_all_active_elements(elem_active, mesh) {
    //inspect an active element if it is larger than 1x1
    double elem_active_width = elem_active->vn[2]->x - elem_active->vn[0]->x; 
    double elem_active_height = elem_active->vn[2]->y - elem_active->vn[0]->y;
    if (elem_active_width > 1 || elem_active_height > 1) { //element is examined as well but it is done just once and it is excluded by later tests
      
      //inspect all edges of an active element
      for(int i = 0; i < EDGE_CNT; i++) {
        //make edge of an active element
        int inx = oriented_edge_pairs[i][0], inx_next = oriented_edge_pairs[i][1];
        EdgeRange elem_active_edge(Range<double>(elem_active->vn[inx]->x, elem_active->vn[inx_next]->x), Range<double>(elem_active->vn[inx]->y, elem_active->vn[inx_next]->y));
        double elem_active_edge_len = std::max(elem_active_edge.x.upper() - elem_active_edge.x.lower(), elem_active_edge.y.upper() - elem_active_edge.y.lower());
        
        //check an edge of an active element if a length of the edge is greater then 1 (it can be refined and the mesh can be improved)
        if (elem_active_edge_len > 1) {
          for(int k = 0; k < EDGE_CNT; k++) {
            if (elem_input_edge_unit[k] && elem_active_edge.is_in_closed(elem_input_edges[k])) {
              priority_queue.push(ElementReference(elem_active->id, 0));
              break;
            }
          }
        }
      }
    }
  }
}

bool H1AdaptImage::should_ignore_element(const int inx_element, const Mesh* mesh, const Element* element) {
  //init 'do-not-adapt' map
  if (inx_element == 0) { //clean a 'do-not-adapt' map
    no_adapt_elem.clear();
  }

  //check if element is in the 'do-not-adapt' map
  int id = element->id;
  if (no_adapt_elem.find(id) != no_adapt_elem.end())
    return true; //element wat found in the map: ignore it for now

  //if element is too small (i.e., it cannot be split), ignore it
  double element_width = element->vn[2]->x - element->vn[0]->x; 
  double element_height = element->vn[2]->y - element->vn[0]->y;
  if (element_width < 2 && element_height < 2) {
    add_neighbours_1px_to_adapt(mesh, element);
    return true;
  }

  //refine element for processing (it not already refined)
  Mesh* rmesh = rsln[0]->get_mesh();
  Element* relem = rmesh->get_element(id);
  if (relem->active)
    rmesh->refine_element(id);

  return false;
}

bool H1AdaptImage::can_refine_element(Mesh* mesh, Element* e, bool refined, ElementToRefine& elem_ref) {
  //if the element was not refined even though there was a high error: this is due to a hanging node (probably)
  if (!refined) {
    add_neighbours_1px_to_adapt(mesh, e);
  } 
  else { //check whether the refinement deas not lead an element below a pixel size
    //prevent creation of elements smaller than 1 pixel
    double element_width = e->vn[2]->x - e->vn[0]->x; 
    double element_height = e->vn[2]->y - e->vn[0]->y;

    //check if refinement would not break the element below size of a pixel
    if ((element_width < 2 && (elem_ref.split == H2D_REFINEMENT_H || elem_ref.split == H2D_REFINEMENT_ANISO_V))
      || (element_height < 2 && (elem_ref.split == H2D_REFINEMENT_H || elem_ref.split == H2D_REFINEMENT_ANISO_H))) {
        return false; /** \todo Replace by removing inappropriate candidates. */
    }

    //element will be refined: add neighbours to a do-not-refine hash table
    //switch (elem_ref.split) {
    //case H2D_REFINEMENT_H:
    //  {
    //    for(int i = 0; i < 4; i++) {
    //      Node *n0 = e->vn[i], *n1 = e->vn[(i+1) % 4];
    //      add_neighbours_nonadapt(mesh, e->id, n0->id, n1->id, Range<double>(n0->x, n1->x), Range<double>(n0->y, n1->y));
    //    }
    //  }
    //  break;

    //case H2D_REFINEMENT_ANISO_H:
    //  {
    //    Node *n0 = e->vn[0], *n1 = e->vn[3];
    //    add_neighbours_nonadapt(mesh, e->id, n0->id, n1->id, Range<double>(n0->x, n1->x), Range<double>(n0->y, n1->y));
    //    n0 = e->vn[1], n1 = e->vn[2];
    //    add_neighbours_nonadapt(mesh, e->id, n0->id, n1->id, Range<double>(n0->x, n1->x), Range<double>(n0->y, n1->y));
    //  }
    //  break;

    //case H2D_REFINEMENT_ANISO_V:
    //  {
    //    Node *n0 = e->vn[0], *n1 = e->vn[1];
    //    add_neighbours_nonadapt(mesh, e->id, n0->id, n1->id, Range<double>(n0->x, n1->x), Range<double>(n0->y, n1->y));
    //    n0 = e->vn[3], n1 = e->vn[2];
    //    add_neighbours_nonadapt(mesh, e->id, n0->id, n1->id, Range<double>(n0->x, n1->x), Range<double>(n0->y, n1->y));
    //  }
    //  break;
    //}
  }

  //add element to a list of non-processed elements
  no_adapt_elem[e->id] = 1;

  return refined;
}

void H1AdaptImage::apply_refinements(std::vector<ElementToRefine>& elems_to_refine) {
  //sort refinements
  std::sort(elems_to_refine.begin(), elems_to_refine.end(), is_elem_ref_ordered);

  //apply refinements
  H1Adapt::apply_refinements(elems_to_refine);
}

bool H1AdaptImage::is_elem_ref_ordered(const ElementToRefine& a, const ElementToRefine& b) {
  if (a.split != b.split)
    return a.split < b.split;
  else { //equal split
    int num_sons = a.get_num_sons();
    for(int i = 0; i < num_sons; i++) {
      if (H2D_GET_H_ORDER(a.p[i]) != H2D_GET_H_ORDER(b.p[i]))
        return H2D_GET_H_ORDER(a.p[i]) < H2D_GET_H_ORDER(b.p[i]);
      else
        return H2D_GET_V_ORDER(a.p[i]) < H2D_GET_V_ORDER(b.p[i]);
    }
    //equal orders
    return a.id < b.id;
  }
}
