#include "himg.h"
#include "pix_by_pix_integrator.h"
#include "h1_image_proj_based_selector.h"

void H1ImageProjBasedSelector::create_candidates(Element* e, int quad_order, int max_ha_quad_order, int max_p_quad_order) {
  H1ProjBasedSelector::create_candidates(e, quad_order, max_ha_quad_order, max_p_quad_order);

  //prevent creation of elements smaller than 1 pixel
  double element_width = e->vn[2]->x - e->vn[0]->x; 
  double element_height = e->vn[2]->y - e->vn[0]->y;

  //remove inappropriate candidates
  if (element_width < 2 || element_height < 2) {
    std::vector<Cand> new_candidates;
    new_candidates.reserve(candidates.size());
    for(unsigned i = 0; i < candidates.size(); i++) {
      const Cand& cand = candidates[i];
      if ((element_width > 1 || (cand.split != H2D_REFINEMENT_H && cand.split != H2D_REFINEMENT_ANISO_V))
        && (element_height > 1 || (cand.split != H2D_REFINEMENT_H && cand.split != H2D_REFINEMENT_ANISO_H))) {
          new_candidates.push_back(cand);
      }
    }

    //switch lists
    candidates.swap(new_candidates);
  }
}


scalar H1ImageProjBasedSelector::evaluate_rsh_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, int shape_inx) {
  const int quad_order = shapeset->get_order(shape_inx);
  EvalRHSParams rhs_params = { sub_trf, shapeset, shape_inx };
  double rsh = pix_by_pix_int->integrate_element(sub_elem, H2D_GET_H_ORDER(quad_order), H2D_GET_V_ORDER(quad_order), eval_right_side, &rhs_params, false);
  return rsh;
}

double H1ImageProjBasedSelector::eval_right_side(const Vector2D& sub_ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr) {
  EvalRHSParams& pars = *(EvalRHSParams*)pars_ptr;
  
  //calculate coordinates in the reference domain of the element
  Vector2D ref(sub_ref.x * pars.sub_trf.trf->m[0] + pars.sub_trf.trf->t[0], sub_ref.y * pars.sub_trf.trf->m[1] + pars.sub_trf.trf->t[1]);

  //get data from image
  double img_value, img_dx, img_dy;
  img_value = image->get_sample(phys.x, phys.y, img_dx, img_dy);
  img_dx *= diff_rescale.x * pars.sub_trf.coef_mx; //convert from physical to reference domain of a sub-element and from that domain to reference domain of an element
  img_dy *= diff_rescale.y * pars.sub_trf.coef_my;

  //get value of a shape function
  double shape_value = pars.shapeset->get_value(H2D_FEI_VALUE, pars.shape_inx, ref.x, ref.y, 0);
  double shape_dx = pars.shapeset->get_value(H2D_FEI_DX, pars.shape_inx, ref.x, ref.y, 0);
  double shape_dy = pars.shapeset->get_value(H2D_FEI_DY, pars.shape_inx, ref.x, ref.y, 0);

  //calculate
  return (img_value * shape_value) + (img_dx * shape_dx) + (img_dy * shape_dy);
}

double H1ImageProjBasedSelector::evaluate_error_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj) {
  EvalProjErrorParams proj_error_params = { sub_trf, shapeset, elem_proj };
  double error = pix_by_pix_int->integrate_element(sub_elem, H2D_GET_H_ORDER(elem_proj.max_quad_order), H2D_GET_V_ORDER(elem_proj.max_quad_order), eval_proj_error, &proj_error_params, false);
  return error;
}

double H1ImageProjBasedSelector::eval_proj_error(const Vector2D& sub_ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr) {
  EvalProjErrorParams& pars = *(EvalProjErrorParams*)pars_ptr;
  
  //calculate coordinates in the reference domain of the element
  Vector2D ref(sub_ref.x * pars.sub_trf.trf->m[0] + pars.sub_trf.trf->t[0], sub_ref.y * pars.sub_trf.trf->m[1] + pars.sub_trf.trf->t[1]);

  //get data from image
  double img_value, img_dx, img_dy;
  img_value = image->get_sample(phys.x, phys.y, img_dx, img_dy);
  img_dx *= diff_rescale.x * pars.sub_trf.coef_mx; //convert from physical to reference domain of a sub-element and from that domain to reference domain of an element
  img_dy *= diff_rescale.y * pars.sub_trf.coef_my;

  //get value of a projection
  double proj_value = 0, proj_dx = 0, proj_dy = 0;
  for(int i = 0; i < pars.elem_proj.num_shapes; i++) {
    const double coef = pars.elem_proj.shape_coefs[i];
    const int shape_inx = pars.elem_proj.shape_inxs[i];
    proj_value += coef * pars.shapeset->get_value(H2D_FEI_VALUE, shape_inx, ref.x, ref.y, 0);
    proj_dx += coef * pars.shapeset->get_value(H2D_FEI_DX, shape_inx, ref.x, ref.y, 0);
    proj_dy += coef * pars.shapeset->get_value(H2D_FEI_DY, shape_inx, ref.x, ref.y, 0);
  }

  //calculate
  return sqr(img_value - proj_value)
    + sqr(img_dx - proj_dx)
    + sqr(img_dy - proj_dy);
}
