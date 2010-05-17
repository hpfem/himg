#include "himg.h"
#include "parameters.h"
#include "adapt_storage.h"

using namespace std;

/* filenames */
#define FILENAME_RASTERISED_SLN "%s.sln.%03d.png" ///< Format string of a rasterised solution: first parameter = input filenane, second paramter = iteration */
#define FILENAME_MESH "%s.%03d.mesh" ///< Format string of a mesh file: first parameter = input filenane, second paramter = iteration */
#define FILENAME_SLN_VECTOR "%s.sln_vec.txt" ///< Format string of solution vectors: first parameter = input filenane */

/* macros */
#define fwrite_value(__value, __type, __file) { __type tmp = (__type)(__value); fwrite(&tmp, sizeof(__type), 1, file); }
#define fread_value(__value, __type, __file, __err_action) { __type tmp; if (fread(&tmp, sizeof(__type), 1, file) != 1) { debug_log("Reading after EOF."); __err_action; } else __value = tmp; }

/// Stores solution vector.
static void store_solution_vector(const char* input_name, const int iteration, const std::vector<double>& sln_vec) {
  //create a filename
  FILE* file = NULL;
  char buffer[1024];
  sprintf(buffer, FILENAME_SLN_VECTOR, input_name);

  //determine open mode
  ios_base::open_mode mode = 0;
  if (iteration == HIMG_FIRST_ITERATION)
    mode = ios_base::out;
  else
    mode = ios_base::app;

  //open
  ofstream fout(buffer, mode);
  error_if(!fout, "Failure to open file \"%s\".", buffer);

  //write mark of little/big endianness
  if (iteration == HIMG_FIRST_ITERATION) {
    fout << string("# vector_size vector_comp_1 vector_comp_2 ...") << endl;
  }

  //write number of results
  fout << (int32_t)sln_vec.size() << ' ';

  //write results
  fout << scientific;
  for(unsigned i = 0; i < sln_vec.size(); i++)
    fout << '\t' << sln_vec[i];
  fout << endl;

  //close
  fout.close();
}

/// Stores solution vector.
static void store_refinements(const char* input_name, const int iteration, const std::vector<ElementToRefine>& refinements) {
  //create a filename
  stringstream out_fname;
  out_fname << input_name << ".ers";

  //determine open mode
  ios_base::open_mode mode = 0;
  if (iteration == HIMG_FIRST_ITERATION)
    mode = ios_base::out;
  else
    mode = ios_base::app;

  //open
  ElementToRefineStream fout(out_fname.str().c_str(), mode | ios_base::binary);
  error_if(!fout.is_open(), "Unable to access refinement stream \"%s\".", out_fname.str().c_str());

  //write
  fout << refinements;

  //close
  fout.close();
}

/// Rasterizes and stores solution.
static void store_solution(const char* input_name, const int iteration, const Image2d& image_sln) {
  //store
  char buffer[1024];
  sprintf(buffer, FILENAME_RASTERISED_SLN, input_name, iteration);
  image_sln.store_to_png_clamped(buffer, 0, 1);
}

/* stores mesh */
void store_mesh(const char* input_name, const int iteration, Solution* sln) {
  //open file
  char buffer[1024];
  sprintf(buffer, FILENAME_MESH, input_name, iteration);
  FILE* file = fopen(buffer, "wb");
  error_if(file == NULL, "Unable to open file \"%s\".", buffer);

  //write number of elements
  Mesh* mesh = sln->get_mesh();
  fwrite_value(mesh->get_num_active_elements(), uint32_t, file);

  //write elements
  Element* e;
  for_all_active_elements(e, mesh) {
    //store ID
    fwrite_value(e->id, uint32_t, file);

    //store vertices
    Node *bl = e->vn[0];				// bottom left node
    Node *tr = e->vn[2];				// top right node
    for(int i = 0; i < 4; i++) {
      Node* vertex = e->vn[i];
      fwrite_value(vertex->x, float, file);
      fwrite_value(vertex->y, float, file);
      fwrite_value(sln->get_ref_value(e, 
        ((vertex->x - bl->x) / (tr->x - bl->x)) * 2 - 1,
        ((vertex->y - bl->y) / (tr->y - bl->y)) * 2 - 1)
        , float, file);
    }
  }

  //finish
  fclose(file);
}

/* order pallette: a copy from hermes2d/order_view.cpp */
static int order_palette[] =
{
  0x7f7f7f,
  0x7f2aff,
  0x2a2aff,
  0x2a7fff,
  0x00d4aa,
  0x00aa44,
  0xabc837,
  0xffd42a,
  0xc87137,
  0xc83737,
  0xff0000
};

/* stores mesh and orders of elements to SVG */
void store_mesh_orders(const char* input_name, const int iteration, const int width, const int height, Space* space) {
  //open file
  stringstream out_fname;
  out_fname << input_name << '.' << setfill('0') << setw(3) << iteration << ".mesh.svg";
  FILE* fout = fopen(out_fname.str().c_str(), "wt");
  error_if(fout == NULL, "Unable to open file \"%s\"", out_fname.str().c_str());

  //write number of elements
  Mesh* mesh = space->get_mesh();

  //write header
#define SVG_HEADER "<?xml version=\"1.0\" standalone=\"no\"?>\n"\
  "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "\
  "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
  fprintf(fout, SVG_HEADER);
  fprintf(fout, "<svg width=\"%d\" height=\"%d\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n", width, height);
  fprintf(fout, "<desc>%s</desc>\n", input_name);
  fprintf(fout, "<g transform=\"scale(1,-1)\"><g transform=\"translate(0,%d)\">\n", -height); //move bottom of mesh flipped over Y to bottom of image

  //write elements
  Element* e;
  for_all_active_elements(e, mesh) {
    int quad_order = space->get_element_order(e->id);
    fprintf(fout, "<g id=\"%d\">", e->id);
    fprintf(fout, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"#%06x\"/>",
      (float)e->vn[0]->x, (float)e->vn[0]->y, (float)(e->vn[2]->x - e->vn[0]->x), (float)(e->vn[2]->y - e->vn[0]->y),
      order_palette[H2D_GET_H_ORDER(quad_order)]);
    if (H2D_GET_H_ORDER(quad_order) != H2D_GET_V_ORDER(quad_order))
      fprintf(fout, "<path d=\"M %.2f %.2f L %.2f %.2f L %.2f %.2f L %.2f %.2f z\" fill=\"#%06x\"/>",
        (float)e->vn[0]->x, (float)e->vn[0]->y, 
        (float)e->vn[2]->x, (float)e->vn[2]->y,
        (float)e->vn[1]->x, (float)e->vn[1]->y,
        (float)e->vn[3]->x, (float)e->vn[3]->y,
        order_palette[H2D_GET_V_ORDER(quad_order)]);
    fprintf(fout, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" stroke=\"black\" stroke-width=\"0.3\" fill=\"none\"/>",
      (float)e->vn[0]->x, (float)e->vn[0]->y, (float)(e->vn[2]->x - e->vn[0]->x), (float)(e->vn[2]->y - e->vn[0]->y));
    fprintf(fout, "</g>\n");
  }

  //finish
  fprintf(fout, "</g></g>\n");
  fprintf(fout, "</svg>\n");
  fclose(fout);
}

void store_adaptivity_step(const char* input_name, const int iteration,
                           Space* space, const Image2d& image_sln, const int image_width, const int image_height,
                           const std::vector<double>& sln_vec, 
                           const std::vector<ElementToRefine>& refinements,
                           TimeStatistics& timing) {
  //store refinements of a mesh
  store_refinements(input_name, iteration, refinements);

  //store DOFs
  store_solution_vector(input_name, iteration, sln_vec);

  //store mesh
  //store_mesh(input_name, iteration, sln);
  store_mesh_orders(input_name, iteration, image_width, image_height, space);

  //store result (store only every second, save time)
  store_solution(input_name, iteration, image_sln);
}

void store_error_estimate(int iteration, const char* input_name, double err_h1, double err_mse, double err_mad, double err_sigma, int dofs) {
  stringstream out_fname;
  out_fname << input_name << ".conv.csv";
  FILE* fconv = fopen(out_fname.str().c_str(), "at");
  if (iteration == HIMG_FIRST_ITERATION)
    fprintf(fconv, "#iteration\tMSE\tMAD\tSIGMA\th1_error\tdofs\n");
  fprintf(fconv, "%d\t%g\t%g\t%g\t%g\t%d\n", iteration, err_mse, err_mad, err_sigma, err_h1, dofs);
  fclose(fconv);
}

void store_timing(int iteration, const char* input_name, const TimeStatistics& timing) {
  stringstream out_fname;
  out_fname << input_name << ".time.csv";
  FILE* fconv = fopen(out_fname.str().c_str(), "at");
  if (iteration == HIMG_FIRST_ITERATION) {
    fprintf(fconv, "#all times in seconds\n");
    fprintf(fconv, "#iteration\tsolving\tmesh_refining\tadaptivity\trasterization\n");
  }
  fprintf(fconv, "%d\t%.01f\t%.01f\t%.01f\t%.01f\n", iteration, (float)timing.solving, (float)timing.mesh_refining, (float)timing.adaptivity, (float)timing.rasterization);
  fclose(fconv);
}
