//
// Adaptive hp-FEM for images
//
#include "himg.h"

#include "shared/bitmap.h"
#include "shared/input_params.h"

#include "parameters.h"
#include "pix_by_pix_integrator.h"
#include "adapt_storage.h"
#include "h1_adapt_image.h"
#include "h1_image_proj_based_selector.h"

using namespace std;
using namespace RefinementSelectors;

/***** Global variables *****/ 
Parameters params; ///< HIMG input parameters.
Image2dExtrapolated image;  ///< An input image. Stored as a global variable due to the definition of the exact solution.
PixByPixInt pix_by_pix_int(&image); ///< A pixel-by-pixel integrator used for projections and error evaluation

OrderView oview("PolyOrders", 420, 0, 400, 400); ///< An order view.
ScalarView sview("Solution"); ///< A scalar view.

// function returning values and derivatives of the continuous image representation
double get_image_sample(double x, double y, double &dx, double &dy)
{
  return image.get_sample(x, y, dx, dy);
}

/***** Boundary conditions for Hermes *****/
BCType bc_mesh_1(int marker) {
  return BC_NONE;
}

/***** Forms *****/
template<typename T>
T bilinear_form(int n, double *wt, Func<T> *u, Func<T> *v, Geom<T> *e, ExtData<T> *data) {
  return int_u_v<T, T>(n, wt, u, v) + int_grad_u_grad_v<T, T>(n, wt, u, v);
}

struct EvalParams {
  Shapeset* shapeset; ///< Shapeset.
  const int shape_inx; ///< Shape index.
  EvalParams(Shapeset* shapeset, int shape_inx) : shapeset(shapeset), shape_inx(shape_inx) {};
};

static double eval_H1_lin_form(const Vector2D& ref, const Vector2D& phys, const Vector2D& diff_coef, const Image2d* image, Element* element, void* pars_ptr) {
  EvalParams& pars = *((EvalParams*)pars_ptr);

  //get value of a shape function
  double v_val_phys = pars.shapeset->get_fn_value(pars.shape_inx, ref.x, ref.y, 0);
  double v_dx_ref = pars.shapeset->get_dx_value(pars.shape_inx, ref.x, ref.y, 0);
  double v_dy_ref = pars.shapeset->get_dy_value(pars.shape_inx, ref.x, ref.y, 0);
  
  //transform derivates to the physical domain
  double v_dx_phys = v_dx_ref * diff_coef.x;
  double v_dy_phys = v_dy_ref * diff_coef.y;

  //get image values
  double img_val, img_dx, img_dy;
  img_val = image->get_sample(phys.x, phys.y, img_dx, img_dy);

  //calculate result
  return (img_val * v_val_phys +
    img_dx * v_dx_phys +
    img_dy * v_dy_phys);
}

scalar linear_form(int point_cnt, double *weights, Func<double> *values_v, Geom<double> *geometry, ExtData<scalar> *values_fnc_ext, Element* element, Shapeset* shape_set, int shape_inx) {
  //obtain order of the element
  int quad_order = shape_set->get_order(shape_inx);
  EvalParams params(shape_set, shape_inx);
  double value = pix_by_pix_int.integrate_element(element, H2D_GET_H_ORDER(quad_order), H2D_GET_V_ORDER(quad_order), eval_H1_lin_form, &params, true);
  return value;
}

Ord linear_form(int point_cnt, double *weights, Func<Ord> *values_v, Geom<Ord> *geometry, ExtData<Ord> *values_fnc_ext, Element* element, Shapeset* shape_set, int shape_inx) {
  return Ord(image.get_interpolation_h_order() + image.get_interpolation_v_order());
}

/****** Initialization ******/
/// Initializes an image.
bool init_image(const char* filename, Image2d& image) {
  //load image
  trace("Loading image \"%s\".", filename);
  Bitmap bmp;
  try { bmp.LoadFromFile(filename); }
  catch (runtime_error& err) { error("Loading of image \"%s\" failed: %s", filename, err.what()); }
  verbose("Loaded image of %d x %d pixels", bmp.GetWidth(), bmp.GetHeight());

  // convert to gray scale
  image.create_from_intensity(&bmp);

  return true;
}

/// Initializes a mesh.
bool init_mesh(Image2d& image, Mesh& mesh) {
  int mesh_width = image.get_width(), mesh_height = image.get_height(); 
  int vertex_num = 4, tria_num = 0, quad_num = 1, marker_num = 4;
  double2 vertex_array[4] = {{0, 0}, {mesh_width, 0}, {mesh_width, mesh_height}, {0, mesh_height}};
  int4 *tria_array = NULL;
  int5 quad_array[1] = {{0, 1, 2, 3, 0}};
  int3 marker_array[4] = {{0, 1, 3}, {1, 2, 2}, {2, 3, 4}, {3, 0, 1}};
  mesh.create(vertex_num, vertex_array, tria_num, tria_array, quad_num, quad_array, marker_num, marker_array);
  verbose("Mesh geometry: (0,%d)x(0,%d).", mesh_width, mesh_height);
  return true;
}

/// Prints info.
void print_info() {
  params.report_settings();
}

/// Initialized viewers.
bool init_viewers() {
  sview.set_palette(H2DV_PT_GRAYSCALE); 
  return true;
}

/// Visualizes the result.
void visualize(int iteration, Space& space, Solution& sln_coarse) {
  trace("Visualization and data backup");
  stringstream sview_str;
  sview_str << "Solution after iter " << iteration;
  sview.set_title(sview_str.str().c_str());
  sview.show(&sln_coarse);

  stringstream oview_str;
  oview_str << "Poly orders after iter " << iteration;
  oview.set_title(oview_str.str().c_str());
  oview.show(&space);
}

/***** Major functions *****/

/// Initializes computation.
bool do_init(Mesh** mesh_out, H1Space** space_out, H1Shapeset** shapeset_out, WeakForm** wf_out) {
  srand(0);

  //initialize mesh
  Mesh* mesh = new Mesh();
  if (!init_image(params.image_filename.c_str(), image))
    return false;
  if (!init_viewers())
    return false;
  init_mesh(image, *mesh);

  //prepare shapeset and space
  H1Shapeset* shapeset = new H1Shapeset();
  H1Space* space = new H1Space(mesh, shapeset);

  //prepare space
  space->set_bc_types(bc_mesh_1);
  space->set_uniform_order(params.initial_poly_degree);
  space->assign_dofs();  // enumerate basis functions

  // initialize the weak formulation
  WeakForm* wf = new WeakForm(1);
  wf->add_biform(0, 0, bilinear_form, bilinear_form, H2D_SYM);
  wf->add_liform(0, linear_form, linear_form, 0, H2D_ANY);

  // return data
  *mesh_out = mesh;
  *space_out = space;
  *shapeset_out = shapeset;
  *wf_out = wf;
  return true;
}

/// Cleans resources.
void do_cleanup(Mesh* mesh, H1Space* space, H1Shapeset* shapeset, WeakForm* wf) {
  if (wf != NULL)
    delete wf;
  if (space != NULL)
    delete space;
  if (shapeset != NULL)
    delete shapeset;
  if (mesh != NULL)
    delete mesh;
}

/// Processes an image. Main processing loop.
void do_adaptivity(Mesh& mesh, H1Space& space, H1Shapeset& shapeset, WeakForm& wf) {
  //prepare structures
  PrecalcShapeset pss(&shapeset);
  UmfpackSolver solver;
  Solution sln_coarse;

  //initialize selector
  H1ImageProjBasedSelector selector(&pix_by_pix_int, params.cand_list, params.conv_exp, H2DRS_DEFAULT_ORDER, &shapeset);
  
  //set weights to equal, otherwise it will try to increase order too much
  selector.set_error_weights(1.0, 1.0, 1.0);
  
  //adaptivity refinement
  TimePeriod cpu_time;
  TimeStatistics timing;
  std::vector<ElementToRefine> refinements;
  int iteration = HIMG_FIRST_ITERATION;
  bool done = false;
  do {
    info("!---- Adaptivity step %d ---------------------------------------------", iteration);

    // build the coarse solution
    trace("Assembling and solving.");
    cpu_time.tick(H2D_SKIP);
    LinSystem sys(&wf, &solver);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(1, &sln_coarse);
    timing.solving = cpu_time.tick().last();

    // build reference solution
    trace("Copying mesh.");
    cpu_time.tick(H2D_SKIP);
    Mesh ref_mesh;
    ref_mesh.copy(space.get_mesh());
    Solution sln_exact;
    sln_exact.set_exact(&ref_mesh, get_image_sample); //mesh will be refined while processing
    timing.mesh_refining = cpu_time.tick().last();

    // rasterize solution
    trace("Rasterizing solution.");
    cpu_time.tick(H2D_SKIP);
    Image2d image_sln;
    image_sln.create_from_solution(image.get_width(), image.get_height(), &sln_coarse);
    timing.rasterization = cpu_time.tick().last();
    
    // calculate element errors and total error estimate
    trace("Error calculation.");
    cpu_time.tick(H2D_SKIP);
    H1AdaptImage hp(&pix_by_pix_int, &space);
    hp.set_solutions(&sln_coarse, &sln_exact);
    double err_h1 = hp.calc_error();
    double err_mse = Image2d::evaluate_mse(image, image_sln);
    double err_mad = Image2d::evaluate_mad(image, image_sln);
    double err_sigma = Image2d::evaluate_sigma(image, image_sln);
    double timing_err_calculation = cpu_time.tick().last(); //time is used in the next iteration because it calculates data for adaptivity
    verbose("H1 error estimation: %g", err_h1);
    verbose("MSE: %g", err_mse);
    verbose("MAD: %g", err_mad);
    verbose("sigma: %g", err_sigma);

    //store result
    trace("Storing results.");
    store_error_estimate(iteration, params.image_filename.c_str(), err_h1, err_mse, err_mad, err_sigma, space.get_num_dofs());
    vector<double> sln_vector;
    sys.get_solution_vector(sln_vector);
    store_adaptivity_step(params.image_filename.c_str(), iteration,
      &space, image_sln, image.get_width(), image.get_height(),
      sln_vector, refinements, timing);

    //store timing since this iteration is done (adaptivity belongs to the next iteration)
    store_timing(iteration, params.image_filename.c_str(), timing);

    //visualize
    if (params.visualize)
      visualize(iteration, space, sln_coarse);

    //adaptivity step
    trace("Adaptivity step.");
    if (err_sigma < params.stop_error) 
    {
      verbose("Done due to an SIGMA below threshold.");
      done = true;
    }
    else {
      cpu_time.tick(H2D_SKIP);
      bool not_refined = hp.adapt(&selector, params.threshold, params.adapt_strategy, params.mesh_regularity);
      if (not_refined) {
        not_refined = true;
      }
      int ndofs = space.assign_dofs();
      timing.adaptivity = timing_err_calculation + cpu_time.tick().last(); //adaptivity counts to the next step (adaptivity that cuased the step)
      refinements = hp.get_last_refinements();

      if (params.max_ndof > 0 && ndofs >= params.max_ndof)
      {
        verbose("Done due to reaching of the maximum DOF.");
        done = true;
      }
    }

    iteration++;
  } while (!done);

  View::wait();
}

/// Main entry point.
int main(int argc, char *argv[]) {
  cout << "-------------------------------------------------" << endl;
  cout << "        We are HIMG. Resistance is futile.       " << endl;
  cout << "  Pavel Solin & David Andrs & Ivo Hanak, 2010.   " << endl;
  cout << "-------------------------------------------------" << endl;

  // parse command line and print settings
  try {
    InputParams input;
    params.set_defaults(input);
    input.parse_command_line(argc, argv);
    params.load(input);
  }
  catch (InputParamsError& e) {
    params.print_legend();
    error("%s", e.what());
    return false;
  }
  print_info();

  //compute
  Mesh* mesh = NULL;
  H1Space* space = NULL;
  H1Shapeset* shapeset = NULL;
  WeakForm* wf = NULL;
  if (do_init(&mesh, &space, &shapeset, &wf)) {
    do_adaptivity(*mesh, *space, *shapeset, *wf);
  }
  do_cleanup(mesh, space, shapeset, wf);

  return 0;
}
