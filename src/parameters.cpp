#include "himg.h"
#include "shared/input_params.h"
#include "parameters.h"

using namespace std;
using namespace RefinementSelectors;

void Parameters::set_defaults(InputParams& input) {
  input.set_default("cl", get_cand_list_str(H2D_HP_ANISO));
  input.set_default("se", 0.002); //a maximum st. deviation is less than 1/255
  input.set_default("t", 0.2);
  input.set_default("mdof", 0);
  input.set_default("strat", 1);
  input.set_default("ipoly", 1);
  input.set_default("mr", -1);
  input.set_default("ce", 1.0);
  input.set_default("vis", false);
}

void Parameters::load(InputParams& input) {
  image_filename = input.get_string("i");
  stop_error = input.get_double("se");
  threshold = input.get_double("t");
  max_ndof = input.get_int32("mdof");
  adapt_strategy = input.get_int32("strat");
  initial_poly_degree = input.get_int32("ipoly");
  mesh_regularity = input.get_int32("mr");
  conv_exp = input.get_double("ce");
  visualize = input.get_bool("vis");

  CandList used_cand_lists[] = { H2D_H_ISO, H2D_H_ANISO, H2D_P_ISO, H2D_P_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO };
  string cand_list_str = input.get_string("cl");
  int inx_found = -1;
  for(int i = 0; i < 8; i++)
    if (cand_list_str.compare(get_cand_list_str(used_cand_lists[i])) == 0)
      inx_found = i;

  if (inx_found < 0)
    throw InputParamsError("Unknown candidate list: " + cand_list_str);

  cand_list = used_cand_lists[inx_found];
}

void Parameters::report_settings()
{
  verbose("!Parameters:");
  verbose(" - Stopping criterion:");
  if (max_ndof == 0)
    verbose("     Max. ndof: infinity");
  else
    verbose("     Max. ndof: %d", max_ndof);
  verbose("     Sigma error: %g", stop_error);

  verbose(" - Adaptivity parameters:");
  verbose("     Convergence exponent: %g", conv_exp);
  verbose("     Threshold: %g", threshold);
  verbose("     Mesh regularity: %d", mesh_regularity);
  verbose("     Adaptivity strategy: %d", adapt_strategy);
  verbose("     Predef. candidate list: %s", get_cand_list_str(cand_list)); 

  verbose(" - Other:");
  verbose("     Input filename: %s", image_filename.c_str());
  verbose("     Initial poly degree: %d", initial_poly_degree);
}

#define ERROR_CHECK(__cond, __action) if (__cond) { __action; return false; }  
bool Parameters::is_valid()
{
  error_if(image_filename.empty(), "No input filename given.");
  error_if(threshold <= 0 || threshold > 1, "Threshold is out of range.");
  return true;
}
#undef ERROR_CHECK

void Parameters::print_legend() {
  cout << "Quick start:" << endl;
  cout << "  Type './run name' where name = 'squares', 'diag', 'lena', or 'sat'." << endl;
  cout << "Basic usage:" << endl;
  cout << "  himg -i image.ppm" << endl;
  cout << "Optional parameters:" << endl;
  cout << "  -i X               Input image filename (PPM required)." << endl;
  cout << "  -se X              Tolerance for standard deviation of image values (default 0.002)." << endl;
  cout << "  -mdof X            Maximum allowed number of DOF (default infinity)." << endl;
  cout << "  -strat X           H2D adaptive strategy (default 1)." << endl;
  cout << "  -t X               H2D adaptivity error threshold (default 0.2)." << endl;
  cout << "  -ipoly X           H2D initial polynomial degree (default 1)." << endl;
  cout << "  -mr X              H2D mesh regularity (default -1)." << endl;
  cout << "  -ce X              H2D convergence exponent (default 1.0)." << endl;
  cout << "  -cl X              H2D refinement candidate list (default HP_ANISO)." << endl;
  cout << "  -vis true          Enable visualization (default false)." << endl;
}
