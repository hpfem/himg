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
  cout << "Basic use:" << endl;
  cout << "  himg -i image.ppm" << endl;
  cout << "Supported parameters:" << endl;
  cout << "  -i X               an input image filename. PPM required." << endl;
  cout << "  -se X              a stop standard deviation of image values" << endl;
  cout << "  -t X               an error threshold for adaptivity" << endl;
  cout << "  -mdof X            a maximum number of DOFs" << endl;
  cout << "  -strat X           an adaptivity strategy" << endl;
  cout << "  -ipoly X           an initial polynomial degree" << endl;
  cout << "  -mr X              a mesh regularity" << endl;
  cout << "  -ce X              a convergence exponent" << endl;
  cout << "  -cl X              a candidate list (HP_ANISO is the default)" << endl;
  cout << "  -vis true          enable visualization" << endl;
}
