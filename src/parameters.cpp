#include "himg.h"
#include "shared/input_params.h"
#include "parameters.h"

using namespace std;

void Parameters::set_defaults(InputParams& input) {
  input.set_default("se", 0.002); //a maximum st. deviation is less than 1/255
  input.set_default("t", 0.2);
  input.set_default("mdof", 0);
  input.set_default("strat", 1);
  input.set_default("ipoly", 1);
  input.set_default("mr", -1);
  input.set_default("ce", 1.0);
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
  cand_list = RefinementSelectors::H2D_HP_ANISO;
}

void Parameters::report_settings()
{
  verbose("!parameters");
  verbose(" finalization conditions");
  if (max_ndof == 0)
    verbose("   max. DOFs: infinity");
  else
    verbose("   max. DOFs: %d", max_ndof);
  verbose("   stop error: %g", stop_error);

  verbose(" adaptivity parameters");
  verbose("   convergence exponent: %g", conv_exp);
  verbose("   threshold: %g", threshold);
  verbose("   mesh regularity: %d", mesh_regularity);
  verbose("   adaptivity strategy: %d", adapt_strategy);
  verbose("   predef. candidate list: %s", get_cand_list_str(cand_list)); 

  verbose(" others");
  verbose("   input filename: %s", image_filename.c_str());
  verbose("   init. poly. degree on elements: %d", initial_poly_degree);
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
  cout << "  -i X               an image filename" << endl;
  cout << "  -se X              a stop standard deviation of image values" << endl;
  cout << "  -t X               an error threshold for adaptivity" << endl;
  cout << "  -mdof X            a maximum number of DOFs" << endl;
  cout << "  -strat X           an adaptivity strategy" << endl;
  cout << "  -ipoly X           an initial polynomial degree" << endl;
  cout << "  -mr X              a mesh regularity" << endl;
  cout << "  -ce X              a conversion exponent" << endl;
}
