#ifndef __PARAMETERS_H
#define __PARAMETERS_H

/* Adaptivity types */
#define ADAPT_TYPE_HP 0
#define ADAPT_TYPE_H 1
#define ADAPT_TYPE_P 2

class InputParams;

struct Parameters
{
public:
  std::string image_filename;	///< file name of the image to process

  int initial_poly_degree; ///< initial polynomial degree of elements

  double stop_error;		///< adaptivity process stops when error wrt. exact solution in H1 norm is less than this number
  double threshold; ///< a quantitative parameter of the adapt(...) function and it has different meanings for various adaptive strategies
  int max_ndof; ///< a maximum number of NDOFs allowed. 0 for infinity
  int adapt_strategy; ///< an adaptive strategy
  RefinementSelectors::CandList cand_list; ///< an type of adaptation
  double conv_exp; ///< A convergence exponent.

  int mesh_regularity; ///< mesh regularity

  Parameters() {};
  void set_defaults(InputParams& input); ///< Sets default values.
  void load(InputParams& input); ///< Loads values. Throws an exception on failure.

  bool is_valid(); ///< Returns true if valid.

  void report_settings(); ///< Report stored settings.

  void print_legend(); ///< Prints legent about parameters.
};


#endif
