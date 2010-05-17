#ifndef __HIMG_ADAPT_STORAGE_H
#define __HIMG_ADAPT_STORAGE_H

class Mesh;

struct TimeStatistics { ///< Time statistics.
  double solving;
  double mesh_refining;
  double adaptivity;
  double rasterization;
  TimeStatistics() : solving(0), mesh_refining(0), adaptivity(0), rasterization(0) {};
};

extern void store_adaptivity_step(const char* input_name, const int iteration, Space* space, const Image2d& image_sln, const int image_width, const int image_height, const std::vector<double>& sln_vec, const std::vector<ElementToRefine>& refinements, TimeStatistics& timing); ///< Stores adaptivity step into a Raw files
extern void store_error_estimate(int iteration, const char* input_name, double err_h1, double err_mse, double err_mad, double err_sigma, int dofs); ///< Stores error estimation.
extern void store_timing(int iteration, const char* input_name, const TimeStatistics& timing); ///< Stores timing.

#endif
