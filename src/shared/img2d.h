//
// 2d image
//

#ifndef _IMG2D_H_
#define _IMG2D_H_

#include "bitmap.h"

// Image2d class ////////////////////////////////////////////////////////////////////////
class Solution;

/// Image.
/// Image sampled using bilinear interpolation. Samples shifted to lower-left corner, edges are repeated.
class Image2d {
public:
  Image2d();
  virtual ~Image2d();

  bool create_from_intensity(Bitmap *bmp);
  bool create_from_png(const char* filename); ///< Creates an image from PNG
  bool create_from_image(const Image2d& img); ///< Creates a copy of an image
  bool create_from_image(const Image2d& img, const int scale); ///< Creates a copy a scaled copy of an image
  bool create_from_solution(int width, int height, Solution* sln, int scale_coef = 1); ///< creates the image from solution

  void store_to_png(const char* filename) const; ///< stores the image into a file. Range is given automatically.
  void store_to_png_clamped(const char* filename, float range_min, float range_max) const; ///< stores the image into a file. All values are clamped to specified range.

  int get_width() const { return width; };
  int get_height() const { return height; };
  float *get_data() { return data; };
  
  void calculate_range(float* range_min, float* range_max) const; ///< calculates range of stored values

  virtual double get_sample(double x, double y, double& dx, double& dy) const; ///< Returns image sample.
  virtual int get_interpolation_h_order() const { return 1; }; ///< Returns order of interpolation along the horizontal axis.
  virtual int get_interpolation_v_order() const { return 1; }; ///< Returns order of interpolation along the vertical axis.

  void damage_by_jpeg(int quality, int* output_size); ///< Damage the image with JPEG.
  void encode_jpeg(const char* out_filename, int quality); ///< Encode with JPEG.

public: //error evaluation
  static double evaluate_mse(const Image2d& accurate, const Image2d& approximate); ///< Evaluate MSE of approximate from accurate
  static double evaluate_mad(const Image2d& accurate, const Image2d& approximate); ///< Evaluate MAD of approximate from accurate
  static double evaluate_snr(const Image2d& accurate, const Image2d& approximate); ///< Evaluate SNR of approximate from accurate. SNR = 20*log_10 (signal/sqrt(MSE))
  static double evaluate_sigma(const Image2d& accurate, const Image2d& approximate); ///< Evaluate standard deviation.
  static double evaluate_inv_snr(const Image2d& accurate, const Image2d& approximate); ///< Evaluate inverted SNR of approximate from accurate. SNR = 20*log_10 (signal/sqrt(MSE)). Used when error should be decreasing as quality increases.

protected:
  int width;	 ///< width of the bitmap (in pixels)
  int height;	 ///< height of the bitmap (in pixels)
  float *data; ///< 2D array of values. First value is left-lower corner.

  void decode_jpeg(const char* in_filename); ///< Decode with JPEG

  void flip_y(); ///< flips image along the Y-axis
  bool is_little_endian() const; ///< returns true if machine is little endian
  float get_pixel(int x, int y) const; ///< returns value of a pixel, x must be in (0 .. Width - 1); y must be in (0 .. Height - 1)
};

/// Image.
/// Image sampled using Catmull-Rom interpolation. Samples shifted to lower-left corner, edges are repeated.
class Image2dCatmullRom : public Image2d {
public:
  Image2dCatmullRom() : Image2d() {};

  virtual double get_sample(double x, double y, double& dx, double& dy) const; ///< Returns image sample.
  virtual int get_interpolation_h_order() const { return 3; }; ///< Returns order of interpolation along the horizontal axis.
  virtual int get_interpolation_v_order() const { return 3; }; ///< Returns order of interpolation along the vertical axis.
};

/// Image.
/// Image sampled using bilinear interpolation. Samples shifted to lower-left corner, edges are extrapolated using linear extrapolation.
class Image2dExtrapolated : public Image2d {
public:
  Image2dExtrapolated() : Image2d() {};

  virtual double get_sample(double x, double y, double& dx, double& dy) const; ///< Returns image sample.
protected:
  double get_pixel_extrapolated(int inx_x, int inx_y) const; ///< Returns a value of a pixel. Extrapolates if neccessary.
};

/// Image.
/// Image not interpolated. Samples shifted to lower-left corner, edges are extrapolated using linear extrapolation.
/// Derivation estimated using side-wise interpolations.
class Image2dExtrPoint : public Image2dExtrapolated {
public:
  Image2dExtrPoint() : Image2dExtrapolated() {};

  virtual double get_sample(double x, double y, double& dx, double& dy) const; ///< Returns image sample.
};


#endif
