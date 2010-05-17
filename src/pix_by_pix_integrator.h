#ifndef __PIX_BY_PIX_INTEGRATOR_H
#define __PIX_BY_PIX_INTEGRATOR_H

struct Vector2D { ///< 2D coordinate.
  double x, y;
  Vector2D() {};
  Vector2D(double x, double y) : x(x), y(y) {};
};

/// Integration of an element over an image, pixel by pixel.
/// Image in considered as a continuous function.
class PixByPixInt {
protected:
  const Image2d* image;
  const Quad1D& quad; ///< Used 1D quadrature.
  double* gip_upos_h; ///< A storage of positions of GIP in the horizontal direction transformed to unity interval. Valid contents only during call of integrate_element.
  double* gip_upos_v; ///< A storage of positions of GIP in the vertical direction transformed to unity interval. Valid contents only during call of integrate_element.

public:
  /// \brief Value evaluation callback function.
  /// \param diff_rescale A multiplicative coefficient that converts d/dx to either physical domain (physical_domain is true) or to referemce domain (physical_domain is false).
  /// \param element A pointer to element. Cannot be const since get functions of Shapeset require non-const element.
  /// \return Evaluted value.
  typedef double (*EvalValFunc)(const Vector2D& ref, const Vector2D& phys, const Vector2D& diff_rescale,
    const Image2d* image, Element* element, void* pars_ptr);

  PixByPixInt(const Image2d* image);
  virtual ~PixByPixInt();

  const Image2d* get_image() const { return image; }; ///< Returns an image.

  virtual double integrate_element(Element* element, int func_order_h, int func_order_v, EvalValFunc eval_value, void* params, bool physical_domain) const; ///< Integrates an element per pixel. Order of the other function is given by \param order_function.
};

///> Integration of an element over an image, pixel by pixel.
///> Image is considered as a discrete function with values at integer locations.
class PixByPixIntComb : public PixByPixInt {
public:
  PixByPixIntComb(const Image2d* image) : PixByPixInt(image) {};

  virtual double integrate_element(Element* element, int func_order_h, int func_order_v, EvalValFunc eval_value, void* params, bool physical_domain) const; ///< Integrates an element per pixel. Order of the other function is given by \param order_function.
};

#endif

