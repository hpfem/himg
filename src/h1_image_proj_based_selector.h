#ifndef __H1_IMAGE_PROJ_BASED_SELECTOR_H
#define __H1_IMAGE_PROJ_BASED_SELECTOR_H

/// Projection based selector of candidates in H1 space for image.
class H1ImageProjBasedSelector : public RefinementSelectors::H1ProjBasedSelector {
protected: //data
  const PixByPixInt* pix_by_pix_int; ///< Pixel-by-pixel integrator.

public: //API
  H1ImageProjBasedSelector(const PixByPixInt* pix_by_pix_int, RefinementSelectors::CandList cand_list = RefinementSelectors::H2D_HP_ANISO, double conv_exp = 1.0, int max_order = H2DRS_DEFAULT_ORDER, H1Shapeset* user_shapeset = NULL)
  : RefinementSelectors::H1ProjBasedSelector(cand_list, conv_exp, max_order, user_shapeset), pix_by_pix_int(pix_by_pix_int) {};

protected: //evaluators
  /// Parameters for evaluation of the right side.
  struct EvalRHSParams {
    const ElemSubTrf& sub_trf;
    Shapeset* shapeset;
    int shape_inx;
  };

  /// Parameters for evaluation of the projection error.
  struct EvalProjErrorParams {
    const ElemSubTrf& sub_trf;
    Shapeset* shapeset;
    const ElemProj& elem_proj;
  };

  /// Prevents calculation of values of shape function at GIP for all transformations.
  /**  Overriden function. For details, see ProjBasedSelector::precalc_shapes(). */
  virtual void precalc_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const std::vector<ShapeInx>& shapes, const int max_shape_inx, TrfShape& svals) {};

  /// Prevents calculation of values of orthogonalized shape function at GIP for all transformations.
  /**  Overriden function. For details, see ProjBasedSelector::precalc_ortho_shapes(). */
  virtual void precalc_ortho_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const std::vector<ShapeInx>& shapes, const int max_shape_inx, TrfShape& svals) {};

  /// Fill a list of candidates.
  /** Removes candidiate that have elements smaller than 1 px in any direction.
   *  Overriden function. For details, see OptimumSelector::create_candidates(). */
  virtual void create_candidates(Element* e, int quad_order, int max_ha_quad_order, int max_p_quad_order);

  /// Evaluates a value of the right-hande side in a subdomain.
  /**  Overriden function. For details, see ProjBasedSelector::evaluate_rhs_subdomain(). */
  virtual scalar evaluate_rhs_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemSubShapeFunc& sub_shape);

  /// Evaluates an squared error of a projection of an element of a candidate onto subdomains.
  /**  Overriden function. For details, see ProjBasedSelector::evaluate_error_squared_subdomain(). */
  virtual double evaluate_error_squared_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj);

  /// Evaluates a value of an integral at right side. Used by integrator, see H1ImageProjBasedSelector::pix_by_pix_int.
  static double eval_right_side(const Vector2D& sub_ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr);

  /// Evaluate a value of a projection error integral. Used by integrator, see H1ImageProjBasedSelector::pix_by_pix_int.
  static double eval_proj_error(const Vector2D& sub_ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr);
};

#endif
