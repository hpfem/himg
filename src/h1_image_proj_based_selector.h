#ifndef __H1_IMAGE_PROJ_BASED_SELECTOR_H
#define __H1_IMAGE_PROJ_BASED_SELECTOR_H

/// Projection based selector of candidates in H1 space for image.
class H1ImageProjBasedSelector : public RefinementSelectors::H1ProjBasedSelector {
protected: //data
  const PixByPixInt* pix_by_pix_int;

public: //API
  H1ImageProjBasedSelector(const PixByPixInt* pix_by_pix_int, RefinementSelectors::CandList cand_list = RefinementSelectors::H2D_HP_ANISO, double conv_exp = 1.0, int max_order = H2DRS_DEFAULT_ORDER, H1Shapeset* user_shapeset = NULL)
  : RefinementSelectors::H1ProjBasedSelector(cand_list, conv_exp, max_order, user_shapeset), pix_by_pix_int(pix_by_pix_int) {};

protected: //evaluators
  struct EvalRHSParams { ///< Parameters for evaluation of the right side.
    const ElemSubTrf& sub_trf;
    Shapeset* shapeset;
    int shape_inx;
  };
  struct EvalProjErrorParams { ///< Parameters for evaluation of the projection error.
    const ElemSubTrf& sub_trf;
    Shapeset* shapeset;
    const ElemProj& elem_proj;
  };

  /// Fill a list of candidates.
  virtual void create_candidates(Element* e, int quad_order, int max_ha_quad_order, int max_p_quad_order);

  virtual scalar evaluate_rsh_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, int shape_inx); ///> Evaluate a single value of the right side for a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.
  virtual double evaluate_error_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj); ///> Evaluate an error of a projection on a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.

  static double eval_right_side(const Vector2D& sub_ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr); ///> Evaluate a value of an integral at right side.
  static double eval_proj_error(const Vector2D& sub_ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr); ///> Evaluate a value of a projection error integral.
};

#endif
