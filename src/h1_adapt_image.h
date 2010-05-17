#ifndef __HIMG_H1_ADAPT_IMAGE_H
#define __HIMG_H1_ADAPT_IMAGE_H

/// HP adaptivity with candidate projection using evaluation instead of an orthonormal base.
class H1AdaptImage : public H1Adapt {
protected: //data
  const Image2d* image;
  const PixByPixInt* pix_by_pix_int;

public: //constructors
  H1AdaptImage(const PixByPixInt* pix_by_pix_int, H1Space* space);
  virtual ~H1AdaptImage();

  virtual void apply_refinements(std::vector<ElementToRefine>& elems_to_refine); ///< Applies refinements. Sorts refinement prir applying and stores refinements for later use.

protected: //error evaluation
  //virtual double calc_error(); ///< Calculates H1 error, sets solutions used by adaptivity, and prepares attributes 'errors' and 'esort'.

  /// Evaluates a square of an absolute error of an active element among a given pair of components.
  virtual scalar eval_error(biform_val_t bi_fn, biform_ord_t bi_ord,
                    MeshFunction *sln1, MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2,
                    RefMap *rv1,        RefMap *rv2,        RefMap *rrv1,        RefMap *rrv2);

  /// Evaluates a square of a norm of an active element in the reference solution among a given pair of components.
  /** Uses H1 form strictly.
   *  \todo Norm is a function of the image: an integral image can be used instead. */
  virtual scalar eval_norm(biform_val_t bi_fn, biform_ord_t bi_ord,
                   MeshFunction *rsln1, MeshFunction *rsln2, RefMap *rrv1, RefMap *rrv2);

  /// Evaluate error of a solution from the image using H1 metrics. Assumes physical domain.
  static double eval_element_error(const Vector2D& ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr);

  static double eval_element_norm(const Vector2D& ref, const Vector2D& phys, const Vector2D& diff_rescale, const Image2d* image, Element* element, void* pars_ptr); ///< Evaluate norm of the image using H1 metrics. Assumes physical domain.

protected: //adaptivity
  struct EdgeRange { ///< A range of an edge
    Range<double> x, y;
    EdgeRange(const Range<double>& x = Range<double>(), const Range<double>& y = Range<double>()) : x(x), y(y) {};
    bool is_in_closed(const EdgeRange& edge) { return x.is_in_closed(edge.x) && y.is_in_closed(edge.y); };
  };

  std::map<int, int> no_adapt_elem; ///< Elements that should not be adapted due to a neighbourhood.

  /// Returns true if a given element should be ignored and not processed through refinement selection.
  /** If element is too small to be refined (i.e., size 1x1 px), the element is skipped and
   *  neighbours larger than 1x1 px are added to the priority queue */
  virtual bool should_ignore_element(const int inx_element, const Mesh* mesh, const Element* element);

  /// Returns true if a given element can be refined using proposed refinement.
  /** This method prevent an element being refined if this caused an element of a size below 1x1 px. */
  virtual bool can_refine_element(Mesh* mesh, Element* e, bool refined, ElementToRefine& elem_ref);

  static bool is_elem_ref_ordered(const ElementToRefine& a, const ElementToRefine& b); ///< Compares whether the pair is ordered. Used to sort records for adaptivity.

  void add_neighbours_nonadapt(Mesh* mesh, const int element_id, const int start_vertex_id, const int end_vertex_id, const Range<double>& range_x, const Range<double>& range_y); ///< Adds elements that are on an edge of the element.

  ///< Adds neighbours of an element along an edge of a length 1px to priority queue for processing.
  void add_neighbours_1px_to_adapt(const Mesh* mesh, const Element* element);
};

#endif
