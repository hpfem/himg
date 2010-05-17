#include "himg.h"
#include "shared/img2d.h"
#include "pix_by_pix_integrator.h"

PixByPixInt::PixByPixInt(const Image2d* image) : image(image), quad(g_quad_1d_std) {
  //allocate buffer for positions
  const int max_num_gip = quad.get_max_order();
  gip_upos_h = new double[max_num_gip];
  gip_upos_v = new double[max_num_gip];
}

PixByPixInt::~PixByPixInt() {
  delete[] gip_upos_h;
  delete[] gip_upos_v;
}

double PixByPixInt::integrate_element(Element* element, int func_order_h, int func_order_v, EvalValFunc eval_value, void* params, bool physical_domain) const {
  int image_width = image->get_width(), image_height = image->get_height();

  //check dimensions
  Node *bl = element->vn[0];				// bottom left node
  Node *tr = element->vn[2];				// top right node
  if (bl->x < 0 || bl->x > image_width || tr->x < 0 || tr->x > image_width
    || bl->y < 0 || bl->y > image_height || tr->y < 0 || tr->y > image_height) {
    error("Element %d out of range of image.", element->id);
  }

  //prepare dimensions and conversion coefficients
  double element_width = tr->x - bl->x;
  double element_height = tr->y - bl->y;
  double conv_coef_dx, conv_coef_dy;
  if (physical_domain) {
    conv_coef_dx = 2 / element_width;
    conv_coef_dy = 2 / element_height;
  } else {
    conv_coef_dx = element_width / 2;
    conv_coef_dy = element_height / 2;
  }

  //prepare pixel-by-pixel Gauss quadrature
  int gip_order_h = image->get_interpolation_h_order() + func_order_h;
  int gip_order_v = image->get_interpolation_v_order() + func_order_v;
  error_if(gip_order_h > quad.get_max_order() || gip_order_v > quad.get_max_order(), "Order (H:%d; V:%d) requested but maximum sopported order is %d.", gip_order_h, gip_order_v, quad.get_max_order());
  const int num_gip_pts_h = quad.get_num_points(gip_order_h);
  const int num_gip_pts_v = quad.get_num_points(gip_order_v);

  //scale GIP locations from [-1; 1] to [0; 1]; modification of values due to change of the domain will be taken care of after integration
  const double2* gip_pts_h = quad.get_points(gip_order_h);
  for(int i = 0; i < num_gip_pts_h; i++)
    gip_upos_h[i] = (gip_pts_h[i][H2D_GIP1D_X] + 1) / 2;
  const double2* gip_pts_v = quad.get_points(gip_order_v);
  for(int i = 0; i < num_gip_pts_v; i++)
    gip_upos_v[i] = (gip_pts_v[i][H2D_GIP1D_X] + 1) / 2;

  //integrate in physical domain
  //for each pixel in the Y-axis
  double value = 0;
  double start_x = floor(bl->x);
  double y = floor(bl->y);
  while (y < tr->y) {
    //calculate coordinates for the Y-axis (consider partial pixels)
    double shift_y_phys = 0, size_y_phys = 1;
    if (y < bl->y) {
      shift_y_phys = bl->y - y;
      size_y_phys = ceil(bl->y) - bl->y;
      if (size_y_phys > (tr->y - bl->y))
        size_y_phys = tr->y - bl->y;
    } else if ((y+1) > tr->y) {
      shift_y_phys = 0;
      size_y_phys = tr->y - y;
    }
    double y_phys_start = y + shift_y_phys;

    //for each pixel in the X-axis
    double x = start_x;
    while (x < tr->x) {
      //calculate coordinates for the X-axis (consider partial pixels)
      double shift_x_phys = 0, size_x_phys = 1;
      if (x < bl->x) {
        shift_x_phys = bl->x - x;
        size_x_phys = ceil(bl->x) - bl->x;
        if (size_x_phys > (tr->x - bl->x))
          size_x_phys = tr->x - bl->x;
      } else if ((x+1) > tr->x) {
        shift_x_phys = 0;
        size_x_phys = tr->x - x;
      }
      double x_phys_start = x + shift_x_phys;

      //calculate quadrature
      double value_pixel = 0;
      for(int r = 0; r < num_gip_pts_v; r++) {
        double y_phys = y_phys_start + gip_upos_v[r] * size_y_phys;
        double y_ref = ((y_phys - bl->y) / (tr->y - bl->y)) * 2 - 1;
        double y_weight = gip_pts_v[r][H2D_GIP1D_W];

        for(int s = 0; s < num_gip_pts_h; s++) {
          const double2& gip_pt_h = gip_pts_h[s];
          double x_phys = x_phys_start + gip_upos_h[s] * size_x_phys;
          double x_ref = ((x_phys - bl->x) / (tr->x - bl->x)) * 2 - 1;
          double x_weight = gip_pts_h[s][H2D_GIP1D_W];

          //calculate
          value_pixel += x_weight * y_weight * eval_value(Vector2D(x_ref, y_ref),
            Vector2D(x_phys, y_phys),
            Vector2D(conv_coef_dx, conv_coef_dy),
            image, element, params);
        }
      }
      value += value_pixel * size_y_phys * size_x_phys;
      x++;
    }
    y++;
  }

  //convert to reference domain
  //weight of GIPs assumes a reference domain which size is 2x2.
  if (physical_domain) { 
    value /= 2*2;
  }
  else {
    value /= (element_width * element_height);
  }

  return value;
}

double PixByPixIntComb::integrate_element(Element* element, int func_order_h, int func_order_v, EvalValFunc eval_value, void* params, bool physical_domain) const {
  int image_width = image->get_width(), image_height = image->get_height();

  //check dimensions
  Node *bl = element->vn[0];				// bottom left node
  Node *tr = element->vn[2];				// top right node
  if (bl->x < 0 || bl->x > image_width || tr->x < 0 || tr->x > image_width
    || bl->y < 0 || bl->y > image_height || tr->y < 0 || tr->y > image_height) {
    debug_log("E element %d out of range of image\n", element->id);
    assert(false);
  }

  //prepare dimensions and conversion coefficients
  double element_width = tr->x - bl->x;
  double element_height = tr->y - bl->y;
  double conv_coef_dx, conv_coef_dy;
  if (physical_domain) {
    conv_coef_dx = 2 / element_width;
    conv_coef_dy = 2 / element_height;
  } else {
    conv_coef_dx = element_width / 2;
    conv_coef_dy = element_height / 2;
  }

  //integrate in physical domain
  //for each pixel in the Y-axis
  double value = 0;
  double start_x = ceil(bl->x);
  double y = ceil(bl->y);
  while (y < tr-> y) {
    double y_phys = y;

    //for each pixel in the X-axis
    double x = start_x;
    while (x < tr->x) {
      double x_phys = x;

      //calculate reference coordinates
      double y_ref = ((y_phys - bl->y) / element_height) * 2 - 1;
      double x_ref = ((x_phys - bl->x) / element_width) * 2 - 1;

      //evaluate a value
      value += eval_value(Vector2D(x_ref, y_ref),
        Vector2D(x_phys, y_phys),
        Vector2D(conv_coef_dx, conv_coef_dy),
        image, element, params);

      x++;
    }
    y++;
  }

  //convert to reference domain
  if (!physical_domain) {
    value = 2*2 * value / (element_width * element_height);
  }

  return value;
}

