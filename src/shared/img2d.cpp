#include <png.h>
#include <hermes2d.h>
#include <cmath>
#include "img2d.h"

#undef EXTERN
#ifndef HIMG_NO_JPEG
# include <jpeglib.h>
#endif

using namespace std;

// PNG reading utilities ////////////////////////////////////////////////////////////////
static void read_PNG_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
	FILE* file = (FILE*)png_get_io_ptr(png_ptr);
	fread(data, 1, length, file);
}

static void write_PNG_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
	FILE* file = (FILE*)png_get_io_ptr(png_ptr);
	fwrite(data, 1, length, file);
}

static void flush_PNG_data(png_structp png_ptr)
{
	FILE* file = (FILE*)png_get_io_ptr(png_ptr);
	fflush(file);
}

// Image2d class ////////////////////////////////////////////////////////////////////////

Image2d::Image2d() {
	data = NULL;
}

Image2d::~Image2d() {
	delete [] data;
}

bool Image2d::create_from_png(const char* filename) {
  //cleanup
  delete[] data; width = height = -1;
  data = NULL;

  //data  
  bool readed = false;
 	FILE* file = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_byte png_signature[8];
	png_bytep row = NULL;
  png_uint_32 dx, dy;
  int bitDepth, colorType, interlaceMethod, compressionMethod, filterMethod;
  int bytes_per_value = 1;
	bool little_endian = is_little_endian();

	//open file
	file = fopen(filename, "rb");
	if (file == NULL) {
		debug_log("Unable to find file \"%s\".", filename);
		goto quit;
	}

	//read signature
	fread(png_signature, 1, 8, file);
	if (png_sig_cmp(png_signature, 0, 8)) {
		debug_log("File \"%s\" is not PNG, PNG signature not found.", filename);
		goto quit;
	}

	//create read structures
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		debug_log("Unable to allocate PNG read structures.");
		goto quit;
	}
  info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		debug_log("Unable to allocate PNG info read structures.");
		goto quit;
	}

  //prepare read operations
  png_set_read_fn(png_ptr, file, read_PNG_data);

  //connect file to PNG and read info
  png_set_sig_bytes(png_ptr, 8); //skip signature
  png_read_info(png_ptr, info_ptr);
  png_get_IHDR(png_ptr, info_ptr, &dx, &dy,
	  &bitDepth, &colorType, &interlaceMethod, &compressionMethod, &filterMethod);

	//check data format
	if (interlaceMethod != PNG_INTERLACE_NONE) {
		debug_log("Interlaced PNGs not supported");
		goto quit;
	}
  if (colorType != PNG_COLOR_TYPE_GRAY) {
    debug_log("Non-gray PNGs not supported");
  }
  switch (bitDepth) {
    case 8:
      bytes_per_value = 1; break;
    case 16:
      bytes_per_value = 2;
      if (little_endian)
        png_set_swap(png_ptr); //PNG is Big-endian
      break;
    default:
      debug_log("Only 8-bit and 16-bit PNGs supported");
      goto quit;
  }
  
  //allocate space
  width = dx; height = dy;
  try {
    row = new png_byte[dx*dy*bytes_per_value];
    data = new float[dx*dy];
  }
  catch(bad_alloc) {
    debug_log("Unable to allocate memory");
    goto quit;
  }

  //read lines (store mirrored along the Y-axis)
  {
    float* tgt_row = data + (height - 1)*width;
    for(int i = 0; i < (int)dy; i++, tgt_row -= width) {
      //read row
      png_read_row(png_ptr, row, NULL);

      //decode row
      if (bytes_per_value == 1) { //8-bit
        unsigned char* src_row = (unsigned char*)row;
        for(int k = 0; k <(int)dx; k++)
          tgt_row[k] = src_row[k] / 255.0f;
      }
      else { //16-bit
        unsigned short* src_row = (unsigned short*)row;
        for(int k = 0; k <(int)dx; k++)
          tgt_row[k] = src_row[k] / 65535.0f;
      }
    }
  }

  //set readed flag
  readed = true;

quit:
	if (row != NULL)
		delete[] row;
  if (file != NULL)
    fclose(file);
	if (info_ptr != NULL)
		png_destroy_info_struct(png_ptr, &info_ptr);
	if (png_ptr != NULL)
		png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
	if (file != NULL)
		fclose(file);
  if (readed)
    return true;
  else {
    delete[] data;
    width = height = -1;
    return false;
  }
}

bool Image2d::create_from_intensity(Bitmap *bmp) {
  //cleanup
  delete[] data; width = height = -1;
	
  //get source data
	RGB *rgb = bmp->GetData();
	width = bmp->GetWidth();
	height = bmp->GetHeight();

  //allocate
  data = new float [width * height];
	float *sample = data;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			*sample = (float)(0.3f * rgb->r + 0.59f * rgb->g + 0.11f * rgb->b) / 255.0f;

			// move to the next element
			sample++;
			rgb++;
		}
	}

  //invert coordinate system
  flip_y();

  return false;
}

bool Image2d::create_from_solution(int width, int height, Solution* sln, int scale_coef) {
  //cleanup
  delete[] data; this->width = this->height = -1;
	
  //allocate
  double width_unscaled = width, height_unscaled = height;
  width *= scale_coef; height *= scale_coef;
  data = new float [width * height];
  this->width = width; this->height = height;

  //rasterize each active element
  double step = 1.0 / scale_coef;
  Element* e;
  Mesh* mesh = sln->get_mesh();
  for_all_active_elements(e, mesh) {
    //check dimensions
    Node *bl = e->vn[0];				// bottom left node
    Node *tr = e->vn[2];				// top right node
    error_if((bl->x < 0 || bl->x > width_unscaled || tr->x < 0 || tr->x > width_unscaled
      || bl->y < 0 || bl->y > height_unscaled || tr->y < 0 || tr->y > height_unscaled)
      , "Element %d out of range of image.", e->id);

    //for each pixel in element
    double start_x = ceil(bl->x * scale_coef) / scale_coef;
    double y = ceil(bl->y * scale_coef) / scale_coef;
    int inx_row = (int)(y * scale_coef) * width + (int)(start_x * scale_coef);
    while (y < tr->y) {
      double y_ref = ((y - bl->y) / (tr->y - bl->y)) * 2 - 1;
      double x = start_x;
      int inx = inx_row;
      while (x < tr->x) {
        double x_ref = ((x - bl->x) / (tr->x - bl->x)) * 2 - 1;
        data[inx] = (float)sln->get_ref_value(e, x_ref, y_ref);
        inx++;
        x += step;
      }
      y += step;
      inx_row += width;
    }
  }

  return false;
}

bool Image2d::create_from_image(const Image2d& img)
{
  //cleanup
  delete[] data; this->width = this->height = -1;

  //allocate
  data = new float [img.width * img.height];
  this->width = img.width; this->height = img.height;

  //copy data
  memcpy(data, img.data, img.width * img.height * sizeof(float));

  return true;
}

bool Image2d::create_from_image(const Image2d& img, const int scale)
{
  //cleanup
  delete[] data; this->width = this->height = -1;

  //allocate
  data = new float [img.width * img.height * scale*scale];
  this->width = img.width * scale; this->height = img.height * scale;

  //copy data
  int inx = 0;
  for(int i = 0; i < this->height; i++) {
    double y = i / (double)scale;
    for(int k = 0; k < this->width; k++, inx++) {
      double x = k / (double)scale;
      double dx, dy;
      data[inx] = (float)img.get_sample(x, y, dx, dy);
    }
  }

  return true;
}

void Image2d::store_to_png(const char* filename) const {
  //get range
  float range_min, range_max;
  calculate_range(&range_min, &range_max);

  //store
  store_to_png_clamped(filename, range_min, range_max);
}

void Image2d::store_to_png_clamped(const char* filename, float range_min, float range_max) const {
  //variables
  unsigned short* image_row;
  FILE* file = NULL;
  png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep* pRowPtrs = NULL;

  //allocate space for a row
  try { image_row = new unsigned short[width]; }
  catch(std::bad_alloc &) { error("Unable to allocate space for image row"); goto quit; }

  //open target file
	file = fopen(filename, "wb");
  error_if(file == NULL, "Cannot open file \"%s\" for writing", filename);

	//prepare structures
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  error_if(png_ptr == NULL, "Unable to allocate PNG write structures");
  info_ptr = png_create_info_struct(png_ptr);
  error_if(info_ptr == NULL, "Unable to allocate PNG info write structures");

	//set compression
	png_set_compression_level(png_ptr, Z_BEST_COMPRESSION);
	png_set_compression_mem_level(png_ptr, 8);
	png_set_compression_strategy(png_ptr, Z_DEFAULT_STRATEGY);
	png_set_compression_window_bits(png_ptr, 15);
	png_set_compression_method(png_ptr, 8);
	png_set_compression_buffer_size(png_ptr, 8192);

	//connect png to file and write info
	png_set_write_fn(png_ptr, file, write_PNG_data, flush_PNG_data);
  png_set_IHDR(png_ptr, info_ptr, width, height, 16, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_ptr, info_ptr);

  //PNG is Big-endian
	if (is_little_endian())
		png_set_swap(png_ptr); 

  //store rows
  {
    float range_length = range_max - range_min;
    if (range_length == 0)
      range_length = 0.0f;
    float* row = data + (height-1) * width;
    for(int i = 0; i < height; i++, row -= width) {
      //convert
      for(int k = 0; k < width; k++) {
        float value = row[k];
        if (value > range_max)
          image_row[k] = 0xFFFF;
        else if (value < range_min)
          image_row[k] = 0;
        else 
          image_row[k] = (unsigned short)floor(65535.0f * (value - range_min) / range_length);
      }

      //store
      png_write_row(png_ptr, (png_bytep)image_row);
    }
  }

  //finish
	png_write_flush(png_ptr);
	png_write_end(png_ptr, info_ptr);
  
quit:
  delete[] image_row;
	if (png_ptr != NULL || info_ptr != NULL)
		png_destroy_write_struct(&png_ptr, &info_ptr);
  if (file != NULL)
    fclose(file);
}

bool Image2d::is_little_endian() const {
	unsigned short value = 0xABCD;
	unsigned char* pValue = (unsigned char*)&value;
	if (*pValue == 0xAB)
		return false;
	else
		return true;
}

void Image2d::calculate_range(float* range_min, float* range_max) const {
  float* sample = data;
  *range_min = *range_max = *sample;
  int cnt = width * height;
  for(int i = 0; i < cnt; i++, sample++) {
    if (*sample > *range_max)
      *range_max = *sample;
    if (*sample < *range_min)
      *range_min = *sample;
  }
}

void Image2d::flip_y() {
  int flip_cnt = height / 2;
  float* row_a = data;
  float* row_b = data + (height-1)*width;

  for(int i = 0; i < flip_cnt; i++, row_a += width, row_b -= height) {
    for(int k = 0; k < width; k++) {
      float tmp = row_a[k];
      row_a[k] = row_b[k];
      row_b[k] = tmp;
    }
  } 
}

float Image2d::get_pixel(int x, int y) const  {
  error_if(x < 0 || x >= width || y < 0 || y >= height, "Pixel (%d, %d) is out of range.", x, y);
  return data[y * width + x]; 
} 

double Image2d::get_sample(double x, double y, double& dx, double& dy) const 
{
  int pixel_x0 = (int)floor(x), pixel_y0 = (int)floor(y);
  double x_diff = (double)(x - pixel_x0), y_diff = (double)(y - pixel_y0);
  int pixel_x1 = pixel_x0 + 1, pixel_y1 = pixel_y0 + 1;

  //create coordinates
  if (pixel_x0 < 0)
    pixel_x0 = 0;
  else if (pixel_x0 >= width)
    pixel_x0 = width - 1;
  if (pixel_x1 < 0)
    pixel_x1 = 0;
  else if (pixel_x1 >= width)
    pixel_x1 = width - 1;

  if (pixel_y0 < 0)
    pixel_y0 = 0;
  else if (pixel_y0 >= height)
    pixel_y0 = height - 1;
  if (pixel_y1 < 0)
    pixel_y1 = 0;
  else if (pixel_y1 >= height)
    pixel_y1 = height - 1;
 
  //calculate
  double value_A_0 = get_pixel(pixel_x0, pixel_y0);
  double value_A_1 = get_pixel(pixel_x1, pixel_y0);
  double value_A = value_A_0 * (1-x_diff) + value_A_1 * x_diff;
  double value_B_0 = get_pixel(pixel_x0, pixel_y1);
  double value_B_1 = get_pixel(pixel_x1, pixel_y1);
  double value_B = value_B_0 * (1-x_diff) + value_B_1 * x_diff;
  double value = value_A * (1-y_diff) + value_B * y_diff;
  double value_0 = value_A_0 * (1-y_diff) + value_B_0 * y_diff;
  double value_1 = value_A_1 * (1-y_diff) + value_B_1 * y_diff;
  dx = value_1 - value_0;
  dy = value_B - value_A;

  return value;
}

double Image2d::evaluate_mse(const Image2d& accurate, const Image2d& approximate)
{
  assert_msg(accurate.width == approximate.width && accurate.height == approximate.height, "Image dimensions do not match.");
  double mse_sum = 0.0;
  int px_cnt = accurate.width * accurate.height;
  for(int i = 0; i < px_cnt; i++)
  {
    float diff = accurate.data[i] - approximate.data[i];
    mse_sum += diff*diff;
  }
  return mse_sum / px_cnt;
}

double Image2d::evaluate_sigma(const Image2d& accurate, const Image2d& approximate)
{
  assert_msg(accurate.width == approximate.width && accurate.height == approximate.height, "Image dimensions do not match.");
  double sum = 0.0, mse_sum = 0.0;
  int px_cnt = accurate.width * accurate.height;
  for(int i = 0; i < px_cnt; i++)
  {
    float diff = accurate.data[i] - approximate.data[i];
    sum += diff;
    mse_sum += diff*diff;
  }
  return sqrt((mse_sum / px_cnt) - sqr(sum / px_cnt));;
}

double Image2d::evaluate_mad(const Image2d& accurate, const Image2d& approximate)
{
  assert_msg(accurate.width == approximate.width && accurate.height == approximate.height, "Image dimensions do not match.");
  float max_diff = 0.0f;
  int px_cnt = accurate.width * accurate.height;
  for(int i = 0; i < px_cnt; i++)
  {
    float diff = abs(accurate.data[i] - approximate.data[i]);
    if (diff > max_diff)
      max_diff = diff;
  }
  return max_diff;
}

double Image2d::evaluate_snr(const Image2d& accurate, const Image2d& approximate) {
  assert_msg(accurate.width == approximate.width && accurate.height == approximate.height, "Image dimensions do not match.");
  double mse_sum = 0.0;
  double sum = 0.0;
  int px_cnt = accurate.width * accurate.height;
  for(int i = 0; i < px_cnt; i++)
  {
    sum += accurate.data[i];
    float diff = abs(accurate.data[i] - approximate.data[i]);
    mse_sum += diff*diff;
  }
  double snr = (sum / px_cnt) / sqrt(mse_sum / px_cnt);
  return 20 * log10(snr);
}

double Image2d::evaluate_inv_snr(const Image2d& accurate, const Image2d& approximate) {
  return 1.0 / evaluate_snr(accurate, approximate);
}

#ifndef HIMG_NO_JPEG
void Image2d::damage_by_jpeg(int quality, int* output_size)
{
  //get temporary filename
  std::stringstream str;
  str << "HIMG_TEMP" << rand() << ".jpg";

  //encode
  encode_jpeg(str.str().c_str(), quality);

  //get size of the file
  FILE* file = fopen(str.str().c_str(), "rb");
  fseek(file, 0, SEEK_END);
  *output_size = ftell(file);
  fclose(file);

  //decode
  decode_jpeg(str.str().c_str());

  //delete file
  remove(str.str().c_str());
}

void Image2d::decode_jpeg(const char* in_filename)
{
  //variables
  FILE* file = NULL;
  JSAMPLE* image_row;

  //allocate temporary space
  try { image_row = new JSAMPLE[width]; }
  catch(std::bad_alloc &) { error("Unable to allocate space for image row"); }

  //initialize structures
 	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);

  //open input file
  file = fopen(in_filename, "rb");
  error_if(file != NULL, "Unable to open temporary file for reading");
  jpeg_stdio_src(&cinfo, file);

  //read header
  jpeg_read_header(&cinfo, TRUE);

  //start decompressing
  jpeg_start_decompress(&cinfo);

  //decompress
  float* row = data + (height-1) * width;
  for(int i = 0; i < height; i++, row -= width) {
    //read
    JSAMPROW row_ptr[1];
    row_ptr[0] = image_row;
    jpeg_read_scanlines(&cinfo, row_ptr, 1);

    //convert
    for(int k = 0; k < width; k++)
      row[k] = image_row[k] / (float)MAXJSAMPLE;
  }

  //finish
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(file);
  delete[] image_row;
}

void Image2d::encode_jpeg(const char* out_filename, int quality)
{
  error_if(width > 0 && height > 0, "Invalid size of image or image not initialized");

  //variables
  FILE* file = NULL;
  JSAMPLE* image_row;

  //allocate temporary space
  try { image_row = new JSAMPLE[width]; }
  catch(std::bad_alloc &) { error("Unable to allocate space for image row"); }

  //initialize structures
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

  //create output file
  file = fopen(out_filename, "wb");
  error_if(file != NULL, "Unable to open temporary file for writing");
  jpeg_stdio_dest(&cinfo, file);

  //set compression parameters
  cinfo.image_width = width; 	/* image width and height, in pixels */
	cinfo.image_height = height;
	cinfo.input_components = 1;	/* # of color components per pixel */
  cinfo.in_color_space = JCS_GRAYSCALE; /* colorspace of input image */
	jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, quality, false);

  //start compressing
  jpeg_start_compress(&cinfo, TRUE);

  //write lines
  float range_min = 0.0f, range_max = 1.0f;
  float range_length = 1.0f;
  float* row = data + (height-1) * width;
  for(int i = 0; i < height; i++, row -= width) {
    //convert
    for(int k = 0; k < width; k++) {
      float value = row[k];
      if (value > range_max)
        image_row[k] = MAXJSAMPLE;
      else if (value < range_min)
        image_row[k] = 0;
      else 
        image_row[k] = (JSAMPLE)floor( MAXJSAMPLE * (value - range_min) / range_length);
    }

    //store
    JSAMPROW row_ptr[1];
    row_ptr[0] = image_row;
    jpeg_write_scanlines(&cinfo, row_ptr, 1);
  }

  //finish
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
  fclose(file);
  delete[] image_row;
}
#else
void Image2d::damage_by_jpeg(int quality, int* output_size) { error("JPEG not enabled by CMAKE."); };
void Image2d::decode_jpeg(const char* in_filename) { error("JPEG not enabled by CMAKE."); };
void Image2d::encode_jpeg(const char* out_filename, int quality) { error("JPEG not enabled by CMAKE."); };
#endif

/* Catmull-Rom interoplation */
template<typename T>
static inline void catmul_rom_f(const T t, const T p0, const T p1, const T p2, const T p3, T& value, T& value_dt) {
  double coefs[4] = {
    (2*p1),
    (-p0 + p2),
    (2*p0 - 5*p1 + 4*p2 - p3),
    (-p0 + 3*p1 - 3*p2 + p3) };

  value = (((coefs[3]
    * t + coefs[2])
    * t + coefs[1])
    * t + coefs[0]) / 2;
  value_dt = ((coefs[3]
    * 1.5 * t + coefs[2])
    * 2 * t + coefs[1]) /2;
}

double Image2dCatmullRom::get_sample(double x, double y, double& dx, double& dy) const {
  int pixel_x0 = (int)floor(x), pixel_y0 = (int)floor(y);
  double x_diff = (double)(x - pixel_x0), y_diff = (double)(y - pixel_y0);

  //create coordinates
  int pixel_x = pixel_x0 - 1, pixel_y = pixel_y0 - 1;
  int pixel_xs[4], pixel_ys[4];
  for(int i = 0; i < 4; i++, pixel_x++, pixel_y++) {
    if (pixel_x < 0)
      pixel_xs[i] = 0;
    else if (pixel_x >= width)
      pixel_xs[i] = width-1;
    else
      pixel_xs[i] = pixel_x;

    if (pixel_y < 0)
      pixel_ys[i] = 0;
    else if (pixel_y >= height)
      pixel_ys[i] = height-1;
    else
      pixel_ys[i] = pixel_y;
  }

  //calculate X direction
  double values[4];
  double values_dx[4];
  for(int i = 0; i < 4; i++)
  {
    int pixel_y = pixel_ys[i];
    double p0 = get_pixel(pixel_xs[0], pixel_y);
    double p1 = get_pixel(pixel_xs[1], pixel_y);
    double p2 = get_pixel(pixel_xs[2], pixel_y);
    double p3 = get_pixel(pixel_xs[3], pixel_y);
    catmul_rom_f(x_diff, p0, p1, p2, p3, values[i], values_dx[i]);
  }

  //calculate Y direction
  double value, tmp;
  catmul_rom_f(y_diff, values_dx[0], values_dx[1], values_dx[2], values_dx[3], dx, tmp);
  catmul_rom_f(y_diff, values[0], values[1], values[2], values[3], value, dy);
  return value;
}

double Image2dExtrapolated::get_pixel_extrapolated(int inx_x, int inx_y) const {
  //make safe indices
  int safe_inx_x = inx_x, safe_inx_y = inx_y;
  if (safe_inx_x < 0)
    safe_inx_x = 0;
  else if (safe_inx_x >= width)
    safe_inx_x = width - 1;
  if (safe_inx_y < 0)
    safe_inx_y = 0;
  else if (safe_inx_y >= height)
    safe_inx_y = height - 1;

  //handle extrapolation
  double value_shift = 0;
  if (inx_x < 0) {
    double value_a = get_pixel(0, safe_inx_y);
    double value_b = get_pixel(1, safe_inx_y);
    value_shift += -inx_x * (value_a - value_b);
  }
  else if (inx_x >= width) {
    double value_a = get_pixel(width-1, safe_inx_y);
    double value_b = get_pixel(width-2, safe_inx_y);
    value_shift += (inx_x - width + 1) * (value_a - value_b);
  }

  if (inx_y < 0) {
    double value_a = get_pixel(safe_inx_x, 0);
    double value_b = get_pixel(safe_inx_x, 1);
    value_shift += -inx_y * (value_a - value_b);
  }
  else if (inx_y >= height) {
    double value_a = get_pixel(safe_inx_x, height-1);
    double value_b = get_pixel(safe_inx_x, height-2);
    value_shift += (inx_y - height + 1) * (value_a - value_b);
  }

  return get_pixel(safe_inx_x, safe_inx_y) + value_shift;
}

double Image2dExtrapolated::get_sample(double x, double y, double& dx, double& dy) const 
{
  int pixel_x0 = (int)floor(x), pixel_y0 = (int)floor(y);
  double x_diff = (double)(x - pixel_x0), y_diff = (double)(y - pixel_y0);

  //create coordinates
  error_if(pixel_x0 < 0 || pixel_y0 < 0, "Accessing negative coordinates (%g, %g).", x, y);

  //get values
  double value_A_0 = get_pixel_extrapolated(pixel_x0, pixel_y0);
  double value_A_1 = get_pixel_extrapolated(pixel_x0+1, pixel_y0);
  double value_B_0 = get_pixel_extrapolated(pixel_x0, pixel_y0+1);
  double value_B_1 = get_pixel_extrapolated(pixel_x0+1, pixel_y0+1);
 
  //calculate
  double value_A = value_A_0 * (1-x_diff) + value_A_1 * x_diff;
  double value_B = value_B_0 * (1-x_diff) + value_B_1 * x_diff;
  double value = value_A * (1-y_diff) + value_B * y_diff;
  double value_0 = value_A_0 * (1-y_diff) + value_B_0 * y_diff;
  double value_1 = value_A_1 * (1-y_diff) + value_B_1 * y_diff;
  dx = value_1 - value_0;
  dy = value_B - value_A;

  return value;
}

double Image2dExtrPoint::get_sample(double x, double y, double& dx, double& dy) const 
{
  int pixel_x0 = (int)floor(x), pixel_y0 = (int)floor(y);

  //create coordinates
  error_if(pixel_x0 < 0 || pixel_y0 < 0, "Accessing negative coordinates (%g, %g).", x, y);

  //get values
  double value = get_pixel_extrapolated(pixel_x0, pixel_y0);
  double value_x_0 = get_pixel_extrapolated(pixel_x0-1, pixel_y0);
  double value_x_1 = get_pixel_extrapolated(pixel_x0+1, pixel_y0);
  double value_y_0 = get_pixel_extrapolated(pixel_x0, pixel_y0-1);
  double value_y_1 = get_pixel_extrapolated(pixel_x0, pixel_y0+1);
 
  //calculate d/dx and d/dy
  dx = (value_x_1 - value_x_0) / 2;
  dy = (value_y_1 - value_y_0) / 2;

  return value;
}
