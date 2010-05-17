//
// bitmap.h
//
// by David Andrs, 2007
//
// Interface for loading bitmap files (PPM, binary versions)
//

#ifndef _BITMAP_H_
#define _BITMAP_H_

// RGB triplet //////////////////////////////////////////////////////////////////////////

struct RGB {
  uint8_t r;
	uint8_t g;
	uint8_t b;
  void set_gray(uint8_t value) { r = value; g = value; b = value; };
  void set_gray_clamped(float value)
  {
    if (value < 0.0f)
    {
      r = 255; g = 0; b = 0;
    }
    else if (value > 1.0f)
    {
      r = 0; g = 0; b = 255;
    }
    else {
      uint8_t val = (uint8_t)floor(value * 255 + 0.5);
      r = g = b = val;
    }
  };
};


// Bitmap class /////////////////////////////////////////////////////////////////////////

class Bitmap {
public:
  Bitmap() : Data(NULL) {};
	virtual ~Bitmap();

	void LoadFromFile(const char *fileName);
  bool SaveToFile(const char *fileName); ///< Saves bitmap to a file.
  bool Create(int width, int height); ///< Creates a bitmap.

	int GetWidth() { return Width; }
	int GetHeight() { return Height; }
	RGB *GetData() { return Data; }

protected:

	int Width;			// width of the bitmap (in pixels)
	int Height;			// height of the bitmap (in pixels)
	RGB *Data;			// pixels (array - width x height elements)

protected: //functions
  bool readRow(FILE *file, unsigned char *row);
};

#endif
