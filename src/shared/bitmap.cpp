//
// bitmap.cpp
//
// by David Andrs, 2007
//
// - load a bitmap files (PPM)
//

#include <inttypes.h>
#include <stdio.h>
#include <new>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include "bitmap.h"

using namespace std;

bool Bitmap::readRow(FILE *file, unsigned char *row) {
	int idx = 0;
	int rd;
	unsigned char ch;
	do {
		if ((rd = fread(&ch, 1, 1, file)) == 1) {
			row[idx] = ch;
			idx++;
		}
	} while (rd == 1 && ch != 0x0A);
	row[idx] = '\0';

  return false;
}

// Bitmap class /////////////////////////////////////////////////////////////////////////

Bitmap::~Bitmap() {
	delete [] Data;
}

bool Bitmap::Create(int width, int height)
{
  try 
  {
    RGB* new_data = new RGB[width * height];
    delete[] Data;
    Data = new_data;
    this->Width = width;
    this->Height = height;
    return true;
  }
  catch (std::bad_alloc)
  {
    return false;
  }
}

void Bitmap::LoadFromFile(const char *fileName)
{
	FILE *img = NULL;
  RGB* data = NULL;
  try {
	  int wd, ht;				// image dimensions
	  if ((img = fopen(fileName, "rb")) == NULL)
      throw std::runtime_error("Cannot open file");

	  size_t rd;
	  unsigned char hdr[3];
	  rd = fread(&hdr, 3, 1, img);
	  if (hdr[0] == 0x50 &&
		  hdr[1] == 0x36 &&
		  hdr[2] == 0x0A)
	  {
		  // do nothing
	  }
	  else
      throw std::runtime_error("Invalid PPM header");

    //read dimensions
  #define BUFFER_SIZE 1024
	  char buffer[BUFFER_SIZE];
    wd = -1; ht = -1;
    do
    {
      char* row = fgets(buffer, BUFFER_SIZE, img);
      if (row[0] != '#')
      {
        if (sscanf(row, "%d %d", &wd, &ht) != 2 || wd <= 0 || ht <= 0)
          throw std::runtime_error("Dimensions are corrupted.");
      }
    } while(buffer[0] == '#' && wd <= 0 && ht <= 0);

    //read maximum range
    int maxRange;
    fgets(buffer, BUFFER_SIZE, img);
    if (sscanf(buffer, "%d", &maxRange) != 1) 
      throw std::runtime_error("Invalid maximum value.");
    if (maxRange > 255)
      throw std::runtime_error("16 bps images not supported.");

	  // read the image data
    int pxs = wd*ht;
    data = new RGB[pxs];
    fread(data, sizeof(RGB), pxs, img);

    // normalize data to be in a range [0, 255]
    int vls = pxs * 3;
    unsigned char* values = (unsigned char*)data;
    for(int i = 0; i < vls; i++)
    {
      int val = values[i];
      values[i] = (unsigned char)((val * 255) / maxRange);
    }

    //close file
    fclose(img); img = NULL;

    // clean the class
    delete[] this->Data;

    // fill the class
	  Width = wd;
	  Height = ht;
	  Data = data;
  }
  catch (std::runtime_error& err) {
    if (img != NULL)
      fclose(img);
    if (data != NULL)
      delete[] data;
    throw err;
  }
}

bool Bitmap::SaveToFile(const char *fileName)
{
  FILE *file = NULL;
  std::string err_message;

  //check the contents
  if (Data == NULL) { err_message = "E bitmap is not initialized\n"; goto error; };

  //open the file
  file = fopen(fileName, "wb");
  if (file == NULL) { err_message = "E unable to open the output bitmap file\n"; goto error; };

  //write header
  fprintf(file, "P6\n"); 
  fprintf(file, "%d %d\n", Width, Height);
  fprintf(file, "255\n");

  //write contents
  fwrite(Data, sizeof(RGB), Width * Height, file);

  //close and finish
  fclose(file);

  return true;

error:
  if (file != NULL)
    fclose(file);
  if (!err_message.empty())
    fprintf(stderr, "%s", err_message.c_str());
  return false;
}


