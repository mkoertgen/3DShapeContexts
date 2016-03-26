#include "glbmp.h"

#define GL_CLAMP_TO_EDGE 0x812f
#define GL_BGRA 0x80E1

CGLBMP::CGLBMP()
{
  m_width = 0;
  m_height = 0;
  m_colorDepth = 0;
  m_pData = 0;
  texID = 0;

  m_texWrapS = GL_CLAMP;
  m_texWrapT = GL_CLAMP;
  m_minFilter = GL_LINEAR;
  m_magFilter = GL_LINEAR;

  m_iFormat = RGB;
  
  jpgQuality = 75;
  jpgSupported = InitJPG();
}

CGLBMP::CGLBMP(char szFileName[])
{
  CGLBMP();
  LoadImage(szFileName); 
}

bool CGLBMP::LoadImage(char szFileName[])
{
  char *ext;

  ClearMem();

  ext = &szFileName[strlen(szFileName)-3];

  if (!_stricmp(ext, "bmp"))
	return LoadBMP(szFileName);

  if (!_stricmp(ext, "jpg"))
	return LoadJPG(szFileName);

  if (!_stricmp(ext, "png"))
	return LoadPNG(szFileName);

  if (!_stricmp(ext, "tga"))
	return LoadTGA(szFileName);

  return false;
}

bool CGLBMP::SaveImage(char szFileName[])
{
  char *ext;

  ext = &szFileName[strlen(szFileName)-3];

  if (!_stricmp(ext, "bmp"))
	return SaveBMP(szFileName);

  if (!_stricmp(ext, "jpg"))
	return SaveJPG(szFileName);

  if (!_stricmp(ext, "png"))
	return SavePNG(szFileName);

  if (!_stricmp(ext, "tga"))
	return SaveTGA(szFileName);
  return false;
}

bool CGLBMP::LoadBMP(const char *filename)
{
  BITMAPFILEHEADER fHeader;
  BITMAPINFOHEADER iHeader;
  FILE *bmpFile;
  int palSize;

  if (!(bmpFile = fopen(filename, "rb")))
	return false;
  
  fread(&fHeader, 1, sizeof(fHeader), bmpFile);

  if (fHeader.bfType != BM)
  {
	fclose(bmpFile);
	return false;
  }

  fread(&iHeader, 1, sizeof(iHeader), bmpFile);

  if (iHeader.biSize == sizeof(BITMAPCOREHEADER))
  {
	fclose(bmpFile);
	return false;
  }

  m_width = iHeader.biWidth;
  m_height = iHeader.biHeight;
  m_colorDepth = iHeader.biBitCount;

  if (iHeader.biCompression != BI_RGB)
  {
	if (iHeader.biCompression == BI_BITFIELDS)
	{
	  if (m_colorDepth == 16)
	  {
		fread(&m_redMask, 1, sizeof(DWORD), bmpFile);
		fread(&m_greenMask, 1, sizeof(DWORD), bmpFile);
		fread(&m_blueMask, 1, sizeof(DWORD), bmpFile);
	  }
	}
	else 
	{
	  fclose(bmpFile);
	  return false;
	}
  }

  if (m_colorDepth < 16)
  {
	palSize = iHeader.biClrUsed*4;

	fread(&m_colors, 1, palSize, bmpFile);
  }

  fseek(bmpFile, fHeader.bfOffBits, SEEK_SET); 

  m_imageSize = iHeader.biSizeImage;

  if (m_imageSize == 0)
	m_imageSize = m_height * GetLineWidth(m_colorDepth);

  switch (m_colorDepth)
  {
	case 32 : 
	  Load32(bmpFile); 
	  break;
	case 24 : 
	  Load24(bmpFile);
	  break;
	case 16 :
	  Load16(bmpFile);
	  break;
	case 8:
	  Load8(bmpFile);
	case 4:
	  Load4(bmpFile);
	default :
	  break;
	
  }

  fclose(bmpFile);

  if (m_pData) 
	return true;
  else
	return false;
}

bool CGLBMP::LoadJPG(const char *filename)
{
  RGBTRIPLE *image;
  int jErr;
  int pad;
  JPEG_CORE_PROPERTIES jpgProps;

  if (!jpgSupported) 
	return false;

  jErr = ijlInit(&jpgProps);
  if (jErr != IJL_OK)
	return false;

  jpgProps.JPGFile = filename;

  jErr = ijlRead(&jpgProps, IJL_JFILE_READPARAMS);
  if (jErr != IJL_OK)
  {
	ijlFree(&jpgProps);
	return false;
  }

  m_width = jpgProps.JPGWidth;
  m_height = jpgProps.JPGHeight;
  m_colorDepth = 24;

  pad = GetLineWidth(24) - m_width*3;

  image = (RGBTRIPLE*) new unsigned char[GetLineWidth(24)*m_height];
  m_pData = new GLRGBQUAD[m_width*m_height];

  jpgProps.DIBWidth = m_width;
  jpgProps.DIBHeight = -m_height;
  jpgProps.DIBChannels = 3;
  jpgProps.DIBColor = IJL_RGB;
  jpgProps.DIBPadBytes = pad;
  jpgProps.DIBBytes = (unsigned char *)image;

  switch (jpgProps.JPGChannels) 
  {
	case 1: 
	  jpgProps.JPGColor = IJL_G;
	  break;
	case 3:
	  jpgProps.JPGColor = IJL_YCBCR;
	  break;
	default:
	  jpgProps.DIBColor = IJL_COLOR(IJL_OTHER);
	  jpgProps.JPGColor = IJL_COLOR(IJL_OTHER);
	  break;
  }

  jErr = ijlRead(&jpgProps, IJL_JFILE_READWHOLEIMAGE);
  if (jErr != IJL_OK) 
  {
	delete [] image;
	delete [] m_pData;
	m_pData = 0;
	return false;
  }

  ExpandTo32(image, pad);
  m_colorDepth = 32;

  delete [] image;

  ijlFree(&jpgProps);

  m_iFormat = RGB;

  return true;
}

bool CGLBMP::LoadPNG(const char *filename)
{
  png_structp png = 0;
  png_infop pngInfo = 0;
  FILE *pngFile;
  int rowBytes;
  int bpp;
  int colorType;
  png_byte **rowPtrs;
  png_byte sig[8];

  RGBTRIPLE *image;
  int i;
  double gamma;

  if (!filename) return false;

  if (!(pngFile = fopen(filename, "rb")))
	return false;

  fread(sig, 1, 8, pngFile);
  if (!png_check_sig(sig, 8))
  {
	fclose(pngFile);
	return false;
  }

  png = png_create_read_struct("1.0.10", 0, 0, 0);
  if (!png)
  {
	fclose(pngFile);
	return false;
  }

  pngInfo = png_create_info_struct(png);
  if (!pngInfo)
  {
	png_destroy_read_struct(&png, 0, 0);
	fclose(pngFile);
	return false;
  }

  png_init_io(png, pngFile);
  png_set_sig_bytes(png, 8);

  png_read_info(png, pngInfo);
  unsigned long w, h;
  png_get_IHDR(png, pngInfo, &w, &h, &m_colorDepth, &colorType, NULL, NULL, NULL);

  m_width = w;
  m_height = h;

  if (m_colorDepth == 16)
	png_set_swap(png);

  if (colorType == PNG_COLOR_TYPE_PALETTE)
	png_set_expand(png);

  if (m_colorDepth < 8)
	png_set_packing(png);

  if (png_get_valid(png, pngInfo, PNG_INFO_tRNS) == PNG_INFO_tRNS)
	png_set_expand(png);

  if ((colorType == PNG_COLOR_TYPE_GRAY) || (colorType == PNG_COLOR_TYPE_GRAY_ALPHA)) 
	png_set_gray_to_rgb(png);

  if ((png_get_gAMA(png, pngInfo, &gamma) == 1))
	png_set_gamma(png, 2.2, gamma);
  else
	png_set_gamma(png, 2.2, 0.45455);

  png_read_update_info(png, pngInfo);

  rowBytes = png_get_rowbytes(png, pngInfo);
  bpp = rowBytes / m_width;

  switch (bpp)
  {
	case 4:
	  m_pData = (GLRGBQUAD*) new unsigned char[rowBytes*m_height];
	  if (!m_pData)
	  {
		fclose(pngFile);
		png_destroy_read_struct(&png, &pngInfo, 0);
		return false;
	  }

	  rowPtrs = new png_byte * [m_height];
	  if (!rowPtrs)
	  {
		fclose(pngFile);
		delete [] m_pData;
		png_destroy_read_struct(&png, &pngInfo, 0);
		return false;
	  }

		for (i = 0; i < m_height; i++)
			rowPtrs[i] = (unsigned char *)m_pData + i * rowBytes;

	  png_read_image(png, rowPtrs);

	  delete [] rowPtrs;
	  rowPtrs = 0;

	  m_colorDepth = 32;
	  break;
	case 3:
	  image = (RGBTRIPLE*) new unsigned char[rowBytes*m_height];
	  if (!image)
	  {
		fclose(pngFile);
		png_destroy_read_struct(&png, &pngInfo, 0);
		return false;
	  }

	  rowPtrs = new png_byte * [m_height];
	  if (!rowPtrs)
	  {
		fclose(pngFile);
		png_destroy_read_struct(&png, &pngInfo, 0);
		delete [] image;
		return false;
	  }

		for (i = 0; i < m_height; i++)
			rowPtrs[i] = (unsigned char *)image + i * rowBytes;

	  png_read_image(png, rowPtrs);

	  delete [] rowPtrs;
	  rowPtrs = 0;

	  m_colorDepth = 32;

	  m_pData = new GLRGBQUAD[m_width*m_height];
	  if (!m_pData)
	  {
		fclose(pngFile);
		png_destroy_read_struct(&png, &pngInfo, 0);
		delete [] image;
		return false;
	  }

	  ExpandTo32(image, 0);
	  delete [] image;
	  SetAlpha(255);

	  break;
  }
  
  m_iFormat = RGB;
  
  png_read_end(png, pngInfo);
  png_destroy_read_struct(&png, &pngInfo, 0);
  
  fclose(pngFile);
  
  return true;
}

bool CGLBMP::LoadTGA(const char *filename)
{
  FILE *tgaFile;
  TGAHEADER tHeader;
  int bytesRead;
  RGBTRIPLE *image;

  if (!(tgaFile = fopen(filename, "rb")))
	return false;
  
  fread(&tHeader, 1, sizeof(tHeader), tgaFile);

  if (tHeader.tfImageType != TGA_RGB)
  {
	fclose(tgaFile);
	return false;
  }

  if (tHeader.tfColorMapType != 0)
  {
	fclose(tgaFile);
	return false;
  }

  m_width = tHeader.tfWidth[0] + tHeader.tfWidth[1] * 256;
  m_height = tHeader.tfHeight[0] + tHeader.tfHeight[1] * 256;
  m_colorDepth = tHeader.tfBpp;
  m_imageSize = m_width*m_height*(m_colorDepth/8);

  if ((m_colorDepth != 24) && (m_colorDepth != 32))
  {
	fclose(tgaFile);
	return false;
  }

  switch (m_colorDepth)
  {
	case 32:
	  m_pData = new GLRGBQUAD[m_width*m_height];

	  if (!m_pData)
		return false;

	  bytesRead = fread(m_pData, 1, m_imageSize, tgaFile);
	  if (bytesRead != m_imageSize)
	  {
		fclose(tgaFile);
		delete [] m_pData;
		m_pData = 0;
		return false;
	  }

	  break;
	case 24:
	  image = (RGBTRIPLE*) new unsigned char[m_imageSize]; 
	  if (!image)
		return false;

	  bytesRead = fread(image, 1, m_imageSize, tgaFile);
	  if (bytesRead != m_imageSize)
	  {
		fclose(tgaFile);
		delete [] image;
		image = 0;
		return false;
	  }

	  m_pData = new GLRGBQUAD[m_width*m_height];
	  if (!m_pData)
	  {
		delete [] image;
		image = 0;
		return false;
	  }

	  ExpandTo32(image, 0);
	  delete [] image;
	  image = 0;

	  break;
	default : 
	  fclose(tgaFile);
	  return false;
	  break;
  }

  m_iFormat = BGR;

  return true;
}

bool CGLBMP::SaveBMP(const char *filename)
{
  BITMAPFILEHEADER fHeader;
  BITMAPINFOHEADER iHeader;
  FILE *bmpFile;
  bool colorsSwapped;

  int bytesWritten;
  int length = m_width*m_height*4;

  if (!m_pData)
	return false;

  if (!(bmpFile = fopen(filename, "wb")))
	return false;

  fHeader.bfType = BM;
  fHeader.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(iHeader) + length;
  fHeader.bfReserved1 = 0; 
  fHeader.bfReserved2 = 0; 

  // This sets the distance from the start of the file to the start
  // of the bitmaps color data
  fHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(iHeader);
  
  bytesWritten = fwrite(&fHeader, 1, sizeof(fHeader), bmpFile);
  if (bytesWritten != sizeof(fHeader))
  {
	fclose(bmpFile);
	return false;
  }

  memset(&iHeader, 0, sizeof(iHeader));

  // Set up the information header
  iHeader.biSize = sizeof(BITMAPINFOHEADER);
  iHeader.biWidth = m_width;         // Current width
  iHeader.biHeight = m_height;       // Current height
  iHeader.biPlanes = 1;              // Number of planes, must be set to 1
  iHeader.biBitCount = 32;           // Current color depth
  iHeader.biCompression = BI_RGB;    // No compression
  iHeader.biSizeImage = length;      // Number of bytes in bitmap

  bytesWritten = fwrite(&iHeader, 1, sizeof(iHeader), bmpFile);
  if (bytesWritten != sizeof(iHeader))
  {
	fclose(bmpFile);
	return false;
  }

  colorsSwapped = SetInternalFormat(BGR);

  bytesWritten = fwrite(m_pData, 1, length, bmpFile);
  if (bytesWritten != length)
  {
	fclose(bmpFile);
	return false;
  }

  if (colorsSwapped)
	SetInternalFormat(RGB);

  fclose(bmpFile);
  return true;
}

bool CGLBMP::SaveJPG(const char *filename)
{
  int jErr;
  JPEG_CORE_PROPERTIES jpgProps;
  bool colorsSwapped;
  
  if (!jpgSupported)  
	return false;

  jErr = ijlInit(&jpgProps);
  if (jErr != IJL_OK)
	return false;

  jpgProps.DIBWidth = m_width;
  jpgProps.DIBHeight = -m_height;
  jpgProps.DIBBytes = (unsigned char *)m_pData;
  jpgProps.DIBPadBytes = 0 ;     
  jpgProps.DIBChannels = 4;
  jpgProps.DIBColor = IJL_RGB;
  
  jpgProps.JPGFile = filename;
  jpgProps.JPGWidth = m_width;
  jpgProps.JPGHeight = m_height;
  jpgProps.JPGChannels = 3;
  jpgProps.JPGColor = IJL_YCBCR;
  jpgProps.JPGSubsampling = IJL_411;
  jpgProps.jquality = jpgQuality;    

  colorsSwapped = SetInternalFormat(RGB);

  jErr = ijlWrite(&jpgProps, IJL_JFILE_WRITEWHOLEIMAGE);
  
  if (colorsSwapped)
	SetInternalFormat(BGR);
  
  if (jErr != IJL_OK)
  {
	ijlFree(&jpgProps);
	return false;
  }

  ijlFree(&jpgProps);

  return true;
}

bool CGLBMP::SavePNG(const char *filename)
{
  png_structp png = 0;
  png_infop pngInfo = 0;
  FILE *pngFile;
  png_byte **rowPtrs;
  int i;
  bool colorsSwapped;

  if (!m_pData)
	return false;

  if (!(pngFile = fopen(filename, "wb")))
	return false;

  png = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
  if (!png)
  {
	fclose(pngFile);
	return false;
  }

  pngInfo = png_create_info_struct(png);
  if (!pngInfo)
  {
	png_destroy_read_struct(&png, 0, 0);
	fclose(pngFile);
  }

  png_init_io(png, pngFile);

  png_set_IHDR(png, pngInfo, m_width, m_height, 8, PNG_COLOR_TYPE_RGB_ALPHA, 
	PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, 
	PNG_FILTER_TYPE_DEFAULT);

  png_write_info(png, pngInfo);

  rowPtrs = new png_byte * [m_height];
  if (!rowPtrs)
  {
	fclose(pngFile);
	png_destroy_read_struct(&png, &pngInfo, 0);
	return false;
  }

  colorsSwapped = SetInternalFormat(RGB);
  
  for (i = 0; i < m_height; i++)
	rowPtrs[(m_height-1) - i] = (unsigned char *)m_pData + (i * m_width * 4);

  png_write_image(png, rowPtrs);
  png_write_end(png, pngInfo);
  png_write_flush(png);

  delete [] rowPtrs;

  png_destroy_read_struct(&png, &pngInfo, 0);

  fclose(pngFile);

  if (colorsSwapped)
	SetInternalFormat(BGR);

  return true;
}

bool CGLBMP::SaveTGA(const char *filename)
{
  TGAHEADER tHeader;
  FILE *tgaFile;
  bool colorsSwapped;
  
  int bytesWritten;
  int length = m_width*m_height*4;

  
  if (!m_pData)
	return false;

  if (!(tgaFile = fopen(filename, "wb")))
	return false;

  memset(&tHeader, 0, sizeof(tHeader));

  tHeader.tfImageType = TGA_RGB;
  tHeader.tfWidth[0] = m_width % 256;
  tHeader.tfWidth[1] = m_width / 256;
  tHeader.tfHeight[0] = m_height % 256;
  tHeader.tfHeight[1] = m_height / 256;
  tHeader.tfBpp = 32;

  bytesWritten = fwrite(&tHeader, 1, sizeof(tHeader), tgaFile);
  if (bytesWritten != sizeof(tHeader))
  {
	fclose(tgaFile);
	return false;
  }

  colorsSwapped = SetInternalFormat(BGR);

  bytesWritten = fwrite(m_pData, 1, length, tgaFile);
  if (bytesWritten != length)
  {
	fclose(tgaFile);
	return false;
  }
  fclose(tgaFile);

  if (colorsSwapped)
	SetInternalFormat(RGB);

  return true;
}

bool CGLBMP::SaveScreen(char szFileName[])
{
  return SaveGLBuffer(szFileName);
}

bool CGLBMP::SaveScreen()
{
  return SaveGLBuffer("");
}

bool CGLBMP::SaveGLBuffer(char szFileName[])
{
  GLint viewport[4];

  glGetIntegerv(GL_VIEWPORT, viewport);

  m_width = viewport[2];
  m_height = viewport[3];
  m_colorDepth = 32;

  ClearMem();

  m_pData = new GLRGBQUAD[m_width*m_height];

  if (!m_pData)
	return false;

  glFinish();

  glPixelStorei(GL_PACK_ALIGNMENT, 4);
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);

  glReadPixels(0, 0, viewport[2], viewport[3], GL_RGBA, GL_UNSIGNED_BYTE, m_pData);
  m_iFormat = RGB;

  if (strcmp(szFileName, ""))
	return SaveImage(szFileName);

  return true;
}

bool CGLBMP::AddAlpha(char szFileName[])
{
  CGLBMP alphaImg;
  GLRGBQUAD *pixel;
  GLRGBQUAD *alpha;
  int i;

  if ((m_pData) && (alphaImg.LoadImage(szFileName)))
  {
	if ((alphaImg.m_width != m_width) || (alphaImg.m_height != m_height))
	{
	  alphaImg.~CGLBMP();
	  return false;
	}

	pixel = m_pData;
	alpha = alphaImg.m_pData;
	for (i = 0; i < (m_width*m_height); i++)
	{
	  pixel->alpha = ((alpha->red + alpha->green + alpha->blue) / 765) * 255;

	  pixel++;
	  alpha++;
	}
  }
  else 
	return false;
  
  return true;
}

bool CGLBMP::ColorKey(BYTE red, BYTE green, BYTE blue, BYTE tolerance)
{
  GLRGBQUAD *pixel;
  int i;

  if (m_pData)
  {
	red += tolerance;
	green += tolerance;
	blue += tolerance;

	pixel = m_pData;
	for (i = 0; i < (m_width*m_height); i++)
	{
	  if ((pixel->red <= red) && (pixel->green <= green) && (pixel->blue <= blue))
		pixel->alpha = 0;
	  else
		pixel->alpha = 255;

	  pixel++;
	}
  }
  else 
	return false;
  
  return true;
}

bool CGLBMP::Stencil(char szFileName[], BYTE red, BYTE green, BYTE blue)
{
  CGLBMP stencil;
  GLRGBQUAD *pixel;
  GLRGBQUAD *alpha;
  int i;

  if ((m_pData) && (stencil.LoadImage(szFileName)))
  {
	if ((stencil.m_width != m_width) || (stencil.m_height != m_height))
	{
	  stencil.~CGLBMP();
	  return false;
	}

	pixel = m_pData;
	alpha = stencil.m_pData;
	for (i = 0; i < (m_width*m_height); i++)
	{
	  if ((alpha->red == red) && (alpha->green == green) && (alpha->blue == blue))
		pixel->alpha = 0;
	  else
		pixel->alpha = 255;

	  pixel++;
	  alpha++;
	}
  }
  else 
	return false;
  
  return true;
}

void CGLBMP::Load32(FILE *imgFile)
{
  int bytesRead;

  m_pData = new GLRGBQUAD[m_width*m_height];

  if (!m_pData)
	return;

  bytesRead = fread(m_pData, 1, m_imageSize, imgFile);
  if (bytesRead != m_imageSize)
  {
	fclose(imgFile);
	delete [] m_pData;
	m_pData = 0;
	return;
  }

  m_iFormat = BGR;

  SetAlpha(255);
}

void CGLBMP::Load24(FILE *imgFile)
{
  RGBTRIPLE *bgrBits;
  int bytesRead;
  int pad;

  pad = GetLineWidth(24) - (m_width*3);

  bgrBits = (RGBTRIPLE*) new unsigned char[m_imageSize];

  if (!bgrBits)
	  return;

  bytesRead = fread(bgrBits, 1, m_imageSize, imgFile);
  if (bytesRead != m_imageSize)
  {
	fclose(imgFile);
	delete [] bgrBits;
	bgrBits = 0;
	return;
  }

  
  m_pData = new GLRGBQUAD[m_width*m_height];

  if (!m_pData)
	return;

  ExpandTo32(bgrBits, pad);

  delete [] bgrBits;
  bgrBits = 0;
  m_iFormat = BGR;
  m_colorDepth = 32;
}

void CGLBMP::Load16(FILE *imgFile)
{
  WORD *bgrBits;
  int bytesRead;
  int pad;
  
  int redBits, greenBits, blueBits;
  int redShr, greenShr;

  int x, y;

  WORD *source16;
  GLRGBQUAD *dest;

  bgrBits = (WORD*) new unsigned char[m_imageSize];

  if (!bgrBits)
	  return;

  bytesRead = fread(bgrBits, 1, m_imageSize, imgFile);
  if (bytesRead != m_imageSize)
  {
	fclose(imgFile);
	delete [] bgrBits;
	bgrBits = 0;
	return;
  }

  m_pData = new GLRGBQUAD[m_width*m_height];

  if (!m_pData)
	return;

  pad = GetLineWidth(16) - (m_width*2);
  
  redShr = GetRightShiftCount(m_redMask);
  greenShr = GetRightShiftCount(m_greenMask);

  redBits = GetNumSet(m_redMask);
  greenBits = GetNumSet(m_greenMask);
  blueBits = GetNumSet(m_blueMask);

  source16 = bgrBits;
  dest = m_pData;

  for (x = 0; x < m_height; x++)
  {
	for (y = 0; y < m_width; y++)
	{
	  dest->red = Scale8(*source16 & (int)m_blueMask, blueBits);
	  dest->green = Scale8((*source16 & (int)m_greenMask) >> greenShr, greenBits);
	  dest->blue = Scale8((*source16 & (int)m_redMask) >> redShr, redBits);
	  dest->alpha = 255;
	  source16++;
	  dest++;
	}
	source16 = (WORD*)((BYTE*)source16 + pad);
  }

  delete [] bgrBits;
  bgrBits = 0;
  m_iFormat = BGR;
  m_colorDepth = 32;
}

void CGLBMP::Load8(FILE *imgFile)
{
  BYTE *bgrBits;
  int bytesRead;
  BYTE *source8;
  GLRGBQUAD *dest;
  
  int x, y;
  int pad;

  bgrBits = new BYTE[m_imageSize];

  if (!bgrBits) 
	return;

  bytesRead = fread(bgrBits, 1, m_imageSize, imgFile);
  if (bytesRead != m_imageSize)
  {
	fclose(imgFile);
	delete [] bgrBits;
	bgrBits = 0;
	return;
  }

  m_pData = new GLRGBQUAD[m_width*m_height];

  if (!m_pData)
	return;

  pad = GetLineWidth(8) - m_width;

  source8 = bgrBits;
  dest = m_pData;

  for (x = 0; x < m_height; x++)
  {
	for (y = 0; y < m_width; y++)
	{
	  dest->red = m_colors[*source8].rgbBlue;
	  dest->green = m_colors[*source8].rgbGreen;
	  dest->blue = m_colors[*source8].rgbRed;
	  dest->alpha = 255;
	  dest++;
	  source8++;
	}
	source8 += pad;
  }

  delete [] bgrBits;
  bgrBits = 0;

  m_iFormat = BGR;
  m_colorDepth = 32;
}

void CGLBMP::Load4(FILE *imgFile)
{
  BYTE *bgrBits;
  int bytesRead;
  BYTE *source4;
  GLRGBQUAD *dest;
  int palEntry;
  
  int x, y;
  int pad;

  bgrBits = new BYTE[m_imageSize];

  if (!bgrBits) 
	return;

  bytesRead = fread(bgrBits, 1, m_imageSize, imgFile);
  if (bytesRead != m_imageSize)
  {
	fclose(imgFile);
	delete [] bgrBits;
	bgrBits = 0;
	return;
  }

  m_pData = new GLRGBQUAD[m_width*m_height];

  if (!m_pData)
	return;

  pad = GetLineWidth(4) - (m_width >> 1);

  source4 = bgrBits;
  dest = m_pData;

  for (x = 0; x < m_height; x++)
  {
	for (y = 0; y < m_width; y++)
	{
	  if (y % 2)
		palEntry = (int)(*source4 & 0xF);
	  else
		palEntry = (int)((*source4 >> 4) & 0xF);

	  dest->red = m_colors[palEntry].rgbBlue;
	  dest->green = m_colors[palEntry].rgbGreen;
	  dest->blue = m_colors[palEntry].rgbRed;
	  dest->alpha = 255;
	  dest++;

	  if (y % 2)
		source4++;
	}
	source4 += pad;
  }

  delete [] bgrBits;
  bgrBits = 0;
  m_iFormat = BGR;
  m_colorDepth = 32;
}

bool CGLBMP::Invert()
{
  int i;
  BYTE *pixel;

  if (!m_pData)
	return false;

  pixel = (BYTE*)m_pData;
  for (i = 0; i < (m_height*m_width*4); i++)
  {
	*pixel = 255 - *pixel;
	pixel++;
  }

  return true;
}


bool CGLBMP::FlipVert()
{
  int i;
  GLRGBQUAD *top;
  GLRGBQUAD *bottom;
  GLRGBQUAD *tmpBits;
  int lineWidth = m_width * 4;

  if (!m_pData)
	return false;

  tmpBits = new GLRGBQUAD[m_width];

  if (!tmpBits)
	return false;

  top = m_pData;
  bottom = (GLRGBQUAD*)((BYTE*)m_pData + lineWidth*(m_height-1));

  for (i = 0; i < (m_height >> 1); i++)
  {
	memcpy(tmpBits, top, lineWidth); 
	memcpy(top, bottom, lineWidth);
	memcpy(bottom, tmpBits, lineWidth);

	top = (GLRGBQUAD*)((BYTE*)top + lineWidth);
	bottom = (GLRGBQUAD*)((BYTE*)bottom - lineWidth);
  }
  
  delete [] tmpBits;
  tmpBits = 0;

  return true;
}

bool CGLBMP::FlipHorz()
{
  int x, y;
  GLRGBQUAD *front;
  GLRGBQUAD *back;
  GLRGBQUAD tmp;

  int line = m_width * 4;

  if (!m_pData)
	return false;

  for (y = 0; y < m_height; y++)
  {
	front = (GLRGBQUAD*)((int)m_pData + line*y);
	back = (GLRGBQUAD*)((int)m_pData + (line*y) + (line-4));

	for (x = 0; x < (m_width >> 1)-1; x++)
	{
	  tmp = *front;
	  *front = *back;
	  *back = tmp;

	  front++;
	  back--;
	}
  }

  return true;
}

bool CGLBMP::Rotate180()
{
  if (!FlipHorz()) 
	return false;

  if (!FlipVert())
	return false;

  return true;
}

void CGLBMP::GenTexture()
{
  if (texID)
	glDeleteTextures(1, &texID);

  glGenTextures(1, &texID);
  glBindTexture(GL_TEXTURE_2D, texID);
  
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, m_texWrapS);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, m_texWrapT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, m_magFilter);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, m_minFilter);
  
  if (m_iFormat == RGB)
	gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, m_width, m_height, GL_RGBA, GL_UNSIGNED_BYTE, m_pData);
  else
	gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, m_width, m_height, GL_BGRA, GL_UNSIGNED_BYTE, m_pData);
}

void CGLBMP::SetFilter(GLint min, GLint mag)
{
  switch (min)
  {
	case GL_NEAREST : 
	  m_minFilter = GL_NEAREST;
	  break;
	case GL_LINEAR : 
	  m_minFilter = GL_LINEAR;
	  break;
	case GL_NEAREST_MIPMAP_NEAREST : 
	  m_minFilter = GL_NEAREST_MIPMAP_NEAREST;
	  break;
	case GL_LINEAR_MIPMAP_NEAREST : 
	  m_minFilter = GL_LINEAR_MIPMAP_NEAREST;
	  break;
	case GL_NEAREST_MIPMAP_LINEAR : 
	  m_minFilter = GL_NEAREST_MIPMAP_LINEAR;
	  break;
	case GL_LINEAR_MIPMAP_LINEAR : 
	  m_minFilter = GL_LINEAR_MIPMAP_LINEAR;
	  break;
	default : 
	  m_minFilter = GL_LINEAR;
	  break;
  }

  switch (mag) 
  {
	case GL_NEAREST : 
	  m_magFilter = GL_NEAREST;
	  break;
	case GL_LINEAR : 
	  m_magFilter = GL_LINEAR;
	  break;
	default : 
	  m_magFilter = GL_LINEAR;
	  break;
  }
}

void CGLBMP::SetTextureWrap(GLint s, GLint t)
{
  switch (s)
  {
	case GL_CLAMP : 
	  m_texWrapS = GL_CLAMP;
	  break;
	case GL_REPEAT : 
	  m_texWrapS = GL_REPEAT;
	  break;
	case GL_CLAMP_TO_EDGE :
	  m_texWrapS = GL_CLAMP_TO_EDGE;
	  break;
	default:
	  m_texWrapS = GL_CLAMP;
	  break;
  }

  switch (t)
  {
	case GL_CLAMP : 
	  m_texWrapT = GL_CLAMP;
	  break;
	case GL_REPEAT : 
	  m_texWrapT = GL_REPEAT;
	  break;
	case GL_CLAMP_TO_EDGE : 
	  m_texWrapT = GL_CLAMP_TO_EDGE;
	  break;
	default:
	  m_texWrapT = GL_CLAMP;
	  break;
  }
}

void CGLBMP::Bind()
{
  if (texID)
	glBindTexture(GL_TEXTURE_2D, texID);
}

bool CGLBMP::SetInternalFormat(INTERNALFORMAT format)
{
  int i;
  GLRGBQUAD *pixel;
  unsigned char tmp;

  if (m_iFormat == format)
	return false;

  pixel = m_pData;

  for (i = 0; i < (m_width*m_height); i++)
  {
	tmp = pixel->blue;
	pixel->blue = pixel->red;
	pixel->red = tmp;

	pixel++;
  }
  m_iFormat = format;

  return true;
}

void CGLBMP::SetAlpha(BYTE alpha)
{
  int i;
  GLRGBQUAD *pixel;

  pixel = m_pData;
  for (i = 0; i < (m_width*m_height); i++)
  {
	pixel->alpha = alpha;
	pixel++;
  }
}

void CGLBMP::ExpandTo32(RGBTRIPLE *image, int pad)
{
  int x, y;
  RGBTRIPLE *pixel24;
  GLRGBQUAD *pixel32;

  pixel24 = image;
  pixel32 = m_pData;

  for (x = 0; x < m_height; x++)
  {
	for (y = 0; y < m_width; y++)
	{
	  pixel32->red = pixel24->rgbtBlue;
	  pixel32->green = pixel24->rgbtGreen;
	  pixel32->blue = pixel24->rgbtRed;
	  pixel32->alpha = 255;

	  pixel24++;
	  pixel32++;
	}
	pixel24 = (RGBTRIPLE*)((BYTE*)pixel24 + pad);
  }
}

bool CGLBMP::InitJPG()
{
  jpgLib = LoadLibrary("ijl15.dll");
  if (!jpgLib)
	return false;

  ijlInit = (PFNIJLINIT) GetProcAddress(jpgLib, "ijlInit");
  ijlFree = (PFNIJLFREE) GetProcAddress(jpgLib, "ijlFree");
  ijlRead = (PFNIJLREAD) GetProcAddress(jpgLib, "ijlRead");
  ijlWrite = (PFNIJLWRITE) GetProcAddress(jpgLib, "ijlWrite");

  if (ijlInit &&
	  ijlFree &&
	  ijlRead &&
	  ijlWrite)
	return true;
  else
	return false;
}

CGLBMP::~CGLBMP()
{
  ClearMem();

  FreeLibrary(jpgLib);

  if (texID)
  {
	glDeleteTextures(1, &texID);
	texID = 0;
  }
}
