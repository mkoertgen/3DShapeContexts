#ifndef __CGLBMP__
#define __CGLBMP__

#define WIN32_LEAN_AND_MEAN	

#include <windows.h>
#include <gl/gl.h>
#include <gl/glu.h>
#include <iostream>

#include "ijl.h"
#include "png\png.h"

#include <memory.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

const int BM = 19778;
const int TGA_RGB = 2;

typedef enum INTERNALFORMAT {BGR, RGB};

typedef struct
{
  BYTE red;
  BYTE green;
  BYTE blue;
} GLRGBTRIPLE;

typedef struct
{
  BYTE red;
  BYTE green;
  BYTE blue;
  BYTE alpha;
} GLRGBQUAD;

typedef struct
{
  BYTE tfType;
  BYTE tfColorMapType;
  BYTE tfImageType;
  BYTE tfColorMapSpec[5];
  BYTE tfOrigX[2];
  BYTE tfOrigY[2];
  BYTE tfWidth[2];
  BYTE tfHeight[2];
  BYTE tfBpp;
  BYTE tfImageDes;
} TGAHEADER;

class CGLBMP
{
  public:
    GLuint texID;

    CGLBMP();   
    CGLBMP(char szFileName[]);

    virtual ~CGLBMP();  

    bool LoadImage(char szFileName[]);
    bool SaveImage(char szFileName[]);
    
    bool SaveScreen(char szFileName[]);
    bool SaveScreen();

    bool AddAlpha(char szFileName[]);
    bool ColorKey(BYTE red, BYTE green, BYTE blue, BYTE tolerance);
    bool Stencil(char szFileName[], BYTE red, BYTE green, BYTE blue);

    bool Invert();
    bool FlipVert();
    bool FlipHorz();
    bool Rotate180();

    void GenTexture();
    void SetTextureWrap(GLint s, GLint y);
    void SetFilter(GLint min, GLint mag);
    void Bind();

    inline int GetWidth()
    { return m_width; }

    inline int GetHeight()
    { return m_height; }

    inline int GetColorDepth()
    { return m_colorDepth; }
    
    inline const void* GetData()
    { 
      if (!m_pData) 
        return 0;

      return (const void*)m_pData;
    }

    inline INTERNALFORMAT GetInternalFormat()
    {
      return m_iFormat;
    }

    inline void SetJPGQuality(int quality)
    {
      if (quality > 100) quality = 100;
      if (quality < 1) quality = 1;

      jpgQuality = quality;
    }

    inline void InitEmpty(int w, int h)
    {
      ClearMem();

      m_pData = new GLRGBQUAD[w*h];
      m_width = w;
      m_height = h;
      m_colorDepth = 32;
    }

  private:
    int m_width;
    int m_height;
    int m_colorDepth;
    INTERNALFORMAT m_iFormat;
    
    GLRGBQUAD *m_pData;
    RGBQUAD m_colors[255];
    
    DWORD m_redMask;
    DWORD m_greenMask;
    DWORD m_blueMask;

    GLint m_texWrapS;
    GLint m_texWrapT;
    GLint m_minFilter;
    GLint m_magFilter;

    int m_imageSize;

    bool pngSupported;
    bool jpgSupported;
    int jpgQuality;

    bool LoadBMP(const char *filename);
    bool LoadJPG(const char *filename);
    bool LoadPNG(const char *filename);
    bool LoadTGA(const char *filename);

    bool SaveBMP(const char *filename);
    bool SaveJPG(const char *filename);
    bool SavePNG(const char *filename);
    bool SaveTGA(const char *filename);

    bool SaveGLBuffer(char szFileName[]);

    void Load32(FILE *imgFile);
    void Load24(FILE *imgFile);
    void Load16(FILE *imgFile);
    void Load8(FILE *imgFile);
    void Load4(FILE *imgFile);

    bool SetInternalFormat(INTERNALFORMAT format);
    void SetAlpha(BYTE alpha);
    void ExpandTo32(RGBTRIPLE *image, int pad);

    bool InitJPG();
    PFNIJLINIT ijlInit;
    PFNIJLFREE ijlFree;
    PFNIJLREAD ijlRead;
    PFNIJLWRITE ijlWrite;
    HINSTANCE jpgLib;

    /*bool InitPNG();
    PFNPNGCREATEREADSTRUCT pngCreateReadStruct;
    PFNPNGCREATEINFOSTRUCT pngCreateInfoStruct;
    PFNPNGINITIO pngInitIO;
    PFNPNGSETREADSTATUSFN pngSetReadStatusFN;
    PFNPNGREADINFO pngReadInfo;
    PFNPNGCHECKSIG pngCheckSig;
    PFNPNGSETSIGBYTES pngSetSigBytes;
    PFNPNGGETIHDR pngGetIHDR;
    PFNPNGSETSWAP pngSetSwap;
    PFNPNGSETEXPAND pngSetExpand;
    PFNPNGSETPACKING pngSetPacking;
    PFNPNGGETVALID pngGetValid;
    PFNPNGSETGRAYTORGB pngSetGrayToRGB;
    PFNPNGGETGAMMA pngGetGamma;
    PFNPNGSETGAMMA pngSetGamma;
    PFNPNGREADUPDATEINFO pngReadUpdateInfo;
    PFNPNGGETROWBYTES pngGetRowBytes;
    PFNPNGREADIMAGE pngReadImage;
    PFNPNGREADEND pngReadEnd;
    PFNPNGDESTROYREADSTRUCT pngDestroyReadStruct;
    PFNPNGCREATEWRITESTRUCT pngCreateWriteStruct;
    PFNPNGSETWRITESTATUSFN pngSetWriteStatusFN;
    PFNPNGSETIHDR pngSetIHDR;
    PFNPNGWRITEINFO pngWriteInfo;
    PFNPNGWRITEIMAGE pngWriteImage;
    PFNPNGWRITEEND pngWriteEnd;
    PFNPNGWRITEFLUSH pngWriteFlush;
    PFNPNGDESTROYWRITESTRUCT pngDestroyWriteStruct;
    HINSTANCE pngLib;*/

    inline int GetLineWidth(int bpp)
    { 
      return ((m_width * bpp + 31) & -32) >> 3;
    } 
    
    inline void ClearMem()
    {
      if (m_pData) 
      {
        delete [] m_pData;
        m_pData = 0;
      }
    }

    inline int GetRightShiftCount(DWORD dwVal)
    {
      for (int i = 0; i < sizeof(DWORD)*8; i++)
      {
        if (dwVal & 1)
          return i;
        dwVal = dwVal >> 1;
      }
      return -1;
    }
    
    inline int GetNumSet(DWORD dwVal)
    {
      int nCount = 0;
      while (dwVal != 0) 
      {
        nCount += dwVal & 1;
        dwVal = dwVal >> 1;
      }

      return nCount;
    }

    inline int Scale8(int colorValue, int colorBits)
    {
      switch (colorBits)
      {
        case 1 : return ((colorValue) ? 255 : 0 );
        case 2 : return (colorValue << 6) | (colorValue << 4) | (colorValue << 2) | colorValue;
        case 3 : return (colorValue << 5) | (colorValue << 2) | (colorValue >> 1);
        case 4 : return (colorValue << 4) | colorValue;
        case 5 : return (colorValue << 3) | (colorValue >> 2);
        case 6 : return (colorValue << 2) | (colorValue >> 4);
        case 7 : return (colorValue << 1) | (colorValue >> 6);
        default : return colorValue;
      }
    }
};

#endif