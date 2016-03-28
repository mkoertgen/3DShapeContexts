#ifndef MATHDEF_H
#define MATHDEF_H

//-------------- constants -------------------------------

/// 2*PI
#define M_2PI			6.28318530718f

/// M_PI is not defined in Microsoft headers
#ifndef M_PI
  #define M_PI			3.14159265358979323846
#endif

/// PI/2
#define PI2				1.57079632679f

/// PI/4
#define PI4				0.7853981634f

/// PI/180
#define PIOVER180		1.74532925199433E-002f

/// 180/PI
#define PIUNDER180		5.72957795130823E+001f

/// 1/PI
#define ONEOVERPI       0.318309886182f

/// 1/2PI
#define ONEOVER2PI      0.159154943091f

/// lower bound for underflow detection with floats
#define EPSF 1.0e-4f
/// lower bound for underflow detection with doubles
#define EPSD 1.0e-8f

/// index for X-value
#define X 0
/// index for Y-value
#define Y 1
/// index for Z-value
#define Z 2

//-------------- macros ----------------------------------
/// \macro return squared value
#define SQR(x) ((x)*(x))

/// \macro return sign of given value: -1,0,1
#define SIGN(x) (((x)>0)?1:-1)

/// \macro zero-test for floats
#define ISZEROF(x) ( fabs(x) < EPSF ? true : false )
/// \macro zero-test for doubles
#define ISZEROD(x) ( fabs(x) < EPSD ? true : false )
/// \macro non-zero-test for floats
#define NONZEROF(x) ( fabs(x) > EPSF ? true : false )
/// \macro non-zero-test for doubles
#define NONZEROD(x) ( fabs(x) > EPSD ? true : false )

#define RAD2DEG(x) ((x)*PIUNDER180)

#endif //MATHDEF_H
