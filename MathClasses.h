/** 
* \author Marcel Koertgen (marcel@koertgen.de)
* \date 02-02-2003
* \version 1.0
*/
#ifndef MATHCLASSSES_H
#define MATHCLASSSES_H

#include <math.h>
#include <ostream>
#include "MathDef.h"

//-------------- classes ---------------------------------

/**
* \class quaternion
* \brief template quaternion class
*
* This class is a template quaternion class that works for any floating point class.
* It is mainly used for doing 3d vector rotations.
*/
template<class Type> class quaternion {
public:
	 /// (default) Constructor - Sets w=1, x=y=z=0.
	quaternion() { w=1; x=y=z=0; }
	 /// Constructor - Sets quaternion to given values
	quaternion(const Type sw, const Type sx, const Type sy, const Type sz)
	{  w=sw; x=sx; y=sy; z=sz; }
	 /// Destruktor
	~quaternion() {}
	/// Copy-Constructor
	quaternion(const quaternion &q) 
	{  w=q.w; x=q.x; y=q.y; z=q.z; }
	/// assign operator
	quaternion& operator = (const quaternion &q) 
	{  w=q.w; x=q.x; y=q.y; z=q.z; return *this; }
	/// add operator
	quaternion operator +  (const quaternion &q) 
	{
		quaternion q2(w+q.w, x+q.x, y+q.y, z+q.z);
		return q2;
	}
	/// scalar mult operator
	quaternion operator *  (const Type &fac)   
	{
		quaternion q2(fac*w, fac*x, fac*y, fac*z);
		return q2;
	}
	/// mult operator
	quaternion operator *  (const quaternion &q) 
	{
		Type s[9], t;
		s[1] = (z-y) * (q.y-q.z); 
		s[2] = (w+x) * (q.w+q.x);
		s[3] = (w-x) * (q.y+q.z); 
		s[4] = (z+y) * (q.w-q.x);
		s[5] = (z-x) * (q.x-q.y); 
		s[6] = (z+x) * (q.x+q.y);
		s[7] = (w+y) * (q.w-q.z);
		s[8] = (w-y) * (q.w+q.z) ;
		s[9] = s[6] + s[7] + s[8] ;
		t = (s[5] + s[9])/2;
		quaternion q2(s[1]+t-s[6], s[2]+t-s[9], s[3]+t-s[8], s[4]+t-s[7] );
		return q2;
	}
	/// accum. add operator
	quaternion& operator += (const quaternion &q) 
	{  return *this = *this + q; }
	/// accum. scalar mult operator
	quaternion& operator *= (const Type &fac) 
	{  return *this = *this * fac; }
	/// accum. mult operator
	quaternion& operator *= (const quaternion &q) 
	{  return  *this = *this * q; }
	/// equal operator
	bool operator == (const quaternion &q) 
	{  return (w==q.w && x==q.x && y==q.y && z==q.z); }
	/// not-equal operator
	bool operator != (const quaternion &q) 
	{  return (w!=q.w || x!=q.x || y!=q.y || z!=q.z); }

	/// printf
	inline void print() const
	{
		printf("[ %.3g, (%.3g, %.3g, %.3g) ]",w,x, y, z);
	}
	/// print to ostream
	inline void Print(std::ostream &output) const 
	{  output << "[ " << w << ", (" << x << ", " << y << ", " << z << ") ]"; }
	/// ostream "<<" operator
	friend inline std::ostream & operator << (std::ostream& output, const quaternion & q)
	{
	  q.Print(output);
	  return output; 
	}

	/// return quaternion magnitude: w*w + <v,v>
	Type magnitude() {  return (w*w + x*x + y*y + z*z); }
	/// return quaternion norm
	Type norm() {  return sqrt(w*w + x*x + y*y + z*z); }
	/// return scaled quaternion
	quaternion scale(const Type fac)
	{
		quaternion q2(fac*w, fac*x, fac*y, fac*z);
		return q2;
	}
	/// return normalized or "union" quaternion
	quaternion normalize()
	{
		quaternion q2(w,x,y,z);
		q2 = q2.scale(1.0/sqrt(w*w + x*x + y*y + z*z));
		return q2;
	}
	/// return conjugated quaternion
	quaternion conjugate()
	{
		quaternion q2(w, -x, -y, -z);
		return q2;
	}
	/// return inverse quaternion
	quaternion inverse()
	{
		quaternion q2(w,x,y,z);
		Type m = 1/q2.magnitude();
		q2 = q2.conjugate();
		q2 = q2.scale(m);
		return q2;
	}
	/// set quaternion from an angle and an axis -> used for vector rotation
	void from_axis(const Type rad, const Type ax, const Type ay, const Type az)
	{
		Type a = rad/2;
		Type sina = sin(a);
		Type v_len = sqrt(ax*ax + ay*ay + az*az); // Rotationsachse normieren
		w = cos(a);
		x = ax * sina / v_len;
		y = ay * sina / v_len;
		z = az * sina / v_len;
	}
	/// set quaternion from euler angles: roll, pitch, yaw
	void EulerToQuat(const Type roll, const Type pitch, const Type yaw)
	{
		Type cr, cp, cy, sr, sp, sy, cpcy, spsy;
		cr = cos(roll/2);
		cp = cos(pitch/2);
		cy = cos(yaw/2);
		sr = sin(roll/2);
		sp = sin(pitch/2);
		sy = sin(yaw/2);
		cpcy = cp * cy;
		spsy = sp * sy;
		w = cr * cpcy + sr * spsy;
		x = sr * cpcy - cr * spsy;
		y = cr * sp * cy + sr * cp * sy;
		z = cr * cp * sy - sr * sp * cy;
	}

	/// the four quaternion components: w - Scalar, (x,y,z) - Vector
	Type w,x,y,z;
	 /// Zero quaternion (w:0,x:0,y:0,z:0)-> for comparison
	static quaternion ZERO;
};

/// \typedef quatf - a float quaternion class
typedef quaternion<float> quatf;
/// \typedef quatv - a double quaternion class
typedef quaternion<double> quatd;

template<class Type> quaternion<Type> quaternion<Type>::ZERO = quaternion<Type>(0,0,0,0);


//---------------------------------------------------------------------------------------
/**
* \class Vec3
* \brief 3d vector template class.
*
* 3d vector class template. Works for any integer or real type.
*/
template<class Type> class Vec3 {
public:
	/// x- coordinate
	Type x;
	/// y- coordinate
	Type y;
	/// z- coordinate
	Type z;

	/// default constructor
	Vec3 (void) {};
	/// constructor - set vector to given (x,y,z) values
	Vec3 (const Type vx, const Type vy, const Type vz) { x=vx; y=vy; z=vz; };
	/// constructor - set vector to given 3d vector
	Vec3 (const Vec3& v) { x=v.x; y=v.y; z=v.z; };
	/// constructor - set vector to given array (0:x, 1:y, 2:z)
	Vec3 (const Type v[3]) { x=v[0]; y=v[1]; z=v[2]; };

	/// Set vector coordinates
	void Set (const Type vx, const Type vy, const Type vz) { x=vx; y=vy; z=vz; }
	/// Set vector coordinates from array
	void Set (const Type v[3]) { x=v[0]; y=v[1]; z=v[2]; };

	/// Type * CONVERSION
	operator Type*()
	{ return (Type *)&x; }
	/// CONST Type * CONVERSION
	operator const Type*() const
	{ return &x; }

	/// COMPARISON (==)
	int operator == (const Vec3& A) const        
	{ return (x==A.x && y==A.y && z==A.z); }
	/// COMPARISON (!=)
	int operator != (const Vec3& A) const        
	{ return (x!=A.x || y!=A.y || z!=A.z); }

	/// ASSIGNMENT (=)
	inline Vec3& operator = (const Vec3& A)            
	{ x=A.x; y=A.y; z=A.z; return(*this);  }
	/// ADDITION (+)
	Vec3 operator + (const Vec3& A) const       
	{ Vec3 Sum(x+A.x, y+A.y, z+A.z); return(Sum); }
	/// SUBTRACTION (-)
	Vec3 operator - (const Vec3& A) const       
	{ Vec3 Diff(x-A.x, y-A.y, z-A.z); return(Diff); }
	/// DOT-PRODUCT (*)
	Type operator * (const Vec3& A) const       
	{ Type DotProd = x*A.x+y*A.y+z*A.z; return(DotProd); }
	/// CROSS-PRODUCT (/)
	Vec3 operator / (const Vec3& A) const       
	{ Vec3 CrossProd(y*A.z-z*A.y, z*A.x-x*A.z, x*A.y-y*A.x); return(CrossProd); }
	/// ALSO CROSS-PRODUCT (^)
	Vec3 operator ^ (const Vec3& A) const       
	{ Vec3 CrossProd(y*A.z-z*A.y, z*A.x-x*A.z, x*A.y-y*A.x); return(CrossProd); }
	/// MULTIPLY BY SCALAR V*s (*)
	Vec3 operator * (const Type s) const        
	{ Vec3 Scaled(x*s, y*s, z*s); return(Scaled); }
	/// DIVIDE BY SCALAR (/)
	Vec3 operator / (const Type s) const        
	{ Vec3 Scaled(x/s, y/s, z/s); return(Scaled); }
	/// COMPONENT MULTIPLY (&)
	Vec3 operator & (const Vec3& A) const       
	{ Vec3 CompMult(x*A.x, y*A.y, z*A.z); return(CompMult); }
	/// SCALAR MULT s*V
	friend inline Vec3 operator *(Type s, const Vec3& v)  
	{ return Vec3(v.x*s, v.y*s, v.z*s); }
	/// ACCUMULATED VECTOR ADDITION (+=)
	Vec3& operator += (const Vec3& A)    
	{ x+=A.x; y+=A.y; z+=A.z; return *this;}
	/// ACCUMULATED VECTOR SUBTRACTION (-=)
	Vec3& operator -= (const Vec3& A)    
	{ x-=A.x; y-=A.y; z-=A.z; return *this; }
	/// ACCUMULATED SCALAR MULT (*=)
	Vec3& operator *= (const Type s)     
	{ x*=s; y*=s; z*=s; return *this; }
	/// ACCUMULATED SCALAR DIV (/=)
	Vec3& operator /= (const Type s)     
	{ x/=s; y/=s; z/=s; return *this; }
	/// ACCUMULATED COMPONENT MULTIPLY (&=)
	Vec3& operator &= (const Vec3& A)  
	{ x*=A.x; y*=A.y; z*=A.z; return *this; }
	/// NEGATION (-)
	Vec3 operator - (void) const        
	{ Vec3 Negated(-x, -y, -z); return(Negated); }

	/// ALLOWS VECTOR ACCESS AS AN ARRAY.
	const Type& operator [] (const int i) const 
	{ return( (i==0)?x:((i==1)?y:z) ); }
//    /// ALLOWS VECTOR ACCESS AS AN ARRAY.
//    Type & operator [] (const int i) { return( (i==0)?x:((i==1)?y:z) ); };

	/// return vector length
	Type Length (void) const
	{ return ((Type)sqrt(x*x+y*y+z*z)); }
	/// return vector magnitude
	Type LengthSqr (void) const
	{ return (x*x+y*y+z*z); }

	/// normalize vector
	inline void Normalize (void)      
	{
		Type L = Length();        
		if (L>0) { x/=L; y/=L; z/=L; } 
	}                           
	/// Updates the given Max/Min vectors if necessary
	inline void UpdateMinMax(Vec3 &Min, Vec3 &Max) {
	  if (x<Min.x) Min.x=x; else if (x>Max.x) Max.x=x;
	  if (y<Min.y) Min.y=y; else if (y>Max.y) Max.y=y;
	  if (z<Min.z) Min.z=z; else if (z>Max.z) Max.z=z;
	}

	/// returns rotated Vector about given radian angle around (normalized) axis
	inline Vec3 Rotated(const Type rad, const Vec3& a) // in radians
	{
		quaternion<Type> p(0, x,y,z), q, u;
		q.from_axis(rad, a.x, a.y, a.z);
		u = q*p*q.inverse();
		return Vec3(u.x, u.y, u.z);
	}

	/// rotates Vector about given radian angle around (normalized) axis
	inline void Rotate(const Type rad, const Vec3& a) // in radians
	{
		quaternion<Type> p(0, x,y,z), q, u;
		q.from_axis(rad, a.x, a.y, a.z);
		u = q*p*q.inverse();
		x=u.x; y=u.y; z=u.z;
	}

	/// swaps two vectors
	inline void swap(Vec3& a)
	{
		Type tx=x, ty=y, tz=z;
		x=a.x; y=a.y; z=a.z;
		a.x=tx; a.y=ty; a.z=tz;
	}

	/// Prints vector coordinates (uses printf)
	inline void print() const
	{
		printf("(%.3g, %.3g, %.3g)",x, y, z);
	}
	/// print to ostream
	inline void Print(std::ostream &output) const 
	{  output << "(" << x << ", " << y << ", " << z << ")"; }
	/// ostream "<<" operator
	friend inline std::ostream & operator << (std::ostream& output, const Vec3 &v)
	{
	  v.Print(output);
	  return output; 
	}

	/// static ZERO-Vector
	static Vec3 VEC0;
	/// static X-Vector (1,0,0)
	static Vec3 VECX;
	/// static Y-Vector (0,1,0)
	static Vec3 VECY;
	/// static Z-Vector (0,0,1)
	static Vec3 VECZ;
};

template<class Type> Vec3<Type> Vec3<Type>::VEC0 = Vec3<Type>(0,0,0);
template<class Type> Vec3<Type> Vec3<Type>::VECX = Vec3<Type>(1,0,0);
template<class Type> Vec3<Type> Vec3<Type>::VECY = Vec3<Type>(0,1,0);
template<class Type> Vec3<Type> Vec3<Type>::VECZ = Vec3<Type>(0,0,1);

/// \typedef Vec3f - A float 3d Vector
typedef Vec3<float> Vec3f;
/// \typedef Vec3d - A double 3d Vector
typedef Vec3<double> Vec3d;

//---------------------------------------------------------------------------------------
/**
* \class Matrix
* \brief 3x3 Matrix template class.
*
* 3x3 Matrix class template. Works for any integer or real type.
* The matrix can represent any Rotation transformation and is used
* for rotation and scaling, etc...
*/
template<class Type> class Matrix
{
	public:
		Type m[9]; /// matrix elements in row-order (one-dimensional arrays are faster)

	/// Default constructor, sets the identity matrix
	Matrix() 
	{
		m[0]=m[4]=m[8] = (Type)1;
		m[1]=m[2]=m[3] = m[5]=m[6]=m[7] = (Type)0;
	}

	/// Copy-constructor
	Matrix(const Matrix& in) 
	{ 
		m[0]=in.m[0]; m[1]=in.m[1]; m[2]=in.m[2];
		m[3]=in.m[3]; m[4]=in.m[4]; m[5]=in.m[5];
		m[6]=in.m[6]; m[7]=in.m[7]; m[8]=in.m[8];
	}

	/// Atribuition operator
	void operator=(const Matrix& in) 
	{ 
		m[0]=in.m[0]; m[1]=in.m[1]; m[2]=in.m[2];
		m[3]=in.m[3]; m[4]=in.m[4]; m[5]=in.m[5];
		m[6]=in.m[6]; m[7]=in.m[7]; m[8]=in.m[8];
	}

	/// Allows Matrix access as a 2d-array
	const Type* operator [] (const int i) const 
	{ return &m[i*3]; }

	/// Nullify all elements
	inline void null(void)
	{
		m[0]=m[1]=m[2]=
		m[3]=m[4]=m[5]=
		m[6]=m[7]=m[8]= (Type)0;
	}

	/// Load the identity matrix
	inline void load_identity(void)
	{
		m[0]=m[4]=m[8] = (Type)1;
		m[1]=m[2]=m[3] = m[5]=m[6]=m[7] = (Type)0;
	}

	/** Set the matrix as the rotation matrix of local system <lx,ly,lz>
	* ex -> lx
	* ey -> ly
	* ez -> lz
	*/
	inline void set_rotation( const Vec3<Type> &lx, const Vec3<Type> &ly, const Vec3<Type> &lz )
	{
		m[0]=lx.x; m[1]=ly.x; m[2]=lz.x;
		m[3]=lx.y; m[4]=ly.y; m[5]=lz.y;
		m[6]=lx.z; m[7]=ly.z; m[8]=lz.z;
	}

	/** Set the matrix as the inverse rotation matrix of local system <lx,ly,lz>
	* lx -> ex
	* ly -> ey
	* lz -> ez
	*/
	inline void set_inverse_rotation( const Vec3<Type> &lx, const Vec3<Type> &ly, const Vec3<Type> &lz )
	{
		m[0]=lx.x; m[1]=lx.y; m[2]=lx.z;
		m[3]=ly.x; m[4]=ly.y; m[5]=ly.z;
		m[6]=lz.x; m[7]=lz.y; m[8]=lz.z;
	}

	/** return axis vectors from matrix
	* lx = (m[0],m[3],m[6])
	* ly = (m[1],m[4],m[7])
	* lz = (m[2],m[5],m[8])
	*/
	inline void get_rotation_axes(Vec3<Type> &lx, Vec3<Type> &ly, Vec3<Type> &lz)
	{
	  lx = Vec3<Type>(m[0],m[3],m[6]);
	  ly = Vec3<Type>(m[1],m[4],m[7]);
	  lz = Vec3<Type>(m[2],m[5],m[8]);
	}

	/// transpose the matrix, i.e. for rotations: invert it
	inline void transpose()
	{
		Type tmp;
		tmp=m[3]; m[3]=m[1]; m[1]=tmp; // swap(1,3)
		tmp=m[6]; m[6]=m[2]; m[2]=tmp; // swap(2,6)
		tmp=m[7]; m[7]=m[5]; m[5]=tmp; // swap(5,7)
	}

	/// compute determinant of Matrix
	Type det() {
		Type d  = m[0] * ( m[4]*m[8] - m[7]*m[5] );
			 d += m[1] * ( m[3]*m[8] - m[6]*m[5] );
			 d += m[2] * ( m[3]*m[7] - m[6]*m[4] );
		return d;
	}

	/// Multiplication operator
	inline Matrix operator*(const Matrix& m1) const
	{
	  Matrix m2;
	  int i,j;
	  for(i=0; i<3; i++)
		for(j=0; j<3; j++)
		{
			m2.m[i][j] = m[i][0]*m1.m[0][j] + m[i][1]*m1.m[1][j] +
						 m[i][2]*m1.m[2][j] + m[i][3]*m1.m[3][j];
		}
	  return m2;
	}

	inline void print()
	{
		printf("[%.3g %.3g %.3g]\n",m[0], m[1], m[2]);
		printf("[%.3g %.3g %.3g]\n",m[3], m[4], m[5]);
		printf("[%.3g %.3g %.3g]\n",m[6], m[7], m[8]);
	}
	/// print to ostream
	inline void Print(std::ostream &output) const 
	{  
		output << "[" << m[0] << " " << m[1] << " " << m[2] << "]" << endl;
		output << "[" << m[3] << " " << m[4] << " " << m[5] << "]" << endl;
		output << "[" << m[6] << " " << m[7] << " " << m[8] << "]" << endl;
	}
	/// ostream "<<" operator
	friend inline std::ostream & operator << (std::ostream& output, const Matrix &M)
	{
	  M.Print(output);
	  return output; 
	}
};

/// Multiplies a vector by a matrix
template<class Type>inline Vec3<Type> operator*(const Vec3<Type>& v,const Matrix<Type>& m)
{
	Vec3<Type> r;
	r.x = v.x*m.m[0] + v.y*m.m[1] + v.z*m.m[2];
	r.y = v.x*m.m[3] + v.y*m.m[4] + v.z*m.m[5];
	r.z = v.x*m.m[6] + v.y*m.m[7] + v.z*m.m[8];
	return r;
}

/// Multiplies a matrix by a vector 
template<class Type>inline Vec3<Type> operator*(const Matrix<Type>& m, const Vec3<Type>& v)
{
	Vec3<Type> r;
	r.x = v.x*m.m[0] + v.y*m.m[1] + v.z*m.m[2];
	r.y = v.x*m.m[3] + v.y*m.m[4] + v.z*m.m[5];
	r.z = v.x*m.m[6] + v.y*m.m[7] + v.z*m.m[8];
	return r;
}

/// \typedef Matrixf - A float Matrix
typedef Matrix<float> Matrixf;
/// \typedef Matrixd - A double Matrix
typedef Matrix<double> Matrixd;

#endif //MATHCLASSSES_H
