#ifndef RAY_H
#define RAY_H

#include "MathClasses.h"

/**
* \class Ray Ray.h
* \brief A template Ray class used for intersections
*
* This class represents a Ray. It is mainly used for intersection tests
* but it can be also used in the implementation of a Raytracer.
*/
template<class Type>class Ray {
public:
	/// ray position
	Vec3<Type> pos;
	/// ray direction
	Vec3<Type> dir;

    /// constructor - sets position and direction of ray to (0,0,0)
	Ray() { pos=Vec3<Type>::VEC0; dir=Vec3<Type>::VEC0; }
    /// constructor - sets position and direction of ray to given vectors
	Ray(const Vec3<Type> p, const Vec3<Type> d) { pos=p; dir=d; }

	void Print() { pos.Print(); dir.Print(); }
};

/// \typedef Rayf - A float Ray
typedef Ray<float> Rayf;
/// \typedef Rayd - A double Ray
typedef Ray<double> Rayd;

#endif //RAY_H