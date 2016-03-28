#ifndef FACE_H
#define FACE_H

#include "GVertexArray.h"
#include "Ray.h"

/**
* \class Face Face.h
* \brief This class represents a facette in a TriMesh.
*
* A Face consists of 3 vertex indices, a normal vector and it´s area.
*/
template<class Type>class Face {
public:
	/// Vertex indices
	int v[3];
	/// Normal Vector
	Vec3<Type> n;    
	/// (unsigned) area of the Triangle
	Type area;

	/// Constructor
	Face() {
		v[0]=v[1]=v[2] = 0;
		n = Vec3<Type>::VEC0;
		area = (Type)0;
	}

	/// Copy Constructor
	Face(const Face<Type> &face) {
		v[0] = face.v[0];
		v[1] = face.v[1];
		v[2] = face.v[2];
		n = face.n;
		area = face.area;
	}

	/** Constructor
	* \param va - pointer to VertexArray where vertices are stored
	* \param i0 - index to vertex 0
	* \param i1 - index to vertex 1
	* \param i2 - index to vertex 2
	*/
	Face(Vec3<Type> v0, Vec3<Type> v1, Vec3<Type> v2, int i0, int i1, int i2) {
		v[0] = i0;
		v[1] = i1;
		v[2] = i2;
		n = (v1 - v0) ^ (v2 - v0);
		area = (Type)0.5 * n.Length();
		n.Normalize();
	}

	Face(std::vector<Vec3<Type>>& vertices, int i0, int i1, int i2)
	: Face(vertices[i0], vertices[i1], vertices[i2], i0,i1,i2)
	{
	}

	/// Destructor
	~Face() {}

	/** Finds an intersection between the Face´s triangle and a given Ray.
	* \param va - vertices for normal computation
	* \param ray - The incoming ray.
	* \param point - return intersection here, if any
	* \return int -\n 
	* -3 = face is degenerate, i.e. a segment or point (normal could not be computed)\n
	* -2 = ray disjoint from triangle\n
	* -1 = ray goes away from triangle\n
	*  0 = intersect with plane but not inside triangle\n
	*  1 = intersect in unique point\n
	*  2 = ray and Face are in the same plane
	*/
	int findIntersection(GVertexArray<Type>* va, Ray<Type> ray, Vec3<Type> &point) {
		if (!va) return -3;

		Vec3<Type>  U = va->getVec3(v[1]) - va->getVec3(v[0]),
			  V = va->getVec3(v[2]) - va->getVec3(v[0]);

		Vec3<Type> w0, w;       // ray vectors
		Type r, a, b;    // params to calc ray-plane intersect
		
		if (n == Vec3<Type>::VEC0) return -3;   // triangle is degenerate

		w0 = ray.pos - va->getVec3(v[0]);
		a = n * w0;
		b = n * ray.dir;

		if (ISZERO(b))      // ray is parallel to triangle plane
		{
			if (a < EPS)    // ray lies in triangle plane
				return 2;
			else            // ray disjoint from plane
				return -2;
		}

		r = -a / b;     // get intersect point of ray with triangle plane
		if (r < 0.0) return -1;   // ray goes away from triangle
		// for a segment, also test if (r > 1.0) => no intersect
		point = ray.pos + r * ray.dir; // intersect point of ray and plane
		
		// is point inside Triangle?
		Type uu, uv, vv, wu, wv, D;
		uu = U*U;
		uv = U*V;
		vv = V*V;
		w = point - va->getVec3(v[0]);
		wu = w*U;
		wv = w*V;
		D = uv*uv - uu*vv;
		
		// get and test parametric coords
		Type s, t;
		s = (uv * wv - vv * wu) / D;
		if (s < 0.0 || s > 1.0) return 0; // I is outside T
		t = (uv * wu - uu * wv) / D;
		if (t < 0.0 || (s + t) > 1.0) return 0; // I is outside T
		
		return 1;                      // I is in T
	}
};
#endif //FACE_H