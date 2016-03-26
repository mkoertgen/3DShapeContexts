#ifndef G_VERTEX_ARRAY_H
#define G_VERTEX_ARRAY_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <memory.h>
#include "MathClasses.h"

/**
* \class GVertexArray GVertexArray.h
* \brief Used for traversing vertices in a TriMesh.
*/
template<class Type> class GVertexArray {
public:
	/// Constructor (default)
	GVertexArray() {
		in=max=0;
		package = 10;
		vert_list = NULL;
		*this = GVertexArray(0);
	}

	/** Constructor
	* \param anzahlVertices - number of GVertexes to reserve space for
	*/
	GVertexArray(int anzahlVertices) {
		in = anzahlVertices;
		max = 0;
		package = 10;

		if (in > 0)	{	// allocate [3*in]-Array of <Type>
			vert_list = new Type[3*in];
			if (vert_list == NULL) {
				std::cerr << "Memory could not be allocated"; 
				exit(0);
			} else {
				max = in;
				memset(vert_list, 0, in * 3 * sizeof(Type));
			}
		}
		else vert_list = NULL; // No Vertices; allocate empty list
	}

	/// Copy Constructor
	GVertexArray(const GVertexArray& Arr) {
		package = Arr.package;
		in = Arr.in;
		max = Arr.max;

		if (max > 0) {
			vert_list = new Type[3*max];
			if (vert_list == NULL) {
				std::cerr << "Memory could not be allocated"; 
				exit(0);
			} else memcpy(vert_list, Arr.vert_list, max * 3 * sizeof(Type));
		} else vert_list = NULL;
	}

	/// Destructor
	~GVertexArray() {
			if (vert_list != NULL) free(vert_list);
	}

	/** Add a Vertex to this array
	* \param p - position of this GVertex
	* \param ind - return index of the added vector here
	* \return int - 1 if successful, 0 if not (out of memory, for example)
	*/
// Add the <Type> components of a new vertex to the list of vertices. 
// Return TRUE on success, otherwise FALSE and nothing has been 
// changed. 
// Argument:
//   vert: points to the vertex components.
//         It is presupposed, that the components of "vert" and the
//         components, used in the actual GVertexArray-object are equal.
//         Order of components: 
//         coords, intensity or color, normal, texture map, alpha, MD.
// Reason for this specific format: HP-Starbase uses this (-> vert)
// format. IRIX-Gl uses no double array.
	int addVertex(const std::vector<Type>& p, int& ind) {
		Type dummy[3];
		dummy[0] = p[0];
		dummy[1] = p[1];
		dummy[2] = p[2];
		return addVert(dummy, ind);
	}

	/** Set a Vertex coordinates
	* \param ind - index of Vertex to change
	* \param P - <Type> vector of coordinates (0:x, 1:y, 2:z)
	* \return int - 1 if successful, 0 if not (invalid index, for example)
	*/
	int setVertex(int ind, const std::vector<Type>& P) {
		int indexvalid = 1;
		if (-1 < index &&
			index < in)
		{
			vert_list[3*index] = p[0];
			vert_list[3*index+1] = p[1];
			vert_list[3*index+2] = p[2];
		}
		else
			indexvalid = 0;
		return indexvalid;
	}

	/** Add a Vertex to this array
	* \param p - new Vector for this GVertex
	* \param ind - resulting index of the added vector
	* \return int - 1 if successful, 0 if not (out of memory, for example)
	*/
	int addVec3(const Vec3<Type> &p, int& ind) {
		Type dummy[3];
		dummy[0] = p.x;
		dummy[1] = p.y;
		dummy[2] = p.z;
		return addVert(dummy, ind);
	}

	/** Set a Vertex coordinates
	* \param ind - index of Vertex to change
	* \param P - a Vector.
	* \return int - 1 if successful, 0 if not (invalid index, for example)
	*/
	int setVec3(int ind, Vec3<Type> &P) {
		int indexvalid = 1;
		if (-1 < ind &&
			ind < in)
		{
			vert_list[3*ind] = P.x;
			vert_list[3*ind+1] = P.y;
			vert_list[3*ind+2] = P.z;
		}
		else
			indexvalid = 0;
		return indexvalid;
	}

	/** Get a Vertex coordinates
	* \param ind - index of Vertex to
	* \return std::vector<Type> - <Type> vector of coordinates (0:x, 1:y, 2:z)
	* will be (0,0,0) for invalid indices
	*/
	std::vector<Type> getVertex(int ind) {
		std::vector<Type> tmp(3);
		if(ind < in) {
			tmp[0] = vert_list[3 * ind];
			tmp[1] = vert_list[3 * ind + 1];
			tmp[2] = vert_list[3 * ind + 2];
		} else {
			tmp[0]=tmp[1]=tmp[2]=0;
		}
		return tmp;
	}

	/** Get a <Type> pointer to a specific Vertex
	* \param i - index of Vertex
	* \return <Type>* - pointer
	* will be NULL for invalid indices
	*/
	Type* getVertexPtr(int i) {
		if (i >= 0 &&
			i <= in)
			return vert_list + 3 * i;
		else
			return NULL;
	}

	/** Get a Vertex coordinates
	* \param ind - index of Vertex to
	* \return Vec3<Type> - a Vector
	* will be (0,0,0) for invalid indices
	*/
	Vec3<Type> getVec3(int ind) {
		Vec3<Type> tmp = Vec3<Type>::VEC0;
		if(ind < in) {
			tmp.x = vert_list[3*ind];
			tmp.y = vert_list[3*ind+1];
			tmp.z = vert_list[3*ind+2];
		}
		return tmp;
	}

	/** Get number of Vertices in this array
	* \return int - number of Vertices
	*/
	int getSize() const	{ return in; }


	/** Get a <Type> pointer to the vertex list
	* \return <Type>* - pointer
	*/
	const Type* getVertexList() const { return vert_list; }

private:
	/// add a Vertex and return index
	int addVert(const Type* const, int& ind) {
		int Ok = 1;
		if (in + 1 >= max) {
			Type* tmp;
				if (vert_list != NULL) {
					tmp = (double*)realloc(vert_list, (max+package) * sizeof(Type) * 3);
					memset(tmp + max * 3, 0, package * sizeof(Type) * 3);
				} else {
					tmp = (Type*)malloc((max + package) * sizeof(Type) * 3);
					memset(tmp, 0, (max + package) * sizeof(Type) * 3);
				}
			if (tmp) {
				vert_list = tmp;
				max += package;
			} else {
				std::cerr << "Out of memory";
				Ok = 0;
			}
		}
		if (Ok) {
			memcpy(vert_list + in * 3, vert, 3 * sizeof(Type));
			ind = in;
			in++;
		}
		return Ok;
	}
	/// the vertex list as raw data
	Type* vert_list;
	
	/// Vertices in this array
	int in;
	
	/// maximum number of vertices that is space allocated for
	int max;
	
	/// increment for extension of list
	int package;
};

/// \typedef GVertexArrayf - a float GVertexArray
typedef GVertexArray<float> GVertexArrayf;
/// \typedef GVertexArrayD - a double GVertexArray
typedef GVertexArray<double> GVertexArrayd;

#endif //G_VERTEX_ARRAY_H