/** 
* \author Marcel Koertgen (marcel.koertgen@gmail.com)
* \date 02-02-2003
* \version 1.0
*/
#ifndef SAMPLEARR_H
#define SAMPLEARR_H

//#include <windows.h>
//#include <math.h>

#include <stdlib.h> // random numbers
#include <time.h>   // random number seed

#include <gl/GL.h>
#include <gl/GLU.h>
#include "gl/glut.h"

#include "TriMesh.h"

/**
* \def SCANRANDOM - enhance sampling method
* \brief When this Symbol is defined an enhanced sampling strategy will be used.
*
*  Shooting all rays exactly to a bounding-sphere큦 center may disturb a
*  uniform distribution of samples on some models.
*  Imagine a mesh that holds a plane like model. Shooting rays to a bounding sphere큦
*  center will result in samples located all near the center.
*
*  The idea here is to change the ray큦 direction by using an angle threshold <t> instead.
*  <t> defines a cone around the original direction-vector. Then the new direction-vector
*  will be randomly taken out of that cone by using two random numbers between [-t, t]
*
*  Control the sampling by tuning <t> and the radius of the bounding sphere.
*/
#define SCANRANDOM 

/// \def NORMAL_LEN - factor to multiply normals on drawing
#define NORMAL_LEN 5;

/**
*  \class SampleArr SampleArr.h
*  \brief This class holds an array of samples on a mesh큦 bounding sphere.
*
*  An Instance of this class may hold hermite data of a given input "TriMesh",
*  i.e. points and normals on the mesh큦 surface.
*
*  It is also used to (re)sample a given TriMesh, although the sampling
*  strategy used here is very simple:
*
*     1. divide a given bounding sphere into sectors in x and z
*     2. for each sector do:
*	      - on the sphere surface find the center of the sector.
*		  - from there shoot a ray to the sphere큦 center
*		  - if intersection with mesh add intersection point (and normal) to array
*
*  The idea is that for compact models (or say models of topological genus 0) the
*  sampling will lead to roughly uniform distributed samples.
*/
template<class Type>class SampleArr {

public:
	/** This constructor will intialize member data to 0, although a GVertexArray
	*   for samples (and normals) will be instantiated
	*/
	SampleArr() {
		mesh = NULL;
		numSamples = 0;
		center = Vec3<Type>::VEC0;
		radius = (Type)0;
		samples=normals = NULL;
	}

	/** This constructor will intialize member data with respect to given TriMesh*
	*   GVertexArrays for samples (and normals) will be instantiated.
	*   \param mesh - A pointer to instance of class TriMesh.
	*   - The bounding sphere is defined here by the mesh큦 Center Of Mass (C.O.M)
	*   and the corresponding radius.
	*   - The number of sectors will be set to 10.
	*
	*   Use this constructor if you have to sample the mesh before use.
	*/
	SampleArr(TriMesh<Type>* mesh) {
		this->mesh = mesh;
		numSamples = 0;
		if (mesh) {
			center = mesh->vCOM;
			radius = mesh->CRadius;
		} else {
			center=Vec3<Type>::VEC0;
			radius=(Type)0;
		}
		samples=normals = NULL;
	}

	/// Destructor
	~SampleArr() {
		if (samples) delete samples;
		if (normals) delete normals;
	}
	
	/** This method draws the sample Array as simple GL_POINTS.
	* \param drawNormals - If set to 1 normals will also be drawn as GL_LINES.
	*/
	void draw(int drawNormals) {
		if (samples!=NULL &&
			normals!=NULL)
		{
			int i;
			glPushAttrib(GL_ENABLE_BIT);
			glDisable(GL_LIGHTING); // is lighting on ?
			
			if (drawNormals) { // double the loop for speed
				for (i=0; i<numSamples; i++) {
					Vec3<Type> v = samples->getVec3(i);
					  glBegin(GL_LINES);
					  glVertex3f((GLfloat)v.x, (GLfloat)v.y, (GLfloat)v.z);
					  v += normals->getVec3(i) * NORMAL_LEN;
					  glVertex3f((GLfloat)v.x, (GLfloat)v.y, (GLfloat)v.z);
					glEnd();
				}
			} else {
				glBegin(GL_POINTS);
				  for (i=0; i<numSamples; i++) {
					  Vec3<Type> v = samples->getVec3(i);
					  glVertex3f((GLfloat)v.x, (GLfloat)v.y, (GLfloat)v.z);
				  }
				glEnd();
			}
			glPopAttrib();
		}
	}

	/** This method samples the mesh to get samples.
	* \param sectors - number of sectors to divide sphere into
	*\n
	* The strategy used here is very simple:
	*
	*     1. divide a given bounding sphere into sectors in x and z
	*     2. for each sector do:
	*	      - on the sphere surface find the center of the sector.
	*		  - from there shoot a ray to the sphere큦 center
	*		  - if intersection with mesh add intersection point (and normal) to array
	*
	*  The idea is that for compact models (or say models of topological genus 0) the
	*  sampling will lead to roughly uniform distributed samples.
	*/
	void sampleMesh(int sectors) {
		if (mesh) {
			int j,i;
			Type xAng, zAng;
			Ray<Type> ray;
			Vec3<Type> onSphere, toCenter, intersection, normal;
			
			if (samples) delete samples;
			samples = new GVertexArray<Type>();
			if (normals) delete normals;
			normals = new GVertexArray<Type>();

			for (i=0; i<sectors-1; i++) { // divide spherical circle in "sectors" intervals [i,i+1]
				// calculate angle for that position
				// (2*i+1) / sectors is in [0,2] so multiply with PI instead of 2PI
				xAng = (Type)M_PI * (Type)(2*i+1) / (Type)sectors; // middle of interval [i,i+1]
				for (j=0; j<sectors-1; j++) {
					zAng = (Type)M_PI * (Type)(2*j+1) / (Type)sectors; // middle of interval [i,i+1]
					
					onSphere = Vec3<Type>::VECX.Rotate(zAng, Vec3<Type>::VECZ);
					onSphere = onSphere.Rotate(xAng, Vec3<Type>::VECX );
					toCenter = - onSphere;  // Ray direction |n|==1
					onSphere *= radius; // scale position
					onSphere += center;     // translate position

					ray.pos = onSphere;
					ray.dir = toCenter;
					//	toCenter.Print();
					if (mesh->NIS(ray, intersection, normal)) {
						//	printf("mesh intersection.\n");
						normals->addVec3(normal, numSamples);
						samples->addVec3(intersection, numSamples);
					}
				}
			}
		}
	}

	/** This is the Osada method to samples a mesh. The strategy used here
	* is also simple:
	*
	*     1. given the mesh큦 area <A> as the sum of all faces area <ai>, generate
	*        a random number <N> in [0,A).
	*
	*     2. Find the number j such that N is in [aj,aj+1) (Binary Search -> O(log n) )
	*
	*     3. With j found, generate a random point on face/triangle j (barycentric coords):
	*
	*        Let A,B,C be the 3 vertices of the triangle. Any point P inside can be
	*        expressed uniquely as P = aA + bB + cC, where a+b+c=1 and a,b,c are
	*        each >= 0.\n
	*        Knowing a and b permits you to calculate c=1-a-b. So just generate two
	*        random numbers a and b, each in [0,1], such that their sum <=1.
	*        If a+b>1, replace a by 1-a, b by 1-b.\n
	*        Then P = aA + bB + cC is uniformly distributed in triangle ABC.
	*/
	void sampleOsada(int numS) {
		if (mesh) {
			int i,j;
			Vec3<Type> A,B,C,P; // sample point and normal
			Type a,b,c;

			if (samples) delete samples;
			if (normals) delete normals;
			numSamples = numS;
			samples = new GVertexArray<Type>(numSamples);
			normals = new GVertexArray<Type>(numSamples);

			for (j=0; j<numSamples; j++) {
				// generate random number in [0, mesh->area)
				pivot = mesh->area * (Type)rand()/(Type)RAND_MAX;
				// find Triangle to this number
				i = binSearch(0,mesh->numFaces);

				// generate two random numbers in [0,1]
				a = (Type)rand()/(Type)RAND_MAX;
				b = (Type)rand()/(Type)RAND_MAX;
				if (a+b>1) { a=1-a; b=1-b; }
				c = 1-a-b;
				// get triangle vertices
				A = mesh->vertices->getVec3(mesh->faces[i]->v[0]);
				B = mesh->vertices->getVec3(mesh->faces[i]->v[1]);
				C = mesh->vertices->getVec3(mesh->faces[i]->v[2]);

				P = a*A + b*B + c*C;

				samples->setVec3(j, P);
				normals->setVec3(j, mesh->faces[i]->n);
			}
		}
	}

	/// prints sample coordinates (uses printf)
	inline void print()
	{
		printf("SampleArr:\n");
		printf(" numSamples %d\n", numSamples);
		printf(" Center: "); center.print(); printf("\n");
		printf(" Radius: %.3g", Radius);

		if (samples!=NULL &&
			normals!=NULL)
		{
			printf(" Coordinates:\n [");
			samples->getVec3(0).print();
			int i;
			for (i=1; i<numSamples; i++)
			{
				printf(", ");
				samples->getVec3(i).print();
			}
			printf("]\n");
		}
	}
	/// print to ostream
	inline void Print(std::ostream &output) const 
	{  
		output << " numSamples: " << numSamples << endl;
		output << " Center: " << center << endl;
		output << " Radius: " << Radius << endl;
		if (samples!=NULL &&
			normals!=NULL)
		{
			output << " Coordinates:" << endl << "[" <<	samples->getVec3(0);
			int i;
			for (i=1; i<numSamples; i++)
			{
				output << ", " << samples->getVec3(i);
			}
			output << "]" << endl;
		}
	}
	/// ostream "<<" operator
	friend inline std::ostream & operator << (std::ostream& output, const SampleArr &s)
	{
	  s.Print(output);
	  return output; 
	}



	/// TriMesh to sample
	TriMesh<Type>* mesh;          
	/// List of Samples
	GVertexArray<Type>* samples;  
	/// List of Normals
	GVertexArray<Type>* normals;  
	/// the Center of the bounding sphere used for sampling
	Vec3<Type> center;
	/// the Radius of the bounding sphere used for sampling
	Type radius;
	/** numbers of samples obtained by sampling
	* \warning READ-ONLY
	*/
	int numSamples;

private:
	int binSearch(int l, int r) {
		int m = (l+r)/2;
		Type al = mesh->aSum_lp[m],
			 ar = mesh->aSum_lp[m+1];

		if (ar < pivot) // look in right half
			return binSearch(m,r);
		else if (pivot < al) // look in left half
			return binSearch(l,m);
		else return m; // found
		//	if (mesh->aSum_lp[i] <= v && v < mesh->aSum_lp[i+1]) // found
	}

	Type pivot; // pivot value - random number 
};

/// \typedef SampleArrf - A float SampleArr
typedef SampleArr<float> SampleArrf;
/// \typedef SampleArrd - A double SampleArr
typedef SampleArr<double> SampleArrd;

#endif //SAMPLEARR_H