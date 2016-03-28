/** 
* \author Marcel Koertgen (marcel.koertgen@gmail.com)
* \date 02-02-2003
* \version 1.0
*/

//#define PRINT_DEBUG

#ifndef TRIMESH_H
#define TRIMESH_H

#include <windows.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include "GL/glut.h"
#include <fstream>
#include <sstream>
#include <vector>

#include "MathClasses.h"
#include "Ray.h"
#include "Face.h"
#include "GVertexArray.h"

/**
*  \class TriMesh TriMesh.h
*  \brief holds a Triangle Mesh
*
*  An Instance of this class may hold a Triangle mesh:
*  Features:
*    - Vertices
*	 - Halfedges and Edges
*	 - Faces and their normals  
*    - Bounding Box
*    - Bounding Sphere around Geometric Center
*	 - Bounding Sphere Center Of Mass
*	 - Calculation of C.O.M and Principal Axes (Covariance method)
*    - Principal Axes Transform
*	 - Loads from and saves to Wavefront OBJ files (ASCII format)
*/
template<class Type>class TriMesh {
public:
	/// vertex list pointer
	GVertexArray<Type>* vertices; 
	/// number of vertices
	int numVertices;
	/// face list
	Face<Type> **faces;
	/// number of faces
	int numFaces;           
	/// area of mesh -> sum of all face areas
	Type area;
	/// area sum - lookup table
	Type* aSum_lp;

	/// Center Of Mass
	Vec3<Type> vCOM;

	/// Principal Axes
	Vec3<Type> PA[3];

	/// barycentric weights for the principal axes
	Type alpha[3];

	/// Radius of Sphere around Center Of Mass
	Type CRadius;   

	/// Geometric Center
	Vec3<Type> vGCenter;     
	/// Radius of Sphere around Geometric Center
	Type SRadius;         
	/// Bounding Box Vectors
	Vec3<Type> vBCenter, vBMin, vBMax;

	/// constructor
	TriMesh() {
		vertices = NULL;
		faces = NULL;
		aSum_lp = NULL;
		vCOM=vBCenter=vGCenter=vBMin=vBMax = Vec3<Type>::VEC0;
		PA[X]=Vec3<Type>::VECX;
		PA[Y]=Vec3<Type>::VECY;
		PA[Z]=Vec3<Type>::VECZ;
		SRadius=CRadius=area=(Type)0;
		alpha[0]=alpha[1]=alpha[2]=1.0/3.0;
	}

	/// destructor
	~TriMesh() {
		if (faces) deleteFaces();
		if (vertices) delete vertices;
	}

	/// assignment operator
	TriMesh& operator = (const TriMesh &q) 
	{
		if (faces) deleteFaces();
		if (vertices) delete vertices;

		numFaces = q.numFaces;
		numVertices = q.numVertices;

		faces = allocFaces(numFaces);
		aSum_lp = new Type[numFaces+1];
		aSum_lp[0] = (Type)0;
		vertices = new GVertexArray<Type>(numVertices);				

		int i;

		// create a copy of the vertices
		Vec3<Type> v;
		for (i=0; i<numVertices; i++)
		{
			v = q.vertices->getVec3(i);
			vertices->setVec3(i, v);
		}

		// create a copy of the faces
		Face<Type> *f;
		for (i=0; i<numFaces; i++)
		{
			f = q.faces[i];
			faces[i] = new Face<Type>(vertices, f->v[0], f->v[1], f->v[2]); // add face
			aSum_lp[i+1] = q.aSum_lp[i+1];
		}

		vGCenter = q.vGCenter;
		vBCenter = q.vBCenter;
		vBMin = q.vBMin;
		vBMax = q.vBMax;

		vCOM = q.vCOM;
		PA[X] = q.PA[X];
		PA[Y] = q.PA[Y];
		PA[Z] = q.PA[Z];
		SRadius = q.SRadius;
		CRadius = q.CRadius;
		area = q.area;
		for (i=0; i<3; i++)
			alpha[i]=q.alpha[i];
		return *this;
	}


	/// draws the Trimesh (OpenGl)
	inline void draw() {
		if (vertices!=NULL &&
			faces!=NULL &&
			numFaces>0)
		{
			int i;
			Vec3<Type> n, v;
			glBegin(GL_TRIANGLES);
			  for (i=0; i<numFaces; i++) {
				  n = faces[i]->n;
				  glNormal3f((GLfloat)n.x, (GLfloat)n.y, (GLfloat)n.z);
				  
				  v = vertices->getVec3(faces[i]->v[0]);
				  glVertex3f((GLfloat)v.x, (GLfloat)v.y, (GLfloat)v.z);
				  v = vertices->getVec3(faces[i]->v[1]);
				  glVertex3f((GLfloat)v.x, (GLfloat)v.y, (GLfloat)v.z);
				  v = vertices->getVec3(faces[i]->v[2]);
				  glVertex3f((GLfloat)v.x, (GLfloat)v.y, (GLfloat)v.z);
			  }
			glEnd();
		}
	}

	/// draws the Trimesh´s Bounding Box
	inline void drawBBox() {
		if (vertices) {
			glBegin(GL_QUADS); // draw box-faces
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMax.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMax.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMin.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMin.y, (GLfloat)vBMax.z); // front

			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMax.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMax.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMin.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMin.y, (GLfloat)vBMin.z); // back
			  
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMax.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMax.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMax.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMax.y, (GLfloat)vBMax.z); // up

			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMin.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMin.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMin.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMin.y, (GLfloat)vBMax.z); // down

			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMax.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMax.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMin.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMin.y, (GLfloat)vBMax.z); // left
			  
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMax.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMax.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMin.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMin.y, (GLfloat)vBMax.z); // right
			glEnd();

			glPushAttrib(GL_ENABLE_BIT); // save GL_LIGHTING state
			glDisable(GL_LIGHTING);      // disable lighting for points
			glBegin(GL_POINTS); // draw box-corners
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMin.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMin.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMax.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMin.x, (GLfloat)vBMax.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMin.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMin.y, (GLfloat)vBMax.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMax.y, (GLfloat)vBMin.z);
			  glVertex3f((GLfloat)vBMax.x, (GLfloat)vBMax.y, (GLfloat)vBMax.z);
			  glVertex3f( (GLfloat)vBCenter.x, (GLfloat)vBCenter.y, (GLfloat)vBCenter.z); // draw box center
			glEnd();
			glPopAttrib(); // restore GL_LIGHTING state
		}
	}

	/// draws the Trimesh´s Bounding Sphere (geometric center)
	inline void drawBSphere() {
		if (vertices!=NULL &&
			SRadius>0)
		{
			// draw sphere around geometric center
			glPushMatrix();
			  glTranslatef( (GLfloat)vGCenter.x, (GLfloat)vGCenter.y, (GLfloat)vGCenter.z);
			   glutSolidSphere((GLdouble)SRadius, 100, 100);
			glPopMatrix();

			// draw some sphere points
			glPushAttrib(GL_ENABLE_BIT);
			glDisable(GL_LIGHTING); // is lighting on ?
			glBegin(GL_POINTS);
			  glVertex3f( (GLfloat)(vGCenter.x+SRadius), (GLfloat)vGCenter.y          , (GLfloat)vGCenter.z);
			  glVertex3f( (GLfloat)vGCenter.x          , (GLfloat)(vGCenter.y+SRadius), (GLfloat)vGCenter.z);
			  glVertex3f( (GLfloat)vGCenter.x          , (GLfloat)vGCenter.y          , (GLfloat)(vGCenter.z+SRadius));
			  glVertex3f( (GLfloat)(vGCenter.x-SRadius), (GLfloat)vGCenter.y          , (GLfloat)vGCenter.z);
			  glVertex3f( (GLfloat)vGCenter.x          , (GLfloat)(vGCenter.y-SRadius), (GLfloat)vGCenter.z);
			  glVertex3f( (GLfloat)vGCenter.x          , (GLfloat)vGCenter.y          , (GLfloat)(vGCenter.z-SRadius));
			  glVertex3f( (GLfloat)vGCenter.x          , (GLfloat)vGCenter.y          , (GLfloat)vGCenter.z );
			glEnd();
			glPopAttrib();
		}
	}

	/// draws the Trimesh´s Bounding Sphere (center of mass)
	inline void drawCOM() {
		if (vertices!=NULL &&
			CRadius>0)
		{
			// draw sphere around center of mass
			glPushMatrix();
			  glTranslatef( (GLfloat)vCOM.x, (GLfloat)vCOM.y, (GLfloat)vCOM.z);
//			  glutSolidSphere((GLdouble)CRadius, 100, 100);

			  Vec3<Type> x=CRadius*PA[X],
						 y=CRadius*PA[Y],
						 z=CRadius*PA[Z];
			  // draw PCA axes
			  glPushAttrib(GL_ENABLE_BIT); // save GL_LIGHTING state
			  glDisable(GL_LIGHTING);      // disable lighting for points
			  glBegin(GL_LINES);
				glColor4f(1,0,0,1);
				glVertex3f(0.0f, 0.0f, 0.0f);
				glVertex3f((GLfloat)x.x, (GLfloat)x.y, (GLfloat)x.z);
			  glEnd();
			  glBegin(GL_LINES);
				glColor4f(0,1,0,1);
				glVertex3f(0.0f, 0.0f, 0.0f);
				glVertex3f((GLfloat)y.x, (GLfloat)y.y, (GLfloat)y.z);
			  glEnd();
			  glBegin(GL_LINES);
				glColor4f(0,0,1,1);
				glVertex3f(0.0f, 0.0f, 0.0f);
				glVertex3f((GLfloat)z.x, (GLfloat)z.y, (GLfloat)z.z);
			  glEnd();
			glPopMatrix();

			// draw some sphere points
/*			glBegin(GL_POINTS);
			  glVertex3f( (GLfloat)(vCOM.x+CRadius), (GLfloat)vCOM.y          , (GLfloat)vCOM.z);
			  glVertex3f( (GLfloat)vCOM.x          , (GLfloat)(vCOM.y+CRadius), (GLfloat)vCOM.z);
			  glVertex3f( (GLfloat)vCOM.x          , (GLfloat)vCOM.y          , (GLfloat)(vCOM.z+CRadius));
			  glVertex3f( (GLfloat)(vCOM.x-CRadius), (GLfloat)vCOM.y          , (GLfloat)vCOM.z);
			  glVertex3f( (GLfloat)vCOM.x          , (GLfloat)(vCOM.y-CRadius), (GLfloat)vCOM.z);
			  glVertex3f( (GLfloat)vCOM.x          , (GLfloat)vCOM.y          , (GLfloat)(vCOM.z-CRadius));
			  glVertex3f( (GLfloat)vCOM.x          , (GLfloat)vCOM.y          , (GLfloat)vCOM.z );
			glEnd();*/
			glPopAttrib();
		}
	}

	std::string get_file_contents(const char *filename)
	{
		std::ifstream in(filename, std::ios::in | std::ios::binary);
		if (in)
		{
			std::string contents;
			in.seekg(0, std::ios::end);
			contents.resize(in.tellg());
			in.seekg(0, std::ios::beg);
			in.read(&contents[0], contents.size());
			in.close();
			return(contents);
		}
		throw(errno);
	}

	/// loads a mesh from given filename
	bool loadMesh(const char* fileName) {
		std::string contents = get_file_contents(fileName);
		std::istringstream stream(contents);

		std::vector<Vec3<Type>> vVertices;
		std::vector<Face<Type>> vFaces;
		Vec3<Type> v;
		std::vector<int> faceIndices;
		Face<Type> face;

		vBMax = Vec3<Type>(-1000000,-1000000,-1000000); // initialize BBox
		vBMin = Vec3<Type>( 1000000, 1000000, 1000000);

		std::string line;
		int lineNo = 0;
		while (!safeGetline(stream, line).eof())
		{
			lineNo++;
			if (line.empty() || line[0] == '#' || line[1] != ' ') continue;
			switch(line[0])
			{
				case 'v':
					v = parseVector(line);
					v.UpdateMinMax(vBMin,vBMax); // update BBox
					vVertices.push_back(v);
					break;
				case 'f':
					parseFace(line, faceIndices);
					switch (faceIndices.size())
					{
						case 4: // quad
							face = Face<Type>(vVertices, faceIndices[2], faceIndices[3], faceIndices[0]);
							vFaces.push_back(face);
						case 3: // triangle
						default:
							face = Face<Type>(vVertices, faceIndices[0], faceIndices[1], faceIndices[2]);
							vFaces.push_back(face);
							break;
					}
					break;
			}
		}


		if (vertices) delete vertices;
		numVertices = vVertices.size();
		vertices = new GVertexArray<Type>(numVertices);
		for (int i = 0; i < numVertices; i++)
			vertices->setVec3(i, vVertices[i]);

		vBCenter = (Type)0.5 * (vBMax + vBMin);

		if (faces) deleteFaces();
		numFaces = vFaces.size();
		faces = allocFaces(numFaces);
		area = 0;
		aSum_lp = new Type[numFaces];
		aSum_lp[0] = (Type)0;
		for (int i = 0; i < numFaces; i++)
		{
			face = vFaces[i];
			faces[i] = new Face<Type>(face);
			area += face.area;
			aSum_lp[i + 1] = area;
		}

		return true;
	}

	/// saves current mesh to given filename
	bool saveMesh(const char* fileName) {
		std::ofstream file(fileName);
		if (!file) {
			std::cerr << "error at file create!" << fileName << "\n"; 
			return false;
		}

		file << "# Created by class <TriMesh.h>" << std::endl;
		file << "# Author: M.Körtgen" << std::endl;
		file << "# Feel free to use." << std::endl;
		file << std::endl;

		int i;
		Vec3<Type> v;

		// write vertices
		for (i=0; i<numVertices; i++)
		{
			v = vertices->getVec3(i);
			file << "v " << v.x << " " << v.y << " " << v.z << std::endl;
		}
		file << "# " << numVertices << "Vertices" << std::endl;

		// write faces
		Face<Type> *f;
		for (i=0; i<numFaces; i++)
		{
			f = faces[i];
			file << "f " << f->v[0]+1 << " " << f->v[1]+1 << " " << f->v[2]+1 << std::endl;
		}
		file << "# " << numFaces << "Faces" << std::endl;

		return true;
	}


	/// returns true if there is any intersection with a given Ray
	int IS(Ray<Type> ray) {
		if (vertices!=NULL &&
			numFaces>0)
		{
			Vec3<Type> p;
			int i;
			for (i=0; i<numFaces; i++)
				if (faces[i]->findIntersection(vertices, ray, p) == 1) return TRUE;
		}
		return FALSE;
	}

	/** returns nearest intersection with a given Ray
	* \param ray - incoming ray
	* \param point - return intersection point here, if any
	* \param normal - return intersection normal here, if any
	* \return int - 1 if intersection found, 0 otherwise
	*/
	int NIS(Ray<Type> ray, Vec3<Type> &point, Vec3<Type> &normal) {
		if (vertices!=NULL &&
			numFaces>0)
		{
			double dmin = 1000000, d;
			Vec3<Type> p;
			int fis;
			int i;

			point = ray.pos + dmin * ray.dir.Normalize(); // maximum distance
			normal = Vec3<Type>::VEC0;

			for (i=0; i<numFaces; i++) {
				fis = faces[i]->findIntersection(vertices, ray, p);
				if (fis==1) {
					d = (ray.pos - p).Length();
					if (d < dmin) {
						dmin = d;
						point = p;
						normal = faces[i]->n;
					}
				}
			}
			if (dmin<1000000) return TRUE;
		}
		return FALSE;
	}

	/// computes Max,Min and Center bectors for Bounding Box
	void computeBBox(){ 	
		if (vertices!=NULL &&
			numVertices>0)
		{
			int i;
			vBMax = Vec3<Type>(-1000000,-1000000,-1000000);
			vBMin = Vec3<Type>( 1000000, 1000000, 1000000);
			for (i=0; i<numVertices; i++) vertices->getVec3(i).UpdateMinMax(vBMin, vBMax);
			vBCenter = (Type)0.5 * (vBMax + vBMin);
		}
	}

	/// computes Geometric Center
	void computeGCenter() { 	
		if (vertices!=NULL &&
			numVertices>0)
		{
			int i;
			vGCenter = Vec3<Type>::VEC0;
			for (i=0; i<numVertices; i++) vGCenter += vertices->getVec3(i);
			vGCenter /= (Type)numVertices;
		}
	}

	/// computes Minimum Radius to Geometric Center
	void computeSRadius() { 	
		if (vertices!=NULL &&
			numVertices>0)
		{
			Type d;
			int i;
			SRadius = (Type)0;
			
			for (i=0; i<numVertices; i++) {
				d = (vGCenter - vertices->getVec3(i)).Length();
				if (d > SRadius) SRadius = d;
			}
		}
	}

	/// computes Center Of Mass and Principal Axes
	void computeCOM() {  	
		if (vertices!=NULL &&
			numVertices>0)
		{
			int i;
			Vec3<Type> v(0,0,0);
			Matrix<Type> CoVar;
			
			CoVar.null();
			vCOM = v;

			// save triangle moments for speed reasons
			GVertexArray<Type> *tmp = new GVertexArray<Type>(numFaces);

			// compute C.O.M and save triangle centroids for speed reasons
			for (i=0; i<numFaces; i++)
			{
				v = faces[i]->area *
					(vertices->getVec3(faces[i]->v[0])+
					 vertices->getVec3(faces[i]->v[1])+
					 vertices->getVec3(faces[i]->v[2]));

				tmp->setVec3(i, v);
				vCOM += v;
			}
			vCOM /= 3 * this->area;

			// compute CoVariance-Matrix - (note: symmetric)
			for (i=0; i<numFaces; i++)
			{
				v = tmp->getVec3(i)-vCOM;
				CoVar.m[0] += SQR(v.x);
				CoVar.m[1] += v.x*v.y;
				CoVar.m[2] += v.x*v.z;
				CoVar.m[4] += SQR(v.y);
				CoVar.m[5] += v.y*v.z;
				CoVar.m[8] += SQR(v.z);

				tmp->setVec3(i, v); // save substraction (refer to vCOM)
			}
			// use matrix symmetry
			CoVar.m[3] = CoVar.m[1];
			CoVar.m[6] = CoVar.m[2];
			CoVar.m[7] = CoVar.m[5];

			// normalize matrix (not necessary -> difference in scale is wiped by normalization)
			for (i=0; i<9; i++)	CoVar.m[i] /= (9*area*area);

			// diagoalize matrix and find eigenvectors
#ifdef PRINT_DEBUG
			printf("\nC.O.M. : "); vCOM.print(); printf("\n");
			printf("CoVar  :\n"); CoVar.print(); printf("\n");
#endif
			// diagonalize CoVar to get PA-Axes
			Type EVal[3];
			i = Jacobi(CoVar, EVal, PA);
#ifdef PRINT_DEBUG
			printf("\n Jacobi solved in %i Rotations\n",i);
#endif
#ifdef PRINT_DEBUG
			printf("The eigenvectors:\n");
			printf("\t [%f]  [%f]  [%f]\n", PA[0].x, PA[1].x, PA[2].x);
			printf("\t [%f], [%f], [%f]\n", PA[0].y, PA[1].y, PA[2].y);
			printf("\t [%f]  [%f]  [%f]\n", PA[0].z, PA[1].z, PA[2].z);
			printf("The eigenvalues: (%f) (%f) (%f)\n", EVal[0],EVal[1],EVal[2]);
#endif

			for (i=0; i<3; i++)	PA[i].Normalize(); // normalize eigenvectors
//--- find heaviest axis ---------------
// !!!odd number of flips -> left-handed coordinate system!!!
//
// Each axis, regarded as a normal, defines a plane origined at vCOM.
// For each axis sum positive positive dot-products, i.e. vertices
// being on the "right-handed" side of the respective plane.
			Type c[6],  // sum of positive and negative dot-products for each axis
				 dot;   // one dot-product
			int	j;

			for (j=0; j<6; j++) c[j]=0; // initialize sums with 0

			// sum dot-products
			for (i=0; i<numFaces; i++)
				for (j=0; j<3; j++)
				{
					dot = PA[j]*tmp->getVec3(i);
					if (dot>0)
						c[2*j] += dot;
					else
						c[2*j+1] -= dot;
				}
			for (i=0; i<6; i++) c[i] /= 3*area;

			delete tmp;	// no longer needed

#ifdef PRINT_DEBUG
			printf("\n dot-sums:\n");
			printf("(%f, %f)\n", c[0],c[1]); 
			printf("(%f, %f)\n", c[2],c[3]); 
			printf("(%f, %f)\n", c[4],c[5]); 
#endif
			// flip axes if necessary
			for (j=0; j<3; j++)
			{
				if (c[2*j] < c[2*j+1])
				{
					dot=c[2*j]; c[2*j]=c[2*j+1]; c[2*j+1]=dot;
					PA[j] = -PA[j];
				}
			}

			// sort axes
			if (c[0] < c[2]) // x-axis weaker than y-axis
			{
				dot=c[0]; c[0]=c[2]; c[2]=dot;
				PA[0].swap(PA[1]);
			}
			if (c[0] < c[4]) // x-axis weaker than z-axis
			{
				dot=c[0]; c[0]=c[4]; c[4]=dot;
				PA[0].swap(PA[2]);
			}
			if (c[2] < c[4]) // y-axis weaker than z-axis
			{
				dot=c[2]; c[2]=c[4]; c[4]=dot;
				PA[1].swap(PA[2]);
			}
			// cross-product guarants a right-handed system
			PA[2] = PA[0]^PA[1];

			/// compute weights alpha[i]
			Type total = c[0]+c[2]+ c[4];
			alpha[0] = c[0] / total;
			alpha[1] = c[2] / total;
			alpha[2] = c[4] / total;

#ifdef PRINT_DEBUG
			printf("\nThe PA-Axes:\n");
			printf("\t [%f]  [%f]  [%f]\n", PA[0].x, PA[1].x, PA[2].x);
			printf("\t [%f], [%f], [%f]\n", PA[0].y, PA[1].y, PA[2].y);
			printf("\t [%f]  [%f]  [%f]\n", PA[0].z, PA[1].z, PA[2].z);
#endif
//--------------------------------------------
			computeCRadius();
		}
	}

	/// computes Minimum Radius to Center Of Mass
	void computeCRadius() {
		if (vertices!=NULL &&
			numVertices>0)
		{
			Type d;
			int i;
			CRadius = (Type)0;

			for (i=0; i<numVertices; i++) {
				d = (vCOM - vertices->getVec3(i)).Length();
				if (d > CRadius) CRadius = d;
			}
		}
	}

	/// translates and rotates all vertices to PA-system
	void PATransform() {
		if (vertices!=NULL &&
			numVertices>0)
		{
			Matrix<Type> mat;
			// invert matrix, i.e. transform PA to (ex,ey,ez)
			mat.set_inverse_rotation(PA[X], PA[Y], PA[Z]);

			// transform normals with (transposed of inverted) Matrix, i.e. same matrix
			int i;
			for (i=0; i<numFaces; i++)
				faces[i]->n = faces[i]->n*mat;

			// transform vertices
			Vec3<Type> v;
			for (i=0; i<numVertices; i++) {
				v = (vertices->getVec3(i)-vCOM)*mat;
				vertices->setVec3(i, v);
			}

			// bmin/max needs to be recomputed because they are axis-aligned
			computeBBox();

			// bsphere
			vGCenter = (vGCenter-vCOM)*mat;

			// clear values
			vCOM = Vec3<Type>::VEC0;  // COM is now <0,0,0>
			PA[X] = Vec3<Type>::VECX; // and PA is now <X,Y,Z>
			PA[Y] = Vec3<Type>::VECY;
			PA[Z] = Vec3<Type>::VECZ;
		}
	}

	void uniform_scale(Type f)
	{
		if (vertices!=NULL &&
			numVertices>0)
		{
			int i;
			Vec3<Type> v;

			// scale all vertices
			for (i=0; i<numVertices; i++) {
				v = vertices->getVec3(i);
				vertices->setVec3(i, v*f);
			}

			// scale bounding volumes
			// bbox
			vBMin *= f;
			vBMax *= f;
			vBCenter = 0.5*(vBMax+vBMin);
			vCOM *= f;

			CRadius *= f;
			SRadius *= f;
		}
	}

private:
	std::istream& safeGetline(std::istream& is, std::string& t)
	{
		t.clear();

		// The characters in the stream are read one-by-one using a std::streambuf.
		// That is faster than reading them one-by-one using the std::istream.
		// Code that uses streambuf this way must be guarded by a sentry object.
		// The sentry object performs various tasks,
		// such as thread synchronization and updating the stream state.

		std::istream::sentry se(is, true);
		std::streambuf* sb = is.rdbuf();

		for (;;) {
			int c = sb->sbumpc();
			switch (c) {
			case '\n':
				return is;
			case '\r':
				if (sb->sgetc() == '\n')
					sb->sbumpc();
				return is;
			case EOF:
				// Also handle the case when the last line has no line ending
				if (t.empty())
					is.setstate(std::ios::eofbit);
				return is;
			default:
				t += (char)c;
			}
		}
	}

	Vec3<Type> parseVector(std::string& line)
	{
		std::vector<std::string> tokens;
		split(tokens, line, " ");
		Type vx = (Type)stod(std::string(tokens[1]));
		Type vy = (Type)stod(std::string(tokens[2]));
		Type vz = (Type)stod(std::string(tokens[3]));
		return Vec3<Type>(vx, vy, vz);
	}

	void parseFace(std::string& line, std::vector<int>& indices)
	{
		std::vector<std::string> tokens;
		split(tokens, line, " ");

		indices.clear();
		for (std::vector<int>::size_type i = 1; i != tokens.size(); i++)
		{
			std::vector<std::string> tokens2;
			split(tokens2, tokens[i], "/");
			int index = stoi(std::string(tokens2[0]));
			indices.push_back(--index);
		}
	}

	void split(std::vector<std::string> & theStringVector,  /* Altered/returned value */
		const std::string  & theString,
		const std::string  & theDelimiter)
	{
		size_t  start = 0, end = 0;
		while (end != std::string::npos)
		{
			end = theString.find(theDelimiter, start);

			// If at end, use length=maxLength.  Else use length=end-start.
			theStringVector.push_back(theString.substr(start,
				(end == std::string::npos) ? std::string::npos : end - start));

			// If at end, use start=maxSize.  Else use start=end+delimiter.
			start = ((end > (std::string::npos - theDelimiter.size()))
				? std::string::npos : end + theDelimiter.size());
		}
	}

	/// deletes Facelist, called by destructor
	void deleteFaces() {
		if (faces) {
			int i;
			for (i=0; i<numFaces; i++)
				if (faces[i]) delete faces[i];
			delete faces;
			numFaces=0;
			faces = NULL;
		}
		delete aSum_lp;
	}

	/// allocates Facelist, called by "loadMesh()" 
	Face<Type>** allocFaces(int numFaces) {
		Face<Type>** result;
		int i;
		result = new Face<Type>*[numFaces];
		for (i=0; i < numFaces; i++) result[i] = NULL;
		return result;
	}

	// calculate Eigen-values and -vectors from a symmetrical 3x3 matrix - Jacobi-method
	int Jacobi(Matrix<Type> &Mat, Type *eval, Vec3<Type> *evec)
	{
		Type A[3][3], V[3][3];
		Type b[3];
		Type z[3];
		
		double  c,g,h,s,sm,t,tau,theta,tresh;
		int     i,j,ip,iq, NROT;

		for (ip=0; ip<3; ip++)
		{  //initialize V to identity matrix
			for (iq=0; iq<3; iq++)
			{
				A[ip][iq] = Mat.m[ip*3+iq];
				V[ip][iq]=0;
			}
			V[ip][ip]=1;
		}  
		
		for (ip=0; ip<3; ip++)
		{
			b[ip]=A[ip][ip];	
			eval[ip]=b[ip];
			z[ip]=0;    
		}
		NROT=0;

		for (i=1; i<=50; i++)
		{
			sm=0;
			for (ip=0; ip<2; ip++)    // sum off-diagonal elements
				for (iq=ip+1; iq<3; iq++)
					sm += fabs(A[ip][iq]);
				
			if (sm==0) // sum is zero? then return number of rotations
			{
				for (ip=0; ip<3; ip++)
					evec[ip] = Vec3<Type>(V[ip][0], V[ip][1], V[ip][2]);
				return NROT;
			}

			if (i < 4)
				tresh=0.2*sm*sm;
			else
				tresh=0;
			
			for (ip=0; ip<2; ip++)
			{
				for (iq=ip+1; iq<3; iq++)
				{
					g = 100*fabs(A[ip][iq]);
					// after 4 sweeps, skip the rotation if the off-diagonal element is small
					if ((i > 4) &&
						(fabs(eval[ip])+g == fabs(eval[ip])) &&
						(fabs(eval[iq])+g == fabs(eval[iq])))
						A[ip][iq]=0;
					else
						if (fabs(A[ip][iq]) > tresh)
						{
							h = eval[iq]-eval[ip];
							if (fabs(h)+g == fabs(h))
								t=A[ip][iq]/h;
							else
							{
								theta=0.5*h/A[ip][iq];
								t=1/(fabs(theta)+sqrt(1.0+theta*theta));
								if (theta < 0)
									t=-t;
							}
							c =1.0 / sqrt(1.0+t*t);
							s = t*c;
							tau = s / (1.0+c);
							h = t*A[ip][iq];
							z[ip] -= h;
							z[iq] += h;
							eval[ip] -= h;
							eval[iq] += h;
							A[ip][iq]=0;
							for (j=0; j<ip; j++)
							{
								g=A[j][ip];
								h=A[j][iq];
								A[j][ip] = g-s*(h+g*tau);
								A[j][iq] = h+s*(g-h*tau);
							}
							for (j=ip+1; j<iq; j++)
							{
								g=A[ip][j];
								h=A[j][iq];
								A[ip][j] = g-s*(h+g*tau);
								A[j][iq] = h+s*(g-h*tau);
							}
							for (j=iq+1; j<3; j++)
							{
								g=A[ip][j];
								h=A[iq][j];
								A[ip][j] = g-s*(h+g*tau);
								A[iq][j] = h+s*(g-h*tau);
							}
							for (j=0; j<3; j++)
							{
								g=V[j][ip];
								h=V[j][iq];
								V[j][ip] = g-s*(h+g*tau);
								V[j][iq] = h+s*(g-h*tau);
							}
							NROT++;
						} //end ((i.gt.4)...else if
				} // main iq loop
			} // main ip loop
			
			for (ip=0; ip<3; ip++)
			{
				b[ip] += z[ip];
				eval[ip]=b[ip];
				z[ip]=0;
			}
		} //end of main i loop
#ifdef PRINT_DEBUG
		printf("\n 50 iterations !\n");
#endif

		for (ip=0; ip<3; ip++)
			evec[ip] = Vec3<Type>(V[ip][0], V[ip][1], V[ip][2]);
		return NROT;  // too many iterations
	}

};

/// \typedef TriMeshf - A float TriMesh
typedef TriMesh<float> TriMeshf;
/// \typedef TriMeshd - A double TriMesh
typedef TriMesh<double> TriMeshd;

#endif //TRIMESH_H