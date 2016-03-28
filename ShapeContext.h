/** 
* \author Marcel Koertgen and Gil-Joo Park
* \date 02-02-2003
* \todo bin calculation in sector model, bin-indexing
* \version 1.0
*/

/**
* \mainpage Local and Global Shape Contexts for Similarity Search and Classification in Spatial Databases
*
* \section intro Introduction
*
* The definition of an appropriate distance function is crucial for the
* effectiveness of any nearest neighbbor classifier. A common approach for
* similarity models is based on the paradigm of feature vectors. A <i>feature
* transform</i>  maps a complex object onto a feature vector in multidimensional
* space. The similarity of two object is then defined as the vicinity of their
* feature vectors in the feature space.<br>
*
* \section General Representation of Shape
*
* We follow this approach by using a general representation of shape -
* a set of points sampled from the contours on the object. They need not, and
* tyipcally will not, correspond to key-points such as maxima of curvature.
* We prefer to sample the shape with roughly uniform spacing, though this is
* also not critical.<br>
* Each point is associated with a novel decsriptor, the <i>local shape context</i>,
* which describes the coarse arrangement of the rest of the shape with respect
* to the point. This descriptor will be different for different points on a
* single shape S; however corresponding (homologous) points on similar shapes
* <b>S</b> and <b>S´</b> will tend to have similar local shape contexts.<br>
* Shape Contexts are distributions and can be compared using statistics.
* Correspondences between the point sets <b>S</b> and <b>S´</b> can be found by
* solving a bipartite weighted graph matching problem with edge weights <b>C(i,j)</b>
* defined by the <b>X²</b>-distances of the shape contexts of points <b>i</b> and <b>j</b>.
* Given correspondences we can calculate a similarity measure between the shapes
* <b>S</b> and <b>S´</b>. This search can be used in a nearest-neighbor classifier
* for object recognition.<br><br>
*
* \section Local Shape Contexts (3D Histograms)
*
* We introduce 3D shape histograms as intuitive feature vectors. In general,
* histograms are based on partitioning of bounding volume in which the object resides.<br>
* We suggest a spiderweb model in log-polar space, i.e. the model is emebedded in
* one of it´s bounding spheres (not necessarily the minimum sphere).<br>
* Then the sphere is partitionend into sectors measured in angles (polar) and
* shells of logarithmic width (log r). For the latter definition it is
* assumed that nearby samples are more important than those further away.<br>
* 
*
* \section Contact
*
* This MSVC files are part of a project on a search engine for 3d models.<br>
* For more details on this project go <a href="http://cg.cs.uni-bonn.de>here</a>
*<br>
*<br>
* You can reach us at <a href="mailto:marcel.koertgen@gmail.com">marcel.koertgen@gmail.com</a>
* and <a href="mailto:parkg@cs.uni-bonn.de">parkg@cs.uni-bonn.de</a>
*
*/
#ifndef SHAPECONTEXT_H
#define SHAPECONTEXT_H

#include "MathClasses.h"
#include "SampleArr.h"
#include "ShapeMatrix.h"

//#define DBINS
/**
* \class ShapeContext ShapeContext.h
* \brief This class holds a ShapeContext to a given TriMesh or Samples.
*
*  An Instance of this class may serve as matching criteria for comparing
*  (rigid) 3d models.
*/
template<class Type>class ShapeContext {
public:
	/// The SampleArray
	SampleArr<Type>* sams;
	/// The number of samples in the SampleArray
	int numSams;
	/// The number of shells to use in all local histograms
	int  numShells;

	/** The number of sectors to use in all local histograms; a bounding sphere will be divided
	* into <numSectorY> half-circles rotated around local Z-axis.
	* Each half-circle will then be divided into <numSectorX> sectors.
	* As a rule of thumb <numSectorX> should be half of <numSectorY>
	*/
	int	 numSectorX,
		 numSectorY;
	/// The number of bins to use in all local histograms; xsectors*ysectors*shells
	int  numBins;

	/// the weights "alpha","gamma" used in local matching
	Type alpha[3],gamma[3];
	
	/** The log-base to use for shell-radii; increase that value to weight nearby samples more
	* than those far away.
	* For rlogbase=1 the shells are equi-distanced
	* For rlogbase=2 shell<i+1> is double the width of shell<i>
	*/
	int rlogbase;
	/** <radius> holds double the radius of the meshes bounding sphere, so with
	* a local radius factor of 1 no local shape will have outliers.
	*
	* 2*meshradius might be an unnecessary large value often resulting in too
	* big local spheres and thus weakening the discriminative power of a local shape
	* descriptor.
	*/
	Type radius;

	Type **histos;
	Matrix<Type> **histosSys;

	/** constructor
	* \param samples - A SampleArray
	* \param shells - number of shells to use
	* \param sectors - number of sectors to use
	* \param sh_rlog - specifies the log-base for shell-radii
	* \param a0 - weights for axe weighting
	* \param a1 - weights for axe weighting
	* \param a2 - weights for axe weighting
	* \param b0 - weights for match term weighting
	* \param b1 - weights for match term weighting
	* \param b2 - weights for match term weighting
	*/
	ShapeContext(SampleArr<Type> *samples=0,
				 const int shells=4,
				 const int xsectors=6,
				 const int ysectors=12,
				 const int sh_rlog=1,
				 const Type a0=1.0/3.0,
				 const Type a1=1.0/3.0,
				 const Type a2=1.0/3.0,
				 const Type g0=1.0/3.0,
				 const Type g1=1.0/3.0,
				 const Type g2=1.0/3.0) {
		numShells=shells;
		numSectorX=xsectors;
		numSectorY=ysectors;
		numBins=xsectors*ysectors*shells;
		rlogbase=sh_rlog;
		alpha[0]=a0; alpha[1]=a1; alpha[2]=a2;
		gamma[0]=g0; gamma[1]=g1; gamma[2]=g2;
		sams=samples;
		if (sams) {
			numSams=sams->numSamples;
			// take double the bounding-radius; so we will get all other samples
			radius=2*sams->radius; 
			histos = new Type*[numSams];
			histosSys = new Matrix<Type>*[numSams];
			int i;
			for (i=0; i<numSams; i++) {
				histos[i]=NULL;
				histosSys[i]=NULL;
			}
		} else {
			numSams=0;
			histos=NULL;
			histosSys=NULL;
			radius=0;
		}
	}

	/// destructor
	~ShapeContext() {
		freeHistos();
	}

	/** computes a local system <sx,sy,sz> for histogram/sample <sInd> -
	/* sx -> (center - sample<sInd>) normalized
	/* sy -> projection of PA[0] (or PA[1] if sx || PA[0]) onto plane with normal <sx>,
	/* sz -> sx^sy
	*/
	void computeSystem(int sInd,
					   Vec3<Type> &sx,
					   Vec3<Type> &sy,
					   Vec3<Type> &sz)
	{
		if (sams!=NULL &&
			sams->mesh!=NULL &&
			sInd>-1 && sInd<numSams) {
			// sx = normal of plane for sy,sz
			sx = (sams->mesh->vCOM - sams->samples->getVec3(sInd));
			sx.Normalize();

/*			Vec3<Type> pa[3];
			int i;

			// compute projection of principal axes
			for (i=0; i<3; i++)
			{
				pa[i] = sams->mesh->PA[i] - (sams->mesh->PA[i]*sx) * sx;
			}

			// take longest projection
			sy = pa[0];
			for (i=1; i<3; i++)
			{
				if (pa[i].LengthSqr() > sy.LengthSqr())
					sy = pa[i];
			}*/

			// sx || PA[0] ?
			Type dot = sams->mesh->PA[0]*sx;
			if ( NONZEROF(1-fabs(dot)) ) // take PA[1] for projection
				sy = sams->mesh->PA[1] - (sams->mesh->PA[1]*sx)*sx;
			else // compute projection of PA[0]
				sy = sams->mesh->PA[0] - dot*sx;

			sy.Normalize();
			// cross-product guarants a right-handed COS
			sz = sx^sy;
		}
	}

	/** computes a local histogram to a given sample
	* \param overwrite - if a local histogram already exists, reompute it?
	* \param sInd - Index of Sample to compute local Shape Histogram for
	* \param lx - Local X Axis for this histogram
	* \param ly - Local Y Axis for this histogram
	* \param lz - Local Z Axis for this histogram
	* \param rad - 0 < rad <= 1.0f, sphere radius
	*/
	void computeHistogram(bool overwrite,
						  const int sInd,
						  const Vec3<Type> &lx=Vec3<Type>::VECX,
						  const Vec3<Type> &ly=Vec3<Type>::VECY,
						  const Vec3<Type> &lz=Vec3<Type>::VECZ,
						  const Type rad=(Type)1.0)
	{
		if (sams!=NULL &&
			numSams>1 &&
			numBins>0 &&
			sInd>=0 &&
			sInd<numSams)
		{
			// return if histogram exists and is not to overwrite
			if (!overwrite &&
				histos[sInd]!=NULL)
				return;

			int i,k;
			Vec3<Type> p = sams->samples->getVec3(sInd); // our sample position
			Type r = radius*rad;                         // our local radius
			
			int* absHisto = new int[numBins]; // count numbers of samples in a bin
			for (k=0; k<numBins; k++) absHisto[k]=0;

			Vec3<Type> vi; // vector i; one of the other n-1 samples
			Type ri;       // vector i´s distance from our sample; used for shell-number calc.

			Matrix<Type> *A = new Matrix<Type>(), // rotation matrix to use for our local system
						 *B = new Matrix<Type>();
			A->set_rotation(lx,ly,lz);
			B->set_inverse_rotation(lx,ly,lz);

			// save the local System for this Sample
			if (histosSys[sInd] != NULL) free(histosSys[sInd]);
			histosSys[sInd] = A;

			// loop through samples (i!= sInd)
			for (i=0; i<numSams; i++)
			{
				if (i!=sInd) //	for sample i compute bin number k
				{
					// let vi refer to local system
					vi = (sams->samples->getVec3(i) - p) * (*B);
					// save distance
					ri = vi.Length();
					// normalize for angle/sector calculation
					vi.Normalize();

					k = getBin(vi,ri,r);

					// increment number if valid index returned
					if (k>-1 && k <numBins) (absHisto[k])++;
#if DBINS
					else
					{
						cerr << " Error: (" << sInd << "/" << i << ")";
						switch(k)
						{
						case -1:
							cerr << " -Sample too far " << ri << "/" << r << endl;
							break;
						case -2:
							cerr << " -Sample too near " << ri << endl;
							break;
						default:
							cerr << " -Invalid Bin-Index (" << k << "). Computation Error!" << endl;
							break;
						}
					}
#endif
				}
			}

			// delete unused Matrix
			delete B;

			// compute relative numbers
			if (histos[sInd]!=NULL)	free(histos[sInd]);
			histos[sInd] = new Type[numBins];
			for (k=0; k<numBins; k++)
				histos[sInd][k] = (Type)absHisto[k] / (Type)(numSams-1);
		}
	}

	/** computes similarity between two local histograms
	* \param B - the second ShapeContex
	* \param i - index of local histogram in this ShapeContext
	* \param j - index of local histogram in ShapeContext B
	* \return a similarity value: 0=exact match .. 1=complete different
	*/
	Type localMatch(ShapeContext* B, int i, int j) {

		if ( (!B) ||                            // no ShapeContext to match
			 (B->numBins != this->numBins) ||   // numBins not equal
			 (this->histos[i]==NULL) ||         // i-th histogram in this not calculated
			 (B->histos[j]==NULL)               // j-th histogram in B not calculated
		   )
		   return (Type)1; // worst matching

		Type st,at,pt,tmp; // Shape-, Appearance-, Position-Term
			 st=at=pt=0;
		int k;

		// compute shape term
		if (fabs(gamma[0]) > EPSF)
		{
			for (k=0; k<numBins; k++)
			{
				tmp = histos[i][k]+B->histos[j][k]; 
				if (tmp>EPSF) // is tmp > 0 ?
					st += SQR(histos[i][k] - B->histos[j][k]) / tmp;
			}
			st /= 2*(Type)numBins;
		}

		// compute appearance term
		if (fabs(gamma[1]) > EPSF)
		{
			Vec3<Type> u = Vec3<Type>(histosSys[i]->m[0],histosSys[i]->m[3],histosSys[i]->m[6]),
					   v = Vec3<Type>(B->histosSys[j]->m[0],B->histosSys[j]->m[3],B->histosSys[j]->m[6]);
			for (k=0; k<3; k++)
				at += alpha[k]*B->alpha[k]* SQR(1-(u*v));
			
			at /= 4*(Type)(B->numSams);
		}

		// compute position term
		if (fabs(gamma[2]) > EPSF)
		{
			Vec3<Type> p, q;
			p = sams->samples->getVec3(i) / sams->radius;       // normalize to [0,1]
			q = B->sams->samples->getVec3(j) / B->sams->radius;

			for (k=0; k<3; k++)
				pt += SQR( this->alpha[k]*p[k] - B->alpha[k]*q[k]);
//				pt += SQR ( p[k] - q[k] );
			pt /= (Type)(B->numSams);
		}

		return (gamma[0]*st + gamma[1]*at + gamma[2]*pt);
	}

	/// Prints histograms (uses printf)
	inline void print() {
		printf("ShapeContext:\n");
		printf(" Samples: %d\n", numSams);
		printf(" x-Sectors: %d\n", numSectorX);
		printf(" y-Sectors: %d\n", numSectorY);
		printf(" Shells: %d\n", numShells);
		printf(" Bins: %d\n", numBins);
		if (sams!=NULL &&
			numSams>1 &&
			numBins>0) {
			int i;
			for (i=0; i<numSams; i++) {
				printf(" Histogram %d: [",i);
				if (histos[i] != NULL) {
					printf("%.2f",histos[i][0]);
					int j;
					for (j=1; j<numBins; j++) printf(", %.2f",histos[i][j]);
				} else {
					printf("NULL");
				}
				printf("]\n");
			}
		}
	}

	/// print to ostream
	inline void Print(std::ostream &output) const 
	{
		output << " Samples: " << numSams << std::endl;
		output << " x-Sectors: " << numSectorX << std::endl;
		output << " y-Sectors: " << numSectorY << std::endl;
		output << " Shells: " << numShells << std::endl;
		output << " Bins: " << numBins << std::endl;
		if (sams!=NULL &&
			numSams>1 &&
			numBins>0) {
			int i;
			for (i=0; i<numSams; i++) {
				output << " Histogram "<< i << ": [";
				if (histos[i] != NULL) {
					output << histos[i][0];
					int j;
					for (j=1; j<numBins; j++)
						output << ", " << histos[i][j];
				} else {
					output << "NULL";
				}
				output << "]" << std::endl;
			}
		}
	}
	/// ostream "<<" operator
	friend inline std::ostream & operator << (std::ostream& output, const ShapeContext<Type> &s)
	{
	  s.Print(output);
	  return output; 
	}

	inline void computeHistograms(Type rfac, bool local_orientation)
	{
		// vectors for a local histogram´s local system
		Vec3<Type> sx=Vec3<Type>::VECX, sy=Vec3<Type>::VECY, sz=Vec3<Type>::VECZ;
		int i;

		for (i=0; i<numSams; i++)
		{
			if (local_orientation)
				this->computeSystem(i, sx, sy, sz);
			this->computeHistogram(true, i, sx,sy,sz, rfac);
		}
	}

	ShapeMatrix<Type>* computeShapeMatrix(ShapeContext<Type> *B, Type rfac, bool local_orientation)
	{
		if (!B || numBins != B->numBins) // equal bins ?
			return NULL;

		// a (quadratic) matrix saving values from local matching
		ShapeMatrix<Type> *A = new ShapeMatrix<Type>(this->numSams, B->numSams);

		// vectors for a local histogram´s local system
		Vec3<Type> sxA=Vec3<Type>::VECX, syA=Vec3<Type>::VECY, szA=Vec3<Type>::VECZ,
				   sxB=Vec3<Type>::VECX, syB=Vec3<Type>::VECY, szB=Vec3<Type>::VECZ;

		int i,j;
		for (i=0; i<numSams; i++) // loop for A (this)
		{
			if (local_orientation)
				this->computeSystem(i, sxA, syA, szA);
			this->computeHistogram(false, i, sxA,syA,szA, rfac);

			for (j=0; j<B->numSams; j++) // loop for B
			{
				if (local_orientation) 
					B->computeSystem(j, sxB, syB, szB);
				B->computeHistogram(false, j, sxB,syB,szB, rfac);
				A->m[i][j] = this->localMatch(B,i,j);
			}
		}
		// return shape matrix
		return A;
	}

	void draw(int sInd, float rfac)
	{
		if (sams!=NULL && 
			sInd>-1 && sInd <numSams &&
			histos[sInd]!=NULL &&
			histosSys[sInd]!=NULL)
		{
//			int i,j,k; // secx, secy, shells
			Vec3<Type> p = sams->samples->getVec3(sInd);

			glPushAttrib(GL_ENABLE_BIT); // save lighting state
			glDisable(GL_LIGHTING); // is lighting on ?
			glPushMatrix();
			glTranslatef(p.x, p.y, p.z);

			Vec3<Type> x(histosSys[sInd]->m[0],
						 histosSys[sInd]->m[3],
						 histosSys[sInd]->m[6]);

			Vec3<Type> y(histosSys[sInd]->m[1],
						 histosSys[sInd]->m[4],
						 histosSys[sInd]->m[7]);

			Vec3<Type> z(histosSys[sInd]->m[2],
						 histosSys[sInd]->m[5],
						 histosSys[sInd]->m[8]);

			x *= rfac*radius;
			y *= rfac*radius;
			z *= rfac*radius;

			// draw C.O.S.
			glColor4f(1,0,0,1);
			glBegin(GL_LINES);
			   glVertex3f(0,0,0);
			   glVertex3f(x.x, x.y, x.z);
			glEnd();
			glColor4f(0,1,0,1);
			glBegin(GL_LINES);
			   glVertex3f(0,0,0);
			   glVertex3f(y.x, y.y, y.z);
			glEnd();
			glColor4f(0,0,1,1);
			glBegin(GL_LINES);
			   glVertex3f(0,0,0);
			   glVertex3f(z.x, z.y, z.z);
			glEnd();
			glPopMatrix();
			glPopAttrib();
		}
	}


private:
	void freeHistos() {
		if (histos!=NULL &&
			histosSys!=NULL)
		{
			int i;
			for (i=0; i<numSams; i++)
			{
				if (histos[i] != NULL) delete [] histos[i];
				if (histosSys[i] != NULL) delete [] histosSys[i];
			}
			delete [] histos;
			delete [] histosSys;
		}
	}

	void cartesian_2_polar(Vec3<Type> v, Type &a, Type &b)
	{
		a = acos(v.z/ v.Length());

		if (NONZEROF(v.x)) // x <> 0
		{
			if (v.x<0)        // x<0 -> b += 180°; b in [-90°, 90°) -> b in [90°, 270°)
				b = M_PI + atan(v.y/v.x);
			else
				if (v.y<0)    // x>0 && y<0 -> b += 360°; b in [-90°, 0°) -> b in [270°, 360°)
					b = M_2PI + atan(v.y/v.x);
				else          // x> 0 && y>= 0; b in [0°, 90°)
					b = atan(v.y/v.x);
		} 
		else // x == 0
		{
			if (NONZEROF(v.y))
				if (v.y < 0) b = M_PI;  // y<0  -> b=180°
				else b = PI2;           // y>=0 -> b=90°
			else // x==y==0; z <> 0
				if (v.z<0)
				{
					a = 0;
					b = M_PI;  // x=y=0; z<0  -> a=0, b=180°
				}
				else
					a = b = 0;          // x=y=0; z>=0 -> a = b = 0°
		}
	}

	int getBin(Vec3<Type> &vi, Type ri, Type r) {
		// clamping outliers
		// r==0 || ri>=r -> outlier; not in any bin
		if ((r-ri) < EPSF) return -1;
		if (ri < EPSF) return -2;

		// angle calculation in local oriented planes.
		// 0 <= a < PI
		// 0 <= b <= 2*PI
		// NOTE: vertex needs to refer to histogram´s local system
		Type a,b;
		cartesian_2_polar(vi, a,b);

		// sector x/y calculate linear from angles
		// todo: half_circles (a) should have only 1/2*numsectors
		int secx = (int)( a*ONEOVERPI * (Type)numSectorX ),
			secy = (int)( b*ONEOVER2PI * (Type)numSectorY ),
			shell;

		/* !NOTE!
		* computing a radius r in [0,1) for shell <i> (i from 0..s-1) is straightforward:
		*
		*   r = log(a, pow(a,s)* (i/s)) / s
		*
		* where
		*   a = logbase
		*   s = number of shells
		*   i = shellindex
		*
		*
		* The value (i/s) (from [0,1]) has to be scaled by a^s, for log(1)=0 and log(<1) < 0.
		*
		* computing index out of a radius requires then the inverse:
		*
		*   i = s * pow(a, s*r-1)
		*/
		if (rlogbase>1)  // shellindex ~ 1- r ; inner shells bigger than outer shells
			shell = (int)( numShells * pow(rlogbase, numShells*((ri/r)-1) ) );
		// for log-base=1 it is linear -> shells are equi-distant
		else
			shell = (int)( (ri/r) * numShells );

		// return strict ordered index
		// - <numShells> bins in each sector
		// - <numSectors> in each circle
		// - <secy> is the circle index, i.e. how often to rotate the original half-circle
		return secx + (secy * numSectorX) + shell * (numSectorX*numSectorY);
		/* !!!NOTE on ordering bins!!!
		* The bin ordering above utilizes a fast access of bins in a distinct shell.
		* For example you can access all bins in the innermost shell by iterating
		* over histos[aSample][0..numSectorX*numSectorY-1]
		*/
	}
};

/// \typedef ShapeContextf - A float ShapeContext
typedef ShapeContext<float> ShapeContextf;
/// \typedef ShapeContextd - A double ShapeContext
typedef ShapeContext<double> ShapeContextd;

#endif //SHAPECONTEXT_H
