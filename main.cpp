/** 
* \author Marcel Koertgen and Gil-Joo Park
* \date 02-02-2003
* \version 1.0
* This file is part of the 3D ShapeContexts Project
*/

//************************** INCLUDES *****************************************************
#include "arg_util.h"        // argument utilities

#include "TriMesh.h"         // handling Triangle Meshes
#include "SampleArr.h"       // sampling Meshes and handling samples in an array
#include "ShapeContext.h"    // handling sets of local shape descriptors
#include "GlobalMatching.h"  // finding of correspondences

#include <time.h> // Benchmarking
clock_t start, finish;
double  duration;

//************************** CONSTANTS / SYMBOLS ******************************************
#define DPRINT       // print Shapematrix and Shapecontexts to log-file
#ifdef DPRINT
 #define logfilename "debug.log"
std::ofstream logfile;
#endif

#define RENDER_GL      // undefine that if you don´t want rendering

#define PAMAP          // undefine if you don´t want a Principal Axes transform
#ifndef PAMAP
  #define OBJRADIUS 200.0
#endif

#define PRECISION float     // floating point type to use for template classes

//************************** VARIABLES ****************************************************

int MATCHING = 2;          // 1 = hungarian, 2=linear assignment, 3=soft assignments

int SAMPLES1 = 200;        // number of samples for mesh 1
int SAMPLES2 = 200;        // number of samples for mesh 2

int SECTORSX =   6;        // number of x-sectors to use for each local shape descriptor
int SECTORSY =  12;        // number of y-sectors to use for each local shape descriptor
int SHELLS   =   6;        // number of shells to use for each local shape descriptor
int RLOG_SHELLS = 3;       // radius log-base for shells <-> weight nearby samples

float RFAC_LOCAL = 0.5;    // radius factor for local descriptors (how far can i see?)
bool XYZ_LOCAL = true;     // should local shapes have an own orientation each?
						   // if so, they are oriented with respect to mesh´s C.O.M.

float GAMMA1 = 1.0;        // weight for distribution term
float GAMMA2 = 0.0;        // weight for appearance term
float GAMMA3 = 0.0;        // weight for position term

TriMesh<PRECISION> *m1=NULL, *m2=NULL;        // two triangle meshes to match
SampleArr<PRECISION> *s1=NULL, *s2=NULL;      // two samplearrays
ShapeContext<PRECISION> *sc1=NULL, *sc2=NULL; // two ShapeContexts
ShapeMatrix<PRECISION> *A=NULL,*B=NULL;       // ShapeMatrix holding n x n loc. match. values
GlobalMatching<PRECISION>* M=NULL;            // class to solve correspondences on the ShapeMatrix

PRECISION GMatchResult = 0.0;

//************************** OPENGL RENDERING *********************************************
#ifdef RENDER_GL
  #define OBJDIST 1.5
  #include "render_gl.h"
#endif

//************************** CLEANUP ******************************************************
// Tried to get this as exit_func-callback for GLUT. Doesn´t work right now.
void exitfunc(int exit_code)
{
	printf("\nCleaning up...\n");

	printf("  Deleting Mesh Instances...");
	if (m1 != NULL) delete m1;
	if (m2 != NULL) delete m2;
	printf("Done.\n");

	printf("  Deleting SampleArray Instances...");
	if (s1 != NULL) delete s1;
	if (s2 != NULL) delete s2;
	printf("Done.\n");

	printf("  Deleting ShapeContext Instances...");
	if (sc1 != NULL) delete sc1;
	if (sc2 != NULL) delete sc2;
	printf("Done.\n");

	printf("  Deleting ShapeMatrix Instances...");
	if (A != NULL) delete A;
	if (B != NULL) delete B;
	printf("Done.\n");

	printf("  Deleting GlobalMatching Instances...");
	if (M != NULL) delete M;
	printf("Done.\n");

	printf("Done.\n");
}

//********************************* METHODS ***********************************************
void usage()
{
	printf("Usage: ShapeContexts.exe <mesh1> <mesh2>\n");
	printf("       where <mesh1> and <mesh2> are filenames to Wavefront OBJ files.\n");
	printf("\n");
	printf("Additional parms:\n");
	printf("  Matching:\n");
	printf("    3. <matching> (int) 1=hungarian,2=shortest augmenting path,3=soft assignments\n");
	printf("    4. <gamma1>    (float) 0.0 - 1.0  (weight for distribution term)\n");
	printf("    5. <gamma2>    (float) 0.0 - 1.0  (weight for appearance term)\n");
	printf("    6. <gamma3>    (float) 0.0 - 1.0  (weight for position term)\n");
	printf("\n");
	printf("  Sampling:\n");
	printf("    7. <samplecount1> (int) 10 - 1000\n");
	printf("    8. <samplecount2> (int) 10 - 1000\n");
	printf("\n");
	printf("  Histograms:\n");
	printf("    9. <sectors in x> (int) 1 - 24\n");
	printf("   10. <sectors in y> (int) 1 - 24\n");
	printf("   11. <shells>       (int) 1 - 24\n");
	printf("   12. <logbase>      (int) 1 - 10\n");
	printf("   13. <radius>       (float) 0.1 - 1.0\n");
	printf("   14. <oriented>     (bool) 0=absolute, 1=compute reference frame\n");
	printf("\n");
	exit(1);
}

void check_args(int argc, char** argv)
{
	if (argc<3)    // check arguments
	{
		usage();
	}
	else
	{
		if (argc>3) 
			MATCHING = argtoi(argv[3],1,3);
		if (argc>4)
			GAMMA1 = argtof(argv[4],0,1);
		if (argc>5)
			GAMMA2 = argtof(argv[5],0,1);
		if (argc>6)
			GAMMA3 = argtof(argv[6],0,1);

		if (argc>7)
			SAMPLES1 = argtoi(argv[7],10,1000);
		if (argc>8)
			SAMPLES2 = argtoi(argv[8],10,1000);

		if (MATCHING<3 && SAMPLES1 != SAMPLES2)
		{
			printf("With matching=%i equal samplecounts are needed!\n",MATCHING);
			printf("Aborting program...\n");
			exit(1);
		}
			
		if (argc>9)
			SECTORSX = argtoi(argv[9],1,24);
		if (argc>10)
			SECTORSY = argtoi(argv[10],1,24);
		if (argc>11)
			SHELLS = argtoi(argv[11],1,24);
		if (argc>12)
			RLOG_SHELLS = argtoi(argv[12],1,10);
		if (argc>13)
			RFAC_LOCAL = argtof(argv[13], 0.1 , 1);
		if (argc>14)
			XYZ_LOCAL = argtob(argv[14], true);
	}
}

void prepare_meshes(const char *fm1, const char *fm2)
{
//--- Prepare Meshes --------------------------------------------------
	printf("\nPreparing Meshes:\n");

	// 1st mesh
	printf("  Mesh 1 (%s): ",fm1);
	start = clock();

	m1 = new TriMesh<PRECISION>();
	if (m1->loadMesh(fm1) == false)
	{
		printf("\nError while loading mesh %s\n",fm1);
		printf("Aborting Program...\n");
		exit(2);
	}
	printf("(V:%i, F:%i)",m1->numVertices, m1->numFaces);

#ifdef PAMAP
	m1->computeCOM();    // compute center of mass, bounding sphere and principal axes
	m1->PATransform();   // rotate mesh so that principal axes and (x,y,z) fall together
	m1->uniform_scale(50.0/m1->CRadius);	// scale model to 50 units
#else
	m1->CRadius = OBJRADIUS;
#endif

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" Done. (%2.4f secs)\n", duration );

	// 2nd mesh
	printf("  Mesh 2 (%s): ",fm2);
	start = clock();

	m2 = new TriMesh<PRECISION>();
	if (m2->loadMesh(fm2) == false)
	{
		printf("\nError while loading mesh %s\n",fm2);
		printf("Aborting Program...\n");
		exit(3);
	}
	printf("(V:%i, F:%i)",m2->numVertices, m2->numFaces);

#ifdef PAMAP
	m2->computeCOM();    // compute center of mass, bounding sphere and principal axes
	m2->PATransform();   // rotate mesh so that principal axes and (x,y,z) fall together
	m2->uniform_scale(50.0/m2->CRadius);	// scale model to 50 units
#else
	m2->CRadius = OBJRADIUS;
#endif

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" Done. (%2.4f secs)\n", duration );
}

void sample_meshes()
{
//--- Mesh Sampling --------------------------------------------------
	srand(time(NULL)); // initialize random number generator
	// mesh 1
	printf("  Sampling mesh 1 (%i samples)...",SAMPLES1);
	start = clock();
	
	s1 = new SampleArr<PRECISION>(m1),
	s1->sampleOsada(SAMPLES1);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" Done. (%2.8f secs)\n", duration );

	// mesh 2
	printf("  Sampling mesh 2 (%i samples)...",SAMPLES2);
	start = clock();

	s2 = new SampleArr<PRECISION>(m2);
	s2->sampleOsada(SAMPLES2);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" Done. (%2.8f secs)\n", duration );
}

void compute_histos()
{
//--- Shape Contexts ----------------------------------------------
	sc1 = new ShapeContext<PRECISION>(s1,SHELLS,SECTORSX,SECTORSY,RLOG_SHELLS,
									  m1->alpha[0],m1->alpha[1],m1->alpha[2],
									  GAMMA1, GAMMA2, GAMMA3);

//	printf(" a1 = %.2f, a2 = %.2f, a3 = %.2f\n", m1->alpha[0],m1->alpha[1],m1->alpha[2]);

	sc2 = new ShapeContext<PRECISION>(s2,SHELLS,SECTORSX,SECTORSY,RLOG_SHELLS,
									  m2->alpha[0],m2->alpha[1],m2->alpha[2],
									  GAMMA1, GAMMA2, GAMMA3);
//	printf(" a1 = %.2f, a2 = %.2f, a3 = %.2f\n", m1->alpha[0],m1->alpha[1],m1->alpha[2]);

 // pre-compute histograms to avoid O(k) in inner loop
	// ShapeContext 1
	printf("  ShapeContext 1 : Compute %i Histograms...",SAMPLES1);
	start = clock();
	sc1->computeHistograms(RFAC_LOCAL, XYZ_LOCAL);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" Done. (%2.4f secs)\n", duration );

	// ShapeContext 2
	printf("  ShapeContext 2 : Compute %i Histograms...",SAMPLES2);
	start = clock();
	sc2->computeHistograms(RFAC_LOCAL, XYZ_LOCAL);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" Done. (%2.4f secs)\n", duration );
}

void compute_shapematrix()
{
//--- ShapeMatrix ------------------
	printf("\n\nLocal Matching:\n");
	printf("  %i Samples(1) <-> %i Samples(2)\n",SAMPLES1,SAMPLES2);
	printf("  %i*%i Sectors * %i Shells = %i Bins\n",SECTORSX,SECTORSY,SHELLS,sc1->numBins);
	printf("  Shell Radius log-base: %i\n",RLOG_SHELLS);
	printf("  Local Radius factor  : %2.2f\n", RFAC_LOCAL);
	printf("  Local Orientation: ");
	if (XYZ_LOCAL) printf("YES\n"); else printf("NO\n");
	printf("\n");
	// compute (quadratic) matrix (SAMSxSAMS) with local shape terms(i,j) as elements
	// this takes n1*n2*k time, where ni is sample-count and k is bin-count
	printf("  Computing %i x %i ShapeMatrix ...",SAMPLES1,SAMPLES2);
	start = clock();

	A = sc1->computeShapeMatrix(sc2, RFAC_LOCAL, XYZ_LOCAL);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" Done. (%2.4f secs)\n", duration );

#ifdef DPRINT
	logfile << "ShapeMatrix A: " << std::endl;
	logfile << *A << std::endl;
	logfile << "ShapeContext 1:" << std::endl;
	logfile << *sc1 << std::endl;
	logfile << "ShapeContext 2:" << std::endl;
	logfile << *sc2 << std::endl;
#endif
}

void do_global_matching()
{

//--- Global Matching --------------------------------------------
	printf("\n\nGlobal Matching:\n");
	// save original matrix for it will be destroyed by global matching
	B = new ShapeMatrix<PRECISION>(1,1);
	(*B) = (*A);

	M = new GlobalMatching<PRECISION>();

	start = clock();

	switch(MATCHING)
	{
		case 1: 
			printf("  Hungarian... "); // Hungarian method: O(n*n*n), 1-1
			GMatchResult = M->HungarianMethod(B);
			break;
		case 2:
			printf("  Hard Linear Assignments... ");
			GMatchResult = M->LinearAssignmentMethod(B);
			break;
		case 3:
			printf("  Soft Linear Assignments... ");
			GMatchResult = M->MultiAssignmentMethod(B);
			break;

		default:
			printf("\nError! Not a valid matching method (1,2,3) %d\n",MATCHING);
			printf("Aborting Program...\n");
			exit(4);
	}

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" Done (%2.4f secs). Result = %2.8f\n", duration, GMatchResult );
//-------------------------------------------------------------------------
}

//********************************* MAIN **************************************************
void main(int argc, char** argv)
{
#ifdef DPRINT
	logfile.open(logfilename);
	if (!logfile)
		std::cerr << "Error: Could not create " << logfilename << std::endl;
#endif

	printf("\n\nMatching 3D models with 3D ShapeContexts\n\n");   // say hello

	check_args(argc, argv);
	
	prepare_meshes(argv[1], argv[2]);
	sample_meshes();

	compute_histos();

	compute_shapematrix();

	do_global_matching();

#ifdef RENDER_GL
	printf("\nNow entering GLUT...\n");
	start_rendering(argc, argv);
#else
	exitfunc(0);
#endif
}