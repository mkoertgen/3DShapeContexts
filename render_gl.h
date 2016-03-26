//#define SCREENSHOT
#ifdef SCREENSHOT 
  #include "glbmp\glbmp.h"
  CGLBMP screenBMP;
  int screenNum = 0;
#endif 

char str[256];  // window title, etc.

bool drawMesh = true;        // draw Meshes ?
bool drawBox = false;        // draw Bounding-Box?
bool drawCOM = true;         // draw C.O.M.
bool drawSamples = true;     // draw Samples
bool drawLocalShape = true;  // draw local Histograms
int SampleIndex = 0;         // index of Histogram to look at

int	Width = 800;				// window Width in Pixels
int	Height = 600;				// window Height in Pixels

float rotX=0.0f,    // camera rotation around x-axis
	  rotY=0.0f,    // camera rotation around y-axis
	  camdist=0.0f, // camera distance
	  zoom_fac,     // camera zoom-factor
	  transX=0.0f,  // camera look at X
	  transY=0.0f;  // camera look at Y

bool fullscreen = false;      // fullscreen?
bool lighting = true;         // lighting?
bool blending = true;         // blending?
GLenum polygonMode = GL_FILL; // polygon filling or outlining?
bool drawText = true;         // drawtext?

GLfloat mDiff[4] = {0.7f, 0.7f, 0.7f, 0.5f}; // gray mesh-material
GLfloat bDiff[4] = {0.8f, 0.0f, 0.0f, 0.1f}; // red Bbox-material
GLfloat cDiff[4] = {0.0f, 0.0f, 0.8f, 0.1f}; // blue C.O.M-material

// Mousecontrol
int			  mButton = -1;       // current button
int           mOldY, mOldX;       // old Mouse position
enum { LEFT = 0, RIGHT, MIDDLE }; // button states

// variables to save window state when going to fullscreen
GLint			oldwinpos_x=20;          
GLint			oldwinpos_y=20;          
GLint			oldwindowwidth=Width;	   
GLint			oldwindowheight=Height;  

//---------- GLUT Methods ------------------------------------------------

// setup blending, polygonmode, lighting and lighting model
void InitGlut()
{
  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ONE);
//  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Initialisiert die Blendfunktion

  glDisable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);		// Legt die Art des Tiefentestes fest
  glDisable(GL_CULL_FACE);	// Deaktiviert Backface Culling
  glPolygonMode (GL_FRONT_AND_BACK, polygonMode);
  glShadeModel(GL_SMOOTH);	// Waehlt smooth shading aus

  GLfloat modelAmb[4] = {0.1f, 0.1f, 0.1f,  1.f};
  GLfloat diffuse[4]  = {0.6f, 0.6f, 0.6f,  1.f};
  GLfloat specular[4] = {0.3f, 0.3f, 0.3f,  1.f};
  GLfloat pos[4]      = {0.f, 0.f, 0.f, 1.f};

  // position a light
  glEnable(GL_LIGHTING);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, modelAmb);
  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, pos);
  glEnable(GL_LIGHT0);
}

// draws a text in 3d
void draw_text_3d(float x, float y, float z, float alpha, char *string)
{
	int len, i;

	glPushAttrib(GL_ENABLE_BIT); // save GL_LIGHTING state
	glDisable(GL_LIGHTING);      // disable lighting for points
	glColor4f(alpha,alpha,alpha,1);
	glRasterPos3f(x, y, z);
	len = (int) strlen(string);
	for (i = 0; i < len; i++) {
	  glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
	}
	glPopAttrib(); // restore GL_LIGHTING state
}

// GLUT reshape-callback
static void Reshape(int width, int height)
{
  glViewport(0, 0, (GLint)width, (GLint)height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, 1.0, 1.0, 2000.0);
  glMatrixMode(GL_MODELVIEW);
}

// GLUT key-callback
static void Key(unsigned char key, int x, int y)
{
  switch (key) {
  case 27:     // quit program
	  exit(0);
	  break;
  case 'a':    // toggle sample-drawing
	  drawSamples = !drawSamples;
	  break;
  case 'b':    // toggle bbox-drawing
	  drawBox = !drawBox;
	  break;
  case 'c':    // toggle com-drawing
	  drawCOM = !drawCOM;
	  break;
  case 'm':     // toggle mesh-drawing
	  drawMesh = !drawMesh;
	  break;
  case 's' :    // toggle drawing of local shapes
	  drawLocalShape = !drawLocalShape;
	  break;
  case '+':    // look at next local shape
	  SampleIndex++;
	  if (SampleIndex>sc1->numSams-1) SampleIndex=0;
	  break;
  case '-':    // look at previous local shape
	  SampleIndex--;
	  if (SampleIndex<0) SampleIndex=sc1->numSams-1;
	  break;
//
  case 'l':    // toggle lighting
	  lighting = !lighting;
	  if (lighting)
		  glEnable(GL_LIGHTING);
	  else
		  glDisable(GL_LIGHTING);
	  break;
  case 'n':    // toggle blending
	  blending = !blending;
	  if (blending)
	  {
		  glEnable(GL_BLEND);
		  glDisable(GL_DEPTH_TEST);
	  }
	  else
	  {
		  glDisable(GL_BLEND);
		  glEnable(GL_DEPTH_TEST);
	  }
	  break;
  case 'p':    // toggle polygonmode
	  if (polygonMode==GL_FILL) polygonMode=GL_LINE;
	  else polygonMode=GL_FILL;
	  glPolygonMode (GL_FRONT_AND_BACK, polygonMode);
	  break;
  case 't':    // toggle drawText
	  drawText = !drawText;
	  break;
//
  default: return;
  }
  glutPostRedisplay ();
}

// Print out usable keys in this program
void PrintKeys(void) 
{
	printf("\nKeys:\n");
	printf("----------------------------------------------\n");
	printf("A:           Toggle Samples drawing\n");
	printf("B:           Toggle Bounding Box drawing\n");
	printf("C:           Toggle Center of Mass drawing\n");
	printf("M:           Toggle Mesh drawing\n");
	printf("S:           Toggle Local Shapes drawing\n");
	printf("+:           Next Local ShapeContext\n");
	printf("-:           Previous Local ShapeContext\n");
	printf("\n");
	printf("N:           Toggle Blending\n");
	printf("L:           Toggle Lighting\n");
	printf("P:           Toggle Polygonmode\n");
	printf("T:           Toggle Text Drawing\n");
	printf("\n");
	printf("F1:          Toggle Fullscreen\n");
#ifdef SCREENSHOT
	printf("F2:          Save Screenshot as BMP\n");
#endif
	printf("\n");
	printf("Mouse Left:   Rotate\n");
	printf("Mouse Right:  Zoom\n");
	printf("Cursor Left:  Translate Left\n");
	printf("Cursor Right: Translate Right\n");
	printf("Cursor Up:    Translate Up\n");
	printf("Cursor Down:  Translate Down\n");
	printf("Pos1/Home:    Reset Look Settings\n");
	printf("\n");
	printf("ESC:         Quit Program\n");
}

// GLUT Special-Key callback
static void SpecialKey(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_F1:
		if (fullscreen==false)		// toggle fullscreen
		{
			fullscreen=true;
			oldwinpos_x=glutGet(GLUT_WINDOW_X);
			oldwinpos_y=glutGet(GLUT_WINDOW_Y);
			oldwindowwidth=glutGet(GLUT_WINDOW_WIDTH);
			oldwindowheight=glutGet(GLUT_WINDOW_HEIGHT);
			glutFullScreen();
			break;
		}  
		else 
		{
			glutPositionWindow(oldwinpos_x,oldwinpos_y);
			glutReshapeWindow(oldwindowwidth,oldwindowheight); 
			fullscreen=false;
			break; 
		}
#ifdef SCREENSHOT
	case GLUT_KEY_F2:
		sprintf(str, "screen%d.bmp", screenNum);
		screenBMP.SetJPGQuality(100);
		screenBMP.SaveScreen(str);
		screenNum++;
		break;
#endif
	case GLUT_KEY_LEFT:  transX -= 5; break;
	case GLUT_KEY_RIGHT: transX += 5; break;
	case GLUT_KEY_UP:    transY += 5; break;
	case GLUT_KEY_DOWN : transY -= 5; break;
	case GLUT_KEY_HOME :
		transX=transY=rotX=rotY=0.0f;
		camdist = (m1->CRadius+m2->CRadius)*3.0;
		break;
	default: return;
	}
	glutPostRedisplay();
}

// GLUT Handlemouse-callback: save button state and mouse position
void HandleMouse(int button, int state, int x, int y) 
{
	if (state == GLUT_DOWN)
	{
		mOldX = x;
		mOldY = y;

		switch (button)  
		{
			case GLUT_LEFT_BUTTON: mButton = LEFT; break;
			case GLUT_RIGHT_BUTTON: mButton = RIGHT; break;
			case GLUT_MIDDLE_BUTTON: mButton = MIDDLE; break;
		}
	}
	else if (state == GLUT_UP) mButton = -1;
}

// GLUT Mousemotion-callback: handle mouse motion (cam rotation and zoom)
void glutMouseMotion(int x, int y) 
{
	if (mButton == LEFT) // rotation
	{
		rotY += (x - mOldX);		// cam around Y
		rotX += (y - mOldY);		// cam around X
		glutPostRedisplay();
	}
	else if (mButton == MIDDLE) 
	{
		glutPostRedisplay();			// no middle button right now
	} 
	else if (mButton == RIGHT) // zoom
	{
		camdist -= (y - mOldY)*zoom_fac;	// cam moves on Z
		glutPostRedisplay();
	} 
	// save mouse pos
	mOldX = x; 
	mOldY = y;
}

// GLUT Display-callback: do all rendering
void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  glTranslatef(transX,transY,-camdist);
  glRotatef (rotY, 0,1,0);
  glRotatef (rotX, 1,0,0);

  if (drawBox) {
	  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, bDiff); // red
	  glColor4f(1, 0, 0, 1);
	  glPushMatrix();
		glTranslatef(-m1->CRadius*OBJDIST,0,0);
		m1->drawBBox();
	  glPopMatrix();
	  glPushMatrix();
		glTranslatef(m2->CRadius*OBJDIST,0,0);
		m2->drawBBox();
	  glPopMatrix();
  }

  if (drawCOM) {
	  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, cDiff);
	  glColor4f(0, 0, 1, 1);
	  glPushMatrix();
		glTranslatef(-m1->CRadius*OBJDIST,0,0);
		m1->drawCOM();
	  glPopMatrix();
	  glPushMatrix();
		glTranslatef(m2->CRadius*OBJDIST,0,0);
		m2->drawCOM();
	  glPopMatrix();
  }

  if (drawMesh) {
	  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mDiff);
	  glColor4f(0.3, 0.3, 0.3, 1);
	  glPushMatrix();
		glTranslatef(-m1->CRadius*OBJDIST,0,0);
		m1->draw();
	  glPopMatrix();
	  glPushMatrix();
		glTranslatef(m2->CRadius*OBJDIST,0,0);
		m2->draw();
	  glPopMatrix();
  }

  if (drawSamples) {
	  glColor4f(1, 0, 1, 1); // purple
	  glPushMatrix();
		glTranslatef(-m1->CRadius*OBJDIST,0,0);
		s1->draw(TRUE);
	  glPopMatrix();
	  glPushMatrix();
		glTranslatef(m2->CRadius*OBJDIST,0,0);
		s2->draw(TRUE);
	  glPopMatrix();
  }

  if (drawLocalShape) {
	  glColor4f(1, 1, 1, 1);
	  glPushMatrix();
		glTranslatef(-m1->CRadius*OBJDIST,0,0);

		Vec3<PRECISION> p  = sc1->sams->samples->getVec3(SampleIndex);
/*		, 
						px,py,pz;

		sc1->histosSys[SampleIndex]->get_rotation_axes(px,py,pz);
		px *= RFAC_LOCAL*sc1->radius;
		py *= RFAC_LOCAL*sc1->radius;
		pz *= RFAC_LOCAL*sc1->radius;

 draw yz - plane of histogram(SampleIndex)
		glPushAttrib(GL_ENABLE_BIT); // save lighting state
		glDisable(GL_LIGHTING); // is lighting on ?
		glColor4f(0.3,0.3,0.3,0.1);
		glBegin(GL_QUADS);
		   glVertex3f(p.x - py.x-pz.x, p.y - py.y-pz.y, p.z - py.z-pz.z);
		   glVertex3f(p.x - py.x+pz.x, p.y - py.y+pz.y, p.z - py.z+pz.z);
		   glVertex3f(p.x + py.x+pz.x, p.y + py.y+pz.y, p.z + py.z+pz.z);
		   glVertex3f(p.x + py.x-pz.x, p.y + py.y-pz.y, p.z + py.z-pz.z);
		glEnd();
//		glPopAttrib();
*/
		sc1->draw(SampleIndex, RFAC_LOCAL);
		
		if (drawText)
		{
			sprintf_s(str, "%i",SampleIndex);
			draw_text_3d(p.x,p.y, p.z, 1, str);
		 }
/*
// draw lines to all other points
		glPushAttrib(GL_ENABLE_BIT); // save lighting state
		glDisable(GL_LIGHTING); // is lighting on ?
		glColor4f(0,1,1,0.1);
		glBegin(GL_LINES);
			for (int i=0; i<sc1->numSams; i++)
				if (i != SampleIndex)
				{
					glVertex3f(p.x,p.y,p.z);
					px = sc1->sams->samples->getVec3(i);
					glVertex3f(px.x,px.y,px.z);
				}
		glEnd();
		glPopAttrib();
*/
	  glPopMatrix();

#if MATCHING == 3
	  glPushAttrib(GL_ENABLE_BIT); // save lighting state
	  glDisable(GL_LIGHTING); // is lighting on ?
	  glColor4f(0,1,1,0.1);
	  int i;
/*
	  glPushMatrix();
	  glTranslatef(-m1->CRadius*OBJDIST,0,0);	   
	  Vec3<PRECISION> p1  = sc1->sams->samples->getVec3(SampleIndex), p2;
   
	  glBegin(GL_LINES);  
		  for (i=0; i<M->k[SampleIndex]; i++)
		  {
			glColor4f(M->w[SampleIndex][i],M->w[SampleIndex][i],M->w[SampleIndex][i],1);
			glVertex3f(p1.x,p1.y,p1.z);
			p2  = sc2->sams->samples->getVec3(M->minassign2[SampleIndex][i]);
			glVertex3f(p2.x+(m1->CRadius+m2->CRadius)*OBJDIST,p2.y,p2.z);
		  }
	  glEnd();
	  glPopAttrib();	  
	  glPopMatrix();
*/

	  glPushMatrix();
		glTranslatef(m2->CRadius*OBJDIST,0,0);

		for (i=0; i<M->k[SampleIndex]; i++)
		{
			sc2->draw(M->minassign2[SampleIndex][i], RFAC_LOCAL);
			if (drawText)
			{
				Vec3<PRECISION> p = sc2->sams->samples->getVec3(M->minassign2[SampleIndex][i]);
				sprintf_s(str, "%i",M->minassign2[SampleIndex][i]);
				draw_text_3d(p.x,p.y, p.z,M->w[SampleIndex][i], str);
			}
		}

	  glPopMatrix();
#else
	  glPushMatrix();
		glTranslatef(m2->CRadius*OBJDIST,0,0);
		sc2->draw(M->minassign[SampleIndex], RFAC_LOCAL);
		if (drawText)
		{
			Vec3<PRECISION> p = sc2->sams->samples->getVec3(M->minassign[SampleIndex]);
			sprintf_s(str, "%i",M->minassign[SampleIndex]);
			draw_text_3d(p.x,p.y, p.z, 1, str);
		 }
	  glPopMatrix();
#endif

  }

  glutSwapBuffers();
  glFlush ();
}

void start_rendering(int argc, char** argv)
{
	PrintKeys();

// rendering stuff
	camdist = (m1->CRadius+m2->CRadius)*3.0;
	zoom_fac = -(m1->CRadius+m2->CRadius)/10.0;

	glutInit (&argc, argv);
	glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(Width, Height);

//	if (__glutCreateWindowWithExit("ShapeContexts", &exitfunc) == GL_FALSE)
	if (glutCreateWindow("ShapeContexts") == GL_FALSE)
	{
		printf("\nError while creating GLUT window...\n");
		exit(1);
	}
	sprintf_s(str, "Matching 3D Models with 3D ShapeContexts @ %s",glGetString(GL_RENDERER));
	glutSetWindowTitle(str);

	InitGlut();
	glutDisplayFunc(display);
	glutReshapeFunc (Reshape);

	glutKeyboardFunc(Key);
	glutSpecialFunc(SpecialKey);

	glutMouseFunc(HandleMouse);
	glutMotionFunc(glutMouseMotion);

	glutMainLoop ();
}
