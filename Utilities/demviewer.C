#include <string>
#include <cstdlib>

#include <signal.h>
#include <string.h>
#include <getopt.h>
#include <libgen.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <GL/glut.h>
#include <GL/freeglut_std.h>
#include <GL/freeglut_ext.h>
#include <SDL/SDL.h>
#include <SDL/SDL_video.h>

#include "FileIO/netcdf_io.h"
#include "Core/engine.h"

class Vector {
  public:
    Vector() : x(0.0), y(0.0), z(0.0), w(0.0) { }
    Vector(double _x, double _y, double _z, double _w=0) : x(_x), y(_y), z(_z), w(_w) { }
    Vector(Vector &rhs) { if (this != &rhs) { x=rhs.x; y=rhs.y; z=rhs.z; w=rhs.w; } }
    Vector operator-(Vector &rhs) { Vector v(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w); return v; }
    void   set(double _x, double _y, double _z, double _w=0) { x=_x; y=_y; z=_z; w=_w; }
    double x, y, z, w;
};

class Triangle {
  public:
    Triangle(Vector &_p1, Vector &_p2, Vector &_p3) : p1(_p1), p2(_p2), p3(_p3) { };
    void   set(Vector &_p1, Vector &_p2, Vector &_p3) { p1 = _p1; p2 = _p2; p3 = _p3; }
    Vector p1, p2, p3;
    Vector color1, color2, color3;
};

class FloodView {
  public:
    FloodView() : m_dem(NULL), m_mouseButtonPressed(false), m_rotateAction(false), m_rendered(false),
    m_minElevation(9999999), m_maxElevation(-9999999), m_zscale(1.0), m_screenCols(640), m_screenRows(480), 
    m_demFile(""), m_outputFile("") { }

    ~FloodView() {
      if (m_dem)
        delete m_dem;
    }

    // OpenGL visualization
    bool start();
    void stop();
    void createMainWindow();
    void makeLight();
    void initCamera();
    void createMeshes();
    void draw();
    void drawFilledPolygon(Triangle *t);
    void render();
    void getColor(Vector &color, float z, Vector &point);
    void writeText(string text);
    void takeScreenshot();

    // GLUT stuff
    void registerCallbacks(void);

    // General operation
    bool parseOptions(int argc, char **argv);
    void usage(const char *appname, int errorCode);
    bool run();

    // Public vars
    Grid  *m_dem;
    Vector m_rotation;
    Vector m_translation;
    Vector m_mouseCoordinate;
    GLuint m_display_list_id;
    bool   m_mouseButtonPressed;
    bool   m_rotateAction;
    bool   m_rendered;

  private:
    float     m_minElevation;
    float     m_maxElevation;
    float     m_zscale;
    int       m_screenCols;
    int       m_screenRows;
    string    m_demFile;
    string    m_outputFile;
    pthread_t m_glut_tid;
};

// OpenGL visualization

static FloodView *thiz;
void *gfxdebug_GlutThread(void *data);

bool FloodView::start()
{
  int ret;
  int depth = m_dem->cols < 25 ? 10 : 5;

  float tx = -((m_dem->cols / 2.0 + m_dem->rows / 2.0) / 2.0) * m_dem->cellsize;
  float ty = -tx;
  float tz = -pow((m_dem->cols / 2.0 + m_dem->rows / 2.0), 2) * depth;
  m_translation.set(tx, ty, tz);

  float rx = m_dem->cols < 25 ? -40 : -30.0;
  float ry = 0.0;
  float rz = 0.0;
  m_rotation.set(rx, ry, rz);

  ret = pthread_create(&m_glut_tid, NULL, gfxdebug_GlutThread, NULL);
  return ret == 0 ? true : false;
}

void FloodView::stop()
{
  pthread_join(m_glut_tid, NULL);
}

void *gfxdebug_GlutThread(void *data)
{
  thiz->createMainWindow();
  thiz->makeLight();
  thiz->initCamera();
  thiz->registerCallbacks();
  thiz->createMeshes();
  glutPostRedisplay();
  glutMainLoop();
  pthread_exit(NULL);
}

void FloodView::createMainWindow(void)
{
  const char *argv[] = { "demviewer", NULL };
  int argc = 1;

  glutInit(&argc, (char **) argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(m_screenCols, m_screenRows);
  glutCreateWindow("DEM Viewer");
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

  /* Set default background color */
  glClearColor(1, 1, 1, 0);
  glClearDepth(1.0);

  /* Hidden element removal settings */
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glDepthMask(GL_TRUE);
  glDepthRange(0.0, 1.0);
  glDisable(GL_CULL_FACE);
}

void FloodView::makeLight(void)
{
  /* Illumination */
  GLfloat light_ambient[] = { 0.5, 0.5, 0.5, 1.0 };
  GLfloat light_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  GLfloat light_specular[] = { 0.4, 0.4, 0.4, 1.0 };
  GLfloat light_position[] = { 1.0, 1.0, 0.0, 0.0 };

  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHT0);

  GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_shininess[] = { 100.0 }; // accepted values between 0 and 128

  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  glShadeModel(GL_SMOOTH);
}

void FloodView::initCamera(void)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, 1.0, 0.1, thiz->m_dem->cellsize * thiz->m_dem->cellsize * 1000.0);
  glTranslated(0, 0, -5);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void FloodView::createMeshes()
{
  /* Create the Terrain mesh */
  m_display_list_id = glGenLists(1);

  glNewList(m_display_list_id, GL_COMPILE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);
  glBegin(GL_TRIANGLES);

  int step = 1;
  for (int i=0; i<m_dem->cols-step; i += step) {
    for (int j=0; j<m_dem->rows-step; j += step) {
      float alt0 = GRID(m_dem, i,      j);
      float alt1 = GRID(m_dem, i+step, j);
      float alt2 = GRID(m_dem, i+step, j+step);
      float alt3 = GRID(m_dem, i,      j+step);

      if (alt0==m_dem->nodata || alt1==m_dem->nodata || alt2==m_dem->nodata || alt3==m_dem->nodata)
        continue;

      Vector p1(i,      j,      alt0);
      Vector p2(i+step, j,      alt1);
      Vector p3(i+step, j+step, alt2);
      Vector p4(i,      j+step, alt3);

      Triangle *t1 = new Triangle(p1, p2, p3);
      Triangle *t2 = new Triangle(p1, p3, p4);
      this->getColor(t1->color1, alt0, p1);
      this->getColor(t1->color2, alt1, p2);
      this->getColor(t1->color3, alt2, p3);
      this->getColor(t2->color1, alt0, p1);
      this->getColor(t2->color2, alt2, p3);
      this->getColor(t2->color3, alt3, p4);
      this->drawFilledPolygon((Triangle*) t1);
      this->drawFilledPolygon((Triangle*) t2);
      delete t2;
      delete t1;
    }
  }

  glEnd();
  glDisable(GL_NORMALIZE);
  glDisable(GL_LIGHTING);
  glEndList();
}

void FloodView::draw()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();

  glTranslated(m_translation.x, m_translation.y, m_translation.z);
  glRotated(m_rotation.x, 1.0, 0.0, 0.0);
  glRotated(m_rotation.y, 0.0, 0.0, 1.0);
  glTranslatef(-m_dem->cols/2.0, -m_dem->rows/2.0, 0);
  glScaled(m_dem->cellsize, -m_dem->cellsize, m_zscale);

  this->render();
  glPopMatrix();
  if (! m_outputFile.size())
	  this->writeText("");
  glFlush();
  glutSwapBuffers();

  this->takeScreenshot();
}

void FloodView::drawFilledPolygon(Triangle *t)
{
  Vector u, v, normal;

  u = t->p2 - t->p1;
  v = t->p3 - t->p1;
  normal.set(
      (u.y * v.z) - (u.z * v.y),
      (u.z * v.x) - (u.x * v.z),
      (u.x * v.y) - (u.y * v.x)
      );

  glNormal3f(normal.x, normal.x, normal.z);

  glColor4f(t->color1.x, t->color1.y, t->color1.z, t->color1.w);
  glVertex3f(t->p1.x, t->p1.y, t->p1.z);

  glColor4f(t->color2.x, t->color2.y, t->color2.z, t->color2.w);
  glVertex3f(t->p2.x, t->p2.y, t->p2.z);

  glColor4f(t->color3.x, t->color3.y, t->color3.z, t->color3.w);
  glVertex3f(t->p3.x, t->p3.y, t->p3.z);
}

void FloodView::render(void)
{
  glCallList(m_display_list_id);
}

void FloodView::getColor(Vector &color, float z, Vector &point)
{
  static double step = (m_maxElevation - m_minElevation) / 14;
  z -= m_minElevation;
  
  /*
   * Terrain colors.
   * Low altitudes: dark green
   * Low->Normal altitudes: green / yellow
   * Normal->High altitudes: orange / brown
   * High altitudes: bright brown / white
   */
  if (z < step)         color.set(0.000, 0.651, 0.000);
	else if (z < step*2)  color.set(0.141, 0.702, 0.000);
	else if (z < step*3)  color.set(0.298, 0.749, 0.000);
	else if (z < step*4)  color.set(0.478, 0.800, 0.000);
	else if (z < step*5)  color.set(0.679, 0.850, 0.000);
	else if (z < step*6)  color.set(0.902, 0.902, 0.000);
	else if (z < step*7)  color.set(0.910, 0.780, 0.153);
	else if (z < step*8)  color.set(0.917, 0.714, 0.306);
	else if (z < step*9)  color.set(0.623, 0.506, 0.274);
	else if (z < step*10) color.set(0.925, 0.694, 0.462);
	else if (z < step*11) color.set(0.933, 0.725, 0.623);
	else if (z < step*12) color.set(0.941, 0.812, 0.784);
	else if (z < step*13) color.set(0.949, 0.949, 0.949);
  else                  color.set(0.960, 0.960, 0.960);
}

void FloodView::writeText(string text)
{
  string final = text;

  glColor3f(0.0, 0.0, 0.0);
  glRasterPos3f(-2.0, -1.9, 0.0);

  // First line
  stringstream ss;
  ss << (m_rotateAction ? "" : "[") << "Translation" << (m_rotateAction ? "" : "]");
  ss << ": (" << m_translation.x << ", " << m_translation.y << ", " << m_translation.z << "), ";
  ss << (m_rotateAction ? "[" : "") << "Rotation" << (m_rotateAction ? "]" : "");
  ss << ": (" << m_rotation.x << ", " << m_rotation.y << ")";

  final += ss.str();
  const char *ptr = final.c_str();
  while (*ptr) {
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, *ptr);
    ptr++;
  }
}

void FloodView::takeScreenshot()
{
  if (m_outputFile.size()) {
    int bpp = 3;
    int w = m_screenCols;
    int h = m_screenRows;
    SDL_Surface *image;

    image = SDL_CreateRGBSurface(SDL_SWSURFACE, w, h, bpp*8, 0x000000FF, 0x0000FF00, 0x00FF0000, 0x00000000);

    /* Set pixel alignment. This affects glReadPixels(). */
    glReadBuffer(GL_BACK);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, image->pixels);

    /* glReadPixels reads the given rectangle from bottom-left to top-right, so we must	reverse it */
    unsigned char *pixels = (unsigned char *) image->pixels;
    for(int y = 0; y < h / 2; ++y)	{
      const int swap_y = h - y - 1;
      for(int x = 0; x < w; ++x) {
        const int offset = bpp * (x + y * w);
        const int swap_offset = bpp * (x + swap_y * w);
        /* Swap R, G and B of the 2 pixels */
        std::swap(pixels[offset + 0], pixels[swap_offset + 0]);
        std::swap(pixels[offset + 1], pixels[swap_offset + 1]);
        std::swap(pixels[offset + 2], pixels[swap_offset + 2]);	
      }
    }

    SDL_SaveBMP(image, m_outputFile.c_str());
    SDL_FreeSurface(image);
   
    printf("Leaving loop after screenshot\n"); 
    glutLeaveMainLoop();
  }
}


// GLUT Callbacks

void idleCallback(void);
void displayCallback(void);
void reshapeCallback(int w, int h);
void keyboardCallback(unsigned char key, int x, int y);
void mouseClickCallback(int button, int state, int x, int y);
void mouseMotionCallback(int x, int y);

void FloodView::registerCallbacks(void)
{
  /* Register callback routines */
  glutIdleFunc(idleCallback);
  glutDisplayFunc(displayCallback);
  glutReshapeFunc(reshapeCallback);
  glutKeyboardFunc(keyboardCallback);
  glutMouseFunc(mouseClickCallback);
  glutMotionFunc(mouseMotionCallback);
}

void idleCallback(void)
{
  thiz->draw();
}

void displayCallback(void)
{
  thiz->draw();
}

void reshapeCallback(int w, int h)
{
#if 0
  thiz->m_glut_window_size.x = w;
  thiz->m_glut_window_size.y = h;
  glViewport(0, 0, (w > h) ? w : h, (w > h) ? w : h);
#endif
}

void keyboardCallback(unsigned char key, int x, int y)
{
  switch (key) {
    case ' ':
      thiz->m_rotateAction = !thiz->m_rotateAction;
      break;
    case 'q':
      printf("Leaving loop after keypress\n"); 
      glutLeaveMainLoop();
      break;
    default:
      cout << "Invalid key '" << key << "'" << endl;
      cout << "Valid keys are:" << endl;
      cout << "space: toggle rotation/translation mode" << endl;
      cout << "    q: exit the application" << endl;
      cout << endl;
  }
  glutPostRedisplay();
}

void mouseClickCallback(int button, int state, int x, int y) 
{
  thiz->m_mouseButtonPressed = false;

  switch (button) {
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN) {
        thiz->m_mouseButtonPressed = true;
        thiz->m_mouseCoordinate.x = x;
        thiz->m_mouseCoordinate.y = y;
      }
      break;
    case 3: // scroll up
      // Each wheel event reports as a button click (first a GLUT_DOWN then a GLUT_UP)
      if (state == GLUT_DOWN) {
        thiz->m_translation.z += (thiz->m_dem->cellsize * thiz->m_dem->cellsize) / 2;
        glutPostRedisplay();
      }
      break;
    case 4: // scroll down
      // Each wheel event reports as a button click (first a GLUT_DOWN then a GLUT_UP)
      if (state == GLUT_DOWN) {
        thiz->m_translation.z -= (thiz->m_dem->cellsize * thiz->m_dem->cellsize) / 2;
        glutPostRedisplay();
      }
      break;
    default:
      break;
  }
}

void mouseMotionCallback(int x, int y)
{
  int modifiers;

  if (! thiz->m_mouseButtonPressed)
    return;

  modifiers = glutGetModifiers();
  if (! modifiers && ! thiz->m_rotateAction) {
    thiz->m_translation.x += ((double)(x - thiz->m_mouseCoordinate.x) * thiz->m_dem->cellsize);
    thiz->m_translation.y -= ((double)(y - thiz->m_mouseCoordinate.y) * thiz->m_dem->cellsize);

  } else if ((modifiers & (GLUT_ACTIVE_SHIFT | GLUT_ACTIVE_CTRL)) || thiz->m_rotateAction) {
    thiz->m_rotation.x += (double)(y - thiz->m_mouseCoordinate.y);
    thiz->m_rotation.y += (double)(x - thiz->m_mouseCoordinate.x);

 //   if (thiz->m_rotation.x > 360.0)
 //     thiz->m_rotation.x = 360.0;
 //   else if (thiz->m_rotation.x < 270.0)
 //     thiz->m_rotation.x = 270.0;

    if (thiz->m_rotation.y >= 360.0)
      thiz->m_rotation.y -= 360.0;
    else if (thiz->m_rotation.y < 0)
      thiz->m_rotation.y += 360.0;
  }

  thiz->m_mouseCoordinate.x = x;
  thiz->m_mouseCoordinate.y = y;
  glutPostRedisplay();
}

// General operation

void FloodView::usage(const char *appname, int errorCode)
{
  cout << "Usage: " << basename((char *) appname) << " [options]" << endl << endl;
  cout << "Accepted arguments are:" << endl;
  cout << "   -d, --dem=FILE         Input DEM." << endl;
  cout << "   -h, --help             This help." << endl;
  cout << "   -o, --output=FILE      Save output to FILE (BMP format). If unset, run in interactive mode." << endl;
  cout << "   -z, --z-scale=NUM      Z axis scale (default: 1.0)" << endl;
  cout << endl;
  exit(errorCode);
}

bool FloodView::parseOptions(int argc, char **argv)
{
  int c, option_index = 0;
  const char * const short_options = "d:ho:z:";
  const struct option long_options[] = {
    { "dem", 1, 0, 'd' },
    { "help", 0, 0, 'h' },
    { "output", 1, 0, 'o' },
    { "z-scale", 1, 0, 'z' },
    { NULL, 0, 0, 0 },
  };

  while (true) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
      case 'd':
        m_demFile.assign(optarg);
        break;
      case 'h':
        usage(argv[0], 0);
        break;
      case 'o':
        m_outputFile.assign(optarg);
        break;
      case 'z':
        m_zscale = atof(optarg);
        break;
      case '?':
      default:
        break;
    }
  }

  if (! m_demFile.size()) {
    cerr << "Error: --dem-file is not set" << endl;
    return false;
  }

  return true;
}

bool FloodView::run()
{
  NetCDF netcdf(90);

  thiz = this;

  m_dem = netcdf.parse(m_demFile);
  if (! m_dem) {
    cerr << "Failed to parse the DEM" << endl;
    return false;
  }

  Grid *rotated_dem = new Grid(m_dem->rows, m_dem->cols);
  rotated_dem->nodata = m_dem->nodata;
  rotated_dem->filename = m_dem->filename;
  rotated_dem->cellsize = m_dem->cellsize;
  for (int x=0; x<m_dem->cols; ++x)
    for (int y=0; y<m_dem->rows; ++y) {
      GRID(rotated_dem, y, x) = GRID(m_dem, x, y);
      if (GRID(m_dem, x, y) != m_dem->nodata) {
        if (GRID(m_dem, x, y) < m_minElevation)
          m_minElevation = GRID(m_dem, x, y);
        if (GRID(m_dem, x, y) > m_maxElevation)
          m_maxElevation = GRID(m_dem, x, y);
      }
  }
  delete m_dem;
  m_dem = rotated_dem;

  // GFX settings
  this->start();
  this->stop();

  return true;
}

int main(int argc, char **argv)
{
  FloodView view;
  if (view.parseOptions(argc, argv) == false)
    return 1;
  if (view.run() == false)
    return 1;
  return 0;
}
