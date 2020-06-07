#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdio.h>
#include <stdlib.h>
#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#include  <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "dirent.h"
#include <thread>
#include <sys/time.h>


#define W_WIDTH 1600
#define W_HEIGHT 1300
#define NX 251
#define NY 126
#define NZ 76

#define x0_ -50
#define x1_ 50
#define y0_ -25
#define y1_ 25
#define z0_ 0.01
#define z1_ 30.01

#define dx ((x1_-x0_) * 1.0 / NX)
#define dy ((y1_-y0_) * 1.0 / NY)
#define dz ((z1_-z0_) * 1.0 / NZ)

#define THREADNUM 16
#define PATH "dns_data/"//"DA_ppp_0_025/"//"E:/projects/splineViewer/0.080ppp.txt"//"DA_ppp_0_160/"//"NoisyTracks/"//"ExactTracks/"//"NoisyTracks/"//"ExactTracks/"//"NoisyTracks/"//"ExactTracks/"//"NoisyTracks/"//"E:/projects/splineViewer/ActiveLongTracks50_0.080ppp_best.txt"
//"ExactTracks/"//"E:/projects/splineViewer/ActiveLongTracks50_0.050ppp_best.txt"
//"DA_ppp_0_025/""DA_ppp_0_160/"//"DA_ppp_0_005/""test/""testNoise10x/""testNoise100x/"

using namespace std;

struct posVelAccel
{
    posVelAccel(double ix, double iy, double iz, double iu, double iv, double iw, double iax, double iay, double iaz)
        : x(ix), y(iy), z(iz), u(iu), v(iv), w(iw), ax(iax), ay(iay), az(iaz) {}
    double x, y, z;
    double u, v, w;
    double ax, ay, az;
    vector<int> neighbours;
    bool isBound = false;
};

extern int clear_w;

extern int mx0,my0;
extern int rotate;
extern float rx0;
extern float ry0;
extern float rx;
extern float ry;
extern double mouse_x,mouse_y;
extern int redr;
extern double ck;
extern double scale;
extern double view_x;
extern double view_y;
extern double view_z;
extern double o_x;
extern double o_y;
extern double o_z;

extern int iNum;
extern int jNum;
extern int kNum;
extern int currTime;
extern int itn;
extern int kCur;
extern double sc;
extern double cv;
extern double xmin;
extern double xmax;
extern double ymin;
extern double ymax;
extern double zmin;
extern double zmax;
extern bool drawPoint;
extern bool drawLineSeg;
extern bool drawSpline;
extern bool drawVelocity;
extern bool drawAcceleration;
extern bool drawGridVel;
extern bool drawGridAccel;
extern bool drawGridPressure;

extern double minVel;
extern double maxVel;

extern double minAccel;
extern double maxAccel;


extern double minPressure;
extern double maxPressure;

extern bool rangeCalculated;

extern int isNoOptimized[THREADNUM];


#endif // GLOBALS_H
