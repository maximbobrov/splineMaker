#include "globals.h"

int clear_w = 1.0;

int mx0,my0;
int rotate = 0;
float rx0 = 0.0;
float ry0 = 0.0;
float rx = rx0;
float ry = ry0;
double mouse_x,mouse_y;
int redr=0;
double ck=0.1;
double scale  = 500.0;
double view_x=14.507903;
double view_y=8.300000;
double view_z=24.2;
double o_x=0.0;
double o_y=0.0;
double o_z=0.0;

int iNum = NX-1;
int jNum = NY-1;
int kNum = NZ-1;
int currTime = 1;
int itn=0;
int kCur= 0;
double sc=1;
double cv=0.001;
double xmin = 1e10;
double xmax = -1e10;
double ymin = 1e10;
double ymax = -1e10;
double zmin = 1e10;
double zmax = -1e10;
bool drawPoint = false;
bool drawLineSeg = false;
bool drawSpline = true;
bool drawVelocity = false;
bool drawAcceleration = false;
bool drawGridVel = false;
bool drawGridAccel = false;
bool drawGridPressure = false;

double minVel = 1e100;
double maxVel = -1e100;

double minAccel = 1e100;
double maxAccel = -1e100;

double minPressure = 1e100;
double maxPressure = -1e100;

bool rangeCalculated = false;

int isNoOptimized[THREADNUM];



