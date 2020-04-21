
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

#define THREADNUM 16
using namespace std;

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
double scale  = 1.0;
double view_x=14.507903;
double view_y=8.300000;
double view_z=24.2;
double o_x=0.0;
double o_y=0.0;
double o_z=0.0;
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
bool drawPoint = true;
bool drawLineSeg = false;
bool drawSpline = false;
bool drawVelocity = false;
bool drawAcceleration = false;

int isNoOptimized[THREADNUM];

void threadOptimize(int threadIdx, int startIdx, int endIdx);
void display(void);
void init();

double get_time(void) {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

struct FileLoader
{
    struct vec3
    {
        vec3();
        vec3(double ix, double iy, double iz);
        vector<double> x;
        vector<double> y;
        vector<double> z;
    };

    FileLoader(string iPath);
    void findFileInDir();
    void loadSplines();
    void loadSplinesFromNewFormat();

    string m_path;
    int m_numbers[10000000];
    vector<vec3> m_data;
    vector<vec3> m_vel;
    vector<vec3> m_accel;
    std::vector<char*> m_fileNames;
};

FileLoader::vec3::vec3()
{
}

FileLoader::vec3::vec3(double ix, double iy, double iz)
{
    x.push_back(ix);
    y.push_back(iy);
    z.push_back(iz);
}

FileLoader::FileLoader(string iPath)
{
    m_path = iPath;
    for (size_t i = 0; i < 1000000; i++)
    {
        m_numbers[i] = -1;
    }
}

int strCompare(const void * a, const void * b)
{
    return strcmp(*(char**)a, *(char**)b);
}

void FileLoader::findFileInDir()
{
    std::string nameFile,nameForOpen;
    std::string nameDir(m_path);
    DIR *mydir = opendir(nameDir.data());
    if(mydir == NULL) {
        perror("opendir");
    }
    struct dirent* entry;
    while ((entry = readdir(mydir)))
    {
        int len = strlen (entry->d_name);
        if (len >= 4) {
            if (strcmp (".txt", entry->d_name + len - 4) == 0)//if (strcmp (".dat", entry->d_name + len - 4) == 0)
            {
                nameFile = entry->d_name;
                nameForOpen = nameDir+ nameFile;
                char *str = new char[nameForOpen.length() + 1];
                strcpy(str, nameForOpen.c_str());
                m_fileNames.push_back(str);
            }
        }
    }
    char* sort_char_array[m_fileNames.size()];
    for (size_t i = 0 ;i < m_fileNames.size();i++)
        sort_char_array[i] = m_fileNames[i];
    qsort (sort_char_array, m_fileNames.size(), sizeof(sizeof(char**)), strCompare);
    for (size_t i =0 ; i < m_fileNames.size();i++)
    {
        m_fileNames[i] = sort_char_array[i];
        printf("FilesName= %s\n",m_fileNames[i]);
    }
    closedir(mydir);
}

void FileLoader::loadSplines()
{
    m_data.clear();
    for (size_t i = 0 ; i < m_fileNames.size(); i++)
    {
        FILE *file_data = fopen(m_fileNames[i], "r");
        printf("name=%s\n", m_fileNames[i]);
        char str[128];
        fgets(str, sizeof(str), file_data);
        fgets(str, sizeof(str), file_data);
        fgets(str, sizeof(str), file_data);
        double xx, yy, zz;
        int num;
        while(!feof (file_data))
        {
            fscanf(file_data,"%lf %lf %lf %d", &xx, &yy, &zz, &num);
            if((num >= 0) && (num < 1000000))
            {
                if(m_numbers[num] == -1)
                {
                    m_numbers[num] = m_data.size();
                    m_data.push_back(vec3(xx, yy, zz));
                }
                else
                {
                    m_data[m_numbers[num]].x.push_back(xx);
                    m_data[m_numbers[num]].y.push_back(yy);
                    m_data[m_numbers[num]].z.push_back(zz);
                    xmin = xx < xmin ? xx : xmin;
                    xmax = xx > xmax ? xx : xmax;
                    ymin = yy < ymin ? yy : ymin;
                    ymax = yy > ymax ? yy : ymax;
                    zmin = zz < zmin ? zz : zmin;
                    zmax = zz > zmax ? zz : zmax;
                }
            }
        }
        fclose(file_data);
    }
}

void FileLoader::loadSplinesFromNewFormat()
{
    m_data.clear();
    for (size_t i = 0 ; i < m_fileNames.size(); i++)
    {
        m_data.push_back(vec3());
        m_vel.push_back(vec3());
        m_accel.push_back(vec3());
        FILE *file_data = fopen(m_fileNames[i], "r");
        printf("name=%s\n", m_fileNames[i]);
        double xx, yy, zz, xxf, yyf, zzf, vx, vy, vz, ax, ay, az;
        m_numbers[i] = 0;
        while(!feof (file_data))
        {
            int out=fscanf(file_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &xx, &yy, &zz,  &xxf, &yyf, &zzf, &vx, &vy, &vz, &ax, &ay, &az);
            if (out>0)
            {
                m_data[i].x.push_back(xx);
                m_data[i].y.push_back(yy);
                m_data[i].z.push_back(zz);
                m_vel[i].x.push_back(vx);
                m_vel[i].y.push_back(vy);
                m_vel[i].z.push_back(vz);
                m_accel[i].x.push_back(ax);
                m_accel[i].y.push_back(ay);
                m_accel[i].z.push_back(az);
                m_numbers[i]++;
                xmin = xx < xmin ? xx : xmin;
                xmax = xx > xmax ? xx : xmax;
                ymin = yy < ymin ? yy : ymin;
                ymax = yy > ymax ? yy : ymax;
                zmin = zz < zmin ? zz : zmin;
                zmax = zz > zmax ? zz : zmax;
            }
        }
        fclose(file_data);
    }
}

struct SplineInterpolator
{
    SplineInterpolator(std::vector<double> &vec);
    ~SplineInterpolator();

    double baseSpline(double x, int derivNum = 0);
    void s_sweep();
    void s_sweep_GaussSeidel(size_t itn); //iterative solver
    double getCoord(double n);
    double getVel(double n);
    double getAccel(double n);
    double S_(double x, int derivNum = 0);
    double getF();
    void optimizeByRand(size_t itn);
    int optimizeByGrad(size_t itn);
    void calcGrad(double fBefore);

    double m_tolerance;
    double m_optimizationStep;
    size_t m_size;
    std::vector<double> m_vec;
    double m_dF;
    bool m_moreMinNumber;

    double* m_c;
    double* m_f;
    double* m_x;
    double* m_A;
    double* m_B;
    double* m_C;
    double* m_alpha;
    double* m_betta;
    double* m_cGrad;
};

SplineInterpolator::SplineInterpolator(std::vector<double> &vec)
{
    m_size = vec.size();
    m_moreMinNumber = vec.size() >= 3;
    m_vec = vec;
    m_c     = new double[m_size + 2];
    m_f     = new double[m_size + 2];
    m_x     = new double[m_size + 2];
    m_A     = new double[m_size + 1];
    m_B     = new double[m_size + 1];
    m_C     = new double[m_size + 1];
    m_alpha = new double[m_size + 1];
    m_betta = new double[m_size + 1];
    m_cGrad = new double[m_size + 1];

    for (size_t i = 1 ; i < m_size + 1; i++)
    {
        m_f[i] = vec[i-1];
        m_x[i] = vec[i-1];
    }
    m_dF = 1.0;
    m_tolerance = 1e-6;
    m_optimizationStep = m_tolerance;
    s_sweep();
}

SplineInterpolator::~SplineInterpolator()
{
    delete m_c;
    delete m_f;
    delete m_x;
    delete m_A;
    delete m_B;
    delete m_C;
    delete m_alpha;
    delete m_betta;
    delete m_cGrad;
}

double SplineInterpolator::baseSpline(double x, int derivNum /*= 0*/)
{
    if (derivNum == 0)
    {
        if (x<=-2.0)
            return 0;
        else if ((x>-2)&&(x<=-1))
            return (1.0/6.0) * (2.0 + x) * (2.0 + x) * (2.0 + x);
        else if ((x>-1)&&(x<=0))
            return (4.0/6.0 - (x) * (x) - 0.5 * (x) * (x) * (x));
        else if ((x>0)&&(x<=1))
            return (4.0/6.0 - x * x + 0.5 * x * x * x);
        else if((x>1)&&(x<=2))
            return (1.0/6.0)*(2.0-x)*(2.0-x)*(2.0-x);
        else if (x>2.0)
            return 0;
        else return 0;
    }
    if (derivNum == 1)
    {
        if (x<=-2.0)
            return 0;
        else if ((x>-2)&&(x<=-1))
            return (3.0/6.0)*(2.0+x)*(2.0+x);
        else if ((x>-1)&&(x<=0))
            return ( -2.0 * x - 1.5 * x * x);
        else if ((x>0)&&(x<=1))
            return (- 2.0 * x + 1.5 * x * x);
        else if((x>1)&&(x<=2))
            return (-3.0/6.0)*(2.0-x)*(2.0-x);
        else if (x>2.0)
            return 0;
        else return 0;
    }
    if (derivNum == 2)
    {
        if (x<=-2.0)
            return 0;
        else if ((x>-2)&&(x<=-1))
            return (2.0 + x);
        else if ((x>-1)&&(x<=0))
            return (-2.0 - 3.0 * x);
        else if ((x>0)&&(x<=1))
            return (-2.0 + 3.0 * x);
        else if((x>1)&&(x<=2))
            return (2.0 - x);
        else if (x>2.0)
            return 0;
        else return 0;
    }
}

void SplineInterpolator::s_sweep()
{
    double m;
    for (size_t i = 1 ; i < m_size+1; i++)
    {
        m_C[i] = 4.0/6.0;
        m_A[i] = 1.0/6.0;
        m_B[i] = 1.0/6.0;
    }
    m_A[1] = 0.0;
    m_B[1] = 0.0;
    m_C[1] = 1.0;
    m_A[m_size] = 0;
    m_B[m_size] = 0;
    m_C[m_size] = m_C[1];
    for (size_t i = 2; i < m_size + 1; i++)
    {
        m = m_A[i] / m_C[i - 1];
        m_C[i] = m_C[i] - m * m_B[i - 1];
        m_x[i] = m_x[i] - m * m_x[i - 1];
    }
    m_c[m_size] = m_x[m_size] / m_C[m_size];
    for (size_t i = m_size - 1; i >= 1; i--)
    {
        m_c[i] = (m_x[i] - m_B[i] * m_c[i + 1]) / m_C[i];
    }
    m_c[0] = 2.0 * m_c[1] - m_c[2];
    m_c[m_size + 1] = 2.0 * m_c[m_size] - m_c[m_size - 1];
    s_sweep_GaussSeidel(100);
}

void SplineInterpolator::s_sweep_GaussSeidel(size_t itn)
{
    for (size_t i = 1 ; i < m_size + 1; i++)
    {
        m_C[i] = 4.0/6.0;
        m_A[i] = 1.0/6.0;
        m_B[i] = 1.0/6.0;
    }
    for (size_t n = 0; n < itn; n++)
    {
        m_c[0] = m_c[3] + 3.0 * m_c[1] - 3.0 * m_c[2];
        for (size_t i = 1 ; i < m_size + 1; i++)
        {
            m_c[i] = m_c[i] * 0.7 + 0.3 * (m_f[i] - m_c[i - 1] * m_A[i] - m_c[i + 1] * m_B[i]) / m_C[i];
        }
        m_c[m_size + 1] = m_c[m_size - 2] - 3.0 * m_c[m_size - 1] + 3.0 * m_c[m_size];

        size_t i = m_size;
        m_c[i] = m_c[i] * 0.7 + 0.3 * (m_f[i] - m_c[i - 1] * m_A[i] - m_c[i + 1] * m_B[i]) / m_C[i];
    }
}

double SplineInterpolator::getCoord(double n)
{
    return S_(n, 0);
}

double SplineInterpolator::getVel(double n)
{
    double dt = 0.6;
    return S_(n, 1)/dt;
}

double SplineInterpolator::getAccel(double n)
{
    double dt = 0.6;
    return S_(n, 2)/dt/dt*1000.0;
}

double SplineInterpolator::S_(double x, int derivNum /*= 0*/)
{
    int j;
    double ff;
    ff = x + 1.0;
    j = (int) ff;
    double sum=0;

    m_c[0] = m_c[3] + 3.0 * m_c[1] - 3.0 * m_c[2];
    m_c[m_size+1] = m_c[m_size-2] - 3.0 * m_c[m_size-1] + 3.0 * m_c[m_size];
    //m_c[0] = 2.0 * m_c[1] - m_c[2];
    //m_c[m_size + 1] = 2.0 * m_c[m_size] - m_c[m_size - 1];
    for (int i = j - 1;i <= j + 2;i++)
    {
        {
            sum += baseSpline(ff - i, derivNum) * m_c[i];
        }
    }
    return sum;
}

double SplineInterpolator::getF()
{
    double coeff = 0.6;
    double lambda = (1 / M_PI / coeff) * (1 / M_PI / coeff) * (1 / M_PI / coeff);
    double FRes = 0.0, tmp;
    for (size_t i = 0; i < m_size; i++)
    {
        tmp = S_(i) - m_vec[i];
        FRes += tmp * tmp;
    }
    for (size_t i = 1; i < m_size; i++)
    {
        tmp = m_c[i - 1] - 3 * m_c[i] + 3 * m_c[i + 1] - m_c[i + 2];
        FRes += lambda * tmp * tmp;
    }
    return FRes;
}

void SplineInterpolator::optimizeByRand(size_t itn)
{
    double FMin = getF();
    for(size_t i = 0; i < itn; i++)
    {
        for(size_t j = 0; j <= m_size + 1; j++)
        {
            double cTmp = m_c[j];
            m_c[j] += 0.0001 * (rand() * 1.0 / RAND_MAX - 0.5);
            if(getF() < FMin)
                FMin = getF();
            else
                m_c[j] = cTmp;
        }
    }
}

void SplineInterpolator::calcGrad(double FBefore)
{
    for(size_t j = 0; j <= m_size + 1; j++)
    {
        double dc = 0.1 * m_optimizationStep;
        double cTmp = m_c[j];
        m_c[j] += dc;
        m_cGrad[j] = (getF() - FBefore) / dc;
        m_c[j] = cTmp;
    }
}

int SplineInterpolator::optimizeByGrad(size_t itn)
{
    if(fabs(m_dF) > m_tolerance)
    {
        double FBefore = getF();
        for(size_t i = 0; i < itn; i++)
        {
            calcGrad(FBefore);
            for(size_t j = 0; j <= m_size + 1; j++)
            {
                m_c[j] -= 0.001 * m_cGrad[j];
            }
        }
        m_dF = FBefore - getF();
        return 1;
    }
    return 0;
}

FileLoader* fileLoader;

vector<SplineInterpolator*> splineInterpolatorsX;
vector<SplineInterpolator*> splineInterpolatorsY;
vector<SplineInterpolator*> splineInterpolatorsZ;



void get_color(double gval, double min, double max)
{
    const int nn=4;
    int i;
    double val;
    val=gval;
    if (val>max) val=max;
    if (val<min) val=min;

    typedef struct {
        double x,y,z;
    } XYZ;

    XYZ col_table[5];

    col_table[0].x = 0.0; col_table[0].y = 0.0; col_table[0].z = 1.0;
    col_table[1].x = 0.0; col_table[1].y = 1.0; col_table[1].z = 1.0;
    col_table[2].x = 0.0; col_table[2].y = 1.0; col_table[2].z = 0.0;
    col_table[3].x = 1.0; col_table[3].y = 1.0; col_table[3].z = 0.0;
    col_table[4].x = 1.0; col_table[4].y = 0.0; col_table[4].z = 0.0;

    double alpha;
    if ((max-min) > 1e-35)
    {
        alpha=(val-min)/(max-min)*nn;
        i=(int)(alpha);
        alpha=alpha-i;
    }
    else
    {
        alpha=0.0;
        i=2;
    }
    glColor3f(col_table[i].x * (1 - alpha) + col_table[i+1].x * alpha, col_table[i].y * (1 - alpha) + col_table[i+1].y * alpha, col_table[i].z * (1 - alpha) + col_table[i+1].z * alpha);
}

void display(void)
{
    double orient_x=0.0;
    double orient_y=0.0;
    double orient_z=5.0;
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glScalef(sc,sc,sc);
    o_x=cos(-ry/20)*cos(-rx/20);
    o_y=cos(-ry/20)*sin(-rx/20);
    o_z=sin(-ry/20);
    orient_x=view_x+o_x;
    orient_y=view_y+o_y;
    orient_z=view_z+o_z;
    gluLookAt(view_x,view_y,view_z,orient_x,orient_y,orient_z,0,0,1);

    double minVel = 1e100;
    double maxVel = -1e100;

    double minAccel = 1e100;
    double maxAccel = -1e100;

    size_t numSplinePoints = 100;
    size_t numVelVec = 20;

    for (size_t s = 0; s < fileLoader->m_data.size() && splineInterpolatorsX[s]->m_moreMinNumber; s++)
    {
        for( size_t i=0; i <= numVelVec; i++ )
        {
            double vel = splineInterpolatorsX[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 1) *  splineInterpolatorsX[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 1)
                    +    splineInterpolatorsY[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 1) *  splineInterpolatorsY[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 1)
                    +    splineInterpolatorsZ[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 1) *  splineInterpolatorsZ[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 1);
            if(minVel > vel)
                minVel = vel;
            if(maxVel < vel)
                maxVel = vel;

            double accel = splineInterpolatorsX[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 2) *  splineInterpolatorsX[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 2)
                    +      splineInterpolatorsY[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 2) *  splineInterpolatorsY[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 2)
                    +      splineInterpolatorsZ[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 2) *  splineInterpolatorsZ[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec, 2);
            if(minAccel > accel)
                minAccel = accel;
            if(maxAccel < accel)
                maxAccel= accel;
        }
    }

    for (size_t s = 0; s < fileLoader->m_data.size() && splineInterpolatorsX[s]->m_moreMinNumber; s++)
    {
        if(drawPoint)
        {
            glPointSize(3);
            glBegin(GL_POINTS);
            for( size_t i = 0; i <= splineInterpolatorsX[s]->m_size - 1; i++ )
            {
                glColor3f(1.0,1.0,1.0);
                glVertex3f(splineInterpolatorsX[s]->m_vec[i], splineInterpolatorsY[s]->m_vec[i], splineInterpolatorsZ[s]->m_vec[i]);
            }
            glEnd();
        }

        if(drawLineSeg)
        {
            glLineWidth(2);
            glBegin(GL_LINE_STRIP);
            for( size_t i = 0; i < splineInterpolatorsX[s]->m_size; i++ )
            {
                glColor3f(1.0,1.0,1.0);
                glVertex3f(splineInterpolatorsX[s]->m_vec[i], splineInterpolatorsY[s]->m_vec[i], splineInterpolatorsZ[s]->m_vec[i]);
            }
            glEnd();
        }
        if(drawSpline)
        {
            /*glLineWidth(3);
            glBegin(GL_LINE_STRIP);
            for( size_t i=0; i<=numSplinePoints; i++ )
            {
                double vel = splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints)
                        +    splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints)
                        +    splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints);
                get_color(vel, minVel, maxVel);
                glVertex3f(splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints), splineInterpolatorsY[s]->getCoord(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numSplinePoints), splineInterpolatorsZ[s]->getCoord(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numSplinePoints));
            }
            glEnd();*/
            glPointSize(5);
            glBegin(GL_POINTS);
            glColor3f(0.0,1,0);
            //glVertex3f(splineInterpolatorsX[s]->getCoord(0), splineInterpolatorsY[s]->getCoord(0), splineInterpolatorsZ[s]->getCoord(0));
            glVertex3f(splineInterpolatorsX[s]->getCoord(splineInterpolatorsX[s]->m_size - 1), splineInterpolatorsY[s]->getCoord(splineInterpolatorsX[s]->m_size - 1), splineInterpolatorsZ[s]->getCoord(splineInterpolatorsX[s]->m_size - 1));
            glEnd();
            glLineWidth(3);
            glBegin(GL_LINE_STRIP);
            for( size_t i = 0; i < splineInterpolatorsX[s]->m_size; i++ )
            {
                double vecX = splineInterpolatorsX[s]->getVel(i) - fileLoader->m_vel[s].x[i];
                double vecY = splineInterpolatorsY[s]->getVel(i) - fileLoader->m_vel[s].y[i];
                double vecZ = splineInterpolatorsZ[s]->getVel(i) - fileLoader->m_vel[s].z[i];
                glColor3f(vecX * vecX + vecY * vecY + vecZ * vecZ,0.01,0.01);
                glVertex3f(splineInterpolatorsX[s]->getCoord(i), splineInterpolatorsY[s]->getCoord(i), splineInterpolatorsZ[s]->getCoord(i));
            }
            glEnd();
        }
        if(drawVelocity)
        {
            glLineWidth(1);
            glBegin(GL_LINES);
            for( size_t i=0; i<=numVelVec; i++)
            {
                double vel = splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                        +    splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                        +    splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec);
                get_color(vel, minVel, maxVel);
                glVertex3f(splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec), splineInterpolatorsY[s]->S_(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numVelVec), splineInterpolatorsZ[s]->S_(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numVelVec));
                glVertex3f(splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)   + scale * 5.0 * splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                           , splineInterpolatorsY[s]->getCoord(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numVelVec) + scale * 5.0 * splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                           , splineInterpolatorsZ[s]->getCoord(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numVelVec) + scale * 5.0 * splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec));
            }
            glEnd();
        }
        if(drawAcceleration)
        {
            glLineWidth(1);
            glBegin(GL_LINES);
            for( size_t i=0; i<=numVelVec; i++)
            {
                double accel = splineInterpolatorsX[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsX[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                        +      splineInterpolatorsY[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsY[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                        +      splineInterpolatorsZ[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsZ[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec);
                get_color(accel, minAccel, maxAccel);
                glVertex3f( splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec), splineInterpolatorsY[s]->getCoord(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numVelVec), splineInterpolatorsZ[s]->getCoord(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numVelVec));
                glVertex3f( splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)  + scale * 30.0 * splineInterpolatorsX[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                            ,splineInterpolatorsY[s]->getCoord(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numVelVec) + scale * 30.0 * splineInterpolatorsY[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                            ,splineInterpolatorsZ[s]->getCoord(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numVelVec) + scale * 30.0 * splineInterpolatorsZ[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec));
            }
            glEnd();
        }
    }

    glutSwapBuffers();
    if (redr==1) glutPostRedisplay();
}

void m_m(int x,int y)
{
    if (rotate==1)
    {
        rx=rx0+0.1*(x-mx0);
        ry=ry0+0.1*(y-my0);
    }
    glutPostRedisplay();
}

void m_d(int button, int state,int x, int y)
{
    if (state==GLUT_UP)
    {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN)
    {
        rotate=1;
        mx0=x;
        my0=y;
    }
    mouse_x=(1.0*x)/W_WIDTH;
    mouse_y=(W_HEIGHT-(1.0*y))/W_HEIGHT;
    glutPostRedisplay();
}

void threadOptimize(int threadIdx, int startIdx, int endIdx)
{
    for (int s = startIdx; s < endIdx; s++)
    {
        isNoOptimized[threadIdx] += splineInterpolatorsX[s]->optimizeByGrad(1);
        isNoOptimized[threadIdx] += splineInterpolatorsY[s]->optimizeByGrad(1);
        isNoOptimized[threadIdx] += splineInterpolatorsZ[s]->optimizeByGrad(1);
    }
}

void calcEfficient()
{
    double averPosError = 0.0;
    double averVelError = 0.0;
    double averAccelError = 0.0;
    double meanVel = 0.0;
    int number = 0;
    int numToSkip = 1;
    for (size_t s = 0; s < fileLoader->m_data.size() /*&& splineInterpolatorsX[s]->m_moreMinNumber*/; s++)
    {
        for( size_t i = numToSkip; i < splineInterpolatorsX[s]->m_size - numToSkip; i++ )
        {
            double posX = splineInterpolatorsX[s]->getCoord(i) - fileLoader->m_data[s].x[i];
            double posY = splineInterpolatorsY[s]->getCoord(i) - fileLoader->m_data[s].y[i];
            double posZ = splineInterpolatorsZ[s]->getCoord(i) - fileLoader->m_data[s].z[i];
            averPosError += posX * posX + posY * posY + posZ * posZ;

            double vecX = splineInterpolatorsX[s]->getVel(i) - fileLoader->m_vel[s].x[i];
            double vecY = splineInterpolatorsY[s]->getVel(i) - fileLoader->m_vel[s].y[i];
            double vecZ = splineInterpolatorsZ[s]->getVel(i) - fileLoader->m_vel[s].z[i];
            averVelError += vecX * vecX + vecY * vecY + vecZ * vecZ;
            meanVel+=sqrt(fileLoader->m_vel[s].x[i]*fileLoader->m_vel[s].x[i] + fileLoader->m_vel[s].y[i]*fileLoader->m_vel[s].y[i] + fileLoader->m_vel[s].z[i]*fileLoader->m_vel[s].z[i]);


            double accelX = splineInterpolatorsX[s]->getAccel(i) - fileLoader->m_accel[s].x[i];
            double accelY = splineInterpolatorsY[s]->getAccel(i) - fileLoader->m_accel[s].y[i];
            double accelZ = splineInterpolatorsZ[s]->getAccel(i) - fileLoader->m_accel[s].z[i];
            averAccelError += accelX * accelX + accelY * accelY + accelZ * accelZ;

            number++;
        }
    }
    printf("DiffVel = %f DiffAccel = %f DiffPos = %f \n", sqrt(averVelError/number),sqrt(averAccelError/number),sqrt(averPosError/number));
}

void saveInFile()
{
    for (size_t s = 0; s < fileLoader->m_data.size() /*&& splineInterpolatorsX[s]->m_moreMinNumber*/; s++)
    {
        FILE *file_data=fopen(fileLoader->m_fileNames[s],"w");
        for( size_t i = 0; i < splineInterpolatorsX[s]->m_size; i++ )
        {
                fprintf(file_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", splineInterpolatorsX[s]->getCoord(i), splineInterpolatorsY[s]->getCoord(i), splineInterpolatorsZ[s]->getCoord(i)
                        , splineInterpolatorsX[s]->getCoord(i), splineInterpolatorsY[s]->getCoord(i), splineInterpolatorsZ[s]->getCoord(i)
                        , splineInterpolatorsX[s]->getVel(i), splineInterpolatorsY[s]->getVel(i), splineInterpolatorsZ[s]->getVel(i)
                        , splineInterpolatorsX[s]->getAccel(i), splineInterpolatorsY[s]->getAccel(i), splineInterpolatorsZ[s]->getAccel(i));
        }
            fclose(file_data);
    }
}

void kb(unsigned char key, int x, int y)
{
    if (key=='1')
    {
        drawPoint = !drawPoint;
    }
    if (key=='2')
    {
        drawLineSeg = !drawLineSeg;
    }
    if (key=='3')
    {
        drawSpline = !drawSpline;
    }
    if (key=='4')
    {
        drawVelocity = !drawVelocity;
    }
    if (key=='5')
    {
        drawAcceleration = !drawAcceleration;
    }
    if (key=='.')
    {
        ck*=1.1;
    }
    if (key==',')
    {
        ck/=1.1;
    }
    if (key=='[')
    {
        scale/=1.1;
    }
    if (key==']')
    {
        scale*=1.1;
    }
    if (key=='w')
    {
        view_x+=(o_x)*1.1;
        view_y+=(o_y)*1.1;
        view_z+=(o_z)*1.1;
    }
    if (key=='s')
    {
        view_x-=(o_x)*1.1;
        view_y-=(o_y)*1.1;
        view_z-=(o_z)*1.1;
    }
    if (key=='q')
    {
        view_z+=1.1;
    }
    if (key=='e')
    {
        view_z-=1.1;
    }
    if (key=='a')
    {
        double l2=sqrt(o_y*o_y+o_x*o_x);
        view_y+=(o_x)*1.1/l2;
        view_x+=-(o_y)*1.1/l2;
    }
    if (key=='d')
    {
        double    l2=sqrt(o_y*o_y+o_x*o_x);
        view_x+=(o_y)*1.1/l2;
        view_y+=-(o_x)*1.1/l2;
    }

    if(key == 'c')
    {
        for (size_t s = 0; s < fileLoader->m_data.size();s++)
        {
            splineInterpolatorsX[s]->optimizeByGrad(1);
            splineInterpolatorsY[s]->optimizeByGrad(1);
            splineInterpolatorsZ[s]->optimizeByGrad(1);
        }
        printf("One minimization step is done\n");
    }
    if(key == 'z')
    {
        calcEfficient();
    }
    if(key == 'x')
    {
        saveInFile();
    }
    if(key == 'v')
    {
        int numNoMin = 1;
        int numToStop = 0;
        double time1=get_time();
        while (numNoMin != 0 && numToStop < 10000)
        {
            numToStop++;
            numNoMin = 0;
            std::vector <std::thread> th_vec;
            th_vec.clear();
            for (int i = 0; i < THREADNUM; ++i)
            {
                int numForOneProc = THREADNUM == 1 ? fileLoader->m_data.size() : (int)(fileLoader->m_data.size() / (THREADNUM - 1));
                int startIdx = i * numForOneProc;
                int endIdx =  (i + 1) * numForOneProc;
                if(i == (THREADNUM - 1))
                    endIdx = fileLoader->m_data.size();
                isNoOptimized[i] = 0;
                th_vec.push_back(std::thread(threadOptimize,i, startIdx, endIdx));
            }
            for (int i = 0; i < THREADNUM; ++i)
            {
                th_vec.at(i).join();
                numNoMin += isNoOptimized[i];
            }
            printf("Number no minimazed splines = %d\n", numNoMin);
        }
        double time2=get_time();
        printf("Optimization time = %e \n", time2-time1);
    }
    if (key==' ')
    {
        redr=!redr;
    }

    glutPostRedisplay();
}

void init()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective(45.0f, W_WIDTH*1.0/W_HEIGHT, 0.1f, 250.0f);
    glMatrixMode (GL_MODELVIEW);
    fileLoader = new FileLoader("test/");//"DA_ppp_0_005/"
    fileLoader->findFileInDir();
    fileLoader->loadSplinesFromNewFormat();//fileLoader->loadSplines();
    printf("Points loaded\n");
    for (size_t i = 0; i < fileLoader->m_data.size(); i++)
    {
        splineInterpolatorsX.push_back(new SplineInterpolator(fileLoader->m_data[i].x));
        splineInterpolatorsY.push_back(new SplineInterpolator(fileLoader->m_data[i].y));
        splineInterpolatorsZ.push_back(new SplineInterpolator(fileLoader->m_data[i].z));
    }
    printf("Splines calculated\n");

    view_x=xmin-0.5*(xmax-xmin);//xmin-0.5*(xmax-xmin);
    view_y=0.5*(ymax-ymin);//ymin-0.5*(ymax-ymin);
    view_z=0.5*(zmax-zmin);//zmin-0.5*(zmax-zmin);
}

int main(int argc, char** argv)
{
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(W_HEIGHT*(1-0)/(1-0),W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();

    delete fileLoader;
    for (size_t i = 0; i < splineInterpolatorsX.size(); i++) {
        delete splineInterpolatorsX[i];
        delete splineInterpolatorsY[i];
        delete splineInterpolatorsZ[i];
    }
}
