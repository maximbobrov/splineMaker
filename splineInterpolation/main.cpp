
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

#define W_WIDTH 1600
#define W_HEIGHT 1300
using namespace std;

int clear_w = 1.0;
float rx = 10.0;
float ry = 10.0;
int mx0,my0;
int rotate = 0;
float rx0 = 0.0;
float ry0 = 0.0;
double mouse_x,mouse_y;
int redr=0;
double ck=0.1;
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
bool drawPoint = false;
bool drawSpline = true;
bool drawLineSeg = true;

void display(void);
void init();

struct FileLoader
{
    struct vec3
    {
        vec3(double ix, double iy, double iz);
        vector<double> x;
        vector<double> y;
        vector<double> z;
    };

    FileLoader(string iPath);
    void findFileInDir();
    void loadSplines();

    string m_path;
    int m_numbers[1000000];
    vector<vec3> m_data;
    std::vector<char*> m_fileNames;
};

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
            if (strcmp (".dat", entry->d_name + len - 4) == 0)
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
            if((num >= 0) && (num < 100000))
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

struct SplineInterpolator
{
    SplineInterpolator(std::vector<double> &vec);
    ~SplineInterpolator();

    double s_(double x);
    void s_sweep();
    double S_(double x);
    double getF();
    void optimizeByRand(size_t itn);
    int optimizeByGrad(size_t itn);
    void calcGrad(double fBefore);

    double m_tolerance;
    size_t m_size;
    std::vector<double> m_vec;
    double m_dF;

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
    m_tolerance = 1e-3;
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

double SplineInterpolator::s_(double x)
{
    if (x<=-2.0)
        return 0;
    else if ((x>-2)&&(x<=-1))
        return (1.0/6.0)*(2.0+x)*(2.0+x)*(2.0+x);
    else if ((x>-1)&&(x<=0))
        return (4.0/6.0 - (-x) * (-x) + 0.5 * (-x) * (-x) * (-x));
    else if ((x>0)&&(x<=1))
        return (4.0/6.0 - x * x + 0.5 * x * x * x);
    else if((x>1)&&(x<=2))
        return (1.0/6.0)*(2.0-x)*(2.0-x)*(2.0-x);
    else if (x>2.0)
        return 0;
    else return 0;
}

void SplineInterpolator::s_sweep()
{
    double m;
    for (size_t i = 1 ; i < m_size+1; i++)
    {
        m_C[i] = 4.0/6;
        m_A[i] = 1.0/6;
        m_B[i] = 1.0/6;
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
}

double SplineInterpolator::S_(double x)
{
    int j;
    double ff;
    ff = x + 1.0;
    j = (int) ff;
    double sum=0;

    m_c[0] = 2.0 * m_c[1] - m_c[2];
    m_c[m_size + 1] = 2.0 * m_c[m_size] - m_c[m_size - 1];
    for (int i = j - 1;i <= j + 2;i++)
    {
        {
            sum += s_(ff - i) * m_c[i];
        }
    }
    return sum;
}

double SplineInterpolator::getF()
{
    double lambda = (1 / M_PI / 0.2) * (1 / M_PI / 0.2) * (1 / M_PI / 0.2);
    double FRes = 0.0;
    for (size_t i = 1; i < m_size-1; i++)
    {
        FRes += (S_(i) - m_vec[i]) * (S_(i) - m_vec[i]);
        FRes += lambda * (m_c[i - 1] - 3 * m_c[i] + 3 * m_c[i + 1] - m_c[i + 2]) * (m_c[i - 1] - 3 * m_c[i] + 3 * m_c[i + 1] - m_c[i + 2]);
    }
    return FRes;
}

void SplineInterpolator::optimizeByRand(size_t itn)
{
    double FMin = getF();
    for(int i = 0; i < itn; i++)
    {
        for(int j = 0; j <= m_size + 1; j++)
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
        double dc = 0.001;
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

    for (size_t s = 0; s < fileLoader->m_data.size(); s ++)
    {
        if(drawPoint)
        {
            glLineWidth(1);
            glBegin(GL_POINTS);
            for( size_t i = 0; i < splineInterpolatorsX[s]->m_size; i++ )
            {
                glColor3f(1.0,1.0,1.0);
                glVertex3f(splineInterpolatorsX[s]->m_vec[i], splineInterpolatorsY[s]->m_vec[i], splineInterpolatorsZ[s]->m_vec[i]);
            }
            glEnd();
        }

        if(drawLineSeg)
        {
            glLineWidth(1);
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
            glLineWidth(3);
            glBegin(GL_LINE_STRIP);
            for( size_t i=0; i<=100; i++ )
            {
                glColor3f(0.5, 0.0, 0.0);
                glVertex3f(splineInterpolatorsX[s]->S_(i*(splineInterpolatorsX[s]->m_size-1)*1.0/100), splineInterpolatorsY[s]->S_(i*(splineInterpolatorsY[s]->m_size-1)*1.0/100), splineInterpolatorsZ[s]->S_(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/100));
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

void kb(unsigned char key, int x, int y)
{
    if (key=='1')
    {
        drawPoint = !drawPoint;
    }
    if (key=='2')
    {
        drawSpline = !drawSpline;
    }
    if (key=='3')
    {
        drawLineSeg = !drawLineSeg;
    }
    if (key=='.')
    {
        ck*=1.1;
    }
    if (key==',')
    {
        ck/=1.1;
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
    if(key == 'v')
    {
        int numNoMin = 1;
        while (numNoMin != 0)
        {
            numNoMin= 0;
            for (size_t s = 0; s < fileLoader->m_data.size();s++)
            {
                numNoMin += splineInterpolatorsX[s]->optimizeByGrad(1);
                numNoMin += splineInterpolatorsY[s]->optimizeByGrad(1);
                numNoMin += splineInterpolatorsZ[s]->optimizeByGrad(1);
            }
            printf("Number no minimazed splines = %d\n", numNoMin);
        }
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
    fileLoader = new FileLoader("DA_ppp_0_005/");
    fileLoader->findFileInDir();
    fileLoader->loadSplines();
    printf("Points loaded\n");
    for (size_t i = 0; i < fileLoader->m_data.size(); i++)
    {
        splineInterpolatorsX.push_back(new SplineInterpolator(fileLoader->m_data[i].x));
        splineInterpolatorsY.push_back(new SplineInterpolator(fileLoader->m_data[i].y));
        splineInterpolatorsZ.push_back(new SplineInterpolator(fileLoader->m_data[i].z));
    }
    printf("Splines calculated\n");

    view_x=xmin-0.5*(xmax-xmin);
    view_y=ymin-0.5*(ymax-ymin);
    view_z=zmin-0.5*(zmax-zmin);
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
