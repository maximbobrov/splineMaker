#include "tools.h"

double get_time(void) {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}


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

double ddx_a(double var_[NX][NY][NZ],int i, int j,int k,double rhs)
{
    double res;
    if (i==0)
    {
        res= (-rhs*2.0*dx+ 4.0*var_[i+1][j][k] - var_[i+2][j][k])/3.0;

    }else
    {
        if (i==NX-1)
        {
            res =(rhs*2.0*dx + 4.0*var_[i-1][j][k] - var_[i-2][j][k])/3.0;
        }else
        {
            res=0;
        }
    }
    return res;
}

double ddy_a(double var_[NX][NY][NZ],int i, int j,int k,double rhs)
{
    double res;

    if (j==0)
    {
        res= (-rhs*2.0*dy+ 4.0*var_[i][j+1][k] - var_[i][j+2][k])/3.0;
    }else
    {
        if (j==NY-1)
            res =(rhs*2.0*dy + 4.0*var_[i][j-1][k] - var_[i][j-2][k])/3.0;
        else
            res=0;

    }
    return res;
}

double ddz_a(double var_[NX][NY][NZ],int i, int j,int k,double rhs)
{
    double res;

    if (k==0)
    {
        res= (-rhs*2.0*dz+ 4.0*var_[i][j][k+1] - var_[i][j][k+2])/3.0;
    }else
    {
        if (k==NZ-1)
            res =(rhs*2.0*dz + 4.0*var_[i][j][k-1] - var_[i][j][k-2])/3.0;
        else
            res=0;
    }
    return res;
}

double d2di_a(double var_[NX][NY][NZ],int i, int j,int k)
{
    double res;
    if (i==0)
    {
        res=1.0;
    }else
    {
        if (i==NX-1)
            res=1.0;
        else
            res=-2.0;
    }
    return res;
}

double d2dj_a(double var_[NX][NY][NZ],int i, int j,int k)
{
    double res;
    if (j==0)
    {
        res=1.0;
    }else
    {
        if (j==NY-1)
            res=1.0;
        else
            res=-2.0;
    }
    return res;
}

double d2dk_a(double var_[NX][NY][NZ],int i, int j,int k)
{
    double res;

    if (k==0)
    {
        res=1.0;
    }else
    {
        if (k==NZ-1)
            res=1.0;
        else
            res=-2.0;
    }
    return res;
}

double d2di_res(double var_[NX][NY][NZ],int i, int j,int k)
{
    double res;
    if (i==0)
    {
        res=( - 2.0*var_[i+1][j][k] + var_[i+2][j][k]);
    }else
    {
        if (i==NX-1)
            res=( - 2.0*var_[i-1][j][k] + var_[i-2][j][k]);
        else
            res=(var_[i+1][j][k] + var_[i-1][j][k]);
    }
    return res;
}

double d2dj_res(double var_[NX][NY][NZ],int i, int j,int k)
{
    double res;
    if (j==0)
    {
        res=( - 2.0*var_[i][j+1][k] + var_[i][j+2][k]);
    }else
    {
        if (j==NY-1)
            res=( - 2.0*var_[i][j-1][k] + var_[i][j-2][k]);
        else
            res=(var_[i][j+1][k]  + var_[i][j-1][k]);
    }
    return res;
}

double d2dk_res(double var_[NX][NY][NZ],int i, int j,int k)
{
    double res;
    if (k==0)
    {
        res=( - 2.0*var_[i][j][k+1] + var_[i][j][k+2]);
    }else
    {
        if (k==NZ-1)
            res=( - 2.0*var_[i][j][k-1] + var_[i][j][k-2]);
        else
            res=(var_[i][j][k+1]  + var_[i][j][k-1]);
    }
    return res;
}
