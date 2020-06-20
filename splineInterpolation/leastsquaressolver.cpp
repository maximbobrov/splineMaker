#include "leastsquaressolver.h"
#include "globals.h"

leastSquaresSolver::leastSquaresSolver()
{

}



void leastSquaresSolver::init_for_benchmark()
{
    node3d n;
    deriv3D d,d0;
    n.x=0.0;
    n.y=0.0;
    n.z=0.0;
    get_derivs_bench(n,d0);
    n.f=d0.d[F];

    m_p.clear();
    m_d0.clear();
    m_d.clear();



    int i=0;
    m_p.push_back(n);
    m_d0.push_back(d0);
    m_d.push_back(d);

    double sc=0.5;
    for ( i=1;i<39;i++)
    {
        n.x=sc*(rand()*1.0/RAND_MAX-0.5);
        n.y=sc*(rand()*1.0/RAND_MAX-0.5);
        n.z=sc*(rand()*1.0/RAND_MAX-0.5);
        get_derivs_bench(n,d0);
        n.f=d0.d[F];
        m_p.push_back(n);
        m_d0.push_back(d0);
        m_d.push_back(d);
    }

    for (int i=0;i<39;i++)
    {
        get_derivs(m_p[i],m_d[i],0.5);
        //printf("i=%d x=%f y=%f z=%f f=%f fx=%f fy=%f fz=%f fnew=%f \n",
         //      i,m_p[i].x,m_p[i].y,m_p[i].z,m_p[i].f,m_d[i].d[FX],m_d[i].d[FY],m_d[i].d[FZ],m_d[i].d[F]);
    }
}

void leastSquaresSolver::LU_decompose(void)
{
    int i,j,k,pivotindex;
    double scales[50];
    double normrow,pivot,size,biggest,mult;

    for (i=0;i<VAR_NUM;i++) //заполнение начальными данными
    {
        ps[i]=i;//маппинг изначального порядка на переставленный.
        normrow=0;//максимум в итой строке

        for (j=0;j<VAR_NUM;j++)
        {
            LU[i][j]=M_[i][j];
            if (normrow<fabs(LU[i][j]))
                normrow=fabs(LU[i][j]);
        }
        if (normrow!=0)
            scales[i]=1.0/normrow;//для общих множителей
        else
        {
            scales[i]=0.0;
            //     err_code(DIV_ZERO);
        }
    }
    //метод гаусса с частичным упорядочиванием

    for (k=0;k<VAR_NUM-1;k++)
    {
        biggest=0;
        for (i=k; i<VAR_NUM;i++)
        {
            size=fabs(LU[ps[i]][k])*scales[ps[i]];
            if (biggest<size)
            {
                biggest=size;
                pivotindex=i;
            }
        }

        if (biggest==0)
        {
            //	err_code(1);
            pivotindex=0;
        }

        if (pivotindex!=k)
        {
            j=ps[k];
            ps[k]=ps[pivotindex];
            ps[pivotindex]=j;
        }

        pivot=LU[ps[k]][k];

        for (i=k+1;i<VAR_NUM;i++)
        {
            mult=LU[ps[i]][k]/pivot;
            LU[ps[i]][k]=mult;

            if (mult!=0.0)
            {
                for (j=k+1; j<VAR_NUM;j++)
                    LU[ps[i]][j]-=mult*LU[ps[k]][j];
            }
        }
    }
    //      if (LU[ps[VAR_NUM-1]][VAR_NUM-1]==0.0) err_code(1);
}

void leastSquaresSolver::m_solve(void)
{
    int i,j;
    double dot;

    for (i=0;i<VAR_NUM;i++)
    {
        dot=0;
        for (j=0;j<i;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=b_m[ps[i]]-dot;
    }

    for (i=VAR_NUM-1; i>=0;i--)
    {
        dot=0.0;

        for (j=i+1;j<VAR_NUM;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=(x_m[i]-dot)/LU[ps[i]][i];
    }
}

void leastSquaresSolver::m_invert(void)
{

    int i,j;
    //err_code(1);
    LU_decompose();

    for (j=0;j<VAR_NUM;j++)
    {

        for (i=0;i<VAR_NUM;i++)
        {
            if (i==j)
                b_m[i]=1;
            else
                b_m[i]=0;
        }

        m_solve();

        for (i=0;i<VAR_NUM;i++)
            Inv[i][j]=x_m[i];
    }
}

void leastSquaresSolver::nullify_m(double m[MAX_EQNS][VAR_NUM])
{
    for (int i=0;i<40;i++)
        for(int j=0;j<40;j++)
            m[i][j]=0;
}

void leastSquaresSolver::nullify_v(double v[VAR_NUM])
{
    for (int i=0;i<10;i++)
        v[i]=0;
}
void leastSquaresSolver::get_derivs(node3d &p, deriv3D &res,double delta)
{
    int var_num=VAR_NUM;
    int eq_num=m_p.size();
   // printf("eq_num=%d \n",eq_num);
    //first index is a row number
    //second index is a column number  M[eq_num][var_num]
    // get Mt*W*M
    nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=m_p[i].x-p.x;
        dyl=m_p[i].y-p.y;
        dzl=m_p[i].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=exp(-r2/(delta*delta));//1.0;

        b_m[i]=0.0;
    }

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }
    m_invert();//now mvm inverted is in Inv; and its[var_num][var_num];

    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=m_p[i].f;
    }

    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    nullify_v(res.d); //now to the solution

    nullify_m(MWM);

    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int k=0;k<var_num;k++)
            {
                MWM[i][j]+=M_[i][k]*Inv[k][j];

            }
        //    printf("%f ",MWM[i][j]);

        }
      //  printf("\n ");
    }

    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            res.d[i]+=Inv[i][j]*mwb[j];
        }
    }
}


void leastSquaresSolver::get_derivs_fast(node3d &p, deriv3D &res,double delta)
{
    int var_num=VAR_NUM;
    int eq_num=m_p.size();
   // printf("eq_num=%d \n",eq_num);
    //first index is a row number
    //second index is a column number  M[eq_num][var_num]
    // get Mt*W*M
    nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=m_p[i].x-p.x;
        dyl=m_p[i].y-p.y;
        dzl=m_p[i].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=exp(-r2/(delta*delta));//1.0;

        b_m[i]=0.0;
    }

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }


    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=m_p[i].f;
    }

    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
            b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose();
    m_solve();

    for (int i=0;i<var_num;i++)
    {

            res.d[i]=x_m[i];

    }
}

void leastSquaresSolver::get_poisson_internal(node3d &p, deriv3D &res, double delta)
{
    int var_num=VAR_NUM;
    int eq_num=m_p.size();

    nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=m_p[i].x-p.x;
        dyl=m_p[i].y-p.y;
        dzl=m_p[i].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=exp(-r2/(delta*delta));//or just 1.0;

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;

    b_m[i]=0.0;

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }


    for (int i=0;i<eq_num-1;i++)
    {
        b_m[i]=m_p[i].f;//from prev interation
    }
    b_m[eq_num-1]=p.rhs;//this is the rhs value


    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
            b_m[i]=mwb[i];  //its mvb

    }

    LU_decompose();
    m_solve();

    /*for (int i=0;i<var_num;i++)
    {
            res.d[i]=x_m[i];
    }*/
    p.f=x_m[0];//res.d[F];
}

void leastSquaresSolver::get_poisson_boundary(node3d &p, deriv3D &res, double delta)
{
    //we do not chack if p is a boundary point this should be done before calling this function

    int var_num=VAR_NUM;
    int eq_num=m_p.size();

    nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=m_p[i].x-p.x;
        dyl=m_p[i].y-p.y;
        dzl=m_p[i].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=exp(-r2/(delta*delta));//1.0;  todo: make w anisotrpopic in the direction normal to the wall

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;
    /////////////////////

    i=eq_num; //constraint for gradient at the wall
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=walls[p.is_boundary].nx;             M_0[i][2]=walls[p.is_boundary].ny;       M_0[i][3]=walls[p.is_boundary].nz;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }


    for (int i=0;i<eq_num-2;i++)
    {
        b_m[i]=m_p[i].f;//from prev interation
    }
    b_m[eq_num-2]=p.rhs;//from prev interation
    b_m[eq_num-1]=p.f_bound;//from prev interation



    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
            b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose();
    m_solve();

    /*for (int i=0;i<var_num;i++)
    {

            res.d[i]=x_m[i];

    }*/
    p.f=x_m[0];//res.d[F];
}


void leastSquaresSolver::draw_points(double sc)
{
    glPointSize(4);
    /*glBegin(GL_POINTS);
    for (int i=0;i<m_p.size();i++)
    {
        glColor3f(m_p[i].f,m_p[i].f,-m_p[i].f);
        glVertex3f(m_p[i].x,m_p[i].y,m_p[i].z);
    }
    glEnd();*/

    glPointSize(4);
    glBegin(GL_POINTS);
    for (int i=0;i<m_p.size();i++)
    {
        //glColor3f(10.0*m_d0[i].d[F],10.0*m_d0[i].d[F],-10.0*m_d0[i].d[F]);
        glColor3f(sc*10.0*m_p[i].f,sc*10.0*m_p[i].f,-sc*10.0*m_p[i].f);
        glVertex3f(m_p[i].x,m_p[i].y,m_p[i].z);
    }
    glEnd();

   /* glBegin(GL_LINES);
    for (int i=0;i<m_p.size();i++)
    {
        glColor3f(1,1,1);
        glVertex3f(m_p[i].x,m_p[i].y,m_p[i].z);
        glVertex3f(m_p[i].x + sc*m_d[i].d[FX] ,
                   m_p[i].y + sc*m_d[i].d[FY] ,
                   m_p[i].z + sc*m_d[i].d[FX]);
    }
    glEnd();

    glBegin(GL_LINES);
    for (int i=0;i<m_p.size();i++)
    {
        glColor3f(1,0,0);
        glVertex3f(m_p[i].x,m_p[i].y,m_p[i].z);
        glVertex3f(m_p[i].x + sc*m_d0[i].d[FX] ,
                   m_p[i].y + sc*m_d0[i].d[FY] ,
                   m_p[i].z + sc*m_d0[i].d[FX]);
    }
    glEnd();*/
}

double bench_f(double x,double y,double z)//just a gauss function for benchmark
{
    double r2=x*x+y*y+z*z;
    return exp(-r2);
}
double bench_dfdx(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return -2*x*exp(-r2);
}

double bench_dfdy(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return -2*y*exp(-r2);
}

double bench_dfdz(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return -2*z*exp(-r2);
}

double bench_d2fdx(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return (4*x*x - 2)*exp(-r2);
}

double bench_d2fdy(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return (4*y*y - 2)*exp(-r2);
}

double bench_d2fdz(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return (4*z*z - 2)*exp(-r2);
}

void leastSquaresSolver::get_derivs_bench(node3d &p, deriv3D &res)
{
    double x,y,z;

    x=p.x-0.2;
    y=p.y;
    z=p.z;

    res.d[F]=bench_f(x,y,z);
    res.d[FX]=bench_dfdx(x,y,z);
    res.d[FY]=bench_dfdy(x,y,z);
    res.d[FZ]=bench_dfdz(x,y,z);

    res.d[FXX]=bench_d2fdx(x,y,z);
    res.d[FYY]=bench_d2fdy(x,y,z);
    res.d[FZZ]=bench_d2fdz(x,y,z);
}
