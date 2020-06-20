#ifndef LEASTSQUARESSOLVER_H
#define LEASTSQUARESSOLVER_H
#include <vector>


//number of variables
#define VAR_NUM 10

//maximum number of points
#define MAX_EQNS 60


typedef struct
{double x ,y,z,f;
int is_boundary; //if<0 than it's not the boundary oterwise it's boundary plane number
double f_bound;//value of the gradient at the boundary
double rhs;// value of the rhs for a poisson equation at the boundary
} node3d;


typedef struct
{double d[VAR_NUM];} deriv3D;

typedef struct
{double nx,ny,nz;
 double d;//normal distance from  coordinate origin  i.e.  nx*x+ny*y+nz*z=d
} boundary_plane;





class leastSquaresSolver
{
private:
    double M_[MAX_EQNS][VAR_NUM],
           M_0[MAX_EQNS][VAR_NUM],
           MWM[MAX_EQNS][VAR_NUM],
           x_m[MAX_EQNS],
           b_m[MAX_EQNS],
           w[MAX_EQNS],
           mwb[MAX_EQNS],
           LU[MAX_EQNS][VAR_NUM],
           Inv[MAX_EQNS][VAR_NUM];
    int ps[MAX_EQNS];

public:
    enum deriv_order
    {
        F,
        FX,
        FY,
        FZ,
        FXX,
        FXY,
        FXZ,
        FYY,
        FYZ,
        FZZ
    };

    leastSquaresSolver();

    std::vector<node3d> m_p; //points cloud
    std::vector<deriv3D> m_d,m_d0; //derivs cloud
    std::vector<boundary_plane> walls; //all boundaries are here for now

    void init_with_points();
    void init_for_benchmark();
    void LU_decompose(void);
    void m_solve(void);
    void m_invert(void);

    void nullify_m(double m[MAX_EQNS][VAR_NUM]);
    void nullify_v(double v[VAR_NUM]);
    void get_derivs(node3d &p, deriv3D &res, double delta);
    void get_derivs_fast(node3d &p, deriv3D &res, double delta);

     void get_poisson_internal(node3d &p, deriv3D &res, double delta);
     void get_poisson_boundary(node3d &p, deriv3D &res, double delta);


    void draw_points(double sc);
    void get_derivs_bench(node3d &p, deriv3D &res);
};


#endif // LEASTSQUARESSOLVER_H
