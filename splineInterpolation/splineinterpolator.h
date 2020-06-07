#ifndef SPLINEINTERPOLATOR_H
#define SPLINEINTERPOLATOR_H
#include "globals.h"

struct SplineInterpolator
{
    SplineInterpolator(std::vector<double> &vec, std::vector<int> &t);
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
    std::vector<int> m_time;
    double m_dF;
    double m_dt;
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

#endif // SPLINEINTERPOLATOR_H
