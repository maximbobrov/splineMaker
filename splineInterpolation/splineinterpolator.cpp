#include "splineinterpolator.h"

SplineInterpolator::SplineInterpolator(std::vector<double> &vec, std::vector<int> &t)
{
    m_size = vec.size();
    m_moreMinNumber = vec.size() >= 4;
    m_vec = vec;
    m_time.resize(t.size());
    for(size_t i = 0 ; i < t.size(); i++)
        m_time[i] = t[i];
    m_c     = new double[m_size + 3];
    m_f     = new double[m_size + 3];
    m_x     = new double[m_size + 3];
    m_A     = new double[m_size + 3];
    m_B     = new double[m_size + 3];
    m_C     = new double[m_size + 3];
    m_alpha = new double[m_size + 3];
    m_betta = new double[m_size + 3];
    m_cGrad = new double[m_size + 3];

    for (size_t i = 1 ; i < m_size + 1; i++)
    {
        m_f[i] = vec[i-1];
        m_x[i] = vec[i-1];
    }
    m_dF = 1.0;
    m_tolerance = 1e-6;
    m_dt = 0.0006;//0.02;//0.0006;//0.04;//0.0006;
    m_optimizationStep = m_tolerance;
    s_sweep();
    delete m_f;
    delete m_x;
    delete m_A;
    delete m_B;
    delete m_C;
    delete m_alpha;
    delete m_betta;
}

SplineInterpolator::~SplineInterpolator()
{
    delete m_c;
    delete m_cGrad;
}

double SplineInterpolator::baseSpline(double x, int derivNum /*= 0*/)
{
    if (derivNum == 0)
    {
        if (x<=-2.0)
            return 0;
        else if ((x>-2.0)&&(x<=-1.0))
            return (1.0/6.0) * (2.0 + x) * (2.0 + x) * (2.0 + x);
        else if ((x>-1.0)&&(x<=0.0))
            return (4.0/6.0 - (x) * (x) - 0.5 * (x) * (x) * (x));
        else if ((x>0.0)&&(x<=1.0))
            return (4.0/6.0 - x * x + 0.5 * x * x * x);
        else if((x>1.0)&&(x<=2.0))
            return (1.0/6.0)*(2.0-x)*(2.0-x)*(2.0-x);
        else if (x>2.0)
            return 0;
        else return 0;
    }
    if (derivNum == 1)
    {
        if (x<=-2.0)
            return 0;
        else if ((x>-2.0)&&(x<=-1.0))
            return (3.0/6.0)*(2.0+x)*(2.0+x);
        else if ((x>-1.0)&&(x<=0))
            return ( -2.0 * x - 1.5 * x * x);
        else if ((x>0)&&(x<=1.0))
            return (- 2.0 * x + 1.5 * x * x);
        else if((x>1.0)&&(x<=2.0))
            return (-3.0/6.0)*(2.0-x)*(2.0-x);
        else if (x>2.0)
            return 0;
        else return 0;
    }
    if (derivNum == 2)
    {
        if (x<=-2.0)
            return 0;
        else if ((x>-2.0)&&(x<=-1.0))
            return (2.0 + x);
        else if ((x>-1.0)&&(x<=0.0))
            return (-2.0 - 3.0 * x);
        else if ((x>0.0)&&(x<=1.0))
            return (-2.0 + 3.0 * x);
        else if((x>1.0)&&(x<=2.0))
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
    //s_sweep_GaussSeidel(10);
}

void SplineInterpolator::s_sweep_GaussSeidel(size_t itn)
{
    for (size_t n = 0; n < itn; n++)
    {
        m_c[0] = m_c[3] + 3.0 * m_c[1] - 3.0 * m_c[2];
        for (size_t i = 1 ; i < m_size + 1; i++)
        {
            m_c[i] = m_c[i] * 0.99 + 0.01 * (m_f[i] - m_c[i - 1] * m_A[i] - m_c[i + 1] * m_B[i]) / m_C[i];
        }
        m_c[m_size + 1] = m_c[m_size - 2] - 3.0 * m_c[m_size - 1] + 3.0 * m_c[m_size];

        size_t i = m_size;
        m_c[i] = m_c[i] * 0.99 + 0.01 * (m_f[i] - m_c[i - 1] * m_A[i] - m_c[i + 1] * m_B[i]) / m_C[i];
    }
}

double SplineInterpolator::getCoord(double n)
{
    return S_(n, 0);
}

double SplineInterpolator::getVel(double n)
{
    return S_(n, 1)/1000/m_dt;//S_(n, 1)/*/1000*//m_dt;//S_(n, 1)/1000/m_dt;//S_(n, 1)/*/1000*//m_dt;
}

double SplineInterpolator::getAccel(double n)
{
    return S_(n, 2)/1000/m_dt/m_dt;//S_(n, 2)/*/1000*//m_dt/m_dt;//S_(n, 2)/1000/m_dt/m_dt;//S_(n, 2)/*/1000*//m_dt/m_dt;
}

double SplineInterpolator::S_(double x, int derivNum /*= 0*/)
{
    int j;
    double ff;
    ff = x + 1.0;
    j = (int) ff;
    double sum=0;
    if (j == m_size)
        j = m_size-1;

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
    double coeff = 0.3;
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
    for(size_t j = 1; j <= m_size; j++)
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
            for(size_t j = 1; j <= m_size ; j++)
            {
                m_c[j] -= 0.001 * m_cGrad[j];
            }
        }
        m_c[0] = m_c[3] + 3.0 * m_c[1] - 3.0 * m_c[2];
        m_c[m_size+1] = m_c[m_size-2] - 3.0 * m_c[m_size-1] + 3.0 * m_c[m_size];
        //m_c[0] = 2.0 * m_c[1] - m_c[2];
        //m_c[m_size + 1] = 2.0 * m_c[m_size] - m_c[m_size - 1];
        m_dF = FBefore - getF();
        return 1;
    }
    return 0;
}
