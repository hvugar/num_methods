#include "heatequation.h"
#include "tomasmethod.h"

HeatEquation::HeatEquation()
{
    a = 1.0;
}

HeatEquation::~HeatEquation()
{
}

void HeatEquation::setLengthInterval(double x0, double x1)
{
    this->x0 = x0;
    this->x1 = x1;
}

void HeatEquation::setTimeInterval(double t0, double t1)
{
    this->t0 = t0;
    this->t1 = t1;
}

void HeatEquation::setPartNumber(unsigned int M, unsigned int N)
{
    this->M = M;
    this->N = N;
    this->ht = (t1-t0)/M;
    this->hx = (x1-x0)/N;
    this->C = (M+1)*(N+1);
}

void HeatEquation::calculate_u(DoubleVector &u)
{
    u.clear();
    u.resize(N+1);

    DoubleVector a1;
    DoubleVector b1;
    DoubleVector c1;
    DoubleVector d1;
    DoubleVector x1;

    a1.resize(N-1);
    b1.resize(N-1);
    c1.resize(N-1);
    d1.resize(N-1);
    x1.resize(N-1);

    double alpha = -(a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[i] = fi(i*hx, i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = alpha;
                b1[i-1] = beta;
                c1[i-1] = alpha;
                d1[i-1] = u[i] + ht * f(i*hx, i, j*ht, j);
                //d1[i-1] = u1[i] + ht * f[j*(N+1)+i];
            }

            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            d1[0]   -= alpha * m1(j*ht, j);
            d1[N-2] -= alpha * m2(j*ht, j);

            TomasAlgorithm(a1, b1, c1, d1, x1);

            u[0] = m1(j*ht, j);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[i] = x1[i-1];
            }
            u[N] = m2(j*ht, j);
        }
    }

    a1.clear();
    b1.clear();
    c1.clear();
    d1.clear();
    x1.clear();
}

