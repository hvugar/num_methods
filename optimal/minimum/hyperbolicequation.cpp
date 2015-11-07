#include "hyperbolicequation.h"
#include "tomasmethod.h"
#include "printer.h"

void HyperbolicEquationInterface::calculateU(DoubleVector& u, unsigned int M, unsigned int N, double t0, double t1, double x0, double x1, double a, double lamda) const
{
    double ht = (t1 - t0) / M;
    double hx = (x1 - x0) / N;

    u.clear();
    u.resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

    DoubleVector da;
    DoubleVector db;
    DoubleVector dc;
    DoubleVector rd;
    DoubleVector rx;

    da.resize(N-1);
    db.resize(N-1);
    dc.resize(N-1);
    rd.resize(N-1);
    rx.resize(N-1);

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j=0; j<=M-1; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = fi1(i);
                u1[i] = fi1(i) + ht*fi2(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = alpha2*(u1[i-1]-2.0*u1[i]+u1[i+1]) + 2.0*u1[i] + alpha3*(u0[i+1] - 2.0*u0[i] + u0[i-1]) - u0[i] + (ht*ht)*f(i, j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;
            rd[0]   -= alpha1 * m1(j+1);
            rd[N-2] -= alpha1 * m2(j+1);
            TomasAlgorithm(da, db, dc, rd, rx);

            u[0] = m1(j+1);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[i] = rx[i-1];
            }
            u[N] = m2(j+1);

            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = u1[i];
                u1[i] = u[i];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();

    u1.clear();
    u0.clear();
}

void HyperbolicEquationInterface::calculateU(DoubleMatrix& u, unsigned int M, unsigned int N, double t0, double t1, double x0, double x1, double a, double lamda) const
{}

void HyperbolicEquationInterface::calculateP(DoubleVector &u, unsigned int M, unsigned int N, double t0, double t1, double x0, double x1, double a, double lamda) const
{}

HyperbolicEquation::HyperbolicEquation(unsigned int M, unsigned int N, double t0, double t1, double x0, double x1,  double a)
    : M(M), N(N), t0(t0), t1(t1), x0(x0), x1(x1 ), a(a)
{
    ht = (t1 - t0) / M;
    hx = (x1 - x0) / N;
    lamda = 0.25;
}

HyperbolicEquation::~HyperbolicEquation()
{

}

void HyperbolicEquation::setLengthInterval(double x0, double x1)
{
    this->x0 = x0;
    this->x1 = x1;
    hx = (x1 - x0) / N;
}

void HyperbolicEquation::setTimeInterval(double t0, double t1)
{
    this->t0 = t0;
    this->t1 = t1;
    ht = (t1 - t0) / M;
}

void HyperbolicEquation::setPartNumber(unsigned int M, unsigned int N)
{
    this->M = M;
    this->N = N;
    this->ht = (t1-t0)/M;
    this->hx = (x1-x0)/N;
}

void HyperbolicEquation::calculate(DoubleVector &u)
{
    u.clear();
    u.resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

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

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j=0; j<=M-1; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = fi1(i);
                u1[i] = fi1(i) + ht*fi2(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = alpha1;
                b1[i-1] = beta1;
                c1[i-1] = alpha1;
                d1[i-1] = alpha2*(u1[i-1]-2.0*u1[i]+u1[i+1]) + 2.0*u1[i] + alpha3*(u0[i+1] - 2.0*u0[i] + u0[i-1]) - u0[i] + (ht*ht)*f(i, j);
            }

            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            d1[0]   -= alpha1 * m1(j+1);
            d1[N-2] -= alpha1 * m2(j+1);
            TomasAlgorithm(a1, b1, c1, d1, x1);

            u[0] = m1(j+1);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[i] = x1[i-1];
            }
            u[N] = m2(j+1);

            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = u1[i];
                u1[i] = u[i];
            }
        }
    }

    a1.clear();
    b1.clear();
    c1.clear();
    d1.clear();
    x1.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HyperbolicEquation2D::HyperbolicEquation2D(double t0, double t1, double x10, double x11, double x20, double x21, double a1, double a2, unsigned int M, unsigned int N1, unsigned int N2)
    : t0(t0), t1(t1), x10(x10), x11(x11), x20(x20), x21(x21), a1(a1), a2(a2), M(M), N1(N1), N2(N2)
{
}

HyperbolicEquation2D::~HyperbolicEquation2D()
{}

void HyperbolicEquation2D::calculateImplicitly(DoubleMatrix& m)
{
    DoubleMatrix u0(N2+1, N1+1);
    DoubleMatrix u1(N2+1, N1+1);
    DoubleMatrix u2(N2+1, N1+1);

    for (unsigned int k=0; k<M; k++)
    {
        if (k==0)
        {

        }
        else
        {

        }
    }
}

void HyperbolicEquation2D::calculateExplicitly(DoubleMatrix& m)
{

}


