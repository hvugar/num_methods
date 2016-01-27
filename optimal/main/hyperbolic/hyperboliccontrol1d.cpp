#include "hyperboliccontrol1d.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>

void HyperbolicControl1D::main()
{
    DoubleVector v;
    HyperbolicControl1D hc;

    v.resize(2*(hc.M+1)+1);
//    for (unsigned int j=0; j<=hc.M; j++)
//    {
//        double t = j*hc.ht;
//        v[j] = t*t;
//        v[(hc.M+1)+j] = t*t + 1.0;
//    }

    hc.initialize();

    for (unsigned int j=0; j<=hc.M; j++)
    {
        double t = j*hc.ht;
        v[j] = 2.0*t;
        v[(hc.M+1)+j] = 2.0*t + 2.0;
    }

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0001);
    g2.setEpsilon2(0.0001);
    g2.setR1MinimizeEpsilon(0.1, 0.000001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(false);
    g2.calculate(v);

    DoubleVector v1(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v1[j] = v[j];
    DoubleVector v2(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v2[j] = v[hc.M+1+j];
    IPrinter::printVector(v1);
    IPrinter::printVector(v2);
}

HyperbolicControl1D::HyperbolicControl1D() : RnFunction(), IPrinter()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    a = 1.0;
    M  = 100;
    N  = 100;
    DM = 10;
    ht = (t1-t0)/M;
    hx = (x1-x0)/N;
    dt = DM * ht;
    lamda = 0.25;
}

double HyperbolicControl1D::fx(const DoubleVector& v)
{
    calculateU(v, uT);

    double sum = 0.0;

    double integral = 0.0;
    for (unsigned int i=0; i<=N-1; i++)
    {
        unsigned int j = i+1;
        double f1 = uT[i] - U[i];
        double f2 = uT[j] - U[j];
        integral = integral + f1*f1 + f2*f2;
    }
    integral = integral * 0.5 * hx;

    double norm1 = 0.0;
    for (unsigned int i=0; i<=M-1; i++)
    {
        unsigned int j = i+1;
        double f1 = v[j] - g1(j*ht);
        double f2 = v[i] - g1(j*ht);

        norm1 += f1*f1 + f2*f2;
    }
    norm1 = ht*0.5*norm1;

    double norm2 = 0.0;
    for (unsigned int i=0; i<=M-1; i++)
    {
        unsigned int j = i+1;
        double f1 = v[j+M+1] - g2(j*ht);
        double f2 = v[i+M+1] - g2(i*ht);

        norm2 += f1*f1 + f2*f2;
    }
    norm2 = ht*0.5*norm2;

//    for (unsigned int j=M; j<=M+DM-1; j++)
//    {
//        for (unsigned int i=0; i<=N-1; i++)
//        {
//            double f1 = u[j+0][i+0] - U;
//            double f2 = u[j+0][i+1] - U;
//            double f3 = u[j+1][i+0] - U;
//            double f4 = u[j+1][i+1] - U;

//            integral = integral + (f1*f1 + f2*f2 + f3*f3 + f4*f4);
//        }
//    }
//    integral = integral * 0.25 * ht *hx;

    return sum + integral + norm1 + norm2;
}

void HyperbolicControl1D::gradient(const DoubleVector& v, DoubleVector& g)
{
    calculateU(v, uT);
    calculareP(uT, g);
}

double HyperbolicControl1D::fi1(unsigned int i) const
{
    double x = i*hx;
    return x*x;
}

double HyperbolicControl1D::fi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControl1D::m1(unsigned int j) const
{
    return (*pv)[j];
}

double HyperbolicControl1D::m2(unsigned int j) const
{
    unsigned int i = M+1 + j;
    return (*pv)[i];
}

double HyperbolicControl1D::f(unsigned int i, unsigned int j) const
{
    return 0.0;
}


void HyperbolicControl1D::calculateU(const DoubleVector& v, DoubleVector &u)
{
    pv = &v;

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
                //u[0][i] = u0[i];
                //u[1][i] = u1[i];
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


void HyperbolicControl1D::calculareP(const DoubleVector &u, DoubleVector &g)
{
    DoubleVector p(N+1);
    DoubleVector p0(N+1);
    DoubleVector p1(N+1);

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

    for (unsigned int j1=0; j1<=M-1; j1++)
    {
        unsigned int j = M - j1;

        if (j==M)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = pfi1(i);
                p1[i] = pfi1(i) - ht*pfi2(i);
            }
            calculateG(p0, g, M);
            calculateG(p1, g, M-1);
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = alpha2*(p1[i-1]-2.0*p1[i]+p1[i+1]) + 2.0*p1[i] - p0[i] + alpha3*(p0[i+1] - 2.0*p0[i] + p0[i-1]);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;
            rd[0]   -= alpha1 * pmu1(j-1);
            rd[N-2] -= alpha1 * pmu2(j-1);
            TomasAlgorithm(da, db, dc, rd, rx);

            p[0] = pmu1(j-1);
            for (unsigned int i=1; i<=N-1; i++)
            {
                p[i] = rx[i-1];
            }
            p[N] = pmu2(j-1);

            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = p1[i];
                p1[i] = p[i];
            }
            calculateG(p, g, j-1);
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();

    p0.clear();
    p1.clear();
}

void HyperbolicControl1D::calculateG(const DoubleVector& psi, DoubleVector& g, unsigned int j)
{
    double t = j*ht;
    g[j]     = -(a*a)*(psi[1]-psi[0])/ht   + 2.0*((*pv)[j] - g1(t));
    g[M+1+j] = -(a*a)*(psi[N-1]-psi[N])/ht + 2.0*((*pv)[M+1+j] - g2(t));
}

void HyperbolicControl1D::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double a, RnFunction *fn) const
{
    HyperbolicControl1D *hc = dynamic_cast<HyperbolicControl1D*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(v));
}

void HyperbolicControl1D::initialize()
{
    DoubleVector v;
    v.resize(2*(M+1));
    for (unsigned int j=0; j<=M; j++)
    {
        double t = j*ht;
        v[j] = g1(t);
        v[(M+1)+j] = g2(t);
    }
    calculateU(v, U);
    puts("-------------------------");
    IPrinter::printVector(U);
    puts("-------------------------");

//    DoubleMatrix u;
//    calculateU(u);
//    U = 0.0;

//    for (unsigned int j=M; j<=M+DM; j++)
//    {
//        //printf("%d\n", j);
//        for (unsigned int i=0; i<=N; i++)
//        {
//            U = u[j][i];
//        }
//    }
}
