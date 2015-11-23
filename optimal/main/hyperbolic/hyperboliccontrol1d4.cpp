#include "hyperboliccontrol1d4.h"

void HyperbolicControl1D4::main()
{
    HyperbolicControl1D4 hc;
    hc.calculateSettings();
    //hc.fx(1.0);

    for (double e = 0.1; e < 1.01; e+=0.1)
    {
        hc.e[0] = e;
        for (double t=0.1; t<2.01; t+=0.1) hc.fx(t);
    }

    //for (double t=0.1; t<1.11; t+=0.1) hc.fx(t);
}

HyperbolicControl1D4::HyperbolicControl1D4()
{
    t1 = 1.0;
    //calculateSettings();
}

HyperbolicControl1D4::~HyperbolicControl1D4()
{
}

void HyperbolicControl1D4::calculateSettings()
{
    t0 = 0.0;

    x0 = 0.0;
    x1 = 1.0;

    N = 100;
    hx = (x1-x0)/N;

    ht = 0.01;
    M = (unsigned int) round((t1-t0)/ht);
    D = 10;

    a = 1.0;
    lamda = 0.25;
    R = 1.0;
    U = 0.0;

    L = 1;
    e.resize(L);
}

double HyperbolicControl1D4::fx(const DoubleVector &v)
{
    DoubleMatrix u;
    calculateU(v, u);

    double sum = 0.0;

    double integral = 0.0;
    for (unsigned int j=M; j<=M+D-1; j++)
    {
        for (unsigned int i=0; i<=N-1; i++)
        {
            double f1 = u[j+0][i+0] - U;
            double f2 = u[j+0][i+1] - U;
            double f3 = u[j+1][i+0] - U;
            double f4 = u[j+1][i+1] - U;

            integral = integral + (f1*f1 + f2*f2 + f3*f3 + f4*f4);
        }
    }
    integral = R*integral * 0.25 * ht *hx;

    return sum + integral;
}

void HyperbolicControl1D4::gradient(const DoubleVector &v, DoubleVector &g, double a)
{
    //printf("%u %u\n", v.size(), g.size());
    DoubleMatrix u;
    calculateU(v, u);
    calculareP(u, g);
}

double HyperbolicControl1D4::fx(double t)
{
    t1 = t;
    calculateSettings();
    //e[0] = 0.25;
    //e[1] = 0.66;
    //e[2] = 0.75;

    //puts("------------------------------------------------------------------------------------------------------");
    //printf("T: %.8f\n", t);

    DoubleVector v((M+D+1)*(L+2));
    for (unsigned int j=0; j<=(M+D); j++)
    {
        //double t = j*ht;
        v[0*(M+D+1)+j] = U+1.0;
        v[1*(M+D+1)+j] = U+2.0;
        v[2*(M+D+1)+j] = U+1.5;
        //v[3*(M+D+1)+j] = U;
        //v[4*(M+D+1)+j] = U;
    }

    //printf("%.8f %.8f %.8f %.8f %.8f %.8f %u %u %u\n", t0, t1, x0, x1, ht, hx, M, N, D);

    ConjugateGradient g;
    g.setFunction(this);
    g.setEpsilon1(0.000000001);
    g.setEpsilon2(0.000000001);
    //g.setGradientStep(0.000000001);
    g.setR1MinimizeEpsilon(0.2, 0.000000001);
    g.setPrinter(this);
    g.setProjection(this);
    g.setNormalize(true);
    g.showEndMessage(false);
    g.calculate(v);

    double rf = fx(v);

    //    DoubleVector v1(D+1); for (unsigned j=M; j<=M+D; j++) v1[j] = v[0*(M+D+1)+j];
    //    DoubleVector v2(D+1); for (unsigned j=M; j<=M+D; j++) v2[j] = v[1*(M+D+1)+j];
    //    DoubleVector v3(D+1); for (unsigned j=M; j<=M+D; j++) v3[j] = v[2*(M+D+1)+j];
    //    Printer::printVector(v1, 10, "v1:\t");
    //    Printer::printVector(v2, 10, "v2:\t");
    //    Printer::printVector(v3, 10, "v2:\t");

    DoubleMatrix u;
    calculateU(v, u);

    //FILE* f = fopen("vugar2.txt", "a");
    //fprintf(f, "------------------------------------------------------------\n");
    //fprintf(f, "T: %.8f Integral: %.16f\n", t, rf);
    //for (unsigned int j=M; j<=M+D; j++)
    //{
        //printf("u[%d]:\t", j);
        //Printer::printVector(u[j]);

        //fprintf(f, "u[%d]:\t", j);
        //for (unsigned int i=0; i<=N; i++)
        //{
        //    fprintf(f, "%.8f ", u[j][i]);
        //}
        //fprintf(f, "\n");
    //}
    //fclose(f);
    //Printer::printVector(u[M], 10, "u[M]:\t");
    //Printer::printVector(u[M+D], 10, "u[M+D]:");

    printf("e1: %f T: %.8f Integral: %.16f M: %d\n", e[0], t, rf, M);
    //printf("%.8f %.16f\n", t, rf);

    return rf;
}

void HyperbolicControl1D4::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const
{}

void HyperbolicControl1D4::project(DoubleVector &v, int index)
{}

double HyperbolicControl1D4::fi1(unsigned int i) const
{
    double x = i*hx;
    return x*x;
}

double HyperbolicControl1D4::fi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControl1D4::m1(unsigned int j) const
{
    double v = (*pv)[j];
    //if (M <= j && j <= M+D) v = U;
    return v;
}

double HyperbolicControl1D4::m2(unsigned int j) const
{
    double v = (*pv)[M+D+1 + j];
    //if (M <= j && j <= M+D) v = U;
    return v;
}

double HyperbolicControl1D4::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    //double t = j*ht;
    double sum = 0.0;

    if (fabs(x-e[0]) < (hx+0.000001))
    {
        double v = (*pv)[2*(M+D+1)+j];
        sum += (1.0/hx) * v * ((hx-fabs(x-e[0]))/hx);
        //printf("x: %f e: %f sum: %f %.20f %.20f\n", x, e[0], sum, fabs(x-e[0]), hx);
    }

    //if (fabs(x-e[1]) < (hx+0.000001))
    //{
    //    double v = (*pv)[3*(M+D+1)+j];
    //    sum += (1.0/hx) * v * ((hx-fabs(x-e[1]))/hx);
    //}

    //if (fabs(x-e[2]) < (hx+0.000001))
    //{
    //    double v = (*pv)[4*(M+D+1)+j];
    //    sum += (1.0/hx) * v * ((hx-fabs(x-e[2]))/hx);
    //}

    return sum;
}

void HyperbolicControl1D4::calculateU(const DoubleVector &v, DoubleMatrix &u)
{
    pv = &v;

    u.clear();
    u.resize(M+D+1);
    for (unsigned int j=0; j<=M+D; j++) u[j].resize(N+1);

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

    for (unsigned int j=0; j<=M+D-1; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = fi1(i);
                u1[i] = fi1(i) + ht*fi2(i);
                u[0][i] = u0[i];
                u[1][i] = u1[i];
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

            u[j+1][0] = m1(j+1);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[j+1][i] = rx[i-1];
            }
            u[j+1][N] = m2(j+1);

            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = u1[i];
                u1[i] = u[j+1][i];
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

void HyperbolicControl1D4::calculareP(const DoubleMatrix &u, DoubleVector &g)
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

    for (unsigned int j1=0; j1<=M+D-1; j1++)
    {
        unsigned int j = M+D-j1;

        if (j==M+D)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = pfi1(i);
                p1[i] = pfi1(i) - ht*pfi2(i);
            }
            calculateG(p0, g, M+D);
            calculateG(p1, g, M+D-1);
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = alpha2*(p1[i-1]-2.0*p1[i]+p1[i+1]) + 2.0*p1[i] - p0[i] + alpha3*(p0[i+1] - 2.0*p0[i] + p0[i-1]);

                if (M<=j-1 && j-1<=M+D)
                {
                    rd[i-1] -= (ht*ht)*(2*R*(u[j][i]-U));
                }
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

void HyperbolicControl1D4::calculateG(const DoubleVector& psi, DoubleVector& g, unsigned int j)
{
    g[0*(M+D+1)+j] = -(a*a)*(psi[1]-psi[0])/ht;
    g[1*(M+D+1)+j] = -(a*a)*(psi[N-1]-psi[N])/ht;

    unsigned int i1 = (unsigned int)(round(e[0]/hx));
    g[2*(M+D+1)+j] = -psi[i1];

    //unsigned int i2 = (unsigned int)(round(e[1]/hx));
    //g[3*(M+D+1)+j] = -psi[i2];

    //unsigned int i3 = (unsigned int)(round(e[2]/hx));
    //g[4*(M+D+1)+j] = -psi[i3];
}

