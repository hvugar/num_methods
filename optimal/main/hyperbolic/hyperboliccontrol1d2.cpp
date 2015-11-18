#include "hyperboliccontrol1d2.h"

void HyperbolicControl1D2::main()
{
    DoubleVector v;
    HyperbolicControl1D2 hc(0.0, 1.0);
    v.resize(2*(hc.M+hc.DM+1));
//    for (unsigned int j=0; j<=hc.M; j++)
//    {
//        double t = j*hc.ht;
//        v[j] = t*t;
//        v[(hc.M+1)+j] = t*t + 1.0;
//    }

    hc.initialize();

    for (unsigned int j=0; j<=hc.M+hc.DM; j++)
    {
        double t = j*hc.ht;
        v[j] = 2.0*t;
        v[(hc.M+hc.DM+1)+j] = 2.0*t + 2.0;
    }

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.00001);
    g2.setEpsilon2(0.00001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.1, 0.0001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(false);
    g2.calculate(v);

    DoubleVector v1(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v1[j] = v[j];
    DoubleVector v2(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v2[j] = v[hc.M+hc.DM+1+j];
    Printer::printVector(v1);
    Printer::printVector(v2);
}

HyperbolicControl1D2::HyperbolicControl1D2(double t0, double t1) : RnFunction(), Printer()
{
    this->t0 = t0;
    this->t1 = t1;
    x0 = 0.0;
    x1 = 1.0;
    UT = 4.0;
    a = 1.0;

    ht = 0.001;

    M  = round((t1-t0)/ht);

    N  = 1000;
    DM = 10;
    //ht = (t1-t0)/M;
    //dt = ht*DM;
    hx = (x1-x0)/N;
    lamda = 0.25;
    R = 1.0;
}

HyperbolicControl1D2::~HyperbolicControl1D2()
{
}

double HyperbolicControl1D2::fx(const DoubleVector& v)
{
    DoubleMatrix u;
    calculateU(v, u);

    double sum = 0.0;//t1;

    double integral = 0.0;
    for (unsigned int j=M; j<=M+DM-1; j++)
    {
        for (unsigned int i=0; i<=N-1; i++)
        {
            double f1 = u[j+0][i+0] - UT;//U[M][i+0];
            double f2 = u[j+0][i+1] - UT;//U[M][i+1];
            double f3 = u[j+1][i+0] - UT;//U[M][i+0];
            double f4 = u[j+1][i+1] - UT;//U[M][i+1];

            integral = integral + (f1*f1 + f2*f2 + f3*f3 + f4*f4);
        }
    }
    integral = R*integral * 0.25 * ht *hx;

    return sum + integral;
}

void HyperbolicControl1D2::gradient(const DoubleVector& v, DoubleVector& g, double)
{
    calculateU(v, uT);
    calculareP(uT, g);
    //R *= 2.0;
}

void HyperbolicControl1D2::calculateU(const DoubleVector& v, DoubleMatrix &u)
{
    pv = &v;

    u.clear();
    u.resize(M+DM+1);
    for (unsigned int j=0; j<=M+DM; j++) u[j].resize(N+1);

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

    for (unsigned int j=0; j<=M+DM-1; j++)
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

        //if (j+1==M) Printer::printVector(u[j+1]);
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();

    u1.clear();
    u0.clear();
}

void HyperbolicControl1D2::calculareP(const DoubleMatrix &u, DoubleVector &g)
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

    for (unsigned int j1=0; j1<=M+DM-1; j1++)
    {
        unsigned int j = M+DM - j1;

        if (j==M+DM)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = pfi1(i);
                p1[i] = pfi1(i) - ht*pfi2(i);
            }
            calculateG(p0, g, M+DM);
            calculateG(p1, g, M+DM-1);
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = alpha2*(p1[i-1]-2.0*p1[i]+p1[i+1]) + 2.0*p1[i] - p0[i] + alpha3*(p0[i+1] - 2.0*p0[i] + p0[i-1]);

                if (M<=j-1 && j-1<=M+DM)
                {
                    rd[i-1] -= (ht*ht)*(2*R*(u[j][i]-UT/*U[j][i]*/));
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

void HyperbolicControl1D2::calculateG(const DoubleVector& psi, DoubleVector& g, unsigned int j)
{
    g[j]        = -(a*a)*(psi[1]-psi[0])/ht;
    g[M+DM+1+j] = -(a*a)*(psi[N-1]-psi[N])/ht;
}

void HyperbolicControl1D2::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double a, RnFunction *fn) const
{
    HyperbolicControl1D2 *hc = dynamic_cast<HyperbolicControl1D2*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(v));
//    DoubleVector v1(M+1); for (unsigned j=0; j<=M; j++) v1[j] = v[j];
//    DoubleVector v2(M+1); for (unsigned j=0; j<=M; j++) v2[j] = v[M+DM+1+j];
//    Printer::printVector(v1);
//    Printer::printVector(v2);
}

void HyperbolicControl1D2::initialize()
{
    DoubleVector v;
    v.resize(2*(M+1+DM));
    for (unsigned int j=0; j<=M+DM; j++)
    {
        double t = j*ht;
        v[j] = g1(t);
        v[(M+1+DM)+j] = g2(t);
    }
    calculateU(v, U);
    //puts("-------------------------");
    //Printer::printVector(U[U.size()-1]);
    //Printer::printVector(U[M]);
    //puts("-------------------------");

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
