#include "hyperboliccontrol1d3.h"

void HyperbolicControl1D3::main()
{
    HyperbolicControl1D3 hc;
    //hc.doSettings();
    //hc.initialize();

    DoubleVector v;
    v.resize(2*(hc.M+hc.DM+1)+1);

    for (unsigned int j=0; j<=hc.M+hc.DM; j++)
    {
        double t = j*hc.ht;
        v[j] = t*t+1.5;
        v[(hc.M+hc.DM+1)+j] = t*t+2.5;
    }
    v[v.size()-1]=0.8;

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.000000001);
    g2.setEpsilon2(0.000000001);
    //g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.1, 0.0000001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(v);

//    DoubleVector v1(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v1[j] = v[j];
//    DoubleVector v2(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v2[j] = v[hc.M+hc.DM+1+j];
//    Printer::printVector(v1);
//    Printer::printVector(v2);

    DoubleMatrix u;
    DoubleVector gr(v.size());
    hc.calculateU(v, u);
    hc.calculareP(u, gr);

    puts("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    IPrinter::printVector(u[hc.M]);
    DoubleVector v1(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v1[j] = v[j];
    DoubleVector v2(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v2[j] = v[hc.M+hc.DM+1+j];
    IPrinter::printVector(v1);
    IPrinter::printVector(v2);
    DoubleVector gr1(hc.M+1); for (unsigned j=0; j<=hc.M; j++) gr1[j] = gr[j];
    DoubleVector gr2(hc.M+1); for (unsigned j=0; j<=hc.M; j++) gr2[j] = gr[hc.M+hc.DM+1+j];
    IPrinter::printVector(gr1);
    IPrinter::printVector(gr2);

    double integral = 0.0;
    for (unsigned int i=0; i<=hc.N-1; i++)
    {
        unsigned int j = i+1;
        double f1 = (u[hc.M+hc.DM][i] + u[hc.M][i] - 2*hc.U)*(u[hc.M+hc.DM][i] - u[hc.M][i]);
        double f2 = (u[hc.M+hc.DM][j] + u[hc.M][j] - 2*hc.U)*(u[hc.M+hc.DM][i] - u[hc.M][j]);
        integral += f1+f2;
    }
    integral = 0.5 * hc.hx * integral;
    gr[gr.size()-1] = 1.0 + integral;

    printf("%.10f\n", gr[gr.size()-1]);
}

HyperbolicControl1D3::HyperbolicControl1D3() : RnFunction(), IPrinter()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    U = 4.0;
    R = 1.0;
    doSettings(t1);
}

void HyperbolicControl1D3::doSettings(double t)
{
    t0 = 0.0; t1 = t;
    x0 = 0.0; x1 = 1.0;
    a = 1.0;
    M  = 100;
    N  = 100;
    DM = 10;

    ht = (t1-t0)/M;
    hx = (x1-x0)/N;

    dt = ht*DM;

    lamda = 0.25;


    //printf("%f %f %f %f %f %f %d %d %d\n", t0, t1, x0, x1, ht, hx, M, N, DM);
}

double HyperbolicControl1D3::fx(const DoubleVector& v)
{
    DoubleMatrix u;
    doSettings(v[v.size()-1]);
    calculateU(v, u);

    double sum = v[v.size()-1];

    double integral = 0.0;
    for (unsigned int j=M; j<=M+DM-1; j++)
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
    integral = R*integral * 0.25 * ht * hx;

    return sum + integral;
}

void HyperbolicControl1D3::gradient(const DoubleVector& v, DoubleVector& g, double gs)
{
    puts("--------------------------------------------------------------------");
    if (R < 10000.0) R *= 2.0;
    double t = v[v.size()-1];
    printf("T: %.10f R: %f\n", t, R);
    doSettings(t);

    DoubleMatrix u;
    calculateU(v, u);
    calculateG2(v, g);

    IPrinter::printVector(u[M]);

    DoubleVector v1(M+1); for (unsigned j=0; j<=M; j++) v1[j] = v[j];
    DoubleVector v2(M+1); for (unsigned j=0; j<=M; j++) v2[j] = v[M+DM+1+j];
    IPrinter::printVector(v1);
    IPrinter::printVector(v2);

    DoubleVector g1(M+1); for (unsigned j=0; j<=M; j++) g1[j] = g[j];
    DoubleVector g2(M+1); for (unsigned j=0; j<=M; j++) g2[j] = g[M+DM+1+j];
    IPrinter::printVector(g1);
    IPrinter::printVector(g2);

    printf("TG: %.10f\n", g[g.size()-1]);


//    DoubleMatrix u;
//    calculateU(v, u);
//    calculareP(u, g);

//    double integral = 0.0;
//    for (unsigned int i=0; i<=N-1; i++)
//    {
//        unsigned int j = i+1;
//        double f1 = (u[M+DM][i] + u[M][i] - 2*U)*(u[M+DM][i] - u[M][i]);
//        double f2 = (u[M+DM][j] + u[M][j] - 2*U)*(u[M+DM][i] - u[M][j]);
//        integral += f1+f2;
//    }
//    integral = 0.5 * hx * integral;

//    g[g.size()-1] = 1.0 + integral;
}

void HyperbolicControl1D3::calculateU(const DoubleVector& v, DoubleMatrix &u)
{
    pv = &v;
    double t = v[v.size()-1];
    doSettings(t);

    // Clear
    for (unsigned int j=0; j<u.size(); j++) u[j].size();
    u.clear();
    // Resize
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
                //Printer::printVector(u[0]);
                //Printer::printVector(u[1]);
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

            //Printer::printVector(u[j+1]);
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

void HyperbolicControl1D3::calculareP(const DoubleMatrix &u, DoubleVector &g)
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

void HyperbolicControl1D3::calculateG(const DoubleVector& psi, DoubleVector& g, unsigned int j)
{
    g[j]        = -(a*a)*(psi[1]-psi[0])/ht;
    g[M+DM+1+j] = -(a*a)*(psi[N-1]-psi[N])/ht;
}

void HyperbolicControl1D3::calculateG2(const DoubleVector &v, DoubleVector &g)
{
    double f0 = fx(v);

    DoubleVector vc = v;
    double h = 0.000001;
    for (unsigned int i=0; i<v.size(); i++)
    {
        vc = v;
        vc[i] = vc[i] + h;
        double f1 = fx(vc);
        g[i] = (f1-f0)/h;
    }
}

void HyperbolicControl1D3::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double a, RnFunction *fn) const
{
    HyperbolicControl1D3 *hc = dynamic_cast<HyperbolicControl1D3*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(v));
}

void HyperbolicControl1D3::initialize()
{
    DoubleVector v;
    v.resize(2*(M+DM+1));
    for (unsigned int j=0; j<=M+DM; j++)
    {
        double t = j*ht;
        v[j] = g1(t);
        v[(M+DM+1)+j] = g2(t);
    }
    DoubleMatrix u;
    calculateU(v, u);
    puts("-------------------------");
    //Printer::printVector(u[u.size()-1]);
    IPrinter::printVector(u[M]);
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
