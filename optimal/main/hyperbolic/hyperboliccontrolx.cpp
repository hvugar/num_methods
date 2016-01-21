#include "hyperboliccontrolx.h"

void HyperbolicControlX::main()
{
    HyperbolicControlX hcx;

    //double a = 0.7;
    //double b = 1.1;
    double t = 0.8945911255;

    //    printf("%.8f %.16f\n", a, hcx.fx(a));
    //    printf("%.8f %.16f\n", b, hcx.fx(b));

    //goldenSectionSearch(a, b, t, &hcx, 0.001);
    hcx.fx(t);
    printf("optimal t: %.10f\n", t);
}

HyperbolicControlX::HyperbolicControlX()
{
    U = 0.0;
    t0 = 0.0;
    t1 = 1.0;
    x0 = 0.0;
    x1 = 1.0;
    a = 1.0;
    hx = 0.01;
    lamda = 0.25;
    L = 1;
}

double HyperbolicControlX::fx(double t)
{
    N = 100;
    hx = (x1-x0)/N;

    t1 = t;
    ht = 0.01;
    M = (unsigned int) round((t1-t0)/ht);
    D = 10;
    xi = 0.4;
    Xi = 40;

    DoubleVector v((L+2)*(M+D+1));
    for (unsigned int j=0; j<=(M+D); j++)
    {
        v[0*(M+D+1)+j] = 1.0;
        v[1*(M+D+1)+j] = 1.0;
        v[2*(M+D+1)+j] = 1.0;
    }

    double min_step = 1.0;
    double gold_eps = 0.0001;

    ConjugateGradient cg;
    cg.setFunction(this);
    cg.setGradient(this);
    cg.setEpsilon1(0.00001);
    cg.setEpsilon2(0.00001);
    cg.setEpsilon3(0.00001);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(this);
    cg.setNormalize(true);
    cg.showEndMessage(false);
    cg.calculate(v);

    DoubleVector gr1(v.size());
    gradient1(v, gr1);
    gr1.L2Normalize();

    DoubleVector gr2(v.size());
    gradient(v, gr2);
    gr2.L2Normalize();

    FILE* file = fopen("20160121.txt", "a");
    fprintf(file, "------------------------------------------------------------\n");
    fprintf(file, "t: %f h: %f\n", t, 0.001);
    IPrinter::printVector(v, "v1: ", (M+D+1)/10, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    IPrinter::printVector(v, "v2: ", (M+D+1)/10, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    IPrinter::printVector(v, "v3: ", (M+D+1)/10, 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    fputs("\n", file);
    IPrinter::printVector(gr1, "g11:", (M+D+1)/10, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr1, "g12:", (M+D+1)/10, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr1, "g13:", (M+D+1)/10, 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    fputs("\n", file);
    IPrinter::printVector(gr2, "g21:", (M+D+1)/10, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr2, "g22:", (M+D+1)/10, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr2, "g23:", (M+D+1)/10, 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    fputs("\n", file);

    DoubleMatrix u;
    pv = &v;
    IHyperbolicEquation::calculateU(u, hx, ht, M+D, N);

    for (unsigned int j=M; j<=M+D; j++)
    {
        char buffer[20];
        int n = sprintf(buffer, "u[%d]: ", j);
        buffer[n] = 0;
        IPrinter::printVector(u[j], buffer, u[j].size(), 0, 0, file);
    }
    fclose(file);

    double rf = fx(v);
    printf("%.8f %.16f\n", t, rf);
    return rf;
}

double HyperbolicControlX::fx(const DoubleVector &v)
{
    pv = &v;
    DoubleMatrix u;
    IHyperbolicEquation::calculateU(u, hx, ht, M+D, N);

    double sum = 0.0;
    for (unsigned int j=M; j<=M+D; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==0 || i==N || j==M+D || j==0) alpha = 0.5;
            if (i==0 && j==M+D) alpha = 0.25;
            if (i==N && j==M+D) alpha = 0.25;
            sum += alpha*(u[j][i]-U)*(u[j][i]-U);
        }
    }
    sum = hx*ht*sum;

    return sum;
}

void HyperbolicControlX::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleMatrix u;
    IHyperbolicEquation::calculateU(u, hx, ht, M+D, N);

    pu = &u;
    DoubleMatrix psi;
    IBackwardHyperbolicEquation::calculateU(psi, hx, ht, M+D, N);

    for (unsigned j=0; j<=M+D; j++)
    {
        g[0*(M+D+1)+j] = -(psi[j][1]-psi[j][0])/hx;
        g[1*(M+D+1)+j] = +(psi[j][N]-psi[j][N-1])/hx;
        g[2*(M+D+1)+j] = -psi[j][Xi];
    }

    //    gradient1(v, g);
}

void HyperbolicControlX::gradient1(const DoubleVector &v, DoubleVector &g)
{
    for (unsigned int i=0; i<v.size(); i++)
    {
        double h = 0.001;
        DoubleVector v1 = v;
        DoubleVector v2 = v;
        v1[i] = v[i] - h;
        v2[i] = v[i] + h;
        g[i] = (fx(v2) - fx(v1))/(2.0*h);
    }
}

double HyperbolicControlX::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    const DoubleVector &v = *pv;
    double v3 = v[2*(M+D+1)+j];

    double sum  = 0.0;
    //version 1
    if (i==Xi)
    {
        sum = (1.0/hx) * v3 * ((hx-fabs(x-xi))/hx);
    }

    // version 2
    //    if (fabs(x-xi) < (hx+0.000001))
    //    {
    //        sum = (1.0/hx) * v3 * ((hx-fabs(x-xi))/hx);
    //    }

    //version 3
    //    double sgm = 3.0*hx;
    //    double a = 1.0/(sgm*sqrt(2.0*M_PI));
    //    double b = 2.0*sgm*sgm;
    //    double g = a * exp(-((x-e[0])*(x-e[0]))/b);
    //    sum += v3 * g;


    return sum;
}

//void HyperbolicControlX::calculateP(const DoubleMatrix &u, const DoubleVector &v, DoubleVector &g)
//{
//    DoubleVector p(N+1);
//    DoubleVector p0(N+1);
//    DoubleVector p1(N+1);

//    DoubleVector da;
//    DoubleVector db;
//    DoubleVector dc;
//    DoubleVector rd;
//    DoubleVector rx;

//    da.resize(N-1);
//    db.resize(N-1);
//    dc.resize(N-1);
//    rd.resize(N-1);
//    rx.resize(N-1);

//    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
//    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
//    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
//    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

//    for (unsigned int j1=0; j1<=M+D; j1++)
//    {
//        unsigned int j = M+D-j1;

//        if (j==M+D)
//        {
//            for (unsigned int i=0; i<=N; i++)
//            {
//                p0[i] = pfi1(i);
//            }
//            calculateG(p0, g, M+D);
//        }
//        else if (j==(M+D-1))
//        {
//            for (unsigned int i=0; i<=N; i++)
//            {
//                p1[i] = p0[i] - ht*pfi2(i);
//            }
//            calculateG(p1, g, M+D-1);
//        }
//        else
//        {
//            for (unsigned int i=1; i<=N-1; i++)
//            {
//                da[i-1] = alpha1;
//                db[i-1] = beta1;
//                dc[i-1] = alpha1;
//                rd[i-1] = alpha2*(p1[i+1]-2.0*p1[i]+p1[i-1]) + 2.0*p1[i] + alpha3*(p0[i+1] - 2.0*p0[i] + p0[i-1]) - p0[i];

//                if (M<=j && j<=M+D-2)
//                {
//                    rd[i-1] -= (ht*ht)*(2*(u[j][i]-U));
//                }
//            }

//            da[0]   = 0.0;
//            dc[N-2] = 0.0;
//            rd[0]   -= alpha1 * pm1(j);
//            rd[N-2] -= alpha1 * pm2(j);
//            TomasAlgorithm(da, db, dc, rd, rx);

//            p[0] = pm1(j);
//            for (unsigned int i=1; i<=N-1; i++)
//            {
//                p[i] = rx[i-1];
//            }
//            p[N] = pm2(j);

//            for (unsigned int i=0; i<=N; i++)
//            {
//                p0[i] = p1[i];
//                p1[i] = p[i];
//            }
//            calculateG(p, g, j);
//        }
//    }

//    da.clear();
//    db.clear();
//    dc.clear();
//    rd.clear();
//    rx.clear();

//    p0.clear();
//    p1.clear();
//}

//void HyperbolicControlX::calculateG(const DoubleVector &psi, DoubleVector &g, unsigned int j)
//{
//    g[0*(M+D+1)+j] = -(psi[1]-psi[0])/hx;
//    g[1*(M+D+1)+j] = +(psi[N]-psi[N-1])/hx;
//    g[2*(M+D+1)+j] = -psi[Xi];
//}

double HyperbolicControlX::fi1(unsigned int i) const
{
    double x = i*hx;
    return x*x;
}

double HyperbolicControlX::fi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControlX::m1(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v1 = v[0*(M+D+1)+j];
    return v1;
}

double HyperbolicControlX::m2(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v2 = v[1*(M+D+1)+j];
    return v2;
}

void HyperbolicControlX::print(unsigned int iteration, const DoubleVector &v, const DoubleVector &gradient, double alpha, RnFunction *fn) const
{
    //printf("J[%d]: %.16f\n", iteration, fn->fx(v));
}
