#include "hyperboliccontrolx.h"

void HyperbolicControlX::main()
{
    HyperbolicControlX hcx;

    double a = 0.7;
    double b = 1.1;
    double t = 0.8945911255;
    //0.9544858420

    //    printf("%.8f %.16f\n", a, hcx.fx(a));
    //    printf("%.8f %.16f\n", b, hcx.fx(b));

//    goldenSectionSearch(a, b, t, &hcx, 0.001);
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
    N = 1000;
    hx = (x1-x0)/N;

    t1 = t;
    ht = 0.001;
    M = (unsigned int) round((t1-t0)/ht);
    D = 10;
    xi = 0.4;
    Xi = 400;

    DoubleVector v((L+2)*(M+D+1));
    for (unsigned int j=0; j<=(M+D); j++)
    {
        v[0*(M+D+1)+j] = 1.0;
        v[1*(M+D+1)+j] = 1.0;
        v[2*(M+D+1)+j] = 1.0;
    }

    double min_step = 1.0;
    double gold_eps = 0.001;

    ConjugateGradient cg;
    cg.setFunction(this);
    cg.setGradient(this);
    cg.setEpsilon1(0.001);
    cg.setEpsilon2(0.001);
    cg.setEpsilon3(0.001);
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

    FILE* file = fopen("20160121_2.txt", "a");
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
    DoubleMatrix p;
    IBackwardHyperbolicEquation::calculateU(p, hx, ht, M+D, N);

    for (unsigned j=0; j<=M+D; j++)
    {
        g[0*(M+D+1)+j] = -(p[j][1]-p[j][0])/hx;
        g[1*(M+D+1)+j] = +(p[j][N]-p[j][N-1])/hx;
        g[2*(M+D+1)+j] = -p[j][Xi];
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

double HyperbolicControlX::bf(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u = *pu;
    if (M<=j)
    {
        return -(2.0*(u[j-1][i]-U));
    }
    return 0.0;
}

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

double HyperbolicControlX::bfi1(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControlX::bfi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControlX::bm1(unsigned int j) const
{
    return 0.0;
}

double HyperbolicControlX::bm2(unsigned int j) const
{
    return 0.0;
}

void HyperbolicControlX::print(unsigned int iteration, const DoubleVector &v, const DoubleVector &gradient, double alpha, RnFunction *fn) const
{
    printf("J[%d]: %.16f\n", iteration, fn->fx(v));
}
