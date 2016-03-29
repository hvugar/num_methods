#include "hyperboliccontrolx.h"
#include <gradient_sd.h>

#define USE_NUMERICAL_GRADIENT1
//#define USE_DELTA_FUNCTION

void HyperbolicControlX::main()
{
    HyperbolicControlX hcx;

    //    double a = 0.6;
    //    double b = 1.4;
    double t = 0.92;
    //0.9544858420

    //    printf("%.8f %.16f\n", a, hcx.fx(a));
    //    printf("%.8f %.16f\n", b, hcx.fx(b));

    //    goldenSectionSearch(a, b, t, &hcx, 0.001);
    //    R1Minimize::HalphIntervalMethod(a, b, t, &hcx, 0.01);
    //stranghLineSearch(0.8, 0.5, a, b, &hcx);
    //printf("%.10f %.10f\n", a, b);

    //    for (double t=0.80; t<1.10; t+=0.01)
    //    {
    //        hcx.fx(t);
    //    }

    hcx.fx(t);
    printf("optimal t: %.10f\n", t);
}

HyperbolicControlX::HyperbolicControlX()
{
    U = 0.0;
    lamda = 0.25;
    a = 1.0;
#ifdef USE_DELTA_FUNCTION
    L = 1;
#else
    L = 0;
#endif
}

double HyperbolicControlX::fx(double t)
{
    t0 = 0.0;
    t1 = t;
    x0 = 0.0;
    x1 = 1.0;
    N = 1000;
    hx = 0.001;

    t1 = t;
    ht = 0.001;
    M = (unsigned int) round((t1-t0)/ht);
    D = 10;
#ifdef USE_DELTA_FUNCTION
    xi = 0.4;
    Xi = 400;
#endif

    printf("%d %d %f %f\n", M, N, ht, hx);

    DoubleVector v((L+2)*(M+D-1));
    for (unsigned int j=0; j<=(M+D-2); j++)
    {
        v[0*(M+D-1)+j] = 1.0;
        v[1*(M+D-1)+j] = 1.0;
#ifdef USE_DELTA_FUNCTION
        v[2*(M+D-1)+j] = 1.0;
#endif
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
    //cg.calculate(v);

    //    double h = 0.01;
    //    DoubleVector gr1(v.size());
    //    IGradient::Gradient(this, h, v, gr1);
    //    gr1.L2Normalize();

    //    DoubleVector gr2(v.size());
    //    gradient(v, gr2);
    //    gr2.L2Normalize();

    FILE* file = fopen("20160130.txt", "a");
    fprintf(file, "------------------------------------------------------------------------------------------------------------------------\n");
#ifdef USE_DELTA_FUNCTION
    fprintf(file, "T:%f hx:%f ht:%f M:%d N:%d x:%f X:%d J[v]:%.20f\n", t, hx, ht, M, N, xi, Xi, fx(v));
#else
    fprintf(file, "T:%f hx:%f ht:%f M:%d N:%d J[v]:%.20f\n", t, hx, ht, M, N, fx(v));
#endif
    unsigned int div = (M+D-1);
    fprintf(file, "Controls. Count:%d\n", div);
    IPrinter::printVector(v, "v1: ", div, 0*(M+D-1), 0*(M+D-1)+(M+D-2), file);
    IPrinter::printVector(v, "v2: ", div, 1*(M+D-1), 1*(M+D-1)+(M+D-2), file);
#ifdef USE_DELTA_FUNCTION
    IPrinter::printVector(v, "v3: ", div, 2*(M+D-1), 2*(M+D-1)+(M+D-2), file);
#endif
    //    fprintf(file, "Numerical gradients. Count:%d\n", div);
    //    IPrinter::printVector(gr1, "gr1:", (M+D+1), 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    //    IPrinter::printVector(gr1, "gr2:", (M+D+1), 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    //    IPrinter::printVector(gr1, "gr3:", (M+D+1), 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    //    fprintf(file, "Analytic gradient. Count:%d\n", div);
    //    IPrinter::printVector(gr2, "gr1:", (M+D+1), 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    //    IPrinter::printVector(gr2, "gr2:", (M+D+1), 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    //    IPrinter::printVector(gr2, "gr3:", (M+D+1), 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    fputs("Amplitudes:\n", file);
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
    fputs("------------------------------------------------------------------------------------------------------------------------\n", file);
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
#ifndef USE_NUMERICAL_GRADIENT
    pv = &v;
    DoubleMatrix u;
    IHyperbolicEquation::calculateU(u, hx, ht, M+D, N);

    pu = &u;
    DoubleMatrix p;
    IBackwardHyperbolicEquation::calculateU(p, hx, ht, M+D, N);

    for (unsigned j=2; j<=M+D; j++)
    {
        g[0*(M+D-1)+(j-2)] = -(p[j][1]-p[j][0])/hx;
        g[1*(M+D-1)+(j-2)] = +(p[j][N]-p[j][N-1])/hx;
#ifdef USE_DELTA_FUNCTION
        g[2*(M+D-1)+(j-2)] = -p[j][Xi];
#endif
    }
#else
    IGradient::Gradient(this, 0.01, v, g);
#endif
}

double HyperbolicControlX::f(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);

    double sum  = 0.0;
#ifdef USE_DELTA_FUNCTION
    double x = i*hx;
    const DoubleVector &v = *pv;
    double v3 = v[2*(M+D-1)+(j-2)];
    //version 1
    if (i==Xi)
    {
        sum = (1.0/(hx)) * v3;
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
#endif
    return sum;
}

double HyperbolicControlX::bf(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u = *pu;
    if (M<=j)
    {
        return -(2.0*(u[j][i]-U));
    }
    return 0.0;
}

double HyperbolicControlX::initial1(unsigned int i) const
{
    C_UNUSED(i);
    return 2.0;
}

double HyperbolicControlX::initial2(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}

double HyperbolicControlX::boundary(Boundary type, unsigned int j) const
{
    const DoubleVector &v = *pv;
    if (type==Left) return v[0*(M+D-1)+(j-2)];
    if (type==Right) return v[1*(M+D-1)+(j-2)];
    return 0.0;
}

double HyperbolicControlX::binitial1(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}

double HyperbolicControlX::binitial2(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}


double HyperbolicControlX::bboundary(Boundary type, unsigned int j) const
{
    C_UNUSED(j);
    C_UNUSED(type);
    return 0.0;
}

void HyperbolicControlX::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    printf("J[%d]: %.16f\n", i, fn->fx(v));
}
