#include "hyperboliccontrolh.h"

void HyperbolicControlH::main()
{
    HyperbolicControlH hcx;

    double a = 0.6;
    double b = 1.1;
    double t = 0.92;
    //0.9544858420

    //    printf("%.8f %.16f\n", a, hcx.fx(a));
    //    printf("%.8f %.16f\n", b, hcx.fx(b));

    //    goldenSectionSearch(a, b, t, &hcx, 0.001);
    //    R1Minimize::HalphIntervalMethod(a, b, t, &hcx, 0.01);
    //    stranghLineSearch(0.8, 0.5, a, b, &hcx);
    //    printf("%.10f %.10f\n", a, b);

    //    for (double t=0.80; t<1.10; t+=0.01)
    //    {
    //        hcx.fx(t);
    //    }

    hcx.fx(t);
    printf("optimal t: %.10f\n", t);
}

HyperbolicControlH::HyperbolicControlH()
{
    U = 0.0;
    lamda = 0.25;
    a = 1.0;
    L = 1;
}

double HyperbolicControlH::fx(double t)
{
    t0 = 0.0;
    t1 = t;
    x0 = 0.0;
    x1 = 1.0;
    N = 100;
    hx = 0.01;

    t1 = t;
    ht = 0.01;
    M = (unsigned int) round((t1-t0)/ht);
    D = 10;
    xi = 0.2;
    Xi = 20;

    printf("%d %d %f %f\n", M, N, ht, hx);

    DoubleVector v((L+2)*(M+D-1));
    for (unsigned int j=0; j<=(M+D-2); j++)
    {
        v[0*(M+D-1)+j] = 1.0;
        v[1*(M+D-1)+j] = 1.0;
        v[2*(M+D-1)+j] = 1.0;
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

//    double h = 0.001;
//    DoubleVector gr1(v.size());
//    IGradient::Gradient(this, h, v, gr1);
//    gr1.L2Normalize();

//    DoubleVector gr2(v.size());
//    gradient(v, gr2);
//    gr2.L2Normalize();

    FILE* file = fopen("20160130_4.txt", "a");
    fprintf(file, "------------------------------------------------------------------------------------------------------------------------\n");
    fprintf(file, "T:%f hx:%f ht:%f M:%d N:%d x:%f X:%d J[v]:%.20f\n", t, hx, ht, M, N, xi, Xi, fx(v));
    unsigned int div = (M+D-1);
    fprintf(file, "Controls. Count:%d\n", div);
    IPrinter::printVector(v, "v1: ", div, 0*(M+D-1), 0*(M+D-1)+(M+D-2), file);
    IPrinter::printVector(v, "v2: ", div, 1*(M+D-1), 1*(M+D-1)+(M+D-2), file);
    IPrinter::printVector(v, "v3: ", div, 2*(M+D-1), 2*(M+D-1)+(M+D-2), file);
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

double HyperbolicControlH::fx(const DoubleVector &v)
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

void HyperbolicControlH::gradient(const DoubleVector &v, DoubleVector &g)
{
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

        double sum = 0.0;
        for (unsigned int i=Xi; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==Xi || i==N) alpha = 0.5;
            sum += alpha*p[j][i];
        }
        g[2*(M+D-1)+(j-2)] = -hx*sum;
    }
    //    IGradient::Gradient(this, 0.01, v, g);
}

double HyperbolicControlH::f(unsigned int i, unsigned int j) const
{
    double sum  = 0.0;
    double x = i*hx;
    const DoubleVector &v = *pv;
    double v3 = v[2*(M+D-1)+(j-2)];
    //version 1
    if (i>=Xi)
    {
        sum = v3;
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

double HyperbolicControlH::bf(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u = *pu;
    if (M<=j)
    {
        return -(2.0*(u[j][i]-U));
    }
    return 0.0;
}

double HyperbolicControlH::fi1(unsigned int i) const
{
    return 2.0;
}

double HyperbolicControlH::fi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControlH::m1(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v1 = v[0*(M+D-1)+(j-2)];
    return v1;
}

double HyperbolicControlH::m2(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v2 = v[1*(M+D-1)+(j-2)];
    return v2;
}

double HyperbolicControlH::bfi1(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControlH::bfi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControlH::bm1(unsigned int j) const
{
    return 0.0;
}

double HyperbolicControlH::bm2(unsigned int j) const
{
    return 0.0;
}

void HyperbolicControlH::print(unsigned int iteration, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    printf("J[%d]: %.16f\n", iteration, fn->fx(v));
}
