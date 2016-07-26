#include "example2.h"
#include <math.h>

void Example2::main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    Example2 e2;
    C_UNUSED(e2);
}

Example2::Example2()
{
    eqX.e = this;
    eqP.e = this;

    t0 = 0.0;
    t1 = 1.0;
    h = 0.001;
    M = 1000;
    N = 9;
    x0 = 1.0;

    S = 7;
    xs.resize(S-1);
    xs[0] = +0.000;
    xs[1] = +2.000;
    xs[2] = +12.000;
    xs[3] = +16.000;
    xs[4] = +20.000;
    xs[5] = +32.000;

    DoubleVector k0(S);
    k0[0] = +2.0;
    k0[1] = +20.0;
    k0[2] = +10.0;
    k0[3] = -5.0;
    k0[4] = -3.0;
    k0[5] = +1.0;
    k0[6] = -1.0;

    DoubleVector k1(S);
    k1[0] = +3.0;
    k1[1] = +3.0;
    k1[2] = +3.0;
    k1[3] = +3.0;
    k1[4] = +3.0;
    k1[5] = +3.0;
    k1[6] = +3.0;

    DoubleVector x;
    calculateX(x, k0);
    xT = x[M];
//    printf("x(T): %.8f\n", x[M]);

//    DoubleVector ga(S);
//    this->gradient(k0, ga);
//    ga.L2Normalize();
//    printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", ga[0], ga[1], ga[2], ga[3], ga[4], ga[5], ga[6]);

//    DoubleVector gn;
//    NumericalGradient(k0, gn, 0.001);
//    gn.L2Normalize();
//    printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", gn[0], gn[1], gn[2], gn[3], gn[4], gn[5], gn[6]);

    /* Minimization */
    ConjugateGradient cg;
    cg.setGradient(this);
    cg.setFunction(this);
    cg.setEpsilon1(0.001);
    cg.setEpsilon2(0.001);
    cg.setR1MinimizeEpsilon(0.1, 0.001);
    cg.setNormalize(true);
    cg.setPrinter(this);
    cg.calculate(k1);
}

unsigned int Example2::getKS(double x) const
{
    unsigned int K = 0;
    if (x < xs[0])    { return K = 0; }
    if (x >= xs[S-2]) { return K = S-1; }

    for (unsigned int i=0; i<=S-2; i++)
    {
        if (xs[i] <= x && x < xs[i+1]) { K = i+1; break; }
    }
    return K;
}

double Example2::getK(double x) const
{
    return (*pK)[getKS(x)];
}

void Example2::calculateX(DoubleVector &x, const DoubleVector &K)
{
    pK = &K;

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    DoubleVector t(M+1);
    x.clear();
    x.resize(M+1);
    t[0] = t0;
    x[0] = x0;
    for (unsigned int i=0; i<M; i++)
    {
        if (i%100==0) xi = x[i];

        k1 = eqX.fx(t[i],       x[i]);
        k2 = eqX.fx(t[i]+h/2.0, x[i]+(h/2.0)*k1);
        k3 = eqX.fx(t[i]+h/2.0, x[i]+(h/2.0)*k2);
        k4 = eqX.fx(t[i]+h,     x[i]+h*k3);

        x[i+1] = x[i] + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        t[i+1] = t[i] + h;
    }
}

double Ex2OrdDifEquationX::fx(double t, double x) const
{
    double xi = e->xi;
    double K = e->getK(xi);
    return e->A(t)*x + e->B(t)*K*xi + e->C(t);
}

void Example2::calculateP(DoubleVector &p)
{
    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    DoubleVector t(M+1);

    p.clear();
    p.resize(M+1);

    t[M] = t1;
    p[M] = 2.0*fabs((*px)[M]-xT);

    for (unsigned int i=M; i>=1; i--)
    {
        k1 = eqP.fx(t[i],       p[i], i);
        k2 = eqP.fx(t[i]-h/2.0, p[i]-(h/2.0)*k1, i);
        k3 = eqP.fx(t[i]-h/2.0, p[i]-(h/2.0)*k2, i);
        k4 = eqP.fx(t[i]-h,     p[i]-h*k3, i);
        p[i-1] = p[i] - (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        t[i-1] = t[i] - h;
    }
}

double Ex2OrdDifEquationP::fx(double t, double p, unsigned int j) const
{
    double sum = 0.0;
    for (unsigned int i=0; i<=e->N; i++)
    {
        if (i*100 == j)
        {
            double int_sum = 0.0;
            for (unsigned int k=j; k<(j+100); k++)
            {
                double x = e->px->at(j);
                double P = e->pp->at(k);
                double B = e->B(t);
                double K = e->getK(x);

                int_sum += P*B*K;
            }
            int_sum = int_sum*e->h;
            sum += int_sum;
        }
    }
    return -e->A(t)*p - (1/e->h)*sum;
}

double Example2::fx(const DoubleVector &k)
{
    DoubleVector x;
    calculateX(x, k);
    return ((x[M]-xT)*(x[M]-xT));
}

void Example2::gradient(const DoubleVector &K, DoubleVector &g)
{
    DoubleVector x(M+1);
    px = &x;
    calculateX(x, K);
    //IPrinter::printVector(x);
    //IPrinter::printVector(x.data(), x.size(), "", x.size(), 0, 0, "dataX.txt");

    DoubleVector p(M+1);
    pp = &p;
    calculateP(p);
    //IPrinter::printVector(p);
    //IPrinter::printVector(p.data(), p.size(), "", p.size(), 0, 0, "dataP.txt");

    for (unsigned int s=0; s<g.size(); s++)
    {
        g[s] = 0.0;

        double global_sum = 0.0;
        for (unsigned int i=0; i<=N; i++)
        {
            double local_sum = 0.0;
            double xi = px->at(i*100);

            bool inS = false;

            if (s == 0   && xi < xs[0])    inS = true; else
                if (s == S-1 && xi >= xs[S-2]) inS = true; else
                    if (xs[s-1] <= xi && xi < xs[s]) inS = true;

            if (inS == true)
            {
                for (unsigned int k=i*100; k<i*100+100; k++)
                {
                    double t = k*h;
                    double a1 = p[k+0]*B(t);
                    double a2 = p[k+1]*B(t);

                    local_sum += (h/2.0)*(a1+a2)*xi;
                }
                global_sum += local_sum;
            }
            //printf("%4d %10.6f %10.6f %.4d %10.6f %10.6f\n", s, xi, i*100*h, inS, xs[s-1], xs[s]);
        }
        g[s] = global_sum;
        //puts("");
    }
}

void Example2::print(unsigned int i, const DoubleVector &k, const DoubleVector &g, double a, RnFunction *fn) const
{
    C_UNUSED(i); C_UNUSED(k); C_UNUSED(g); C_UNUSED(a); C_UNUSED(fn);
    DoubleVector x1;
    DoubleVector k1 = k;
    const_cast<Example2*>(this)->calculateX(x1, k1);
    printf("J[%2d]=%10.6f X(T):%10.6f ", i, fn->fx(k), x1[M]);
    printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", k[0], k[1], k[2], k[3], k[4], k[5], k[6]);
    IPrinter::printVector(x1.data(), x1.size(), "", x1.size(), 0, 0, "dataX.txt");
}

void Example2::NumericalGradient(const DoubleVector &k, DoubleVector &g, double h)
{
    g.clear();
    g.resize(k.size());
    for (unsigned int s=0; s<k.size(); s++)
    {
        DoubleVector k1 = k;
        DoubleVector k2 = k;
        k1[s] = k1[s] - h;
        k2[s] = k2[s] + h;
        g[s] = (fx(k2) - fx(k1))/(2.0*h);
    }
}
