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
    k0[0] = +0.0;
    k0[1] = +20.0;
    k0[2] = +10.0;
    k0[3] = -5.0;
    k0[4] = -3.0;
    k0[5] = +1.0;
    k0[6] = -1.0;

    pK = &k0;
    xT = 25.03126661;

    DoubleVector ga(S);

    this->gradient(k0, ga);
    printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", ga[0], ga[1], ga[2], ga[3], ga[4], ga[5], ga[6], ga[7]);

    DoubleVector gn(S);
    for (unsigned int s=0; s<k0.size(); s++)
    {
        double h1 = 0.0001;
        DoubleVector k1 = k0;
        DoubleVector k2 = k0;
        k1[s] = k1[s] - h1;
        k2[s] = k2[s] + h1;
        gn[s] = (fx(k2) - fx(k1))/(2.0*h1);
    }
    gn.L2Normalize();
    printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", gn[0], gn[1], gn[2], gn[3], gn[4], gn[5], gn[6], gn[7]);
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

void Example2::calculateX(DoubleVector &x)
{
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
    return 3.0*t*x + K*xi + (2.0*t - 3.0*t*t*t);
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
    p[M] = -2.0*(*px)[M];

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
            for (unsigned int k=j; k<(j+99); k++)
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
    return -(3.0*t)*p - (1/e->h)*sum;
}

double Example2::fx(const DoubleVector &k)
{
    pK = &k;
    DoubleVector x;
    calculateX(x);
    return x[M];
}

void Example2::gradient(const DoubleVector &K, DoubleVector &g)
{
    DoubleVector x(M+1);
    calculateX(x);
    px = &x;
    IPrinter::printVector(x);
    IPrinter::printVector(x.data(), x.size(), "", x.size(), 0, 0, "dataX.txt");

    DoubleVector p(M+1);
    pp = &p;
    calculateP(p);
    IPrinter::printVector(p);
    IPrinter::printVector(p.data(), p.size(), "", p.size(), 0, 0, "dataP.txt");

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
                for (unsigned int k=i*100; k<i*100+99; k++)
                {
                    double t = k*h;
                    double a1 = p[k+0]*B(t)*x[k+0];
                    double a2 = p[k+1]*B(t)*x[k+1];

                    local_sum += (h/2.0)*(a1+a2)*xi;
                }
                global_sum += local_sum;
            }
            printf("%4d %10.6f %10.6f %.4d %10.6f %10.6f\n", s, xi, i*100*h, inS, xs[s-1], xs[s]);
        }
        g[s] = global_sum;
        puts("");
    }
    g.L2Normalize();
}

void Example2::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double a, RnFunction *fn) const
{
    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(a); C_UNUSED(fn);
}
