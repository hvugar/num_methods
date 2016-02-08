#include "pointcontrol11.h"

void PointControl11::main()
{
    DoubleVector q(3);
    q[0] = 0.0;
    q[1] = 0.0;
    q[2] = 0.0;

    PointControl11 pc;
    pc.x1 = 4.61141362;

    /* Minimization */
    ConjugateGradient g;
    g.setFunction(&pc);
    g.setEpsilon1(0.000001);
    g.setEpsilon2(0.000001);
    g.setR1MinimizeEpsilon(0.1, 0.00001);
    g.setPrinter(&pc);
    g.setNormalize(true);
    g.calculate(q);

    //    printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7], q[8], q[0], q[0]+q[1]+q[2]+q[3]+q[4]+q[5]+q[6]+q[7]+q[8]+q[9]);
    printf("%.8f %.8f %.8f %.8f\n", q[0], q[1], q[2], q[0]+q[1]+q[2]);

    DoubleVector x;
    pc.calculateX(q, x);
    printf("x[N]: %.8f\n", x[pc.N]);

    FILE* f =  fopen("data.txt", "w");
    for (unsigned int i=0; i<=pc.N; i++)
    {
        fprintf(f, "%.16f\n", x[i]);
    }
    //Printer::printVector(x, pc.N+1, "x:", f);
    fclose(f);
}

PointControl11::PointControl11()
{
    t0 = 0.0;
    t1 = 1.0;
    x0 = 0.0;

    N = 10000;
    ht = (t1-t0)/N;
}

PointControl11::~PointControl11()
{}

double PointControl11::fx(const DoubleVector &q)
{
    DoubleVector x;
    calculateX(q, x);
    return (x[N]-x1)*(x[N]-x1);
}

void PointControl11::gradient(const DoubleVector &q, DoubleVector &g, double gradient_step)
{
    DoubleVector x;
    DoubleVector p;
    calculateX(q, x);
    calculateP(q, x, p);

    g[0] = p[(unsigned int)(0.2*N)];
    g[1] = p[(unsigned int)(0.5*N)];
    g[2] = p[(unsigned int)(0.8*N)];
}

void PointControl11::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    printf("J[%d]: %.14f\n", i, fn->fx(x));
}

void PointControl11::calculateX(const DoubleVector &q, DoubleVector &x)
{
    x.clear();
    x.resize(N+1);

//    double sgm = 5.0*ht;
//    double a = 1.0/(sgm*sqrt(2.0*M_PI));
//    double b = 2.0*sgm*sgm;

    x[0] = x0;
    double t = t0;
    for (unsigned int j=0; j<N; j++)
    {
        double k1 = f(t, x[j]);
        double k2 = f(t+ht/2.0, x[j]+(ht/2.0)*k1);
        double k3 = f(t+ht/2.0, x[j]+(ht/2.0)*k2);
        double k4 = f(t+ht, x[j]+ht*k3);
        x[j+1] = x[j] + (ht/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        //x[j+1] = x[j] + f(t, x[j])*ht;

        //x[j] = x[j] + q[0];// * a * exp(-((t-0.2)*(t-0.2))/b);
        //x[j] = x[j] + q[1];// * a * exp(-((t-0.5)*(t-0.5))/b);
        //x[j] = x[j] + q[2];// * a * exp(-((t-0.8)*(t-0.8))/b);

        if (j == (unsigned int)(0.2*N)) x[j+1] = x[j+1] + q[0];
        if (j == (unsigned int)(0.5*N)) x[j+1] = x[j+1] + q[1];
        if (j == (unsigned int)(0.8*N)) x[j+1] = x[j+1] + q[2];

        t = t + ht;
    }
    //x[N] = x[N] + q[0]+q[1]+q[2]+q[3]+q[4]+q[5]+q[6]+q[7]+q[8]+q[9];
    //x[N] = x[N] + (q[0]+q[1]+q[2]);

    //Printer::printVector(x, 10, "x:");
    //for (unsigned int j=0; j<x.size(); j++) printf("%.8f\n", x[j]);
}

double PointControl11::f(double t, double x) const
{
    return x + 2.0*t - t*t;
    //return 2.0*t;
}


void PointControl11::calculateP(const DoubleVector &q, const DoubleVector &x, DoubleVector &p)
{
    p.clear();
    p.resize(N+1);

    //p[N] = -1.0;
    p[N] = 2.0*(x[N]-x1);
    double t = t1;
    for (unsigned int j=N; j>0; j--)
    {
        double k1 = pf(t, p[j], x[j]);
        double k2 = pf(t-ht/2.0, p[j]-(ht/2.0)*k1, x[j]);
        double k3 = pf(t-ht/2.0, p[j]-(ht/2.0)*k2, x[j]);
        double k4 = pf(t-ht,     p[j]-ht*k3, x[j]);
        p[j-1] = p[j] - (ht/6) * (k1 + 2*k2 + 2*k3 + k4);
        //p[j-1] = p[j] - pf(t, p[j], x[j])*ht;
        t = t - ht;
    }

    //Printer::printVector(p, 10, "p:");
}

double PointControl11::pf(double t, double p, double x) const
{
    //return p;
    double h = ht/100.0;
    return -p * (f(t, x+h)-f(t, x-h))/(2.0*h);
}
