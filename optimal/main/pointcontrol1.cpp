#include "pointcontrol1.h"
#include <math.h>
#include <gradient_cjt.h>
#include "utils.h"

PointControl1::PointControl1(double t0, double t1, double x0, double x1, double dt, double dx)
{
    this->t0 = t0;
    this->t1 = t1;
    this->x0 = x0;
    this->x1 = x1;
    this->dt = dt;
    this->dx = dx;

    n = (unsigned int)(ceil(fabs(t1 - t0)/dt))+1;

    x.resize(n);
    psi.resize(n);

    T.resize(3);
    T[0] = 0.2;
    T[1] = 0.5;
    T[2] = 0.8;
}

double PointControl1::fx(const DoubleVector &p)
{
    calculate_x(p);
    return (x[n-1] - x1)*(x[n-1] - x1);
}

void PointControl1::gradient(double step, const DoubleVector& p, DoubleVector& g)
{
    calculate_x(p);
    calculate_psi();

    g[0] = psi[2000];
    g[1] = psi[5000];
    g[2] = psi[8000];

//    printf("p: %f %f %f\n", p[0], p[1], p[2]);
//    printX("x:", x);
//    printX("psi:", psi);
//    printf("g: %f %f %f\n", g[0], g[1], g[2]);
}

void PointControl1::calculate_x(const DoubleVector &p)
{
    x[0] = x0;
    double t = t0;
    double _x0 = x0;

    for (unsigned int i=0; i<n; i++)
    {
        double k1 = dxdt(t,        _x0, p);
        double k2 = dxdt(t+dt/2.0, _x0+(dt/2.0)*k1, p);
        double k3 = dxdt(t+dt/2.0, _x0+(dt/2.0)*k2, p);
        double k4 = dxdt(t+dt,     _x0+dt*k3, p);
        _x0 = _x0 + (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        t = t + dt;
        x[i] = _x0;
    }
//    printX("x:", x);

            FILE *f = fopen("d:/test.txt", "w");
            for (unsigned int i=0; i<x.size(); i++)
            {
                if (i%100==0)
                fprintf(f, "%f\n", x[i]);
            }
            fclose(f);
}

double PointControl1::px(double t, double psi, double x)
{
    return psi;
}

void PointControl1::calculate_psi()
{
    double t = t1;
    psi[n-1] = -1.0;
    double psi1 = psi[n-1];
    double h = -dt;

    for (int i=n-2; i>=0; i--)
    {
        double k1 = px(t,       psi1, x[i]);
        double k2 = px(t+h/2.0, psi1+(h/2.0)*k1, x[i]);
        double k3 = px(t+h/2.0, psi1+(h/2.0)*k2, x[i]);
        double k4 = px(t+h,     psi1+h*k3, x[i]);
        psi1 = psi1 + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        t = t + h;
        psi[i] = psi1;
    }
}

double PointControl1::f(double t, double x)
{
    return 3*t*t + t*t*t - x;
}

double PointControl1::dxdt(double t, double x, const DoubleVector& p)
{
    double sum = f(t, x);

    for (unsigned int i=0; i<T.size(); i++)
    {
        double t_ = T[i];
        double diff = fabs(t - t_);

        if (diff <= (dt/10.0))
        {
            sum += p[i]*1000;
        }
    }

    return sum;
}

double PointControl1::delta(double t)
{
    return 0.0;
}

void PointControl1::main()
{

    DoubleVector p(3, 0.0);
    p[0] = 20.4659868736;//2.3;
    p[1] = 28.5725057266;//3.5;
    p[2] = 22.0457546073;//1.7;

    PointControl1 f(0.0, 1.0, 0.0, +2.17901513, 0.0001, 0.0001);
    PointControl1Printer printer;

    ConjugateGradient g1;
    g1.setFunction(&f);
    g1.setEpsilon1(0.0000001);
    g1.setEpsilon2(0.000000000001);
    g1.setGradientStep(0.0000001);
    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
    g1.setPrinter(&printer);
    g1.calculate(p);
    printf("p: %f %f %f\n", p[0], p[1], p[2]);
}

void PointControl1Printer::print(unsigned int iterationCount, const DoubleVector &p, const DoubleVector &s, double m_alpha, RnFunction *f) const
{
    printf("J[%2d]: %.10f %.10f %.10f %.10f\n", iterationCount, f->fx(p), p[0], p[1], p[2]);
}
