#include "pointcontrol1.h"
#include <math.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>
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
//    printX("          fx: x", x);
    return (x[n-1] - x1)*(x[n-1] - x1);
}

void PointControl1::gradient(double step, const DoubleVector& p, DoubleVector& g)
{
    calculate_x(p);
    printX("x", x);
    calculate_psi();
    printX("psi", psi);

    g[0] = psi[2000];
    g[1] = psi[5000];
    g[2] = psi[8000];
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
    return 2*t - x + t*t;
}

double PointControl1::dxdt(double t, double x, const DoubleVector& p)
{
    double sum = f(t, x);
    if (fabs(t-T[0]) < dt/10.0) sum += p[0];
    if (fabs(t-T[1]) < dt/10.0) sum += p[1];
    if (fabs(t-T[2]) < dt/10.0) sum += p[2];
    return sum;
}

double PointControl1::delta(double t)
{
    return 0.0;
}

void PointControl1::main()
{
    DoubleVector p(3, 0.0);
    p[0] = 0.0;
    p[1] = 0.0;
    p[2] = 0.0;

    PointControl1 f(0.0, 1.0, 0.0, +22.0, 0.0001, 0.0001);
    PointControl1Printer printer;

    ConjugateGradient g1;
    g1.setFunction(&f);
    g1.setEpsilon1(0.0000001);
    g1.setEpsilon2(0.0000001);
    g1.setGradientStep(0.0000001);
    g1.setR1MinimizeEpsilon(1, 0.001);
    g1.setPrinter(&printer);
    g1.calculate(p);
    printf("p: %f %f %f\n", p[0], p[1], p[2]);

    f.write(f.x);
}

void PointControl1Printer::print(unsigned int iterationCount, const DoubleVector &p, const DoubleVector &s, double m_alpha, RnFunction *f) const
{
    printf("J[%2d]: %.10f %.10f %.10f %.10f\n", iterationCount, f->fx(p), p[0], p[1], p[2]);
    puts("*******************************************************************************");
}

void PointControl1::write(DoubleVector &x)
{
    FILE *f = fopen("pointcontrol1.txt", "w");
    for (unsigned int i=0; i<x.size(); i++)
    {
        if (i%10==0)
        fprintf(f, "%.10f\n", x[i]);
    }
    fclose(f);
}
