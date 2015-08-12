#include "pointcontrol.h"
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>
#include <rungekutta.h>
#include <math.h>
#include "utils.h"

PointControl::PointControl(double t0, double t1, double x0, double x1, double dt, double dx)
{
    this->t0 = t0;
    this->t1 = t1;
    this->x0 = x0;
    this->x1 = x1;
    this->dx = dx;
    this->dt = dt;

    n = (unsigned int)(ceil((t1 - t0)/dt))+1;
    epsilon = 0.02;

    T.resize(2);
    T[0] = 0.2;
    T[1] = 0.5;

//    p.resize(2);
//    p[0] = 0.3;
//    p[1] = 0.4;
}

double PointControl::fx(const DoubleVector &p)
{
    DoubleVector x(n, 0.0);
    calculate_x(x, p);
    return (x[n-1] - x1)*(x[n-1] - x1);// + (p[0]*p[0] + p[1]*p[1]);
}

void PointControl::gradient(double step, const DoubleVector &p, DoubleVector &g)
{
    DoubleVector x(n, 0.0);
    DoubleVector psi(n, 0.0);

    calculate_x(x, p);
    calculate_psi(x, psi);

    g[0] = psi[30000];
    g[1] = psi[40000];

    printf("p: %f %f\n", p[0], p[1]);
    printX("x:", x);
    printX("psi:", psi);
    printf("g: %f %f\n", g[0], g[1]);

//    FILE *f = fopen("d:/test.txt", "w");
//    for (unsigned int i=0; i<x.size(); i++)
//    {
//        if (i%100==0)
//        fprintf(f, "%f\n", x[i]);
//    }
//    fclose(f);
}

double PointControl::f(double t, double x)
{
    return 3.0*t*t + t*t*t - x;
}

double PointControl::delta(double t)
{
    for (unsigned int i=0; i<T.size(); i++)
    {
        if (fabs(t - T[i]) < ((epsilon / 2.0) + 0.000001)) return 1.0 / epsilon;
    }
    return 0.0;
}

double PointControl::dxdt(double t, double x)
{
    double sum = f(t, x);
    for (unsigned int i=0; i<p.size(); i++)
    {
        sum += p[i]*delta(t);
    }
    return sum;
}

double PointControl::px(double t, double psi, double x)
{
    return psi;
}

void PointControl::calculate_x(DoubleVector &x, const DoubleVector& p)
{
    this->p = p;
    x[0] = x0;
    double _t0 = t0;
//    double _t1 = t1;
    double _x0 = x0;
//    double _x1 = x1;

    for (unsigned int i=0; i<n; i++)
    {
        double k1 = dxdt(_t0,        _x0);
        double k2 = dxdt(_t0+dt/2.0, _x0+(dt/2.0)*k1);
        double k3 = dxdt(_t0+dt/2.0, _x0+(dt/2.0)*k2);
        double k4 = dxdt(_t0+dt,     _x0+dt*k3);
        _x0 = _x0 + (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        _t0 = _t0 + dt;
        x[i] = _x0;
    }

//    while (_t0 <= _t1)
//    {
//        double k1 = dxdt(_t0,        _x0);
//        double k2 = dxdt(_t0+dt/2.0, _x0+(dt/2.0)*k1);
//        double k3 = dxdt(_t0+dt/2.0, _x0+(dt/2.0)*k2);
//        double k4 = dxdt(_t0+dt,     _x0+dt*k3);

//        _x0 = _x0 + (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
//        _t0 = _t0 + dt;

//        i++;
//        x[i] = _x0;
//    }
}

void PointControl::calculate_psi(DoubleVector &x, DoubleVector &psi)
{
//    int i = n-1;
    psi[n-1] = -1.0;
    double _t0 = t0;
    double _t1 = t1;
    double psi1 = psi[n-1];
//    double _x1 = x1;
    double h = -dt;

    //i = n - 2;

    for (int i=n-2; i>=0; i--)
    {
        double k1 = px(_t1, psi1, x[i]);
        double k2 = px(_t1+h/2.0, psi1+(h/2.0)*k1, x[i]);
        double k3 = px(_t1+h/2.0, psi1+(h/2.0)*k2, x[i]);
        double k4 = px(_t1+h, psi1+h*k3, x[i]);
        psi1 = psi1 + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        _t1 = _t1 + h;
        psi[i] = psi1;
    }

//    while (_t0 <= _t1)
//    {
//        double k1 = px(_t1, psi1, x[i]);
//        double k2 = px(_t1+h/2.0, psi1+(h/2.0)*k1, x[i]);
//        double k3 = px(_t1+h/2.0, psi1+(h/2.0)*k2, x[i]);
//        double k4 = px(_t1+h, psi1+h*k3, x[i]);

//        psi1 = psi1 + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
//        _t1 = _t1 + h;

//        //printf("%f %f %f %d\n", t0, _t1, psi1, i);
//        psi[i+1] = psi1;
//        i--;
//    }
}

void PointControl::main()
{
    DoubleVector p0(2,0.0);
    PointControl pc(0.0, 1.0, 0.0, +1.73925985, 0.00001, 0.00001);

//    /* initial point */
//    for (unsigned int i=0; i<p0.size(); i++) p0[i] = 0.00001;
    p0[0] = +10.0;
    p0[1] = +0.7925423850;

//    p0[0] = +0.3;
//    p0[1] = +0.4;
//    /* Minimization */
    ConjugateGradient g1;
    g1.setFunction(&pc);
    g1.setEpsilon1(0.0000001);
    g1.setEpsilon2(0.000000000001);
    g1.setGradientStep(0.0000001);
    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
//    g1.setPrinter(&cfp);
    g1.calculate(p0);
    printf("p: %f %f\n", p0[0], p0[1]);
}
