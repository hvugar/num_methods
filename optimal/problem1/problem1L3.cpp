#include "problem1L3.h"
#include <time.h>

void Problem1L3::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem1L3 p;

    DoubleMatrix u(p.M+1,p.N+1);

    DoubleVector k; k << -1.1 << -2.5 << -0.1;
    DoubleVector z; z << 10.2 << 12.5 << 20.5;
    DoubleVector e; e << 0.25 << 0.50 << 0.75;

    //clock_t t = clock();
    p.calculateU2(u,k,z,e);
    IPrinter::printVector(14,10,u.row(1));
    p.calculateU3(u,k,z,e);
    IPrinter::printVector(14,10,u.row(1));
    //t = clock() - t;
    //printf("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);


    IPrinter::printSeperatorLine();
}

Problem1L3::Problem1L3()
{
    L = 3;

    N = 100;
    M = 100;
    hx = 0.01;
    ht = 0.01;
    //h  = 0.001;

    lambda0 = 0.001;
    lambda1 = 1000.0;
    lambda2 = 1.0;
    a = 1.0;

    fi = 2.0;
    tt = 3.0;
    return;
    vfi << fi;
    vtt << tt;

    alpha0 = 1.0;
    alpha1 = 0.00;
    alpha2 = 0.00;
    alpha3 = 0.00;

    k0 << 0.0 << 0.0 << 0.0;
    z0 << 0.0 << 0.0 << 0.0;
    e0 << 0.0 << 0.0 << 0.0;

    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    V.resize(N+1);
    R = 0.0;
    vmin = -0;
    vmax = +0;
    d0 = (vmax+vmin)/2.0;
    d1 = (vmax-vmin)/2.0;

    /* init V */
    for (unsigned int n=0; n<=N; n++)
    {
        //double h1 = 0.4/N;
        V[n] = 4.0;//4.2 - n*h1;
    }
    IPrinter::printVector(14, 10, V,"V: ");

    DoubleVector x0;
    x0 << -3.5000 << -3.7000 << -3.4000; //k
    x0 << +4.1000 << +3.9000 << +3.7000; //z
    x0 << +0.2500 << +0.5000 << +0.7500; //e

    //k << 3.50 << 3.70;
    //z << 4.20 << 4.20;
    //e << 0.40 << 0.90;

    ConjugateGradient g;
    g.setFunction(this);
    g.setGradient(this);
    g.setPrinter(this);
    g.setProjection(this);
    g.setOptimalityTolerance(0.00000001);
    g.setStepTolerance(0.00000001);
    g.setFunctionTolerance(0.00000001);
    g.setR1MinimizeEpsilon(2.0, 0.00000001);
    g.setNormalize(true);
    g.calculate(x0);
}

void Problem1L3::project(DoubleVector &x, size_t i)
{
    /* z lower/upper limits */
    if (x.at(3) < 3.80) x.at(3) = 3.80;
    if (x.at(3) > 4.20) x.at(3) = 4.20;

    if (x.at(4) < 3.80) x.at(4) = 3.80;
    if (x.at(4) > 4.20) x.at(4) = 4.20;

    if (x.at(5) < 3.80) x.at(5) = 3.80;
    if (x.at(5) > 4.20) x.at(5) = 4.20;

    /* e lower/upper limits */
    if (x.at(6) < 5*hx) x.at(6) = 5*hx;
    if (x.at(6) > 0.30) x.at(6) = 0.30;

    if (x.at(7) < 0.35) x.at(7) = 0.35;
    if (x.at(7) > 0.65) x.at(7) = 0.65;

    if (x.at(8) < 0.70)  x.at(8) = 0.70;
    if (x.at(8) > (N-5)*hx) x.at(8) = (N-5)*hx;
}
