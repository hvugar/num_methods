#include "problem1L2.h"

void Problem1L2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem1L2 p;
    p.initialize();
    p.startOptimize();
}

Problem1L2::Problem1L2() {}

void Problem1L2::initialize()
{
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;
    withError = false;

    L = 2;

    N = 100;
    hx = 0.001;

    M = 1000;
    ht = 0.1;

    // initial temperatures
    vfi << +0.0;// << -5.0 << +5.0;
    // environment temperatures
    vtt << +0.0;// << -4.0 << +4.0;

    // initial temperature
    fi = 0.0;
    // environment temperature
    tt = 0.0;

    // золото
    {
        const double r0 = 19320.0; // kg/m^3   // плотность
        const double c0 = 130.0;   // C/(kg*S) // удельная теплоемкость
        const double k0 = 312.0;   //          // коэффициент теплопроводности

        const double h1 = 1000.0;      // коэффициент теплообмена ->
        const double h2 = 10.0;        // коэффициент теплообмена ->

        a = sqrt((k0/(c0*r0)));        // коэффициент температуропроворности
        lambda0 = h2/(c0*r0);          // коэффициент теплообмена ->
        lambda1 = h1/k0;               // коэффициент теплообмена ->
        lambda2 = h2/k0;               // коэффициент теплообмена ->
    }

    // мед
    {
        const double r0 = 8900.0;  // kg/m^3     // плотность
        const double c0 = 400.0;   // C/(kg*S)   // удельная теплоемкость
        const double k0 = 380.0;   // Vt/(m*S)   // коэффициент теплопроводности

        const double h1 = 1000.0;      // коэффициент теплообмена ->
        const double h2 = 10.0;        // коэффициент теплообмена ->

        a = sqrt((k0/(c0*r0)));        // коэффициент температуропроворности
        lambda0 = h2/(c0*r0);          // коэффициент теплообмена ->
        lambda1 = h1/k0;               // коэффициент теплообмена ->
        lambda2 = h2/k0;               // коэффициент теплообмена ->

        //printf("%.10f %.10f %.10f %.10f\n", lambda0, lambda1, lambda2, a);
    }

    /* коэффициенты регуляризации */
    alpha0 = 1.0;
    alpha1 = 0.0001;
    alpha2 = 0.0001;
    alpha3 = 0.0001;

    if (!optimizeK) alpha1 = 0.0;
    if (!optimizeZ) alpha2 = 0.0;
    if (!optimizeE) alpha3 = 0.0;

    k0 << 0.0 << 0.0;//-20.57571017653454 << -30.63314593795166;
    z0 << 0.0 << 0.0;//+10.33818417154749 << +10.47968970008047;
    e0 << 0.0 << 0.0;// +0.04500000000000 <<  +0.09500000000000;

    /* шаги числовых производных */
    hk = 0.001;
    hz = 0.001;
    he = 0.001;

    R = 0.0;

    /* пределы z параметров */
    vmin = -100.0;
    vmax = +100.0;
    d0 = (vmax+vmin)/2.0;
    d1 = (vmax-vmin)/2.0;

    /* температура стержня */
    V.resize(N+1);
    for (unsigned int n=0; n<=N; n++) V[n] = 10.0;
}

void Problem1L2::startOptimize()
{
    DoubleVector x0;
    if (optimizeK)
    {
        x0 << -8.5000 << -2.7000; //k
        //x0 << +1.0000 << +1.0000; //k
    }
    else
    {
        K.clear();
        K << -8.5000 << -2.7000; //k
    }

    if (optimizeZ)
    {
        x0 << +2.1000 << +4.9000; //z
        //x0 << +1.0000 << +1.0000; //z
    }
    else
    {
        Z.clear();
        Z << +2.1000 << +4.9000; //z
    }

    if (optimizeE)
    {
        x0 << +0.02000 << +0.08000; //e
    }
    else
    {
        E.clear();
        E << +0.02000 << +0.08000; //e
    }

    R = 1.0;
    optimize(x0);

        while (R < 10000000000.0)
        {
            IPrinter::printSeperatorLine();
            optimize(x0);
            R *= 10.0;
        }
}

void Problem1L2::optimize(DoubleVector &x0) const
{
    Problem1L2* p = const_cast<Problem1L2*>(this);

    ConjugateGradient g;
    g.setFunction(p);
    g.setGradient(p);
    g.setPrinter(p);
    g.setProjection(p);
    g.setOptimalityTolerance(0.0001);//0.00000001
    g.setStepTolerance(0.0001);//0.00000001
    g.setFunctionTolerance(0.0001);//0.00000001
    g.setR1MinimizeEpsilon(1.0, 0.00001); //0.00000001
    g.setNormalize(true);
    g.showExitMessage(true);
    g.setResetIteration(false);
    g.calculate(x0);

    DoubleMatrix u;
    DoubleVector k,z,e;
    getParameters(k,z,e,x0);

    IPrinter::printSeperatorLine();
    printf("Optimal k: %20.14f %20.14f\n", k[0], k[1]);
    printf("Optimal z: %20.14f %20.14f\n", z[0], z[1]);
    printf("Optimal e: %20.14f %20.14f\n", e[0], e[1]);
    IPrinter::printSeperatorLine();
}

void Problem1L2::project(DoubleVector &x UNUSED_PARAM, unsigned int i UNUSED_PARAM)
{
    unsigned int p = 0;
    if (optimizeK) p+=2;

    if (optimizeZ) p+=2;

    /* e lower/upper limits */
    if (optimizeE)
    {
        if (x.at(p) < 5*hx)  x.at(p) = 5*hx;
        if (x.at(p) > 0.045) x.at(p) = 0.045;

        if (x.at(p+1) < 0.055)     x.at(p+1) = 0.055;
        if (x.at(p+1) > (N-5)*hx)  x.at(p+1) = (N-5)*hx;
    }
}
