#include "heatcontrol.h"
#include <math.h>

HeatControl::HeatControl() : gradient(NULL)
{
    t0 = 0.0;
    t1 = 1.0;

    x0 = 0.0;
    x1 = 1.0;

    dx = 0.01;
    dt = 0.01;

    n = (int)(ceil((x1-x0)/dx)) + 1;
    m = (int)(ceil((t1-t0)/dt)) + 1;

    t.resize(m, 0.0);
    x.resize(n, 0.0);

    for (int i=0; i<n; i++) x[i] = dx * i;
    for (int j=0; j<m; j++) t[j] = dt * j;

    u.resize(m, std::vector<double>(n, 0.0));
    //p->u = (double **) malloc( sizeof(double*) * p->m );
    //p->f = (double **) malloc( sizeof(double*) * p->m );
    //p->p = (double **) malloc( sizeof(double*) * p->m );
    //p->u1 = (double **) malloc( sizeof(double*) * p->m );

}

HeatControl::~HeatControl()
{

}

void HeatControl::calculateGradient()
{
}

double HeatControl::_y(double x, double t)
{
    return x*x + 2.0*x + 1.0;
}

double HeatControl::_f(double x, double t)
{
   return 2.0*t - 2.0;
}

double HeatControl::_u(double x, double t)
{
    return 0.0;
}

double HeatControl::_fi(double x)
{
    return x*x + 2.0*x;
}

double HeatControl::_m1(double t)
{
    return t*t;
}

double HeatControl::_m2(double t)
{
    return t*t + 3.0;
}

double HeatControl::JSum()
{
    return 0.0;
}


