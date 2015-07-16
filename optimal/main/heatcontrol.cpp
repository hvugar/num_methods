#include "heatcontrol.h"
#include <math.h>

void HeatGradientMethod::calculateGradient()
{
    unsigned int n = heatControl->mg.size();

    if (m_g.size() != n)
        m_g.resize(n, 0.0);

    for (unsigned int i=0; i<m_g.size(); i++) m_g[i] = heatControl->mg[i];
}

HeatControl::HeatControl()
{
    t0 = 0.0;
    t1 = 1.0;

    x0 = 0.0;
    x1 = 1.0;

    dx = 0.01;
    dt = 0.01;

    n = (int)(ceil((x1-x0)/dx)) + 1;
    m = (int)(ceil((t1-t0)/dt)) + 1;

    mt.resize(m, 0.0);
    mx.resize(n, 0.0);
    mu.resize(m*n, 0.0);
    mf.resize(m*n, 0.0);
    mg.resize(m*n, 0.0);

    for (int i=0; i<n; i++) mx[i] = dx * i;
    for (int j=0; j<m; j++) mt[j] = dt * j;

//    for (int j=0; j<m; j++) mu.push_back(DoubleVector(n, 0.0));
//    for (int j=0; j<m; j++) mf.push_back(DoubleVector(n, 0.0));
//    for (int j=0; j<m; j++) mg.push_back(DoubleVector(n, 0.0));

    gradient.heatControl = this;
}

HeatControl::~HeatControl()
{
    mt.clear();
    mx.clear();
    mu.clear();
    mf.clear();
    mg.clear();
}

double HeatControl::U(double x)
{
    return x*x + 2.0*x + 1.0;
}

double HeatControl::f(double x, double t)
{
   return 2.0*t - 2.0;
}

double HeatControl::fxt1(double x, double t)
{
    unsigned int j = (unsigned int)(ceil(t/dt));
    unsigned int i = (unsigned int)(ceil(x/dx));
    return mf[j*m+i];
}

double HeatControl::u(double x, double t)
{
    return 0.0;
}

double HeatControl::fi(double x)
{
    return x*x + 2.0*x;
}

double HeatControl::m1(double t)
{
    return t*t;
}

double HeatControl::m2(double t)
{
    return t*t + 3.0;
}

void HeatControl::calculate_u()
{
    GridMethod gm;
    gm.alpha = 1.0;
    gm.dt = dt;
    gm.dx = dx;
    gm.t0 = t0;
    gm.t1 = t1;
    gm.x0 = x0;
    gm.x1 = x1;

    gm.m = m;
    gm.n = n;

    gm.setF(f);
    gm.m1 = m1;
    gm.m2 = m2;
    gm.fi = fi;

    gm.implicitDifferenceScheme();

//    int j;
//    implicit_difference_scheme1(fxt1, fi, m1, m2, p->alpha, p->dx, p->dt, p->x0, p->x1, p->t0, p->t1, &g);
//    for (j=0; j<g.m; j++)
//    {
//        memcpy(p->u[j], g.u[j], sizeof(double) * g.n);
//        free(g.u[j]);
//    }
//    free(g.u);
}

double HeatControl::fx( const std::vector<double>& f )
{
    double sum = 0.0;
    for (int i=0; i<(n-1); i++)
    {
        int i0=i;
        int i1=i+1;
        double fj = mu[m-1*m+i1] - U(dx*i1);
        double fi = mu[m-1*m+i0] - U(dx*i0);
        double f0 = fj*fj + fi*fi;
        sum = sum + 0.5 * f0 * (dx);
    }

    double F = 0.0;
    for (int j=0; j<m; j++)
    {
        for (int i=0; i<n; i++)
        {
            F += (f[j*m+i]-(2.0*(dt*j)-2.0)) * (f[j*m+i]-(2.0*(dt*j)-2.0));
        }
    }
    sum = sum + F;

    return sum;
}


