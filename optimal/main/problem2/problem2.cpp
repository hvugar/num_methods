#include "problem2.h"

IProblem2::IProblem2()
{
    a = 1.0;
    lambda0 = 1.0;
    theta = 5.0;

    Lc = 2;
    Lo = 3;
    k.resize(Lc, Lo);
    z.resize(Lc, Lo);
    eta.resize(Lc);
    eta[0] = 0.33;
    eta[1] = 0.66;
    xi.resize(Lo);
    xi[0] = 0.25;
    xi[1] = 0.50;
    xi[2] = 0.75;
}

double IProblem2::initial(const SpaceNodePDE &sn) const
{
    return 0.0;
}

double IProblem2::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary) const
{

}

double IProblem2::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double t = tn.t;
    double x = sn.x;

    double res = x*x - 2.0 * a*a * t + lambda0 * (x*x*t - theta);

    double W = 0.0;
    for (unsigned int i=0; i<Lc; i++)
    {
        if ( x == eta[i] )
        {
            double vi = 0.0;
            for (unsigned int j=0; j<Lo; j++)
            {
                vi = k.at(i,j) * (U(xi[j], t) - z.at(j,i));
            }
            W += vi;
        }
    }
    res -= W;

    return res;
}

double IProblem2::U(double x, double t) const
{
    return x*x*t;
}


