#include "iproblem2backward.h"

IProblem2Backward::IProblem2Backward()
{
//    a = 1.0;
//    lambda0 = 1.0;
//    lambda1 = 2.0;
//    lambda2 = 1.5;
//    theta = 5.0;

//    Lc = 2;
//    Lo = 3;
//    k.resize(Lc, Lo);
//    z.resize(Lc, Lo);

//    eta.resize(Lc);
//    eta[0] = 0.33;
//    eta[1] = 0.66;

//    xi.resize(Lo);
//    xi[0] = 0.25;
//    xi[1] = 0.50;
//    xi[2] = 0.75;

//    for (unsigned int i=0; i<Lc; i++)
//    {
//        for (unsigned int j=0; j<Lo; j++)
//        {
//            k[i][j] = 5.0;
//            z[i][j] = 8.0;
//        }
//    }
}

void IProblem2Backward::gridMethod(DoubleMatrix &p) const
{
    Dimension sd = spaceDimension(Dimension::DimensionX);
    Dimension td = timeDimension();

    unsigned int N = sd.sizeN();
    unsigned int M = td.sizeN();
    double hx = sd.step();
    double ht = td.step();

    p.clear();
    p.resize(M+1, N+1, 1.0);

    SpaceNodePDE sn;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = n; sn.x = n*hx;
        p[M][n] = initial(sn);
    }

    TimeNodePDE tn;
    for (unsigned int m=M-1; m!=UINT32_MAX; m--)
    {
        tn.i = m;
        tn.t = m*ht;

        DoubleMatrix w(N+1, N+1, 0.0);
        DoubleVector d(N+1, 0.0);

        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            d[n] = p[m+1][n] - ht*f(sn, tn);

            if (n == 0)
            {
                w[0][0] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*ht*lambda1)/(hx);
                w[0][1] = -(2.0*a*a*ht)/(hx*hx);

                d[0] += -((2.0*a*a*ht)/(hx))*g0(tn);
            }
            else if (n == N)
            {
                w[N][N-1] = -(2.0*a*a*ht)/(hx*hx);
                w[N][N-0] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*ht*lambda2)/(hx);

                d[N] += +((2.0*a*a*ht)/(hx))*g1(tn);
            }
            else
            {
                w[n][n-1] = -(a*a*ht)/(hx*hx);
                w[n][n+0] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                w[n][n+1] = -(a*a*ht)/(hx*hx);
            }


            for (unsigned int j=0; j<Lo; j++)
            {
                double _delta = delta(sn, j);

                for (unsigned int i=0; i<Lc; i++)
                {
                    //d[n] += -ht*k[j][j]*z[j][j] * _delta;

                    //unsigned int jinx = (unsigned int)(xi[j]*N);
                    //w[n][jinx] += -ht*k[i][j];

                    double h3 = hx*hx*hx;
                    double h32 = (1.0/(2.0*h3));
                    double h36 = (1.0/(6.0*h3));

                    for (unsigned int n1=0; n1<=N; n1++)
                    {
                        double dh = fabs(n1*hx - eta[i]);

                        if (dh <= hx)
                        {
                            w[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(hx+dh)) * h32 * _delta;
                        }

                        if (hx < dh && dh <= 2.0*hx)
                        {
                            w[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(3.0*hx-dh)) * h36 * _delta;
                        }
                    }
                }
            }
        }

        DoubleVector x(N+1);
        LinearEquation::GaussianElimination(w,d,x);
        //IPrinter::printVector(x);

        w.clear();
        d.clear();

        for (unsigned int n=0; n<=N; n++) p[m][n] = x[n];
    }
}

double IProblem2Backward::initial(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return -2.0*mu[sn.i]*(uT[sn.i]-U[sn.i]) + h(sn);
}

double IProblem2Backward::h(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    double x = sn.x;
    return x*x + 2.0*mu[sn.i]*(uT[sn.i] - U[sn.i]);
}

//double IProblem2Backward::U(const SpaceNodePDE &sn UNUSED_PARAM) const
//{
//    return 1.0;
//}

//double IProblem2Backward::mu(const SpaceNodePDE &sn UNUSED_PARAM) const
//{
//    return 1.0;
//}

double IProblem2Backward::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double IProblem2Backward::g0(const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 0.0;
}

double IProblem2Backward::g1(const TimeNodePDE &tn UNUSED_PARAM) const
{
    double t = tn.t;
    return 2.0*t + lambda2*t;
}

double IProblem2Backward::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    double x = sn.x;
    double t = tn.t;

    double res = x*x + 2.0*a*a*t - lambda0*x*x*t;

    double W = 0.0;
    for (unsigned int j=0; j<Lo; j++)
    {
        double vi = 0.0;
        for (unsigned int i=0; i<Lc; i++)
        {
            vi += k.at(i,j) * P(eta[i],t);
        }
        W += vi*delta(sn, j);
    }
    res += W;

    return res;
}

double IProblem2Backward::delta(const SpaceNodePDE &sn UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    //Dimension sd = spaceDimension(Dimension::DimensionX);
    //double hx = sd.step();

    // Approximation delta function using normal distribution formula
    //double sigma = 3.0*hx;
    //return (1.0/(sqrt(2.0*M_PI)*sigma)) * exp( -((sn.x-xi[j])*(sn.x-xi[j])) / (2.0*sigma*sigma) );

    // Approximation delta function using L4 Lagrange interpolation
    //double h3 = hx*hx*hx;
    //double h32 = (1.0/(2.0*h3));
    //double h36 = (1.0/(6.0*h3));
    //double dh = fabs(sn.x - xi[j]);
    //if (dh <= hx)                return ((2.0*hx-dh)*(hx-dh)*(hx+dh)) * h32;
    //if (hx < dh && dh <= 2.0*hx) return ((2.0*hx-dh)*(hx-dh)*(3.0*hx-dh)) * h36;
    //return 0.0;

    if ( sn.x == xi[j] ) return 1.0;
    return 0.0;
}

double IProblem2Backward::P(double x, double t) const
{
    return x*x*t;
}
