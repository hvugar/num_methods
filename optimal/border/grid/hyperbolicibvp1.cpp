#include "hyperbolicibvp1.h"

void HyperbolicIBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    HyperbolicIBVP1 p;
    p.setTimeDimension(Dimension(0.1, 10, 0));
    p.addSpaceDimension(Dimension(0.1, 10, 0));
    {
        unsigned int minN = p.spaceDimension(Dimension::DimensionX).min();
        unsigned int maxN = p.spaceDimension(Dimension::DimensionX).max();
        unsigned int N = p.spaceDimension(Dimension::DimensionX).size();

        unsigned int minM = p.timeDimension().min();
        unsigned int maxM = p.timeDimension().max();
        unsigned int M = p.timeDimension().size();

        DoubleMatrix u(M+1, N+1);

        clock_t t = clock();
        for (unsigned int j=minM; j<=maxM; j++)
        {
            for (unsigned int i=minN; i<=maxN; i++)
            {
                u[j-minM][i-minN] = p.U(i,j);
            }
        }
        t = clock() - t;
        IPrinter::printVector(14,10,u.row(0));
        IPrinter::printVector(14,10,u.row(1));
        IPrinter::printVector(14,10,u.row(2));
        IPrinter::printSeperatorLine();
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.gridMethod0(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    //    {
    //        DoubleMatrix u;
    //        clock_t t = clock();
    //        p.gridMethod1(u);
    //        t = clock() - t;
    //        IPrinter::printMatrix(14,10,u);
    //        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    //        IPrinter::printSeperatorLine();

    //        FILE* file = fopen("data.txt", "w");
    //        IPrinter::printMatrix(14, 10, u, u.rows(), u.cols(), NULL, file);
    //        fclose(file);
    //    }
    //    {
    //        DoubleMatrix u;
    //        clock_t t = clock();
    //        p.gridMethod2(u);
    //        t = clock() - t;
    //        IPrinter::printMatrix(14,10,u);
    //        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    //        IPrinter::printSeperatorLine();
    //    }
}

void HyperbolicIBVP2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    HyperbolicIBVP2 hibvp;
    hibvp.setTimeDimension(Dimension(0.005, 0, 200*50));
    hibvp.addSpaceDimension(Dimension(0.01, 0, 100));
    DoubleVector u;
    hibvp.calculateU1(u, 2.0, 0.25);
    IPrinter::printSeperatorLine();
    //    hibvp.calculateU2(u, 1.0);
    //    IPrinter::printVector(u);
    //    IPrinter::printSeperatorLine();
}

HyperbolicIBVP1::HyperbolicIBVP1()
{
}

double HyperbolicIBVP1::initial1(const SpaceNodePDE &sn) const
{
    return sn.x*sn.x;
}

double HyperbolicIBVP1::initial2(const SpaceNodePDE &sn) const
{
    return 2.0*sn.x*sn.x;
}

double HyperbolicIBVP1::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    C_UNUSED(sn);
    double t = tn.t;
    if (sn.i == 0)   return 0.0;
    if (sn.i == 100) return (t+1.0)*(t+1.0);
    return 0.0;
}

double HyperbolicIBVP1::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 2.0*sn.x*sn.x - 2.0 * a(sn,tn)*(tn.t+1.0)*(tn.t+1.0);
}

double HyperbolicIBVP1::a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    return 1.0;
}

double HyperbolicIBVP1::U(unsigned int n, unsigned int m) const
{
    double t = m*timeDimension().step();
    double x = n*spaceDimension(Dimension::DimensionX).step();
    return x*x*(t+1.0)*(t+1.0);
}

double HyperbolicIBVP2::initial1(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    const static double hx = spaceDimension(Dimension::DimensionX).step();
    static double sigma = 3.0 * hx;

    //return sn.x*sn.x;
    //return sn.x*sn.x*sn.x;
    //return sn.x*sn.x*sn.x;
    //return 0.0;
    return 1.0/(sqrt(2.0*M_PI)*sigma)*exp(-((0.5-sn.i*hx)*(0.5-sn.i*hx))/(2.0*sigma*sigma))*0.05;
}

double HyperbolicIBVP2::initial2(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    const static double hx = spaceDimension(Dimension::DimensionX).step();
    static double sigma = 3.0 * hx;

    //return 0.0;
    //return 0.0;
    //return 0.0;
    //return 1.0/(sqrt(2.0*M_PI)*sigma)*exp(-((0.5-sn.i*hx)*(0.5-sn.i*hx))/(2.0*sigma*sigma))*0.05;
    return 0.0;
}

double HyperbolicIBVP2::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    //return sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t*tn.t;
    return 0.0;
}

double HyperbolicIBVP2::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    //return 2.0 - 2.0*a*a;
    //return 2.0 - 6.0*a*a*sn.x;
    //return 6.0*tn.t - 6.0*a*a*sn.x;
    return 0.0;
}

void HyperbolicIBVP2::layerInfo(const DoubleVector &u, unsigned int ln) const
{
    const static double ht = timeDimension().step();
    const static unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );

    DoubleVector &_u0 = const_cast<HyperbolicIBVP2*>(this)->u0;
    DoubleVector &_u1 = const_cast<HyperbolicIBVP2*>(this)->u1;
    DoubleVector &_u2 = const_cast<HyperbolicIBVP2*>(this)->u2;
    DoubleVector &_ut = const_cast<HyperbolicIBVP2*>(this)->ut;

    _ut.resize(N+1, 0.0);

    if (ln == 0) { const_cast<HyperbolicIBVP2*>(this)->u0 = u; }
    if (ln == 1) { const_cast<HyperbolicIBVP2*>(this)->u1 = u; }
    if (ln == 2) { const_cast<HyperbolicIBVP2*>(this)->u2 = u; }

    if (ln >= 3)
    {
        _u0 = _u1; _u1 = _u2; _u2 = u;
    }

    if (ln >= 2)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            _ut[n] = (_u2[n] - _u0[n])/(2.0*ht);
        }
    }

    //if (ln%2000==0) IPrinter::printVector(u, nullptr, u.length());
    if (ln > 1) printf("%d %f %f %f\n", ln, integralU1(u), integralU2(u), integralU1(u)+integralU2(u));
}

double HyperbolicIBVP2::integralU1(const DoubleVector &) const
{
    const static double hx = spaceDimension(Dimension::DimensionX).step();
    const static unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );

    double sum = 0.0;

    sum += 0.50 * u1[0]*u1[0];
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += u1[n]*u1[n];
    }
    sum += 0.50 * u1[N]*u1[N];

    return sum*hx;
}

double HyperbolicIBVP2::integralU2(const DoubleVector &) const
{
    const double hx = spaceDimension(Dimension::DimensionX).step();
    const unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );

    double sum = 0.0;

    sum += 0.50 * ut[0]*ut[0];
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += ut[n]*ut[n];
    }
    sum += 0.50 * ut[N]*ut[N];

    return sum*hx;
}
