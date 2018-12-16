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

void CCIHyperbolicIBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    QGuiApplication app(argc, argv);

    CCIHyperbolicIBVP1 hibvp;
    hibvp.setTimeDimension(Dimension(0.01, 0, 100));
    hibvp.addSpaceDimension(Dimension(0.01, 0, 100));
    hibvp.addSpaceDimension(Dimension(0.01, 0, 100));

    //DoubleVector u;
    //hibvp.calculateD11(u, 2.0, 0.25);
    //IPrinter::printSeperatorLine();

    DoubleMatrix u;
    hibvp.calculateD2V1(u, 1.0);
    IPrinter::printSeperatorLine();
}

void CC1IHyperbolicIBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    QGuiApplication app(argc, argv);

    CC1IHyperbolicIBVP1 hibvp;
    hibvp.setTimeDimension(Dimension(0.01, 0, 100));
    hibvp.addSpaceDimension(Dimension(0.01, 0, 100));
    hibvp.addSpaceDimension(Dimension(0.01, 0, 100));

    DoubleMatrix u;
    hibvp.calculateD2V1(u, hibvp.a, hibvp.sigma);
    IPrinter::printSeperatorLine();
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

double CCIHyperbolicIBVP1::initial1(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    //const static double hx = spaceDimension(Dimension::DimensionX).step();
    //const static double hy = spaceDimension(Dimension::DimensionY).step();
    //const static double sigma = 3.0 * hx;

    //return sn.x*sn.x;
    //return sn.x*sn.x*sn.x;
    //return sn.x*sn.x*sn.x;
    //return 0.0;
    //return 1.0/(sqrt(2.0*M_PI)*sigma)*exp(-((0.5-sn.i*hx)*(0.5-sn.i*hx))/(2.0*sigma*sigma));//*0.05;

    // 2D
    //return sn.x*sn.x + sn.y*sn.y;
    //return sn.x*sn.x + sn.y*sn.y;
    return 0.0;
}

double distribute_blows(const Dimension &dimX, const Dimension &dimY, const SpaceNodePDE &sn)
{
    std::vector<SpacePoint> theta(2);
    DoubleVector q(2);
    q[0] = +0.145; theta[0].x = 0.2500; theta[0].y = 0.7200;
    q[1] = +0.147; theta[1].x = 0.6400; theta[1].y = 0.2700;

    double hx = dimX.step();
    double hy = dimY.step();

    const double sigmaX = 12.0*hx;
    const double sigmaY = 12.0*hy;
    double sum = 0.0;
    for (unsigned int s=0; s<2; s++)
    {
        double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);
        sum += q[s]*factor*exp(-0.5*(((sn.x-theta[s].x)*(sn.x-theta[s].x))/(sigmaX*sigmaX)+((sn.y-theta[s].y)*(sn.y-theta[s].y))/(sigmaY*sigmaY)));
    }
    return sum;
}

double distribute_control(const Dimension &dimX, const Dimension &dimY, const SpaceNodePDE &sn, const DoubleMatrix &u10)
{
    const static unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const static unsigned int M = static_cast<unsigned int> ( dimY.size() );

    DoubleMatrix k(2, 2);
    DoubleMatrix z(2, 2);
    std::vector<SpacePoint> xi(2);
    std::vector<SpacePoint> eta(2);

    k[0][0] = -0.1820; k[0][1] = -0.1250; k[1][0] = -0.1550; k[1][1] = -0.1310;
    z[0][0] = -0.0262; z[0][1] = -0.0773; z[1][0] = -0.0570; z[1][1] = +0.0653;
    xi[0].x  = +0.3849; xi[0].y  = +0.5442; xi[1].x  = +0.7661; xi[1].y  = +0.6785;
    eta[0].x = +0.6656; eta[0].y = +0.7909; eta[1].x = +0.4856; eta[1].y = +0.3810;

    double hx = dimX.step();
    double hy = dimY.step();

    const double sigmaX = hx;
    const double sigmaY = hy;

    static double u[2] = {0.0,0.0};
    if (sn.i == 1 && sn.j == 1)
    {
        for (unsigned int j=0; j<2; j++)
        {
            SpaceNodePDE sn0;
            for (unsigned int m=0; m<=M; m++)
            {
                sn0.j = m; sn0.y = sn0.j*hy;
                for (unsigned int n=0; n<=N; n++)
                {
                    sn0.i = n; sn0.x = sn0.i*hx;
                    double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);
                    double w = factor*exp(-0.5*(((sn0.x-xi[j].x)*(sn0.x-xi[j].x))/(sigmaX*sigmaX)+
                                                ((sn0.y-xi[j].y)*(sn0.y-xi[j].y))/(sigmaY*sigmaY)));
                    u[j] += u10[m][n] * hx * hy * w;
                }
            }
        }
    }

    double v[2] = {0.0, 0.0};
    double sum = 0.0;
    for (unsigned int i=0; i<2; i++)
    {
        for (unsigned int j=0; j<2; j++)
        {
            v[i] += k[i][j]*(u[j]-z[i][j]);
        }
        double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);
        double w = factor*exp(-0.5*(((sn.x-eta[i].x)*(sn.x-eta[i].x))/(sigmaX*sigmaX)+
                                    ((sn.y-eta[i].y)*(sn.y-eta[i].y))/(sigmaY*sigmaY)));
        sum += v[i] * w;
    }
    return sum;
}

double CCIHyperbolicIBVP1::initial2(const SpaceNodePDE &sn) const
{
    const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const static Dimension &dimY = spaceDimension(Dimension::DimensionY);
    //const static double hx = spaceDimension(Dimension::DimensionX).step();
    //const static double hy = spaceDimension(Dimension::DimensionY).step();
    //const static double sigma = 3.0 * hx;

    //return 0.0;
    //return 0.0;
    //return 0.0;
    //return 1.0/(sqrt(2.0*M_PI)*sigma)*exp(-((0.5-sn.i*hx)*(0.5-sn.i*hx))/(2.0*sigma*sigma))*0.05;
    //return 0.0;

    // 2D
    //return 1.0;
    //return 0.0;
    return distribute_blows(dimX, dimY, sn);
}

double CCIHyperbolicIBVP1::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    //const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    //const static Dimension &dimY = spaceDimension(Dimension::DimensionY);

    //return sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t*tn.t;
    //return 0.0;

    // 2D
    //return sn.x*sn.x + sn.y*sn.y + tn.t;
    //return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
    return 0.0;
}

double CCIHyperbolicIBVP1::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const static Dimension &dimY = spaceDimension(Dimension::DimensionY);

    //return 2.0 - 2.0*a*a;
    //return 2.0 - 6.0*a*a*sn.x;
    //return 6.0*tn.t - 6.0*a*a*sn.x;
    //return 0.0;

    // 2D
    //return 0.0 - (a*a)*(2.0+2.0);
    //return 2.0 - (a*a)*(2.0+2.0);
    return distribute_control(dimX, dimY, sn, this->u10);
}

void CCIHyperbolicIBVP1::layerInfo(const DoubleVector &u, unsigned int ln) const
{
    const static double ht = timeDimension().step();
    const static unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );

    DoubleVector &_u0 = const_cast<CCIHyperbolicIBVP1*>(this)->u0;
    DoubleVector &_u1 = const_cast<CCIHyperbolicIBVP1*>(this)->u1;
    DoubleVector &_u2 = const_cast<CCIHyperbolicIBVP1*>(this)->u2;
    DoubleVector &_ut = const_cast<CCIHyperbolicIBVP1*>(this)->ut;

    _ut.resize(N+1, 0.0);

    if (ln == 0) { const_cast<CCIHyperbolicIBVP1*>(this)->u0 = u; }
    if (ln == 1) { const_cast<CCIHyperbolicIBVP1*>(this)->u1 = u; }
    if (ln == 2) { const_cast<CCIHyperbolicIBVP1*>(this)->u2 = u; }

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

    if (ln%20000==0) IPrinter::printVector(u, nullptr, u.length());
    //if (ln > 1) printf("%d %f %f %f\n", ln, integralU1(u), integralU2(u), integralU1(u)+integralU2(u));
}

void CCIHyperbolicIBVP1::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    if (ln == 0) const_cast<CCIHyperbolicIBVP1*>(this)->u10 = u;
    if (ln == 1) const_cast<CCIHyperbolicIBVP1*>(this)->u10 = u;
    if (ln >= 2) const_cast<CCIHyperbolicIBVP1*>(this)->u10 = u;

    //if (ln == 0 || ln == 1 || ln == 2 || ln%100==0)
    {
        //IPrinter::printMatrix(u);
        //IPrinter::printSeperatorLine();

        char filename1[40];
        int size1 = sprintf(filename1, "e:/data/img1/image%d.png", ln);
        filename1[size1] = 0;

        double min = u10.min();
        double max = u10.max();

        QPixmap pxm;
        visualizeMatrixHeat(u10, min, max, pxm, 0, 0);
        pxm.save(QString(filename1), "PNG");
    }
}

double CCIHyperbolicIBVP1::integralU1(const DoubleVector &) const
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

double CCIHyperbolicIBVP1::integralU2(const DoubleVector &) const
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

CC1IHyperbolicIBVP1::~CC1IHyperbolicIBVP1() {}

double CC1IHyperbolicIBVP1::initial1(const SpaceNodePDE &sn) const
{
    //return sn.x*sn.x + sn.y*sn.y;
    //return sn.x*sn.x + sn.y*sn.y;
    return 0.0;
}

double CC1IHyperbolicIBVP1::initial2(const SpaceNodePDE &sn) const
{
    const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const static Dimension &dimY = spaceDimension(Dimension::DimensionY);

    //return 1.0;
    //return 0.0;
    return distribute_blows(dimX, dimY, sn);
}

double CC1IHyperbolicIBVP1::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    //return sn.x*sn.x + sn.y*sn.y + tn.t;
    //return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
    return 0.0;
}

double CC1IHyperbolicIBVP1::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const static Dimension &dimY = spaceDimension(Dimension::DimensionY);
    //return sigma - 4.0*a*a + 0.0;
    //return sigma - 4.0*a*a + 2.0;
    return distribute_control(dimX, dimY, sn, this->u10);
}

void CC1IHyperbolicIBVP1::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    if (ln == 0) const_cast<CC1IHyperbolicIBVP1*>(this)->u10 = u;
    if (ln == 1) const_cast<CC1IHyperbolicIBVP1*>(this)->u10 = u;
    if (ln >= 2) const_cast<CC1IHyperbolicIBVP1*>(this)->u10 = u;

    //if (ln == 0 || ln == 1 || ln == 2 || ln%100==0)
    {
        //IPrinter::printMatrix(u);
        //IPrinter::printSeperatorLine();

        char filename1[40];
        int size1 = sprintf(filename1, "e:/data/img2/image%d.png", ln);
        filename1[size1] = 0;

        double min = u10.min();
        double max = u10.max();

        QPixmap pxm;
        visualizeMatrixHeat(u10, min, max, pxm, 0, 0);
        pxm.save(QString(filename1), "PNG");
    }
}
