#include "problem4ex2.h"
#include <math.h>
#include <utils/random.h>

void Problem4Ex2::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem4Ex2 p2;
    p2.grid = ODEGrid(Dimension(0.001, 1000, 0));
    p2.initialize();

    DoubleVector x0;
    x0 << Random::value(-3,3,5);
    x0 << Random::value(-3,3,5);
    x0 << Random::value(-3,3,5);
    x0 << Random::value(-3,3,5);
    x0 << Random::value(-3,3,5);
    x0 << Random::value(-3,3,5);
    IPrinter::print(x0);

    DoubleVector x;

    NonLinearEquationEx2 nle(p2);
    nle.calculateNewtonMethodMod(x0, x, 0.001, 0.00001);

    printf("0 %14.10f\n",nle.fx(x,0));
    printf("1 %14.10f\n",nle.fx(x,1));
    printf("2 %14.10f\n",nle.fx(x,2));
    printf("3 %14.10f\n",nle.fx(x,3));
    printf("4 %14.10f\n",nle.fx(x,4));
    printf("5 %14.10f\n",nle.fx(x,5));

    IPrinter::print(x,x.length());
    IPrinter::printSeperatorLine();
    p2.printResult1(x);

    Dimension dim = p2.grid.dimension();
    unsigned int N = dim.sizeN();
    double h = dim.step();

    IPrinter::printSeperatorLine();
    DoubleVector x1(N+1); for (unsigned int i=0; i<=N; i++) x1[i] = p2.X(i*h, 0); IPrinter::printVector(10,6,x1);
    DoubleVector x2(N+1); for (unsigned int i=0; i<=N; i++) x2[i] = p2.X(i*h, 1); IPrinter::printVector(10,6,x2);
    DoubleVector x3(N+1); for (unsigned int i=0; i<=N; i++) x3[i] = p2.X(i*h, 2); IPrinter::printVector(10,6,x3);

//    double aa=0.0;
//    unsigned int row = 2;
//    DoubleVector xx;
//    xx <<p2.X(0.3,0) << p2.X(0.3,1) << p2.X(0.3,2) << p2.X(0.6,0) << p2.X(0.6,1) << p2.X(0.6,2) ;
//    for (unsigned int col=0; col<3; col++) aa += p2.C(0.0, 0, 0, row, col) * p2.g(xx,0,col);
//    for (unsigned int col=0; col<3; col++) aa += p2.C(0.0, 0, 1, row, col) * p2.g(xx,1,col);
//    printf("%.10f\n", aa);
}


Problem4Ex2::Problem4Ex2()
{}

Problem4Ex2::~Problem4Ex2()
{}

void Problem4Ex2::initialize()
{
    unsigned int n = 3;
    Dimension dim = grid.dimension();
    unsigned int N = dim.sizeN();

    LinearODE1stOrder::Condition c0;
    c0.time = 0.0;
    c0.nmbr = 0;
    c0.mtrx.resize(n, n);
//    Random::fillMatrix(c0.mtrx, -3, 3, 5);
    c0.mtrx[0][0] = 1.0000000000; c0.mtrx[0][1] = 7.0000000000; c0.mtrx[0][2] = 4.0000000000;
    c0.mtrx[1][0] = 0.0000000000; c0.mtrx[1][1] = 9.0000000000; c0.mtrx[1][2] = 4.0000000000;
    c0.mtrx[2][0] = 8.0000000000; c0.mtrx[2][1] = 8.0000000000; c0.mtrx[2][2] = 2.0000000000;

    LinearODE1stOrder::Condition c1;
    c1.time = 0.2;
    c1.nmbr = N/5;
    c1.mtrx.resize(n, n);
//    Random::fillMatrix(c1.mtrx, -3, 3, 5);
    c1.mtrx[0][0] = 4.0000000000; c1.mtrx[0][1] = 5.0000000000; c1.mtrx[0][2] = 5.0000000000;
    c1.mtrx[1][0] = 1.0000000000; c1.mtrx[1][1] = 7.0000000000; c1.mtrx[1][2] = 2.0000000000;
    c1.mtrx[2][0] = 1.0000000000; c1.mtrx[2][1] = 5.0000000000; c1.mtrx[2][2] = 2.0000000000;

    LinearODE1stOrder::Condition c2;
    c2.time = 0.5;
    c2.nmbr = N/2;
    c2.mtrx.resize(n, n);
//    Random::fillMatrix(c2.mtrx, -3, 3, 5);
    c2.mtrx[0][0] = 7.0000000000; c2.mtrx[0][1] = 6.0000000000; c2.mtrx[0][2] = 1.0000000000;
    c2.mtrx[1][0] = 4.0000000000; c2.mtrx[1][1] = 2.0000000000; c2.mtrx[1][2] = 3.0000000000;
    c2.mtrx[2][0] = 2.0000000000; c2.mtrx[2][1] = 2.0000000000; c2.mtrx[2][2] = 1.0000000000;

    LinearODE1stOrder::Condition c3;
    c3.time = 0.8;
    c3.nmbr = 4*(N/5);
    c3.mtrx.resize(n, n);
//    Random::fillMatrix(c3.mtrx, -3, 3, 5);
    c3.mtrx[0][0] = 6.0000000000; c3.mtrx[0][1] = 8.0000000000; c3.mtrx[0][2] = 5.0000000000;
    c3.mtrx[1][0] = 7.0000000000; c3.mtrx[1][1] = 6.0000000000; c3.mtrx[1][2] = 1.0000000000;
    c3.mtrx[2][0] = 8.0000000000; c3.mtrx[2][1] = 9.0000000000; c3.mtrx[2][2] = 2.0000000000;

    LinearODE1stOrder::Condition c4;
    c4.time = 1.0;
    c4.nmbr = N;
    c4.mtrx.resize(n, n);
//    Random::fillMatrix(c4.mtrx, -3, 3, 5);
    c4.mtrx[0][0] = 7.0000000000; c4.mtrx[0][1] = 9.0000000000; c4.mtrx[0][2] = 5.0000000000;
    c4.mtrx[1][0] = 4.0000000000; c4.mtrx[1][1] = 3.0000000000; c4.mtrx[1][2] = 1.0000000000;
    c4.mtrx[2][0] = 2.0000000000; c4.mtrx[2][1] = 3.0000000000; c4.mtrx[2][2] = 3.0000000000;

    std::vector<LinearODE1stOrder::Condition> cs;
    cs.push_back(c0);
    cs.push_back(c1);
    cs.push_back(c2);
    cs.push_back(c3);
    cs.push_back(c4);

    DoubleVector betta(n);
    unsigned int L0 = cs.size();

    for (unsigned int row=0; row<n; row++)
    {
        betta[row] = 0.0;
        for (unsigned int s=0; s<L0; s++)
        {
            const LinearODE1stOrder::Condition &c = cs.at(s);
            for (unsigned int i=0; i<n; i++) betta[row] += c.mtrx[row][i] * X(c.time, i);
        }
    }

   // IPrinter::print(betta);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    Problem4Ex2Zetta0 zett01(*this);
    zett01.setGrid(grid);
    zett01.calculate(cs, betta, zm0);

    DoubleMatrix betta0(n, n, 0.0);

    Problem4Ex2Zettai zett1(*this, 0);
    zett1.setGrid(grid);
    zett1.calculateM(cs, betta0, zm1);

    Problem4Ex2Zettai zett2(*this, 1);
    zett2.setGrid(grid);
    zett2.calculateM(cs, betta0, zm2);
}

double Problem4Ex2::A(double t UNUSED_PARAM, unsigned int k UNUSED_PARAM, unsigned int row UNUSED_PARAM, unsigned int col UNUSED_PARAM) const
{
#ifdef SAMPLE_1
    if (row == 0) { if (col == 0) { return +2.0; } if (col == 1) { return t; }      if (col == 2) { return -3.0; } }
    if (row == 1) { if (col == 0) { return +3.0; } if (col == 1) { return -4.0*t; } if (col == 2) { return -8.0; } }
    if (row == 2) { if (col == 0) { return +t; }   if (col == 1) { return +1.0; }   if (col == 2) { return -1.0; } }
#endif
#ifdef SAMPLE_2
    if (row == 0) { if (col == 0) return +1.0; if (col == 1) return +3.0; if (col == 2) return +5.0; }
    if (row == 1) { if (col == 0) return +8.0; if (col == 1) return +7.0; if (col == 2) return +4.0; }
    if (row == 2) { if (col == 0) return +3.0; if (col == 1) return +1.0; if (col == 2) return +9.0; }
#endif
    return NAN;
}

double Problem4Ex2::B(double t UNUSED_PARAM, unsigned int k UNUSED_PARAM, unsigned int row UNUSED_PARAM) const
{
    double K = 1.0;
#ifdef SAMPLE_1
    if (row == 0) return 3.0*t*t*t - 4.0*t*t + 6.0*t - 3.0 - 0.0414023610*K;
    if (row == 1) return 8.0*t*t*t + 5.0*t*t - 7.0*t - 4.0 - 0.0735341450*K;
    if (row == 2) return 2.0*t*t - 3.0*t             + 4.0 - 0.0639390740*K;
#endif
#ifdef SAMPLE_2
    if (row == 0) return +6.0*cos(6.0*t)    - (sin(6.0*t)     + 3.0*cos(4.0*t) + 5.0*(t*t*t*t-2.0*t*t+1.0)) - 0.1403064559*K;
    if (row == 1) return -4.0*sin(4.0*t)    - (8.0*sin(6.0*t) + 7.0*cos(4.0*t) + 4.0*(t*t*t*t-2.0*t*t+1.0)) - 0.4107082389*K;
    if (row == 2) return +4.0*t*t*t - 4.0*t - (3.0*sin(6.0*t) + cos(4.0*t)     + 9.0*(t*t*t*t-2.0*t*t+1.0)) - 0.2880777651*K;
#endif
    return NAN;
}

double Problem4Ex2::C(double t UNUSED_PARAM, unsigned int k UNUSED_PARAM, unsigned int num UNUSED_PARAM, unsigned int row UNUSED_PARAM, unsigned int col UNUSED_PARAM) const
{
    double K = 1.0;
#ifdef SAMPLE_1

    if ( num == 0 )
    {
        if (row == 0) { if (col == 0) { return +0.002*K; } if (col == 1) { return +0.005*K; } if (col == 2) { return +0.003*K; } }
        if (row == 1) { if (col == 0) { return +0.004*K; } if (col == 1) { return +0.008*K; } if (col == 2) { return +0.001*K; } }
        if (row == 2) { if (col == 0) { return +0.001*K; } if (col == 1) { return +0.003*K; } if (col == 2) { return +0.004*K; } }
    }

    if ( num == 1 )
    {
        if (row == 0) { if (col == 0) { return +0.001*K; } if (col == 1) { return +0.003*K; } if (col == 2) { return +0.004*K; } }
        if (row == 1) { if (col == 0) { return +0.002*K; } if (col == 1) { return +0.005*K; } if (col == 2) { return +0.001*K; } }
        if (row == 2) { if (col == 0) { return +0.005*K; } if (col == 1) { return +0.002*K; } if (col == 2) { return +0.008*K; } }
    }

#endif
#ifdef SAMPLE_2
    if ( num == 0 )
    {
        if (row == 0) { if (col == 0) { return +0.01*K; } if (col == 1) { return +0.03*K; } if (col == 2) { return +0.09*K; } }
        if (row == 1) { if (col == 0) { return +0.05*K; } if (col == 1) { return +0.08*K; } if (col == 2) { return +0.01*K; } }
        if (row == 2) { if (col == 0) { return +0.03*K; } if (col == 1) { return +0.05*K; } if (col == 2) { return +0.04*K; } }
    }
    if ( num == 1 )
    {
        if (row == 0) { if (col == 0) { return +0.08*K; } if (col == 1) { return +0.06*K; } if (col == 2) { return +0.05*K; } }
        if (row == 1) { if (col == 0) { return +0.01*K; } if (col == 1) { return +0.03*K; } if (col == 2) { return +0.04*K; } }
        if (row == 2) { if (col == 0) { return +0.02*K; } if (col == 1) { return +0.05*K; } if (col == 2) { return +0.02*K; } }
    }
#endif
    return NAN;
}

double Problem4Ex2::X(double t, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0) { return t*t+t+2.0; }
    if (i == 1) { return 2.0*t-3.0; }
    if (i == 2) { return t*t*t+t; }
#endif
#ifdef SAMPLE_2
    if (i == 0) { return sin(6.0*t); }
    if (i == 1) { return cos(4.0*t); }
    if (i == 2) { return t*t*t*t - 2.0*t*t + 1.0; }
#endif
    return NAN;
}

double Problem4Ex2::g(const DoubleVector &x, unsigned int num, unsigned int row) const
{
#ifdef SAMPLE_1
    double x1_03 = x[0]; double x1_06 = x[3];
    double x2_03 = x[1]; double x2_06 = x[4];
    double x3_03 = x[2]; double x3_06 = x[5];

    if (num == 0)
    {
        if (row == 0) return x1_03*x1_03 + x2_03*x2_03 - x3_03;
        if (row == 1) return x2_03*x2_03 - x1_03*x3_03;
        if (row == 2) return x1_03 + x2_03 + x3_03*x3_03;
    }
    if (num == 1)
    {
        if (row == 0) { return  x1_06*x1_06; }
        if (row == 1) return x2_06*x2_06*x2_06;
        if (row == 2) return x3_06*x3_06;
    }
#endif
#ifdef SAMPLE_2
    double x1_03 = x[0]; double x1_06 = x[3];
    double x2_03 = x[1]; double x2_06 = x[4];
    double x3_03 = x[2]; double x3_06 = x[5];

    if (num == 0)
    {
        if (row == 0) return x1_03*x1_03*x1_03 + x2_03*x2_03*x2_03 + x3_03*x3_03*x3_03;
        if (row == 1) return x1_03*x1_03*x1_03 + exp(x2_03) + exp(x3_03);
        if (row == 2) return x1_03*x1_03*x1_03 + x2_03*x2_03*x2_03 + x3_03;
    }
    if (num == 1)
    {
        if (row == 0) return 3.0*x1_06*x1_06*x1_06 + 4.0*x2_06*x2_06*x2_06 + x3_06;
        if (row == 1) return x1_06 + x2_06 + 2.0*x3_06;
        if (row == 2) return x1_06 + x2_06 + x3_06;
    }
#endif

    return NAN;
}

void Problem4Ex2::printResult1(const DoubleVector &x)
{
    Dimension dim = grid.dimension();
    unsigned int N = dim.sizeN();

    DoubleVector x1(N+1);
    DoubleVector x2(N+1);
    DoubleVector x3(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x1[i] =  zm0[0][i] + (zm1[0][0][i]*g(x,0,0) + zm1[0][1][i]*g(x,0,1) + zm1[0][2][i]*g(x,0,2)) + (zm2[0][0][i]*g(x,1,0) + zm2[0][1][i]*g(x,1,1) + zm2[0][2][i]*g(x,1,2));
        x2[i] =  zm0[1][i] + (zm1[1][0][i]*g(x,0,0) + zm1[1][1][i]*g(x,0,1) + zm1[1][2][i]*g(x,0,2)) + (zm2[1][0][i]*g(x,1,0) + zm2[1][1][i]*g(x,1,1) + zm2[1][2][i]*g(x,1,2));
        x3[i] =  zm0[2][i] + (zm1[2][0][i]*g(x,0,0) + zm1[2][1][i]*g(x,0,1) + zm1[2][2][i]*g(x,0,2)) + (zm2[2][0][i]*g(x,1,0) + zm2[2][1][i]*g(x,1,1) + zm2[2][2][i]*g(x,1,2));
    }
    IPrinter::printVector(10,6,x1);
    IPrinter::printVector(10,6,x2);
    IPrinter::printVector(10,6,x3);
}

////////////////////////////////////////////////////////////////////////////////////////////////

NonLinearEquationEx2::NonLinearEquationEx2(const Problem4Ex2 &p)
    : NonLinearEquation(), p(p) {}

double NonLinearEquationEx2::fx(const DoubleVector &x, unsigned int num) const
{
    Dimension dim = p.grid.dimension();
    unsigned int N = dim.sizeN();

    const std::vector<DoubleVector> &zm0 = p.zm0;
    const std::vector<std::vector<DoubleVector>> &zm1 = p.zm1;
    const std::vector<std::vector<DoubleVector>> &zm2 = p.zm2;

    if (num == 0)
    {
        unsigned int n = 0.3*N;
        return zm0[0][n]
                + (zm1[0][0][n]*p.g(x,0,0) + zm1[0][1][n]*p.g(x,0,1) + zm1[0][2][n]*p.g(x,0,2))
                + (zm2[0][0][n]*p.g(x,1,0) + zm2[0][1][n]*p.g(x,1,1) + zm2[0][2][n]*p.g(x,1,2)) - x[0];
    }
    if (num == 1)
    {
        unsigned int n = 0.3*N;
        return zm0[1][n]
                + (zm1[1][0][n]*p.g(x,0,0) + zm1[1][1][n]*p.g(x,0,1) + zm1[1][2][n]*p.g(x,0,2))
                + (zm2[1][0][n]*p.g(x,1,0) + zm2[1][1][n]*p.g(x,1,1) + zm2[1][2][n]*p.g(x,1,2)) - x[1];
    }
    if (num == 2)
    {
        unsigned int n = 0.3*N;
        return zm0[2][n]
                + (zm1[2][0][n]*p.g(x,0,0) + zm1[2][1][n]*p.g(x,0,1) + zm1[2][2][n]*p.g(x,0,2))
                + (zm2[2][0][n]*p.g(x,1,0) + zm2[2][1][n]*p.g(x,1,1) + zm2[2][2][n]*p.g(x,1,2)) - x[2];
    }

    if (num == 3)
    {
        unsigned int n = 0.6*N;
        return zm0[0][n]
                + (zm1[0][0][n]*p.g(x,0,0) + zm1[0][1][n]*p.g(x,0,1) + zm1[0][2][n]*p.g(x,0,2))
                + (zm2[0][0][n]*p.g(x,1,0) + zm2[0][1][n]*p.g(x,1,1) + zm2[0][2][n]*p.g(x,1,2)) - x[3];
    }
    if (num == 4)
    {
        unsigned int n = 0.6*N;
        //printf("+++++ %f %f %f %f %f %f %f %f\n", zm0[0][n], p.g(x,0,0), p.g(x,0,1), p.g(x,0,2), p.g(x,1,0), p.g(x,1,1), p.g(x,1,2), x[4]);
        return zm0[1][n]
                + (zm1[1][0][n]*p.g(x,0,0) + zm1[1][1][n]*p.g(x,0,1) + zm1[1][2][n]*p.g(x,0,2))
                + (zm2[1][0][n]*p.g(x,1,0) + zm2[1][1][n]*p.g(x,1,1) + zm2[1][2][n]*p.g(x,1,2)) - x[4];
    }
    if (num == 5)
    {
        unsigned int n = 0.6*N;
        //printf("+++++ %f %f %f %f %f %f %f %f\n", zm0[0][n], p.g(x,0,0), p.g(x,0,1), p.g(x,0,2), p.g(x,1,0), p.g(x,1,1), p.g(x,1,2), x[5]);
        return zm0[2][n]
                + (zm1[2][0][n]*p.g(x,0,0) + zm1[2][1][n]*p.g(x,0,1) + zm1[2][2][n]*p.g(x,0,2))
                + (zm2[2][0][n]*p.g(x,1,0) + zm2[2][1][n]*p.g(x,1,1) + zm2[2][2][n]*p.g(x,1,2)) - x[5];
    }
    //printf("double NonLinearEquationEx2::fx(const DoubleVector &x, unsigned int num) const %d\n", num);
    return NAN;
}

////////////////////////////////////////////////////////////////////////////////////////////////

Problem4Ex2Zetta0::Problem4Ex2Zetta0(const Problem4Ex2 &p)
    : LinearODE1stOrder(), p(p)
{}

unsigned int Problem4Ex2Zetta0::equationsNumber() const
{
    return 3;
}

double Problem4Ex2Zetta0::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p.A(t,k,row,col);
}

double Problem4Ex2Zetta0::B(double t, unsigned int k, unsigned int row) const
{
    return p.B(t,k,row);
}

////////////////////////////////////////////////////////////////////////////////////////////////

Problem4Ex2Zettai::Problem4Ex2Zettai(const Problem4Ex2 &p, unsigned int i)
    : LinearODE1stOrder(), p(p), i(i)
{}

double Problem4Ex2Zettai::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p.A(t, k, row, col);
}

double Problem4Ex2Zettai::B(double t, unsigned int k, unsigned int row) const
{
    return p.C(t, k, i, row, cur_col);
}

double Problem4Ex2Zettai::C(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p.C(t, k, i, row, col);
}

unsigned int Problem4Ex2Zettai::equationsNumber() const
{
    return 3;
}

void Problem4Ex2Zettai::calculateM(const std::vector<LinearODE1stOrder::Condition> &cs, const DoubleMatrix &betta, std::vector< std::vector<DoubleVector> > &zmi)
{
    unsigned int n = equationsNumber();

    zmi.resize(n);
    for (unsigned int row=0; row<n; row++)
        zmi[row].resize(n);

    for (unsigned int col=0; col<n; col++)
    {
        cur_col = col;

        std::vector<DoubleVector> z;
        calculate(cs, betta.col(cur_col), z);
        //highOderAccuracy(cs, betta.col(cur_col), z, 4);

        for (unsigned int row=0; row<n; row++)
        {
            zmi[row][cur_col] = z[row];
            z[row].clear();
        }
        z.clear();
    }
}
