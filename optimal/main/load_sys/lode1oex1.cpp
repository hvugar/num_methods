#include "lode1oex1.h"
#include <math.h>
#include <printer.h>

void LinearODE1stOrderEx1::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    LinearODE1stOrderEx1 cpnlcs;
    cpnlcs.setGrid(UniformODEGrid(0.01, 0, 100));
#ifdef EXAMPLE_1
    cpnlcs.example1();
#endif
#ifdef EXAMPLE_2
    cpnlcs.example2();
#endif
#ifdef EXAMPLE_3
    cpnlcs.example3();
#endif
#ifdef EXAMPLE_4
    cpnlcs.example4();
#endif
#ifdef EXAMPLE_5
    cpnlcs.example5();
#endif
#ifdef EXAMPLE_6
    cpnlcs.example6();
#endif
}

void LinearODE1stOrderEx1::example1()
{
    unsigned int N = grid().sizeN();

    unsigned int n = equationsNumber();

    std::vector<Condition> nscs;
    DoubleVector betta;

    Condition nsc0;
    nsc0.time = 0.0;
    nsc0.mtrx.resize(n, n);
    Random::fillMatrix(nsc0.mtrx, -1, +1, 5);

    Condition nsc1;
    nsc1.time = 0.3;
    nsc1.mtrx.resize(n, n);
    Random::fillMatrix(nsc1.mtrx, -3, +3, 5);

    Condition nsc2;
    nsc2.time = 0.73;
    nsc2.mtrx.resize(n, n);
    Random::fillMatrix(nsc2.mtrx, -1, +1, 5);

    Condition nsc3;
    nsc3.time = 0.8;
    nsc3.mtrx.resize(n, n);
    Random::fillMatrix(nsc3.mtrx, -3, +3, 5);

    Condition nsc4;
    nsc4.time = 1.0;
    nsc4.mtrx.resize(n, n);
    Random::fillMatrix(nsc4.mtrx, -1, +1, 5);

    nscs.push_back(nsc0);
    nscs.push_back(nsc1);
    nscs.push_back(nsc2);
    nscs.push_back(nsc3);
    nscs.push_back(nsc4);

    unsigned int L = nscs.size();

    betta.resize(n);
    for (unsigned int row=0; row<n; row++)
    {
        betta[row] = 0.0;
        for (unsigned int s=0; s<L; s++)
        {
            const Condition &c = nscs.at(s);
            for (unsigned int i=0; i<n; i++) betta[row] += c.mtrx[row][i] * X(c.time, i);
        }
    }

    std::vector<DoubleVector> x2;
    highOderAccuracy(nscs, betta, x2, 2);
    IPrinter::printVector(x2[0]);
    IPrinter::printVector(x2[1]);
    IPrinter::printVector(x2[2]);
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x4;
    highOderAccuracy(nscs, betta, x4, 4);
    IPrinter::printVector(x4[0]);
    IPrinter::printVector(x4[1]);
    IPrinter::printVector(x4[2]);
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x6;
    highOderAccuracy(nscs, betta, x6, 6);
    IPrinter::printVector(x6[0]);
    IPrinter::printVector(x6[1]);
    IPrinter::printVector(x6[2]);
}

void LinearODE1stOrderEx1::example2()
{
    unsigned int N = grid().sizeN();
    double h = grid().step();

    unsigned int n = equationsNumber();

    std::vector<Condition> nscs;
    DoubleVector betta;

    Condition nsc0;
    nsc0.time = 0.0;
    nsc0.mtrx.resize(n, n);
    nsc0.mtrx.at(0,0) = 5.0;

    Condition nsc1;
    nsc1.time = 0.4;
    nsc1.mtrx.resize(n, n);
    nsc1.mtrx.at(0,0) = 4.2;

    Condition nsc2;
    nsc2.time = 1.0;
    nsc2.mtrx.resize(n, n);
    nsc2.mtrx.at(0,0) = 10.0;

    nscs.push_back(nsc0);
    nscs.push_back(nsc1);
    nscs.push_back(nsc2);

    betta.resize(n);
    betta[0] = nsc0.mtrx.at(0,0)*X(nsc0.time,0)
             + nsc1.mtrx.at(0,0)*X(nsc1.time,0)
             + nsc2.mtrx.at(0,0)*X(nsc2.time,0);

    DoubleVector x1(N+1);
    for (unsigned int i=0; i<=N; i++) x1[i] = X(i*h, 0);
    IPrinter::printVector(x1);

    std::vector<DoubleVector> x2;
    highOderAccuracy(nscs, betta, x2, 2);
    IPrinter::printVector(x2[0]);

    std::vector<DoubleVector> x4;
    highOderAccuracy(nscs, betta, x4, 4);
    IPrinter::printVector(x4[0]);

    std::vector<DoubleVector> x6;
    highOderAccuracy(nscs, betta, x6, 6);
    IPrinter::printVector(x6[0]);
}

void LinearODE1stOrderEx1::example3()
{
    unsigned int N = grid().sizeN();
    double h = grid().step();

    unsigned int n = equationsNumber();

    std::vector<Condition> nscs;
    DoubleVector betta;

    Condition nsc0;
    nsc0.time = 0.0;
    nsc0.mtrx.resize(n, n);
    nsc0.mtrx.at(0,0) = 5.0;

    Condition nsc1;
    nsc1.time = 0.4;
    nsc1.mtrx.resize(n, n);
    nsc1.mtrx.at(0,0) = 4.2;

    Condition nsc2;
    nsc2.time = 1.0;
    nsc2.mtrx.resize(n, n);
    nsc2.mtrx.at(0,0) = 10.0;

    nscs.push_back(nsc0);
    nscs.push_back(nsc1);
    nscs.push_back(nsc2);

    betta.resize(n);

    betta[0] = nsc0.mtrx.at(0,0)*X(nsc0.time,0)
             + nsc1.mtrx.at(0,0)*X(nsc1.time,0)
             + nsc2.mtrx.at(0,0)*X(nsc2.time,0);

    DoubleVector x1(N+1);
    for (unsigned int i=0; i<=N; i++) x1[i] = X(i*h, 0);
    IPrinter::printVector(x1);

    std::vector<DoubleVector> x2;
    highOderAccuracy(nscs, betta, x2, 2);
    IPrinter::printVector(x2[0]);

    std::vector<DoubleVector> x4;
    highOderAccuracy(nscs, betta, x4, 4);
    IPrinter::printVector(x4[0]);

    std::vector<DoubleVector> x6;
    highOderAccuracy(nscs, betta, x6, 6);
    IPrinter::printVector(x6[0]);
}

void LinearODE1stOrderEx1::example4()
{
    unsigned int N = grid().sizeN();
    double h = grid().step();

    unsigned int n = equationsNumber();

    std::vector<Condition> nscs;
    DoubleVector betta;

    Condition nsc0;
    nsc0.time = 0.0;
    nsc0.mtrx.resize(n, n);
    Random::fillMatrix(nsc0.mtrx, -1, +1, 5);

    Condition nsc1;
    nsc1.time = 0.3;
    nsc1.nmbr = 30;
    nsc1.mtrx.resize(n, n);
    Random::fillMatrix(nsc1.mtrx, -3, +3, 5);

    Condition nsc2;
    nsc2.time = 0.5;
    nsc2.mtrx.resize(n, n);
    Random::fillMatrix(nsc2.mtrx, -1, +1, 5);

    Condition nsc3;
    nsc3.time = 0.8;
    nsc3.mtrx.resize(n, n);
    Random::fillMatrix(nsc3.mtrx, -3, +3, 5);

    Condition nsc4;
    nsc4.time = 1.0;
    nsc4.mtrx.resize(n, n);
    Random::fillMatrix(nsc4.mtrx, -1, +1, 5);

    nscs.push_back(nsc0);
    nscs.push_back(nsc1);
    nscs.push_back(nsc2);
    nscs.push_back(nsc3);
    nscs.push_back(nsc4);

    unsigned int L = nscs.size();

    betta.resize(n);
    for (unsigned int row=0; row<n; row++)
    {
        betta[row] = 0.0;
        for (unsigned int s=0; s<L; s++)
        {
            const Condition &c = nscs.at(s);
            for (unsigned int i=0; i<n; i++) betta[row] += c.mtrx[row][i] * X(c.time, i);
        }
    }

    std::vector<DoubleVector> x0;
    x0.resize(n);
    x0[0].resize(N+1); for (unsigned int n=0; n<=N; n++) x0[0][n] = X(h*n, 0);
    x0[1].resize(N+1); for (unsigned int n=0; n<=N; n++) x0[1][n] = X(h*n, 1);
    x0[2].resize(N+1); for (unsigned int n=0; n<=N; n++) x0[2][n] = X(h*n, 2);
    IPrinter::printVector(x0[0]);
    IPrinter::printVector(x0[1]);
    IPrinter::printVector(x0[2]);
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x2;
    highOderAccuracy(nscs, betta, x2, 2);
    IPrinter::printVector(x2[0]);
    IPrinter::printVector(x2[1]);
    IPrinter::printVector(x2[2]);
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x4;
    highOderAccuracy(nscs, betta, x4, 4);
    IPrinter::printVector(x4[0]);
    IPrinter::printVector(x4[1]);
    IPrinter::printVector(x4[2]);
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x6;
    highOderAccuracy(nscs, betta, x6, 6);
    IPrinter::printVector(x6[0]);
    IPrinter::printVector(x6[1]);
    IPrinter::printVector(x6[2]);
}

void LinearODE1stOrderEx1::example5()
{
    unsigned int N = grid().sizeN();
    double h = grid().step();

    unsigned int n = equationsNumber();

    std::vector<Condition> nscs;
    DoubleVector betta;

    Condition nsc0;
    nsc0.time = 0.0;
    nsc0.mtrx.resize(n, n);
    Random::fillMatrix(nsc0.mtrx, -1, +1, 5);

    Condition nsc1;
    nsc1.time = 0.3;
    nsc1.nmbr = 30;
    nsc1.mtrx.resize(n, n);
    Random::fillMatrix(nsc1.mtrx, -3, +3, 5);

    Condition nsc2;
    nsc2.time = 0.5;
    nsc2.mtrx.resize(n, n);
    Random::fillMatrix(nsc2.mtrx, -1, +1, 5);

    Condition nsc3;
    nsc3.time = 0.8;
    nsc3.mtrx.resize(n, n);
    Random::fillMatrix(nsc3.mtrx, -3, +3, 5);

    Condition nsc4;
    nsc4.time = 1.0;
    nsc4.mtrx.resize(n, n);
    Random::fillMatrix(nsc4.mtrx, -1, +1, 5);

    nscs.push_back(nsc0);
    nscs.push_back(nsc1);
    nscs.push_back(nsc2);
    nscs.push_back(nsc3);
    nscs.push_back(nsc4);

    unsigned int L = nscs.size();

    betta.resize(n);
    for (unsigned int row=0; row<n; row++)
    {
        betta[row] = 0.0;
        for (unsigned int s=0; s<L; s++)
        {
            const Condition &c = nscs.at(s);
            for (unsigned int i=0; i<n; i++) betta[row] += c.mtrx[row][i] * X(c.time, i);
        }
    }

    std::vector<DoubleVector> x0;
    x0.resize(n);
    x0[0].resize(N+1); for (unsigned int n=0; n<=N; n++) x0[0][n] = X(h*n, 0);
    x0[1].resize(N+1); for (unsigned int n=0; n<=N; n++) x0[1][n] = X(h*n, 1);
    x0[2].resize(N+1); for (unsigned int n=0; n<=N; n++) x0[2][n] = X(h*n, 2);
    IPrinter::printVector(x0[0]);
    IPrinter::printVector(x0[1]);
    IPrinter::printVector(x0[2]);
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x2;
    highOderAccuracy(nscs, betta, x2, 2);
    IPrinter::printVector(x2[0]);
    IPrinter::printVector(x2[1]);
    IPrinter::printVector(x2[2]);
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x4;
    highOderAccuracy(nscs, betta, x4, 4);
    IPrinter::printVector(x4[0]);
    IPrinter::printVector(x4[1]);
    IPrinter::printVector(x4[2]);
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x6;
    highOderAccuracy(nscs, betta, x6, 6);
    IPrinter::printVector(x6[0]);
    IPrinter::printVector(x6[1]);
    IPrinter::printVector(x6[2]);
}

void LinearODE1stOrderEx1::example6()
{
    unsigned int N = grid().sizeN();
    double h = grid().step();

    DoubleVector x(N+1); x[0] = 0.0;
    for (unsigned int i=0; i<N; i++) x[i+1] = (1.0+h*A(i*h, i))*x[i] + h*B(i*h, i);
    IPrinter::printVector(x);

    DoubleVector x2(N+1); x2[0] = 0.0;
    for (unsigned int i=0; i<N; i++)
    {
        double m = (4.0-2.0*h*A((i+1)*h, i+1));
        double k1 = (4.0+2.0*h*A(i*h, i)) / m;
        double k2 = (2.0*h*(B(i*h, i)+B((i+1)*h, i+1))) / m;
        x2[i+1] = k1*x2[i] + k2;
    }
    IPrinter::printVector(x2);
}

double LinearODE1stOrderEx1::A(double t UNUSED_PARAM, unsigned int, unsigned int row UNUSED_PARAM, unsigned int col UNUSED_PARAM) const
{
#ifdef EXAMPLE_1
    if (row==0) { if (col==0) { return +2.0; } if (col==1) { return -3.0; } if (col==2) { return +1.0; } }
    if (row==1) { if (col==0) { return +3.0; } if (col==1) { return +1.0; } if (col==2) { return -2.0; } }
    if (row==2) { if (col==0) { return +1.0; } if (col==1) { return -5.0; } if (col==2) { return -3.0; } }
#endif
#ifdef EXAMPLE_2
    return t;
#endif
#ifdef EXAMPLE_3
    return t;
#endif
#ifdef EXAMPLE_4
    if (row==0) { if (col==0) { return +2.0; } if (col==1) { return -3.0; } if (col==2) { return +1.0; } }
    if (row==1) { if (col==0) { return +3.0; } if (col==1) { return +1.0; } if (col==2) { return -2.0; } }
    if (row==2) { if (col==0) { return +1.0; } if (col==1) { return -5.0; } if (col==2) { return -3.0; } }
#endif
#ifdef EXAMPLE_5
    if (row==0) { if (col==0) { return +2.0; } if (col==1) { return -3.0; } if (col==2) { return +1.0; } }
    if (row==1) { if (col==0) { return +3.0; } if (col==1) { return +1.0; } if (col==2) { return -2.0; } }
    if (row==2) { if (col==0) { return +1.0; } if (col==1) { return -5.0; } if (col==2) { return -3.0; } }
#endif
#ifdef EXAMPLE_6
    return 3.0;
#endif

    return NAN;
}

double LinearODE1stOrderEx1::B(double t, unsigned int, unsigned int row UNUSED_PARAM) const
{
#ifdef EXAMPLE_1
    if (row==0) return +11.0*t*t - 7.0*t - 5.0;
    if (row==1) return -2.0*t*t + t - 12.0;
    if (row==2) return +23.0*t*t + 2.0*t - 3.0;
#endif
#ifdef EXAMPLE_2
    return -t*t*t;
#endif
#ifdef EXAMPLE_3
    return 4.0*cos(4.0*t) - t*sin(4.0*t);
#endif
#ifdef EXAMPLE_4
    if (row==0) return +4.0*cos(4.0*t) - 2.0*sin(4.0*t) + 3.0*cos(5.0*t) - 2.0 + exp(t);
    if (row==1) return -5.0*sin(5.0*t) - 3.0*sin(4.0*t) - cos(5.0*t)     + 4.0 - 2.0*exp(t);
    if (row==2) return -sin(4.0*t) + 5.0*cos(5.0*t) + 6.0 - 4.0*exp(t);
#endif
#ifdef EXAMPLE_5
    if (row==0) return +40.0*cos(40.0*t) - 2.0*sin(40.0*t) + 3.0*cos(50.0*t) - 2.0 +     exp(t);
    if (row==1) return -50.0*sin(50.0*t) - 3.0*sin(40.0*t) -     cos(50.0*t) + 4.0 - 2.0*exp(t);
    if (row==2) return                   -     sin(40.0*t) + 5.0*cos(50.0*t) + 6.0 - 4.0*exp(t);
#endif
#ifdef EXAMPLE_6
    return 2.0*t - 3.0*t*t;
#endif
    return NAN;
}

unsigned int LinearODE1stOrderEx1::equationsNumber() const
{
#ifdef EXAMPLE_1
    return 3;
#endif
#ifdef EXAMPLE_2
    return 1;
#endif
#ifdef EXAMPLE_3
    return 1;
#endif
#ifdef EXAMPLE_4
    return 3;
#endif
#ifdef EXAMPLE_5
    return 3;
#endif
#ifdef EXAMPLE_6
    return 1;
#endif
}

double LinearODE1stOrderEx1::X(double t, int row UNUSED_PARAM) const
{
#ifdef EXAMPLE_1
    if (row==0) return 3.0*t+4.0;
    if (row==1) return 4.0*t*t;
    if (row==2) return t*t+t;
#endif
#ifdef EXAMPLE_2
    return t*t+2.0;
#endif
#ifdef EXAMPLE_3
    return sin(4.0*t);
#endif
#ifdef EXAMPLE_4
    if (row==0) return sin(4.0*t);
    if (row==1) return cos(5.0*t);
    if (row==2) return 2.0-exp(t);
#endif
#ifdef EXAMPLE_5
    if (row==0) return sin(40.0*t);
    if (row==1) return cos(50.0*t);
    if (row==2) return 2.0-exp(t);
#endif
#ifdef EXAMPLE_6
    return t*t;
#endif
    return NAN;
}
