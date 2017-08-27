#include "lode1oex1.h"
#include <math.h>
#include <printer.h>

void LinearODE1stOrderEx1::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    LinearODE1stOrderEx1 cpnlcs;
    cpnlcs.setGrid(ODEGrid(Dimension(0.01, 100, 0)));
#ifdef EXAMPLE_1
    cpnlcs.example1();
#endif
#ifdef EXAMPLE_2
    cpnlcs.example2();
#endif
}

void LinearODE1stOrderEx1::example1()
{
    Dimension dim = grid().dimension();
    unsigned int N = dim.sizeN();

    unsigned int n = equationsNumber();

    std::vector<Condition> nscs;
    DoubleVector betta;
    std::vector<DoubleVector> x;

    Condition nsc0;
    nsc0.time = 0.0;
    nsc0.nmbr = 0;
    nsc0.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc0.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc1;
    nsc1.time = 0.25;
    nsc1.nmbr = N/4;
    nsc1.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc1.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc2;
    nsc2.time = 0.5;
    nsc2.nmbr = N/2;
    nsc2.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc2.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc3;
    nsc3.time = 0.75;
    nsc3.nmbr = 3*(N/4);
    nsc3.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc3.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc4;
    nsc4.time = 1.0;
    nsc4.nmbr = N;
    nsc4.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc4.mtrx[row][col] = (rand() % 1000) / 1000.0;

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
    highOder2Accuracy(nscs, betta, x2);
//    IPrinter::printVector(x2[0]);
//    IPrinter::printVector(x2[1]);
//    IPrinter::printVector(x2[2]);
}

void LinearODE1stOrderEx1::example2()
{
    unsigned int n = equationsNumber();
    std::vector<Condition> cs;

    Condition c0;
    c0.time = 0.0;
    c0.mtrx.resize(n, n);
    c0.mtrx.at(0,0) = 5.0;
    c0.index = 0;
    cs.push_back(c0);

    Condition c1;
    c1.time = 0.4;
    c1.mtrx.resize(n, n);
    c1.mtrx.at(0,0) = 4.2;
    c1.index = 0;
    cs.push_back(c1);

    Condition c2;
    c2.time = 1.0;
    c2.mtrx.resize(n, n);
    c2.mtrx.at(0,0) = 10.0;
    c2.index = 0;
    cs.push_back(c2);

    DoubleVector bt(n);

    bt[0] = c0.mtrx.at(0,0)*X(c0.time,0)
            + c1.mtrx.at(0,0)*X(c1.time,0)
            + c2.mtrx.at(0,0)*X(c2.time,0);

    unsigned int N = grid().dimension().sizeN();
    double h = grid().dimension().step();
    DoubleVector x1(N+1);
    for (unsigned int i=0; i<=N; i++) x1[i] = X(i*h, 0);
    IPrinter::printVector(x1);

    std::vector<DoubleVector> x2;
    highOder2Accuracy(cs, bt, x2);
    IPrinter::printVector(x2[0]);

    std::vector<DoubleVector> x4;
    highOder4Accuracy(cs, bt, x4);
    IPrinter::printVector(x4[0]);

    std::vector<DoubleVector> x6;
    highOder6Accuracy(cs, bt, x6);
    IPrinter::printVector(x6[0]);
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
    return NAN;
}

double LinearODE1stOrderEx1::B(double t, unsigned int, unsigned int row UNUSED_PARAM) const
{
#ifdef EXAMPLE_1
    if (row==0) return 3.0     - (2.0*(3.0*t+4.0) - 3.0*(4.0*t*t) + 1.0*(t*t+t));
    if (row==1) return 8.0*t   - (3.0*(3.0*t+4.0) + 1.0*(4.0*t*t) - 2.0*(t*t+t));
    if (row==2) return 2.0*t+1 - (1.0*(3.0*t+4.0) - 5.0*(4.0*t*t) - 3.0*(t*t+t));
#endif
#ifdef EXAMPLE_2
    //return 1.0 - t*t - t;
    return -t*t*t;
    //return 10.0*cos(10.0*t)-t*sin(10.0*t);
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
}

double LinearODE1stOrderEx1::X(double t, int row UNUSED_PARAM) const
{
#ifdef EXAMPLE_1
    if (row==0) return 3.0*t+4.0;
    if (row==1) return 4.0*t*t;
    if (row==2) return t*t+t;
#endif
#ifdef EXAMPLE_2
    //return t+1.0;
    return t*t+2.0;
    //return sin(10.0*t);
#endif
    return NAN;
}
