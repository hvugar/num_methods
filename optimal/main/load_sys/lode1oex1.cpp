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
    Dimension dim = grid().dimension();
    unsigned int N = dim.sizeN();

    unsigned int n = equationsNumber();

    std::vector<Condition> nscs;
    DoubleVector betta;

    Condition nsc0;
    nsc0.time = 0.0;
    nsc0.nmbr = 0;
    nsc0.mtrx.resize(n, n);
    Random::fillMatrix(nsc0.mtrx, -1, +1, 5);

    Condition nsc1;
    nsc1.time = 0.3;
    nsc1.nmbr = 30;
    nsc1.mtrx.resize(n, n);
    Random::fillMatrix(nsc1.mtrx, -3, +3, 5);

    Condition nsc2;
    nsc2.time = 0.73;
    nsc2.nmbr = 70;
    nsc2.mtrx.resize(n, n);
    Random::fillMatrix(nsc2.mtrx, -1, +1, 5);

    Condition nsc3;
    nsc3.time = 0.8;
    nsc3.nmbr = 80;
    nsc3.mtrx.resize(n, n);
    Random::fillMatrix(nsc3.mtrx, -3, +3, 5);

    Condition nsc4;
    nsc4.time = 1.0;
    nsc4.nmbr = 100;
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
    unsigned int N = grid().dimension().sizeN();
    double h = grid().dimension().step();

//    printf("%14.10f %14.10f\n", -3.0*X(h*0) - 10.0*X(h*1) + 18.0*X(h*2) - 6.0*X(h*3) + X(h*4),
//                                12.0*h*(A(h*1,1)*X(h*0)+B(h*1,1)));

//    double m = 25.0-12.0*h*A(h*4,4);
//    printf("%14.10f %14.10f\n", (3.0+3.0*m)*X(h*0)
//                               - (16.0-(10.0+12.0*h*A(h,1))*m)*X(h*1)
//                               + (36.0-18.0*m)*X(h*2)
//                               - (48.0-6.0*m)*X(h*3),
//                                              12.0*h*B(h*4,4)
//                                             -12.0*h*B(h*1,1)*m);

//    double m = 7.0-3.0*h*A(h*3,3);
//    printf("%14.10f %14.10f\n",  (3.0 - h*A(3*h,3))*X(h*0)
//                               + (3.0+3.0*h*A(1*h,1)+m*(3.0+2.0*h*A(1*h,1)))*X(h*1)
//                               - (6.0 +m*(3.0-2.0*h*A(2*h,2)))*X(h*2),
//                                              h*B(h*3,3)
//                                             -3.0*h*B(h*1,1)
//                                             -m*(2.0*h*B(2*h,2)+2.0*h*B(1*h,1)));

    DoubleVector x(N+1);
    x[0] = 2.0;

    for (unsigned int i=0; i<N; i++)
    {
        double A0 = A((i+0)*h, i+0);
        double A1 = A((i+1)*h, i+1);
        double A2 = A((i+2)*h, i+2);
        double A3 = A((i+3)*h, i+3);
        double B0 = B((i+0)*h, i+0);
        double B1 = B((i+1)*h, i+1);
        double B2 = B((i+2)*h, i+2);
        double B3 = B((i+3)*h, i+3);

        double k1 = ((3.0-h*A3)*(3.0-h*A2)-(6.0+(7.0-3.0*h*A3)*(3.0-2.0*h*A2))*(3.0+h*A0));
        double k2 = ((3.0+3.0*h*A1+(7.0-3.0*h*A3)*(3.0+2.0*h*A1))*(3.0-h*A2)-(6.0+(7.0-3.0*h*A3)*(3.0-2.0*h*A2))*4.0*h*A1);
        double k3 = (3.0-h*A2)*(h*B3-3.0*h*B1-(7.0-3.0*h*A3)*(2.0*h*B2+2.0*h*B1))+(6.0+(7.0-3.0*h*A3)*(3.0-2.0*h*A2))*(h*B0+4.0*h*B1+h*B2);

        x[i+1] = -k1/k2*x[i] + k3/k2;

    }
    IPrinter::printVector(x);

    //printf("%14.10f %14.10f\n", k1*X(0*h) + k2*X(1*h), k3);

    return;


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

    DoubleVector x1(N+1);
    for (unsigned int i=0; i<=N; i++) x1[i] = X(i*h, 0);
    IPrinter::printVector(x1);

    std::vector<DoubleVector> x2;
    highOderAccuracy(cs, bt, x2, 2);
    IPrinter::printVector(x2[0]);

    std::vector<DoubleVector> x4;
    highOderAccuracy(cs, bt, x4, 4);
    IPrinter::printVector(x4[0]);

    std::vector<DoubleVector> x6;
    highOderAccuracy(cs, bt, x6, 6);
    IPrinter::printVector(x6[0]);
}

void LinearODE1stOrderEx1::example3()
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
    highOderAccuracy(cs, bt, x2, 2);
    IPrinter::printVector(x2[0]);

    std::vector<DoubleVector> x4;
    highOderAccuracy(cs, bt, x4, 4);
    IPrinter::printVector(x4[0]);

    std::vector<DoubleVector> x6;
    highOderAccuracy(cs, bt, x6, 6);
    IPrinter::printVector(x6[0]);
}

void LinearODE1stOrderEx1::example4()
{
    Dimension dim = grid().dimension();
    unsigned int N = dim.sizeN();
    double h = dim.step();

    unsigned int n = equationsNumber();

    std::vector<Condition> nscs;
    DoubleVector betta;

    Condition nsc0;
    nsc0.time = 0.0;
    nsc0.nmbr = 0;
    nsc0.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc0.mtrx[row][col] = (rand() % 1000) / 1000.0;

//    Condition nsc1;
//    nsc1.time = 0.3;
//    nsc1.nmbr = 30;
//    nsc1.mtrx.resize(n, n);
//    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc1.mtrx[row][col] = (rand() % 1000) / 1000.0;

//    Condition nsc2;
//    nsc2.time = 0.5;
//    nsc2.nmbr = 50;
//    nsc2.mtrx.resize(n, n);
//    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc2.mtrx[row][col] = (rand() % 1000) / 1000.0;

//    Condition nsc3;
//    nsc3.time = 0.8;
//    nsc3.nmbr = 80;
//    nsc3.mtrx.resize(n, n);
//    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc3.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc4;
    nsc4.time = 1.0;
    nsc4.nmbr = N;
    nsc4.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc4.mtrx[row][col] = (rand() % 1000) / 1000.0;

    nscs.push_back(nsc0);
//    nscs.push_back(nsc1);
//    nscs.push_back(nsc2);
//    nscs.push_back(nsc3);
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
    x0.resize(3);
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
    example4();
}

void LinearODE1stOrderEx1::example6()
{
    Dimension dim = grid().dimension();
    unsigned int N = dim.sizeN();
    double h = dim.step();

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
