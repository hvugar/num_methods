#include "loadedsystems.h"
#include <math.h>

struct CauchyProblem : public CauchyProblem2
{
    CauchyProblem(const LoadedSystems &ls) : ls(ls) {}
    virtual ~CauchyProblem() {}

    virtual double f(double t UNUSED_PARAM, const DoubleVector &x) const
    {
        DoubleVector alpha = x.mid(0, ls.n-1);
        std::vector<DoubleVector> betta(ls.k1);
        for (unsigned int s1=0; s1<ls.k1; s1++)
        {
            unsigned int a = ls.n + s*ls.n;
            unsigned int b = a+ls.n-1;
            betta[s1] = x.mid(a, b);
        }
        double q = x.at(ls.n + ls.k1*ls.n);
        double M = x.at(ls.n + ls.k1*ls.n+1);

        // alpha
        if (type == 0)
        {
            double sum = 0.0;
            for (unsigned int row=0; row<ls.n; row++) sum += alpha.at(row) * ls.A(row, i, t);
            return ls.S(x, t) * x.at(i) - sum;
        }

        // betta
        if (type == 1)
        {
            double sum = 0.0;
            for (unsigned int row=0; row<ls.n; row++) sum +=  alpha.at(row) * ls.B(row, i, s, t);
            return ls.S(x, t) * betta[s].at(i) - sum;
        }

        // qamma
        if (type == 2)
        {
            double sum = 0.0;
            for (unsigned int row=0; row<ls.n; row++) sum += ls.C(row, t) * alpha.at(row);
            return ls.S(x, t)*q + sum;
        }

        //M
        if (type == 3)
        {
            return ls.S(x, t)*M;
        }
        return 0.0;
    }
    unsigned int type;
    unsigned int i;
    unsigned int s;
    const LoadedSystems &ls;
};

LoadedSystems::LoadedSystems()
{
    init();

    DoubleMatrix u;
    DoubleVector b;

    u.resize(15, 15);
    b.resize(15);

    ////////////////////////////////////////////////////////////
    puts("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    {
        // j = 1
        DoubleVector alpha1(n);
        alpha1[0] = 0.0; alpha1[1] = 1.0; alpha1[2] = 0.0;
        std::vector<DoubleVector> bettas(2);
        bettas[0].resize(n);
        bettas[1].resize(n);
        bettas[0][0] = 0.0; bettas[0][1] = 0.0; bettas[0][2] = 0.0;
        bettas[1][0] = 0.0; bettas[1][1] = 0.0; bettas[1][2] = 0.0;
        double qmm = +3.68294;

        DoubleVector alpha2(n); alpha2[0] = 0.0; alpha2[1] = 0.0; alpha2[2] = 0.0;
        DoubleVector alpha3(n); alpha3[0] = 2.0; alpha3[1] = 0.0; alpha3[2] = 0.0;

        N = 120;
        DoubleMatrix m;
        calculate(m, N, alpha1, bettas, qmm, n);

        DoubleVector a221(N+1); for (unsigned int i=0; i<=N; i++) a221.at(i) = alpha2[0] * m.row(10).at(i);
        DoubleVector a222(N+1); for (unsigned int i=0; i<=N; i++) a222.at(i) = alpha2[1] * m.row(10).at(i);
        DoubleVector a223(N+1); for (unsigned int i=0; i<=N; i++) a223.at(i) = alpha2[2] * m.row(10).at(i);

        DoubleVector a321(N+1); for (unsigned int i=0; i<=N; i++) a321.at(i) = alpha3[0] * m.row(10).at(i);
        DoubleVector a322(N+1); for (unsigned int i=0; i<=N; i++) a322.at(i) = alpha3[1] * m.row(10).at(i);
        DoubleVector a323(N+1); for (unsigned int i=0; i<=N; i++) a323.at(i) = alpha3[2] * m.row(10).at(i);

        unsigned int N1 = 8;
        IPrinter::printVector(m.row(0),NULL,N1);
        IPrinter::printVector(m.row(1),NULL,N1);
        IPrinter::printVector(m.row(2),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m.row(3),NULL,N1);
        IPrinter::printVector(m.row(4),NULL,N1);
        IPrinter::printVector(m.row(5),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m.row(6),NULL,N1);
        IPrinter::printVector(m.row(7),NULL,N1);
        IPrinter::printVector(m.row(8),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m.row(9),NULL,N1);
        IPrinter::printVector(m.row(10),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(a221,NULL,N1);
        IPrinter::printVector(a222,NULL,N1);
        IPrinter::printVector(a223,NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(a321,NULL,N1);
        IPrinter::printVector(a322,NULL,N1);
        IPrinter::printVector(a323,NULL,N1);

        // t=0.00
        u.at(5,0)  = m.row(0).at(0); u.at(5,1)  = m.row(1).at(0); u.at(5,2)  = m.row(2).at(0);
        u.at(5,3)  = a221.at(0);     u.at(5,4)  = a222.at(0);     u.at(5,5)  = a223.at(0);
        u.at(5,6)  = a321.at(0);     u.at(5,7)  = a322.at(0);     u.at(5,8)  = a323.at(0);
        u.at(5,9)  = m.row(3).at(0); u.at(5,10) = m.row(4).at(0); u.at(5,11) = m.row(5).at(0);
        u.at(5,12) = m.row(6).at(0); u.at(5,13) = m.row(7).at(0); u.at(5,14) = m.row(8).at(0);
        b.at(5)    = m.row(9).at(0);

        // t=0.25
        u.at(6,0)  = m.row(0).at(30); u.at(6,1)  = m.row(1).at(30); u.at(6,2)  = m.row(2).at(30);
        u.at(6,3)  = a221.at(30);     u.at(6,4)  = a222.at(30);     u.at(6,5)  = a223.at(30);
        u.at(6,6)  = a321.at(30);     u.at(6,7)  = a322.at(30);     u.at(6,8)  = a323.at(30);
        u.at(6,9)  = m.row(3).at(30); u.at(6,10) = m.row(4).at(30); u.at(6,11) = m.row(5).at(30);
        u.at(6,12) = m.row(6).at(30); u.at(6,13) = m.row(7).at(30); u.at(6,14) = m.row(8).at(30);
        b.at(6)    = m.row(9).at(30);

        // t=0.50
        u.at(7,0)  = m.row(0).at(60); u.at(7,1)  = m.row(1).at(60); u.at(7,2)  = m.row(2).at(60);
        u.at(7,3)  = a221.at(60);     u.at(7,4)  = a222.at(60);     u.at(7,5)  = a223.at(60);
        u.at(7,6)  = a321.at(60);     u.at(7,7)  = a322.at(60);     u.at(7,8)  = a323.at(60);
        u.at(7,9)  = m.row(3).at(60); u.at(7,10) = m.row(4).at(60); u.at(7,11) = m.row(5).at(60);
        u.at(7,12) = m.row(6).at(60); u.at(7,13) = m.row(7).at(60); u.at(7,14) = m.row(8).at(60);
        b.at(7)    = m.row(9).at(60);

        // t=0.75
        u.at(8,0)  = m.row(0).at(90); u.at(8,1)  = m.row(1).at(90); u.at(8,2)  = m.row(2).at(90);
        u.at(8,3)  = a221.at(90);     u.at(8,4)  = a222.at(90);     u.at(8,5)  = a223.at(90);
        u.at(8,6)  = a321.at(90);     u.at(8,7)  = a322.at(90);     u.at(8,8)  = a323.at(90);
        u.at(8,9)  = m.row(3).at(90); u.at(8,10) = m.row(4).at(90); u.at(8,11) = m.row(5).at(90);
        u.at(8,12) = m.row(6).at(90); u.at(8,13) = m.row(7).at(90); u.at(8,14) = m.row(8).at(90);
        b.at(8)    = m.row(9).at(90);

        // t=1.00
        u.at(9,0)  = m.row(0).at(120); u.at(9,1)  = m.row(1).at(120); u.at(9,2)  = m.row(2).at(120);
        u.at(9,3)  = a221.at(120);     u.at(9,4)  = a222.at(120);     u.at(9,5)  = a223.at(120);
        u.at(9,6)  = a321.at(120);     u.at(9,7)  = a322.at(120);     u.at(9,8)  = a323.at(120);
        u.at(9,9)  = m.row(3).at(120); u.at(9,10) = m.row(4).at(120); u.at(9,11) = m.row(5).at(120);
        u.at(9,12) = m.row(6).at(120); u.at(9,13) = m.row(7).at(120); u.at(9,14) = m.row(8).at(120);
        b.at(9)    = m.row(9).at(120);
    }

    puts("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    {
        // j = 2
        DoubleVector alpha1(n);
        alpha1[0] = +0.0; alpha1[1] = +0.0; alpha1[2] = +3.0;
        std::vector<DoubleVector> bettas(2);
        bettas[0].resize(n);
        bettas[1].resize(n);
        bettas[0][0] = +0.0; bettas[0][1] = +0.0; bettas[0][2] = +0.0;
        bettas[1][0] = +0.0; bettas[1][1] = +0.0; bettas[1][2] = +0.0;
        double qmm = +3.0806;

        DoubleVector alpha2(n); alpha2[0] = +0.0; alpha2[1] = +0.0; alpha2[2] = +0.0;
        DoubleVector alpha3(n); alpha3[0] = +0.0; alpha3[1] = -1.0; alpha3[2] = +0.0;

        N = 120;
        DoubleMatrix m;
        calculate(m, N, alpha1, bettas, qmm, n);

        DoubleVector a231(N+1); for (unsigned int i=0; i<=N; i++) a231.at(i) = alpha2[0] * m.row(10).at(i);
        DoubleVector a232(N+1); for (unsigned int i=0; i<=N; i++) a232.at(i) = alpha2[1] * m.row(10).at(i);
        DoubleVector a233(N+1); for (unsigned int i=0; i<=N; i++) a233.at(i) = alpha2[2] * m.row(10).at(i);

        DoubleVector a331(N+1); for (unsigned int i=0; i<=N; i++) a331.at(i) = alpha3[0] * m.row(10).at(i);
        DoubleVector a332(N+1); for (unsigned int i=0; i<=N; i++) a332.at(i) = alpha3[1] * m.row(10).at(i);
        DoubleVector a333(N+1); for (unsigned int i=0; i<=N; i++) a333.at(i) = alpha3[2] * m.row(10).at(i);

        unsigned int N1 = 8;
        IPrinter::printVector(m.row(0),NULL,N1);
        IPrinter::printVector(m.row(1),NULL,N1);
        IPrinter::printVector(m.row(2),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m.row(3),NULL,N1);
        IPrinter::printVector(m.row(4),NULL,N1);
        IPrinter::printVector(m.row(5),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m.row(6),NULL,N1);
        IPrinter::printVector(m.row(7),NULL,N1);
        IPrinter::printVector(m.row(8),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m.row(9),NULL,N1);
        IPrinter::printVector(m.row(10),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(a231,NULL,N1);
        IPrinter::printVector(a232,NULL,N1);
        IPrinter::printVector(a233,NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(a331,NULL,N1);
        IPrinter::printVector(a332,NULL,N1);
        IPrinter::printVector(a333,NULL,N1);

        // t=0.00
        u.at(10,0)  = m.row(0).at(0); u.at(10,1)  = m.row(1).at(0); u.at(10,2)  = m.row(2).at(0);
        u.at(10,3)  = a231.at(0);     u.at(10,4)  = a232.at(0);     u.at(10,5)  = a233.at(0);
        u.at(10,6)  = a331.at(0);     u.at(10,7)  = a332.at(0);     u.at(10,8)  = a333.at(0);
        u.at(10,9)  = m.row(3).at(0); u.at(10,10) = m.row(4).at(0); u.at(10,11) = m.row(5).at(0);
        u.at(10,12) = m.row(6).at(0); u.at(10,13) = m.row(7).at(0); u.at(10,14) = m.row(8).at(0);
        b.at(10)    = m.row(9).at(0);

        // t=0.25
        u.at(11,0)  = m.row(0).at(30); u.at(11,1)  = m.row(1).at(30); u.at(11,2)  = m.row(2).at(30);
        u.at(11,3)  = a231.at(30);     u.at(11,4)  = a232.at(30);     u.at(11,5)  = a233.at(30);
        u.at(11,6)  = a331.at(30);     u.at(11,7)  = a332.at(30);     u.at(11,8)  = a333.at(30);
        u.at(11,9)  = m.row(3).at(30); u.at(11,10) = m.row(4).at(30); u.at(11,11) = m.row(5).at(30);
        u.at(11,12) = m.row(6).at(30); u.at(11,13) = m.row(7).at(30); u.at(11,14) = m.row(8).at(30);
        b.at(11)    = m.row(9).at(30);

        // t=0.50
        u.at(12,0)  = m.row(0).at(60); u.at(12,1)  = m.row(1).at(60); u.at(12,2)  = m.row(2).at(60);
        u.at(12,3)  = a231.at(60);     u.at(12,4)  = a232.at(60);     u.at(12,5)  = a233.at(60);
        u.at(12,6)  = a331.at(60);     u.at(12,7)  = a332.at(60);     u.at(12,8)  = a333.at(60);
        u.at(12,9)  = m.row(3).at(60); u.at(12,10) = m.row(4).at(60); u.at(12,11) = m.row(5).at(60);
        u.at(12,12) = m.row(6).at(60); u.at(12,13) = m.row(7).at(60); u.at(12,14) = m.row(8).at(60);
        b.at(12)    = m.row(9).at(60);

        // t=0.75
        u.at(13,0)  = m.row(0).at(90); u.at(13,1)  = m.row(1).at(90); u.at(13,2)  = m.row(2).at(90);
        u.at(13,3)  = a231.at(90);     u.at(13,4)  = a232.at(90);     u.at(13,5)  = a233.at(90);
        u.at(13,6)  = a331.at(90);     u.at(13,7)  = a332.at(90);     u.at(13,8)  = a333.at(90);
        u.at(13,9)  = m.row(3).at(90); u.at(13,10) = m.row(4).at(90); u.at(13,11) = m.row(5).at(90);
        u.at(13,12) = m.row(6).at(90); u.at(13,13) = m.row(7).at(90); u.at(13,14) = m.row(8).at(90);
        b.at(13)    = m.row(9).at(90);

        // t=1.00
        u.at(14,0)  = m.row(0).at(120); u.at(14,1)  = m.row(1).at(120); u.at(14,2)  = m.row(2).at(120);
        u.at(14,3)  = a231.at(120);     u.at(14,4)  = a232.at(120);     u.at(14,5)  = a233.at(120);
        u.at(14,6)  = a331.at(120);     u.at(14,7)  = a332.at(120);     u.at(14,8)  = a333.at(120);
        u.at(14,9)  = m.row(3).at(120); u.at(14,10) = m.row(4).at(120); u.at(14,11) = m.row(5).at(120);
        u.at(14,12) = m.row(6).at(120); u.at(14,13) = m.row(7).at(120); u.at(14,14) = m.row(8).at(120);
        b.at(14)    = m.row(9).at(120);
    }

    puts("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    {
        // j = 0
        DoubleVector alpha1(n);
        alpha1[0] = +2.0; alpha1[1] = +0.0; alpha1[2] = +0.0;
        std::vector<DoubleVector> bettas(2);
        bettas[0].resize(n);
        bettas[1].resize(n);
        bettas[0][0] = +0.0; bettas[0][1] = +0.0; bettas[0][2] = +0.0;
        bettas[1][0] = +0.0; bettas[1][1] = +0.0; bettas[1][2] = +0.0;
        double qmm = -1.3569;

        DoubleVector alpha2(n); alpha2[0] = +0.0; alpha2[1] = +3.0; alpha2[2] = +0.0;
        DoubleVector alpha3(n); alpha3[0] = +0.0; alpha3[1] = +0.0; alpha3[2] = +1.0;

        N = 60;
        DoubleMatrix m1;
        calculate(m1, N, alpha1, bettas, qmm, n);

        DoubleVector a211(N+1); for (unsigned int i=0; i<=N; i++) a211.at(i) = alpha2[0] * m1.row(10).at(i);
        DoubleVector a212(N+1); for (unsigned int i=0; i<=N; i++) a212.at(i) = alpha2[1] * m1.row(10).at(i);
        DoubleVector a213(N+1); for (unsigned int i=0; i<=N; i++) a213.at(i) = alpha2[2] * m1.row(10).at(i);

        DoubleVector a311(N+1); for (unsigned int i=0; i<=N; i++) a311.at(i) = alpha3[0] * m1.row(10).at(i);
        DoubleVector a312(N+1); for (unsigned int i=0; i<=N; i++) a312.at(i) = alpha3[1] * m1.row(10).at(i);
        DoubleVector a313(N+1); for (unsigned int i=0; i<=N; i++) a313.at(i) = alpha3[2] * m1.row(10).at(i);

        unsigned int N1 = 4;
        IPrinter::printVector(m1.row(0),NULL,N1);
        IPrinter::printVector(m1.row(1),NULL,N1);
        IPrinter::printVector(m1.row(2),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m1.row(3),NULL,N1);
        IPrinter::printVector(m1.row(4),NULL,N1);
        IPrinter::printVector(m1.row(5),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m1.row(6),NULL,N1);
        IPrinter::printVector(m1.row(7),NULL,N1);
        IPrinter::printVector(m1.row(8),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m1.row(9),NULL,N1);
        IPrinter::printVector(m1.row(10),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(a211,NULL,N1);
        IPrinter::printVector(a212,NULL,N1);
        IPrinter::printVector(a213,NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(a311,NULL,N1);
        IPrinter::printVector(a312,NULL,N1);
        IPrinter::printVector(a313,NULL,N1);

        // t=0.00
        u.at(0,0)  = m1.row(0).at(0); u.at(0,1)  = m1.row(1).at(0); u.at(0,2)  = m1.row(2).at(0);
        u.at(0,3)  = a211.at(0);     u.at(0,4)  = a212.at(0);     u.at(0,5)  = a213.at(0);
        u.at(0,6)  = a311.at(0);     u.at(0,7)  = a312.at(0);     u.at(0,8)  = a313.at(0);
        u.at(0,9)  = m1.row(3).at(0); u.at(0,10) = m1.row(4).at(0); u.at(0,11) = m1.row(5).at(0);
        u.at(0,12) = m1.row(6).at(0); u.at(0,13) = m1.row(7).at(0); u.at(0,14) = m1.row(8).at(0);
        b.at(0)    = m1.row(9).at(0);

        // t=0.25
        u.at(1,0)  = m1.row(0).at(30); u.at(1,1)  = m1.row(1).at(30); u.at(1,2)  = m1.row(2).at(30);
        u.at(1,3)  = a211.at(30);     u.at(1,4)  = a212.at(30);     u.at(1,5)  = a213.at(30);
        u.at(1,6)  = a311.at(30);     u.at(1,7)  = a312.at(30);     u.at(1,8)  = a313.at(30);
        u.at(1,9)  = m1.row(3).at(30); u.at(1,10) = m1.row(4).at(30); u.at(1,11) = m1.row(5).at(30);
        u.at(1,12) = m1.row(6).at(30); u.at(1,13) = m1.row(7).at(30); u.at(1,14) = m1.row(8).at(30);
        b.at(1)    = m1.row(9).at(30);

        // t=0.50
        u.at(2,0)  = m1.row(0).at(60); u.at(2,1)  = m1.row(1).at(60); u.at(2,2)  = m1.row(2).at(60);
        u.at(2,3)  = a211.at(60);     u.at(2,4)  = a212.at(60);     u.at(2,5)  = a213.at(60);
        u.at(2,6)  = a311.at(60);     u.at(2,7)  = a312.at(60);     u.at(2,8)  = a313.at(60);
        u.at(2,9)  = m1.row(3).at(60); u.at(2,10) = m1.row(4).at(60); u.at(2,11) = m1.row(5).at(60);
        u.at(2,12) = m1.row(6).at(60); u.at(2,13) = m1.row(7).at(60); u.at(2,14) = m1.row(8).at(60);
        b.at(2)    = m1.row(9).at(60);

        //---------------------------------------------------------------
        puts("++++++++++++++++++++++++++++++++++");

        alpha1[0] = m1.row(0).at(60) + a211.at(60);
        alpha1[1] = m1.row(1).at(60) + a212.at(60);
        alpha1[2] = m1.row(2).at(60) + a213.at(60);

        bettas[0][0] = m1.row(3).at(60);
        bettas[0][1] = m1.row(4).at(60);
        bettas[0][2] = m1.row(5).at(60);

        bettas[1][0] = m1.row(6).at(60);
        bettas[1][1] = m1.row(7).at(60);
        bettas[1][2] = m1.row(8).at(60);

        alpha3[0] = a311.at(60);
        alpha3[1] = a312.at(60);
        alpha3[2] = a313.at(60);

        qmm = m1.row(9).at(60);

        N = 60;
        DoubleMatrix m2;
        calculate(m2, N, alpha1, bettas, qmm, n);

        for (unsigned int i=0; i<=N; i++) a311.at(i) = alpha3[0] * m1.row(10).at(i);
        for (unsigned int i=0; i<=N; i++) a312.at(i) = alpha3[1] * m1.row(10).at(i);
        for (unsigned int i=0; i<=N; i++) a313.at(i) = alpha3[2] * m1.row(10).at(i);

        IPrinter::printVector(m2.row(0),NULL,N1);
        IPrinter::printVector(m2.row(1),NULL,N1);
        IPrinter::printVector(m2.row(2),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m2.row(3),NULL,N1);
        IPrinter::printVector(m2.row(4),NULL,N1);
        IPrinter::printVector(m2.row(5),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m2.row(6),NULL,N1);
        IPrinter::printVector(m2.row(7),NULL,N1);
        IPrinter::printVector(m2.row(8),NULL,N1);
        puts("------------------------------");
        IPrinter::printVector(m2.row(9),NULL,N1);
        IPrinter::printVector(m2.row(10),NULL,N1);

        puts("------------------------------");
        IPrinter::printVector(a311,NULL,N1);
        IPrinter::printVector(a312,NULL,N1);
        IPrinter::printVector(a313,NULL,N1);

//        // t=0.75
//        u.at(3,0)  = m2.row(0).at(30); u.at(3,1)  = m2.row(1).at(90); u.at(3,2)  = m2.row(2).at(30);
//        u.at(3,3)  = a221.at(90);      u.at(3,4)  = a222.at(30);      u.at(3,5)  = a223.at(30);
//        u.at(3,6)  = a321.at(90);      u.at(3,7)  = a322.at(30);      u.at(3,8)  = a323.at(30);
//        u.at(3,9)  = m2.row(3).at(30); u.at(3,10) = m2.row(4).at(30); u.at(3,11) = m2.row(5).at(30);
//        u.at(3,12) = m2.row(6).at(30); u.at(3,13) = m2.row(7).at(30); u.at(3,14) = m2.row(8).at(30);
//        b.at(3)    = m2.row(9).at(30);

//        // t=1.00
//        u.at(4,0)  = m2.row(0).at(60); u.at(4,1)  = m2.row(1).at(60); u.at(4,2)  = m2.row(2).at(60);
//        u.at(4,3)  = a221.at(60);      u.at(4,4)  = a222.at(60);      u.at(4,5)  = a223.at(60);
//        u.at(4,6)  = a321.at(60);      u.at(4,7)  = a322.at(60);      u.at(4,8)  = a323.at(60);
//        u.at(4,9)  = m2.row(3).at(60); u.at(4,10) = m2.row(4).at(60); u.at(4,11) = m2.row(5).at(60);
//        u.at(4,12) = m2.row(6).at(60); u.at(4,13) = m2.row(7).at(60); u.at(4,14) = m2.row(8).at(60);
//        b.at(4)    = m2.row(9).at(60);
    }

    puts("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");


    IPrinter::printMatrix(u, 15, 15);
}

void LoadedSystems::init()
{
    k1 = 2;
    k2 = 3;
    n  = 3;
    N  = 120;
    ht = 1.0 / N;

    t1.resize(k1);
    t1[0] = 0.25;
    t1[1] = 0.75;

    t2.resize(k2);
    t2[0] = 0.00;
    t2[1] = 0.50;
    t2[2] = 1.00;

    alpha.resize(k2); for (unsigned int i=0; i<k2; i++) alpha[i].resize(n,n);

    alpha[0].at(0,0) = +2.0; alpha[0].at(0,1) = +0.0; alpha[0].at(0,2) = +0.0;
    alpha[0].at(1,0) = +0.0; alpha[0].at(1,1) = +1.0; alpha[0].at(1,2) = +0.0;
    alpha[0].at(2,0) = +0.0; alpha[0].at(2,1) = +0.0; alpha[0].at(2,2) = +3.0;

    alpha[1].at(0,0) = +0.0; alpha[1].at(0,1) = +3.0; alpha[1].at(0,2) = +0.0;
    alpha[1].at(1,0) = +0.0; alpha[1].at(1,1) = +0.0; alpha[1].at(1,2) = +0.0;
    alpha[1].at(2,0) = +0.0; alpha[1].at(2,1) = +0.0; alpha[1].at(2,2) = +0.0;

    alpha[2].at(0,0) = +0.0; alpha[2].at(0,1) = +0.0; alpha[2].at(0,2) = +1.0;
    alpha[2].at(1,0) = +2.0; alpha[2].at(1,1) = +0.0; alpha[2].at(1,2) = +0.0;
    alpha[2].at(2,0) = +0.0; alpha[2].at(2,1) = -1.0; alpha[2].at(2,2) = +0.0;

    betta.resize(k1); for (unsigned int s=0; s<k1; s++) betta[s].resize(n,n);

    betta[0].at(0,0) = +0.0; betta[0].at(0,1) = +0.0; betta[0].at(0,2) = +0.0;
    betta[0].at(1,0) = +0.0; betta[0].at(1,1) = +0.0; betta[0].at(1,2) = +0.0;
    betta[0].at(2,0) = +0.0; betta[0].at(2,1) = +0.0; betta[0].at(2,2) = +0.0;

    betta[1].at(0,0) = +0.0; betta[1].at(0,1) = +0.0; betta[1].at(0,2) = +0.0;
    betta[1].at(1,0) = +0.0; betta[1].at(1,1) = +0.0; betta[1].at(1,2) = +0.0;
    betta[1].at(2,0) = +0.0; betta[1].at(2,1) = +0.0; betta[1].at(2,2) = +0.0;

    qamma.resize(n);

    qamma[0] = -1.3569;
    //(alpha[0].at(0,0)*x1(0.00) + alpha[0].at(0,1)*x2(0.00) + alpha[0].at(0,2)*x3(0.00)) +
    //(alpha[1].at(0,0)*x1(0.50) + alpha[1].at(0,1)*x2(0.50) + alpha[1].at(0,2)*x3(0.50)) +
    //(alpha[2].at(0,0)*x1(1.00) + alpha[2].at(0,1)*x2(1.00) + alpha[2].at(0,2)*x3(1.00)) +
    //(betta[0].at(0,0)*x1(0.25) + betta[0].at(0,1)*x2(0.25) + betta[0].at(0,2)*x3(0.25)) +
    //(betta[1].at(0,0)*x1(0.75) + betta[1].at(0,1)*x2(0.75) + betta[1].at(0,2)*x3(0.75));

    qamma[1] = +3.68294;
    //(alpha[0].at(1,0)*x1(0.00) + alpha[0].at(1,1)*x2(0.00) + alpha[0].at(1,2)*x3(0.00)) +
    //(alpha[1].at(1,0)*x1(0.50) + alpha[1].at(1,1)*x2(0.50) + alpha[1].at(1,2)*x3(0.50)) +
    //(alpha[2].at(1,0)*x1(1.00) + alpha[2].at(1,1)*x2(1.00) + alpha[2].at(1,2)*x3(1.00)) +
    //(betta[0].at(1,0)*x1(0.25) + betta[0].at(1,1)*x2(0.25) + betta[0].at(1,2)*x3(0.25)) +
    //(betta[1].at(1,0)*x1(0.75) + betta[1].at(1,1)*x2(0.75) + betta[1].at(1,2)*x3(0.75));

    qamma[2] = 3.0806;
    //(alpha[0].at(2,0)*x1(0.00) + alpha[0].at(2,1)*x2(0.00) + alpha[0].at(2,2)*x3(0.00)) +
    //(alpha[1].at(2,0)*x1(0.50) + alpha[1].at(2,1)*x2(0.50) + alpha[1].at(2,2)*x3(0.50)) +
    //(alpha[2].at(2,0)*x1(1.00) + alpha[2].at(2,1)*x2(1.00) + alpha[2].at(2,2)*x3(1.00)) +
    //(betta[0].at(2,0)*x1(0.25) + betta[0].at(2,1)*x2(0.25) + betta[0].at(2,2)*x3(0.25)) +
    //(betta[1].at(2,0)*x1(0.75) + betta[1].at(2,1)*x2(0.75) + betta[1].at(2,2)*x3(0.75));
}

void LoadedSystems::calculate(unsigned int j, DoubleMatrix &m, unsigned int N,
                              unsigned int k2 UNUSED_PARAM, const std::vector<DoubleMatrix> &alpha,
                              unsigned int k1, const std::vector<DoubleMatrix> &betta,
                              unsigned int n,  const DoubleVector &qamma)
{
    std::vector<CauchyProblem2*> cs(n+k1*n+2);
    for (unsigned int i=0; i<n; i++)
    {
        CauchyProblem* ca = new CauchyProblem(*this);
        ca->x0 = 0.0;
        ca->y0 = alpha[0].at(j, i);
        ca->i = i;
        ca->s = 0;
        ca->type = 0;
        cs[i] = ca;
    }

    for (unsigned int s=0; s<k1; s++)
    {
        for (unsigned int i=0; i<n; i++)
        {
            CauchyProblem* cb = new CauchyProblem(*this);
            cb->x0 = 0.0;
            cb->y0 = betta[s].at(j, i);
            cb->i = i;
            cb->s = s;
            cb->type = 1;
            cs[n+s*n+i] = cb;
        }
    }

    CauchyProblem *cq = new CauchyProblem(*this);
    cq->x0 = 0.0;
    cq->y0 = qamma[j];
    cq->i = 0;
    cq->s = 0;
    cq->type = 2;
    cs[n+k1*n] = cq;

    CauchyProblem *cm = new CauchyProblem(*this);
    cm->x0 = 0.0;
    cm->y0 = 1.0;
    cm->i = 0;
    cm->s = 0;
    cm->type = 3;
    cs[n+k1*n+1] = cm;

    CauchyProblem2::rungeKutta(cs, 0.0, ht, N, m);
}

void LoadedSystems::calculate(DoubleMatrix &m, unsigned int N,
                              const DoubleVector &alpha, const std::vector<DoubleVector> &bettas,
                              double qamma, unsigned int n)
{
    unsigned int k1 = betta.size();

    // alpha + betta * k1 + qamma + M
    unsigned int size = n + k1*n + 2;

    std::vector<CauchyProblem2*> cps(size);

    // alpha
    for (unsigned int i=0; i<n; i++)
    {
        CauchyProblem* ca = new CauchyProblem(*this);
        ca->x0 = 0.0;
        ca->y0 = alpha.at(i);
        ca->i = i;
        ca->s = 0;
        ca->type = 0;
        cps[i] = ca;
    }

    // betta
    for (unsigned int s=0; s<k1; s++)
    {
        for (unsigned int i=0; i<n; i++)
        {
            CauchyProblem* cb = new CauchyProblem(*this);
            cb->x0 = 0.0;
            cb->y0 = bettas[s].at(i);
            cb->i = i;
            cb->s = s;
            cb->type = 1;
            cps[n+s*n+i] = cb;
        }
    }

    //
    CauchyProblem *cq = new CauchyProblem(*this);
    cq->x0 = 0.0;
    cq->y0 = qamma;
    cq->i = 0;
    cq->s = 0;
    cq->type = 2;
    cps[n+k1*n] = cq;

    CauchyProblem *cm = new CauchyProblem(*this);
    cm->x0 = 0.0;
    cm->y0 = 1.0;
    cm->i = 0;
    cm->s = 0;
    cm->type = 3;
    cps[n+k1*n+1] = cm;

    CauchyProblem2::rungeKutta(cps, 0.0, ht, N, m);

    for (unsigned int i=0; i<cps.size(); i++)
    {
        CauchyProblem *cp = dynamic_cast<CauchyProblem*>(cps.at(i));
        delete cp;
    }
    cps.clear();
}

double LoadedSystems::R(const DoubleVector &args) const
{
    double res = 0.0;

    double a111 = args.at(0);
    double a112 = args.at(1);
    double a113 = args.at(2);
    double b111 = args.at(3);
    double b112 = args.at(4);
    double b113 = args.at(5);
    double b211 = args.at(6);
    double b212 = args.at(7);
    double b213 = args.at(8);
    double qmm1 = args.at(9);
    double M    = args.at(10);

    res += a111*a111+a112*a112+a113*a113;
    //    res += b111*b111+b112*b112+b113*b113;
    //    res += b211*b211+b212*b212+b213*b213;
    //    res += qmm1*qmm1;

    //    DoubleVector alpha1 = args.mid(0, n-1);
    //    for (unsigned int i=0; i<n; i++) res += alpha1.at(i)*alpha1.at(i);

    //    for (unsigned int s=0; s<k1; s++)
    //    {
    //        unsigned int a = k2*n + s*n;
    //        unsigned int b = a+n-1;
    //        DoubleVector bettaS = args.mid(a,b);
    //        for (unsigned int i=0; i<n; i++)
    //        {
    //            res += bettaS.at(i)*bettaS.at(i);
    //        }
    //    }

    //    double q = args.at(n+n*k1);
    //    res += q*q;

    return res;
}

double LoadedSystems::S(const DoubleVector &args, double t) const
{
    double res = 0.0;

    double a111 = args.at(0);
    double a112 = args.at(1);
    double a113 = args.at(2);
    double b111 = args.at(3);
    double b112 = args.at(4);
    double b113 = args.at(5);
    double b211 = args.at(6);
    double b212 = args.at(7);
    double b213 = args.at(8);
    double qmm1 = args.at(9);
    double M    = args.at(10);

    res += a111*(a111*A(0,0,t)+a112*A(1,0,t)+a113*A(2,0,t))
            +a112*(a111*A(0,1,t)+a112*A(1,1,t)+a113*A(2,1,t))
            +a113*(a111*A(0,2,t)+a112*A(1,2,t)+a113*A(2,2,t));

    //    res += a111*(b111*B(0,0,0,t)+b112*B(0,1,0,t)+b113*B(0,2,0,t))
    //          +a112*(b111*B(1,0,0,t)+b112*B(1,1,0,t)+b113*B(1,2,0,t))
    //          +a113*(b111*B(2,1,0,t)+b112*B(2,1,0,t)+b113*B(2,2,0,t));

    //    res += a111*(b211*B(0,0,1,t)+b212*B(0,1,1,t)+b213*B(0,2,1,t))
    //          +a112*(b211*B(1,0,1,t)+b212*B(1,1,1,t)+b213*B(1,2,1,t))
    //          +a113*(b211*B(2,1,1,t)+b212*B(2,1,1,t)+b213*B(2,2,1,t));

    //    res -= qmm1*(a111*C(0,t)+a112*C(1,t)+a113*C(2,t));

    //    DoubleVector alpha1 = args.mid(0, n-1);
    //    for (unsigned int i1=0; i1<n; i1++)
    //    {
    //        double sum = 0.0;
    //        for (unsigned int i2=0; i2<n; i2++)
    //        {
    //            sum += alpha1.at(i2)*A(i2, i1, t);
    //        }
    //        res += alpha1.at(i1)*sum;
    //    }

    //    for (unsigned int s=0; s<k1; s++)
    //    {

    //        unsigned int a = k2*n + s*n;
    //        unsigned int b = a+n-1;
    //        DoubleVector bettaS = args.mid(a,b);

    //        for (unsigned int i1=0; i1<n; i1++)
    //        {
    //            double sum = 0.0;
    //            for (unsigned int i2=0; i2<n; i2++)
    //            {
    //                sum += alpha1.at(i2)*B(i2, i1, s, t);
    //            }
    //            res += bettaS.at(i1)*sum;
    //        }
    //    }

    //    double q = args.at(n+n*k1);
    //    for (unsigned int i2=0; i2<n; i2++)
    //    {
    //        res -= q*alpha1.at(i2)*C(i2, t);
    //    }

    res /= R(args);

    return res;
}

double LoadedSystems::A(unsigned int row, unsigned int col, double t UNUSED_PARAM) const
{
    if (row==0) { if (col==0) return +0.0;   if (col==1) return t;    if (col==2) return +2.0*t; }
    if (row==1) { if (col==0) return +3.0*t; if (col==1) return +0.0; if (col==2) return -1.0;}
    if (row==2) { if (col==0) return +1.0;   if (col==1) return +2.0; if (col==2) return +0.0;}
    return 0.0;
}

double LoadedSystems::B(unsigned int row, unsigned int col, unsigned int s, double t UNUSED_PARAM) const
{
    if (s==0)
    {
        if (row==0) { if (col==0) return +0.0; if (col==1) return +1.0; if (col==2) return +0.0; }
        if (row==0) { if (col==0) return +0.0; if (col==1) return +0.0; if (col==2) return +0.0; }
        if (row==0) { if (col==0) return +0.0; if (col==1) return -1.0; if (col==2) return +0.0; }
    }

    if (s==1)
    {
        if (row==0) { if (col==0) return +0.0; if (col==1) return +0.0; if (col==2) return +0.0; }
        if (row==0) { if (col==0) return +0.0; if (col==1) return +0.0; if (col==2) return +1.0; }
        if (row==0) { if (col==0) return +0.0; if (col==1) return +0.0; if (col==2) return +0.0; }
    }
    return 0.0;
}

double LoadedSystems::C(unsigned int row, double t) const
{
    if (row == 0) return -t*t*t - 6.0*t*t + 2.0*t*(cos(t)+sin(t)-1.0) + cos(t) + 3.87532;//3.8753248434212895682891908989884;
    if (row == 1) return -6.0*t*t + t*(5.0 - 3.0*sin(t)) + sin(t) - 1.56836;//1.5683612399766658332667580472201;
    if (row == 2) return -2.0*t*t - 2.0*t + 3.0*cos(t) - sin(t) + 1.12468;//1.1246751565787104317108091010116;
    return 0.0;
}

//double LoadedSystems::A(unsigned int row, unsigned int col, double t UNUSED_PARAM) const
//{
//    if (row==0)
//    {
//        if (col==0) return +2.0;
//        if (col==1) return -3.0;
//    }
//    if (row==1)
//    {
//        if (col==0) return +1.0;
//        if (col==1) return +4.0;
//    }
//    return 0.0;
//}

//double LoadedSystems::B(unsigned int i, unsigned int j, unsigned int s, double t UNUSED_PARAM) const
//{
//    if (s==0)
//    {
//        if (i==0)
//        {
//            if (j==0) return +1.0;
//            if (j==1) return +2.0;
//        }
//        if (i==1)
//        {
//            if (j==0) return +3.0;
//            if (j==1) return +1.0;
//        }
//    }

//    if (s==1)
//    {
//        if (i==0)
//        {
//            if (j==0) return +3.0;
//            if (j==1) return +1.0;
//        }
//        if (i==1)
//        {
//            if (j==0) return +2.0;
//            if (j==1) return +4.0;
//        }
//    }
//    return 0.0;
//}

//double LoadedSystems::C(unsigned int row, double t) const
//{
//    if (row == 0) return 0.0;
//    if (row == 1) return 0.0;
//    return 0.0;
//}

double LoadedSystems::x1(double t) const { return 2.0*t + sin(t); }
double LoadedSystems::x2(double t) const { return t*t - 2.0*cos(t); }
double LoadedSystems::x3(double t) const { return 3.0*t - sin(t) + 1.0; }
