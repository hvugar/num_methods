#include "bordertest.h"

void BorderTest::Main(int argc, char *argv[])
{
    BorderTest bt;
    DoubleMatrix u;
    bt.calculateU(u, bt.hx, bt.ht, bt.N, bt.M, bt.a);
    //IPrinter::printMatrix(u);
    puts("-------------------------------------------------------------------");
    DoubleMatrix u1;
    bt.calculateU1(u1, bt.hx, bt.ht, bt.N, bt.M, bt.a);
    //IPrinter::printMatrix(u1);
    puts("-------------------------------------------------------------------");
    DoubleMatrix u2;
    bt.calculateU2(u2, bt.hx, bt.ht, bt.N, bt.M, bt.a);
    //IPrinter::printMatrix(u2);
    puts("-------------------------------------------------------------------");
    //IPrinter::printVector(18,14,u.row(bt.M));
    //IPrinter::printVector(18,14,u1.row(bt.M));
    IPrinter::printVector(18,14,u2.row(bt.M));
}

BorderTest::BorderTest()
{}

double BorderTest::initial(unsigned int i) const
{
    double x = i*hx;
    return x*x*x;
}

double BorderTest::boundary(Boundary type, unsigned int j) const
{
    double t = j*ht;
    if (type == Left) return t*t;
    if (type == Right) return t*t+1.0;
    return 0.0;
}

double BorderTest::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    return 2.0*t - 6.0*x*a*a;
}

double BorderTest::U(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    return x*x*x + t*t;
}

void BorderTest::calculateU2(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector dd(N-1);
    DoubleVector rx(N-1);

    double alpha = -(a*a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*a*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[j][i] = initial(i);
            }
        }
        else
        {
            if (j != M)
            {
                for (unsigned int i=0; i<=N; i++) u.at(j,i) = (i*hx)*(i*hx)*(i*hx) + (j*ht)*(j*ht);
            }
            else
            {
                u.at(j, 0) = boundary(Left, j);
                u.at(j, N) = boundary(Right, j);

                DoubleVector b1(N+1);
                DoubleVector b2(N+1);

                double qamma1 = u.at(j-1,1)   + ht * f(1,   j) - alpha * u.at(j, 0);
                double qamma2 = u.at(j-1,N-1) + ht * f(N-1, j) - alpha * u.at(j, N);

                DoubleVector p(N+1);
                DoubleVector q(N+1);
                DoubleVector k(N+1);

                for (unsigned int i=2; i<=N-2; i++)
                {
                    p.at(i) = -beta/alpha;
                    q.at(i) = -1.0;
                    k.at(i) = (u.at(j-1,i)+ht*f(i,j))/alpha;
                }

                b1.at(1) = 1.0 + (2.0*a*a*ht)/(hx*hx);
                b1.at(2) = -(a*a*ht)/(hx*hx);

                b2.at(N-2) = -(a*a*ht)/(hx*hx);
                b2.at(N-1) = 1.0 + (2.0*a*a*ht)/(hx*hx);

                for (unsigned int i=2; i<=N-2; i++)
                {
                    b1.at(i+0) = b1.at(i+0) + b1.at(i-1)*p.at(i);
                    b1.at(i+1) = b1.at(i+1) + b1.at(i-1)*q.at(i);
                    qamma1     = qamma1 - b1.at(i-1)*k.at(i);
                }

                printf("%20.8f %20.8f %20.8f\n", b1.at(N-2), b1.at(N-1), qamma1);
                printf("%20.8f %20.8f %20.8f\n", b2.at(N-2), b2.at(N-1), qamma2);

                DoubleMatrix m(2,2);
                m.at(0,0) = b1.at(N-2);
                m.at(0,1) = b1.at(N-1);
                m.at(1,0) = b2.at(N-2);
                m.at(1,1) = b2.at(N-1);
                DoubleVector b(2);
                b.at(0) = qamma1;
                b.at(1) = qamma2;
                DoubleVector x(2);
                GaussianElimination(m, b, x);
                printf("%20.8f %20.8f\n", x.at(0), x.at(1));

                u.at(j,N-1) = ((N-1)*hx)*((N-1)*hx)*((N-1)*hx) + (j*ht)*(j*ht);
                u.at(j,N-2) = ((N-2)*hx)*((N-2)*hx)*((N-2)*hx) + (j*ht)*(j*ht);
                for (unsigned int i=N-3; i>=1; i--)
                {
                    //u.at(j,i) = p.at(i+1)*U(i+1,j) + q.at(i+1)*U(i+2,j) + k.at(i+1);
                    u.at(j,i) = p.at(i+1)*u.at(j,i+1) + q.at(i+1)*u.at(j,i+2)+ k.at(i+1);
                    printf("%d %10.8f %10.8f\n", i, u.at(j,i), U(i,j));
                }
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

