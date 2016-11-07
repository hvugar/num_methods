#include "bordertest.h"

void BorderTest::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderTest bt;
    DoubleMatrix u0(bt.M+1, bt.N+1);
    for (unsigned int j=0; j<=bt.M; j++)
    {
        for (unsigned int i=0; i<=bt.N; i++)
            u0.at(j,i) = bt.U(i,j);
    }
    IPrinter::printMatrix(u0);
    puts("-------------------------------------------------------------------");
    DoubleMatrix u;
    bt.calculateU(u, bt.hx, bt.ht, bt.N, bt.M, bt.a);
    IPrinter::printMatrix(u);
    puts("-------------------------------------------------------------------");
    //DoubleMatrix u1;
    //bt.calculateU1(u1, bt.hx, bt.ht, bt.N, bt.M, bt.a);
    //IPrinter::printMatrix(u1);
    //puts("-------------------------------------------------------------------");
    //DoubleMatrix u2;
    //bt.calculateUError2(u2, bt.hx, bt.ht, bt.N, bt.M, bt.a);
    //IPrinter::printMatrix(u2);
    //puts("-------------------------------------------------------------------");

    DoubleMatrix u3f;
    bt.calculateU3Forward(u3f, bt.hx, bt.ht, bt.N, bt.M, bt.a);
    IPrinter::printMatrix(u3f);
    puts("-------------------------------------------------------------------");

    DoubleMatrix u3b;
    bt.calculateU3Backward(u3b, bt.hx, bt.ht, bt.N, bt.M, bt.a);
    IPrinter::printMatrix(u3b);
    puts("-------------------------------------------------------------------");

    IPrinter::printVector(18,14,u0.row(1));
    for (unsigned int i=0; i<=bt.N; i++) u0.at(bt.M,i) = bt.U(i,bt.M);
    IPrinter::printVector(18,14,u0.row(bt.M));
    //IPrinter::printVector(14,12,u.row(1));
    //IPrinter::printVector(14,12,u1.row(1));
    //IPrinter::printVector(14,12,u2.row(bt.M));
    puts("---");
    IPrinter::printVector(18,14,u3f.row(1));
    IPrinter::printVector(18,14,u3f.row(bt.M));
    puts("---");
    //IPrinter::printVector(18,14,u3b.row(1));
    //IPrinter::printVector(18,14,u3b.row(bt.M));
}

BorderTest::BorderTest()
{}

double BorderTest::initial(unsigned int i) const
{
    double x = i*hx;
    //return x*x*x;
    return x*x;
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
    //return 2.0*t - 6.0*x*a*a;
    return 2.0*t - 2.0*a*a;
}

double BorderTest::U(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    //return x*x*x + t*t;
    return x*x + t*t;
}

void BorderTest::calculateUError2(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(M+1, N+1);

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
            //if (j != M)
            //{
            //    for (unsigned int i=0; i<=N; i++) u.at(j,i) = U(i,j);
            //}
            //else
            {
                //for (unsigned int i=0; i<=N; i++) u.at(j,i) = U(i,j);
                //continue;

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

                //printf("%20.8f %20.8f %20.8f\n", b1.at(N-2), b1.at(N-1), qamma1);
                //printf("%20.8f %20.8f %20.8f\n", b2.at(N-2), b2.at(N-1), qamma2);

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

                //printf("%20.14f %20.14f\n", x.at(0), x.at(1));
                //printf("%20.14f %20.14f\n", U(N-2,M), U(N-1,M));

                //u.at(j,N-2) = U(N-2,M);
                //u.at(j,N-1) = U(N-1,M);

                u.at(j,N-2) = x.at(0);
                u.at(j,N-1) = x.at(1);

                for (unsigned int i=N-2; i>=2; i--)
                {
                    //u.at(j,i) = p.at(i+1)*U(i+1,j) + q.at(i+1)*U(i+2,j) + k.at(i+1);
                    //u.at(j,i) = p.at(i+1)*u.at(j,i+1) + q.at(i+1)*u.at(j,i+2)+ k.at(i+1);
                    //u.at(j,i-1) = (-beta/alpha)*u.at(j,i) + (-1.0)*u.at(j,i+1) + (u.at(j-1,i)+ht*f(i,j))/alpha;
                    u.at(j,i-1) = p.at(i)*u.at(j,i) + q.at(i)*u.at(j,i+1) + k.at(i);
                    //printf("%d %10.8f %10.8f\n", i, u.at(j,i), U(i,j));
                }
            }
        }
    }
}

void BorderTest::calculateU3Forward(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(M+1, N+1);

    //double ka = -((a*a)*ht)/(hx*hx);
    //double kb = 1.0 - 2.0*ka;

    double ka = -1000.000000;
    double kb = +2001.000000;

    printf("%18.14f %18.14f\n", ka, kb);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u.at(j,i) = initial(i);
            }
        }
        else
        {
            //if (j != M)
            //{
            //    for (unsigned int i=0; i<=N; i++) u.at(j,i) = U(i,j);
            //}
            //else
            {
                u.at(j, 0) = boundary(Left, j);
                u.at(j, N) = boundary(Right, j);

                DoubleVector betta(N-1, 0.0);
                betta.at(0) = kb;
                betta.at(1) = ka;
                double qamma = u.at(j-1,1) + ht * f(1, j) - ka * u.at(j, 0);


                for (unsigned int i=0; i<N-3; i++)
                {
                    betta.at(i+1) = betta.at(i+1) + betta.at(i)*(-kb/ka);
                    betta.at(i+2) = betta.at(i+2) + betta.at(i)*(-1.0);
                    qamma = qamma - betta.at(i)*((u.at(j-1,i+2) + ht * f(i+2, j))/ka);
                }

                DoubleMatrix m(2,2);
                m.at(0,0) = betta.at(N-3);
                m.at(0,1) = betta.at(N-2);
                m.at(1,0) = ka;
                m.at(1,1) = kb;
                DoubleVector b1(2);
                b1.at(0) = qamma;
                b1.at(1) = u.at(j-1,N-1) + ht * f(N-1, j) - ka * u.at(j, N);

                DoubleVector x1(2);
                GaussianElimination(m,b1,x1);

                if (j==1)
                printf("%18.14f %18.14f\n", x1.at(0), x1.at(1));

                u.at(j, N-2) = x1.at(0);
                u.at(j, N-1) = x1.at(1);
                for (unsigned int i=N-2; i!=1; i--)
                {
                    u.at(j, i-1) = (-kb/ka)*u.at(j, i) - u.at(j, i+1) + ((u.at(j-1,i) + ht * f(i, j))/(ka));
                }
            }
        }
    }
}

void BorderTest::calculateU3Backward(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(M+1, N+1);

    double ka = -((a*a)*ht)/(hx*hx);
    double kb = 1.0 - 2.0*ka;

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u.at(j,i) = initial(i);
            }
        }
        else
        {
            u.at(j, 0) = boundary(Left, j);
            u.at(j, N) = boundary(Right, j);

            DoubleVector betta(N-1, 0.0);
            betta.at(N-2) = kb;//-->N-1
            betta.at(N-3) = ka;//-->N-2
            double qamma = u.at(j-1,N-1) + ht * f(N-1, j) - ka * u.at(j, N);

            for (unsigned int i=N-2; i>=2; i--)
            {
                betta.at(i-1) = betta.at(i-1) + betta.at(i)*(-kb/ka);
                betta.at(i-2) = betta.at(i-2) + betta.at(i)*(-1.0);
                qamma = qamma - betta.at(i)*((u.at(j-1,i) + ht * f(i, j))/ka);
            }

            DoubleMatrix m(2,2);
            m.at(0,0) = kb;
            m.at(0,1) = ka;
            m.at(1,0) = betta.at(0);
            m.at(1,1) = betta.at(1);
            DoubleVector b1(2);
            b1.at(0) = u.at(j-1,1) + ht * f(1, j) - ka * u.at(j, 0);
            b1.at(1) = qamma;

            DoubleVector x1(2);
            GaussianElimination(m,b1,x1);

            //if (j==1)
            //    printf("%.8f %.8f\n", x1.at(0), x1.at(1));

            u.at(j, 1) = x1.at(0);
            u.at(j, 2) = x1.at(1);
            for (unsigned int i=2; i<=N-2; i++)
            {
                u.at(j, i+1) = (-kb/ka)*u.at(j, i) + (-1.0)*u.at(j, i-1) + ((u.at(j-1,i) + ht * f(i, j))/(ka));
            }
        }
    }
}
