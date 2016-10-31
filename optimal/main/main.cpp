#include "problem1.h"
#include "problem1k.h"
#include "problem1z.h"
#include "problem1x.h"
#include "problem1x1.h"
#include "problem1kz.h"
#include "problem1kzx.h"

#include "problem3.h"
#include "loadedsystems.h"
#include "example1.h"
#include "example2.h"
#include "example3.h"
#include <cmethods.h>

#include <QtGui/QGuiApplication>
//#include <imaging.h>

#include <float.h>
#include <time.h>

#include <parabolicequation.h>

class A : public IParabolicEquation
{
public:
    double hx = 0.001;
    double ht = 0.001;
    unsigned int N = 1000;
    unsigned int M = 1000;
    double a = 1.0;

    virtual double initial(unsigned int i) const
    {
        double x = i*hx;
        return x*x*x;
    }

    virtual double boundary(Boundary type, unsigned int j UNUSED_PARAM) const
    {
        if (type == Left)  return 0.0;
        if (type == Right) return 3.0;
        return 0.0;
    }
    virtual double f(unsigned int i, unsigned int j) const
    {
        double x = i*hx;
        double t = j*ht;
        return 2.0*t-6.0*x*a*a;
    }

    double U(unsigned int i,unsigned int j)
    {
        double x = i*hx;
        double t = j*ht;
        return x*x*x + t*t;
    }

    void calculate(DoubleMatrix &u)
    {
        u.clear();
        u.resize(M+1, N+1);

        ////////////////////////////////////////
        {
            for (unsigned int j=0; j<=M; j++)
            {
                for (unsigned int i=0; i<=N; i++)
                {
                    u.at(j,i) = U(i,j);
                }
            }

            double alpha = -(a*a*ht)/(hx*hx);
            double betta  = 1.0 + (2.0*a*a*ht)/(hx*hx);

            DoubleMatrix alfa(N-1, 3);

            for (unsigned int i=0; i<=N-2; i++)
            {
                double d = 1.0 - (a*a*ht)/(hx*hx);
                alfa.at(i, 1) = ((-2.0*a*a*ht)/(hx*hx))/d;
                alfa.at(i, 2) = ((a*a*ht)/(hx*hx))/d;
                alfa.at(i, 0) = (u.at(M-1,i)+ht*f(i,M))/d;


//                alfa.at(i,1) = -betta/alpha;
//                alfa.at(i,2) = -alpha/alpha;
//                alfa.at(i,0) = (u.at(M-1,i) + ht * f(i, M))/alpha;
            }

            u.at(M, N-1) = U(N-1, M);
            u.at(M, N-0) = U(N-0, M);
            for (unsigned int i=N-2; i != UINT_MAX; i--)
            {
                u.at(M,i) = alfa.at(i,1)*u.at(M,i+1) + alfa.at(i,2)*u.at(M,i+2) + alfa.at(i,0);
                //u.at(M,i) = alfa.at(i,1)*U(i+1,M) + alfa.at(i,2)*U(i+2,M) + alfa.at(i,0);
            }

            puts("-------------------------");
            IPrinter::printMatrix(u);
            puts("-------------------------");

            return;
        }
        //////////////////////////////////////////

        for (unsigned int j=0; j<=M; j++)
        {
            if (j == 0)
            {
                for (unsigned int i=0; i<=N; i++)
                {
                    u.at(0,i) = initial(i);
                }
            }
            else
            {
                //                if (j==M)
                {
                    double alpha = -(a*a*ht)/(hx*hx);
                    double betta  = 1.0 + (2.0*a*a*ht)/(hx*hx);

                    DoubleMatrix alfa(N-1, 3);

                    for (unsigned int i=1; i<=N-1; i++)
                    {
                        alfa.at(i-1,1) = -betta/alpha;
                        alfa.at(i-1,2) = -alpha/alpha;
                        alfa.at(i-1,0) = (u.at(j-1,i) + ht * f(i, j))/alpha;
                    }

                    DoubleMatrix beta(2, N+1, 0.0);
                    DoubleVector qamma(2);

                    beta.at(0,0) = 1.0 + (a*a*ht)/(hx*hx);
                    beta.at(0,1) = -(a*a*ht)/(hx*hx);
                    qamma.at(0) = u.at(j-1,0) + ht * f(0, j) - (a*a*ht)/(hx)*boundary(Left, j);

                    beta.at(1,N-1) = -(a*a*ht)/(hx*hx);
                    beta.at(1,N-0) = 1.0 + (a*a*ht)/(hx*hx);
                    qamma.at(1) = u.at(j-1,N) + ht * f(N, j) + (a*a*ht)/(hx)*boundary(Right, j);

                    for (unsigned int i=0; i<=N-2; i++)
                    {
                        beta.at(0,i+1) = beta.at(0,i+1) + beta.at(0,i)*alfa.at(i,1);
                        beta.at(0,i+2) = beta.at(0,i+2) + beta.at(0,i)*alfa.at(i,2);
                        qamma.at(0)    = qamma.at(0)    - beta.at(0,i)*alfa.at(i,0);

                        beta.at(1,i+1) = beta.at(1,i+1) + beta.at(1,i)*alfa.at(i,1);
                        beta.at(1,i+2) = beta.at(1,i+2) + beta.at(1,i)*alfa.at(i,2);
                        qamma.at(1)    = qamma.at(1)    - beta.at(1,i)*alfa.at(i,0);
                    }

                    DoubleMatrix m1(2,2);
                    m1.at(0,0) = beta.at(0,N-1);
                    m1.at(0,1) = beta.at(0,N-0);
                    m1.at(1,0) = beta.at(1,N-1);
                    m1.at(1,1) = beta.at(1,N-0);

                    DoubleVector b1(2);
                    b1.at(0) = qamma.at(0);
                    b1.at(1) = qamma.at(1);

                    DoubleVector u1(2);
                    GaussianElimination(m1, b1, u1);

                    m1.clear();
                    b1.clear();

                    u.at(j, N-1) = u1.at(0);
                    u.at(j, N-0) = u1.at(1);

                    if (j==M)
                    {
                        u.at(j, N-1) = 1.9980010000;
                        u.at(j, N-0) = 2.0000000000;
                    }

                    for (unsigned int i=N-2; i != UINT_MAX; i--)
                    {
                        u.at(j,i) = alfa.at(i,1)*u.at(j,i+1) + alfa.at(i,2)*u.at(j,i+2) + alfa.at(i,0);
                    }
                }
                //                else
                //                {
                //                    for (unsigned int i=0; i <= N; i++)
                //                    {
                //                        u.at(j,i) = U(j,i);
                //                    }
                //                }
            }
        }
    }
};

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    srand(time(NULL));

//    A a;

//    DoubleVector u0(a.N+1);
//    for (unsigned int i=0; i<=a.N; i++)
//        u0.at(i) = a.U(i,a.M);
//    IPrinter::printVector(u0);

//    DoubleMatrix u1;
//    a.calculateN1(u1, a.hx, a.ht, a.N, a.M);
//    IPrinter::printVector(u1.row(a.M));

//    //    DoubleMatrix u3;
//    //    a.calculateN(u3, a.hx, a.ht, a.N, a.M);
//    //    IPrinter::printVector(u3.row(a.M));

//    DoubleMatrix u2;
//    a.calculate(u2);
//    IPrinter::printVector(u2.row(a.M));

//    FILE *file1 = fopen("temp1.txt", "w");
//    IPrinter::printVector(u0,NULL,u0.size(),0,0,file1);
//    IPrinter::printVector(u1.row(a.M),NULL,u1.cols(),0,0,file1);
//    IPrinter::printVector(u2.row(a.M),NULL,u2.cols(),0,0,file1);
//    fclose(file1);

    //    Problem3::Main(argc, argv);
    //    Problem1KZX::Main(argc, argv);
        Example3::Main(argc, argv);
    //    Example2::Main(argc, argv);
    //    Problem1K::Main(argc, argv);




    return 0;
}
