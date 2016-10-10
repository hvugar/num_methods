#include "problem1.h"
//#include <imaging.h>

Problem1::Problem1()
{
    //t0 = 0.0;
    //t1 = 2.0;
    //x0 = 0.0;
    //x1 = 1.0;
    a = 1.0;
    lambda = 1.0;
    Te = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = 1000;
    M = 1000;
}

double Problem1::v(unsigned int j) const
{
    C_UNUSED(j);
    return 2.0;
}

void Problem1::calculate1()
{
    DoubleMatrix u(M+1, N+1);

    DoubleVector a1(N+1);
    DoubleVector b1(N+1);
    DoubleVector c1(N+1);
    DoubleVector d1(N+1);
    DoubleVector rx(N+1);

    for (uint32_t j=0; j<=M; j++)
    {
        if (j==0)
        {
            for (uint32_t i=0; i<=N; i++)
                u.at(j,i) = Te;
        }
        else
        {
            a1[0] = 0.0;
            b1[0] = 1.0 + (a*a*ht)/(hx*hx) - (lambda*a*a*ht)/hx + lambda*ht;
            c1[0] = -(a*a*ht)/(hx*hx);
            d1[0] = u.at(j-1,0) - v(j)*((lambda*(a*a)*ht)/(hx)) + lambda*ht*Te;

            for (uint32_t i=1; i<=N-1; i++)
            {
                a1[i] = -(a*a*ht)/(hx*hx);
                b1[i] = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + ht*lambda;
                c1[i] = -(a*a*ht)/(hx*hx);
                d1[i] = u.at(j-1, i) + lambda*ht*Te;
            }

            a1[N] = -(a*a*ht)/(hx*hx);
            b1[N] = 1.0 + (a*a*ht)/(hx*hx) + (a*a*ht*lambda)/hx + lambda*ht;
            c1[N] = 0.0;
            d1[N] = u.at(j-1,N) + ((lambda*(a*a)*ht)/(hx))*Te + lambda*ht*Te;

            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), rx.size());

            for (uint32_t i=0; i<=N; i++) u.at(j,i) = rx[i];
        }
    }

    IPrinter::printVector(u.row(0));
    IPrinter::printVector(u.row(1));
    IPrinter::printVector(u.row(2));
    IPrinter::printMatrix(u);
}

void Problem1::calculate2()
{
    double hx = 0.001;
    double ht = 0.001;
    double N = 1000;
    double M = 50000;
    double a = 1.0;
    double lambda1 = -1.0;
    double lambda2 = 0.0;
    double lambda3 = 0.0;
    double T0 = 1.0;
    double T1 = 2.0;
    double T2 = 2.0;
    double T3 = 2.0;

    DoubleMatrix u(M+1, N+1);

    DoubleVector a1(N+1);
    DoubleVector b1(N+1);
    DoubleVector c1(N+1);
    DoubleVector d1(N+1);
    DoubleVector rx(N+1);

    for (uint32_t j=0; j<=M; j++)
    {
        if (j==0)
        {
            for (uint32_t i=0; i<=N; i++)
                u.at(j,i) = T0;
        }
        else
        {
            a1[0] = 0.0;
            b1[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambda2*a*a*ht)/hx - lambda1*ht;
            c1[0] = -(a*a*ht)/(hx*hx);
            d1[0] = u.at(j-1,0) + ((lambda2*a*a*ht)/(hx))*T2 - lambda1*ht*T1;

            for (uint32_t i=1; i<=N-1; i++)
            {
                a1[i] = -(a*a*ht)/(hx*hx);
                b1[i] = 1.0 + 2.0*((a*a)*ht)/(hx*hx) - ht*lambda1;
                c1[i] = -(a*a*ht)/(hx*hx);
                d1[i] = u.at(j-1, i) - lambda1*ht*T1;
            }

            a1[N] = -(a*a*ht)/(hx*hx);
            b1[N] = 1.0 + (a*a*ht)/(hx*hx) - lambda3*(a*a*ht)/hx - lambda1*ht;
            c1[N] = 0.0;
            d1[N] = u.at(j-1,N) - ((lambda3*a*a*ht)/(hx))*T3 - lambda1*ht*T1;

            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), rx.size());

            for (uint32_t i=0; i<=N; i++) u.at(j,i) = rx[i];
        }
    }

    IPrinter::printVector(u.row(0));
    IPrinter::printVector(u.row(1));
    IPrinter::printVector(u.row(2));
    puts("...");
    IPrinter::printVector(u.row(M-2));
    IPrinter::printVector(u.row(M-1));
    IPrinter::printVector(u.row(M));
    puts("---");
    IPrinter::printMatrix(u);
}

void Problem1::calculate3()
{
    double hx = 0.001;
    double ht = 0.001;
    double N = 1000;
    double M = 5000;
    double a = 1.0;
    double lambda1 = +1.0;
    double lambda2 = +1.0;
    double lambda3 = +1.0;
    double T0 = 1.0;
    double T1 = 3.0;
    double T2 = 3.0;
    double T3 = 3.0;

    DoubleMatrix u(M+1, N+1);

    DoubleVector a1(N+1);
    DoubleVector b1(N+1);
    DoubleVector c1(N+1);
    DoubleVector d1(N+1);
    DoubleVector rx(N+1);

    for (uint32_t j=0; j<=M; j++)
    {
        if (j==0)
        {
            for (uint32_t i=0; i<=N; i++) u.at(j,i) = T0;
        }
        else
        {
            a1[0] = 0.0;
            b1[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambda2*a*a*ht)/hx + lambda1*ht;
            c1[0] = -(a*a*ht)/(hx*hx);
            d1[0] = u.at(j-1,0) + ((lambda2*a*a*ht)/(hx))*T2 + lambda1*ht*T1;

            for (uint32_t i=1; i<=N-1; i++)
            {
                a1[i] = -(a*a*ht)/(hx*hx);
                b1[i] = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + lambda1*ht;
                c1[i] = -(a*a*ht)/(hx*hx);
                d1[i] = u.at(j-1, i) + lambda1*ht*T1;
            }

            a1[N] = -(a*a*ht)/(hx*hx);
            b1[N] = 1.0 + (a*a*ht)/(hx*hx) + lambda3*(a*a*ht)/hx + lambda1*ht;
            c1[N] = 0.0;
            d1[N] = u.at(j-1,N) + ((lambda3*a*a*ht)/(hx))*T3 + lambda1*ht*T1;

            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), rx.size());

            for (uint32_t i=0; i<=N; i++) u.at(j,i) = rx[i];
        }
    }

//    QPixmap img;
//    visualizeMatrixHeat(u, 1.0, 3.0, img);
//    img.save("img.png");

//    for (unsigned int j=0; j<=M; j++)
//    {
//        QPixmap img;
//        visualizeVectorHeat(u[j], 1.0, 3.0, img, N);
//        img.save(QString("img%1.png").arg(j));
//        IPrinter::printVector(u[j]);
//    }
    IPrinter::printVector(u.row(0));
    IPrinter::printVector(u.row(1));
    IPrinter::printVector(u.row(2));
    puts("...");
    IPrinter::printVector(u.row(M-2));
    IPrinter::printVector(u.row(M-1));
    IPrinter::printVector(u.row(M));
    puts("---");
    IPrinter::printMatrix(u);
}

/**
 * @brief Problem1::calculate4
 * u_t = a^2 u_xx
 * u(x,0) = T_e
 * u_x(0,t) = lambda * (v(t) - u(0,t)
 * u_x(l,t) = 0
 */
void Problem1::calculate4()
{
    double hx = 0.001;
    double ht = 0.001;
    double N = 1000;
    double M = 1000;
    double a = 1.0;
    double lambda1 = +1.0;
    double lambda2 = -1.0;
    double lambda3 = +1.0;
    double T0 = 1.0;
    //double T1 = 1.0;
    double T2 = 2.0;
    //double T3 = 1.0;

    DoubleMatrix u(M+1, N+1);

    DoubleVector a1(N+1);
    DoubleVector b1(N+1);
    DoubleVector c1(N+1);
    DoubleVector d1(N+1);
    DoubleVector rx(N+1);

    for (uint32_t j=0; j<=M; j++)
    {
        if (j==0)
        {
            for (uint32_t i=0; i<=N; i++)
                u.at(j,i) = T0;
        }
        else
        {
            a1[0] = 0.0;
            b1[0] = 1.0 + (a*a*ht)/(hx*hx) - (lambda2*a*a*ht)/hx;
            c1[0] = -(a*a*ht)/(hx*hx);
            d1[0] = u.at(j-1,0) - ((lambda2*a*a*ht)/(hx))*T2;

            for (uint32_t i=1; i<=N-1; i++)
            {
                a1[i] = -(a*a*ht)/(hx*hx);
                b1[i] = 1.0 + 2.0*((a*a)*ht)/(hx*hx);
                c1[i] = -(a*a*ht)/(hx*hx);
                d1[i] = u.at(j-1, i);
            }

            a1[N] = -(a*a*ht)/(hx*hx);
            b1[N] = 1.0 + (a*a*ht)/(hx*hx);
            c1[N] = 0.0;
            d1[N] = u.at(j-1,N);

            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), rx.size());

            for (uint32_t i=0; i<=N; i++) u.at(j,i) = rx[i];
        }
    }

    IPrinter::printVector(u.row(0));
    IPrinter::printVector(u.row(1));
    IPrinter::printVector(u.row(2));
    puts("...");
    IPrinter::printVector(u.row(M-2));
    IPrinter::printVector(u.row(M-1));
    IPrinter::printVector(u.row(M));
    puts("---");
    IPrinter::printMatrix(u);
}
