#include "example4.h"

void Example4::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Example4 e;
    e.calculate3();
}

Example4::Example4()
{
}

void Example4::calculate1()
{
    h = 0.00001;
    N = 100000;

    DoubleVector x01(N+1);
    DoubleVector x02(N+1);
    DoubleVector x03(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x01.at(i) = X1(i);
        x02.at(i) = X2(i);
        x03.at(i) = X3(i);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    IPrinter::printVector(14,10,x03);

    DoubleVector x1(N+1);
    DoubleVector x2(N+2);
    DoubleVector x3(N+2);

    x1.at(0) = X1(0); x2.at(0) = X2(0); x3.at(0) = X3(0);
    x1.at(1) = X1(1); x2.at(1) = X2(1); x3.at(1) = X3(1);
    x1.at(2) = X1(2); x2.at(2) = X2(2); x3.at(2) = X3(2);
    x1.at(3) = X1(3); x2.at(3) = X2(3); x3.at(3) = X3(3);

    puts("---");
    for (unsigned int k=4; k<=N; k++)
    {
        unsigned int k1 = k-1;
        double alpha1 = +1.92;//+48.0/25.0;
        double alpha2 = -1.44;//-36.0/25.0;
        double alpha3 = +0.64;//+16.0/25.0;
        double alpha4 = -0.12;//-3.0/25.0;
        double alpha5 = +0.48*h;//12.0/25.0;

        x1.at(k) = alpha1*x1.at(k-1) + alpha2*x1.at(k-2) + alpha3*x1.at(k-3) + alpha4*x1.at(k-4)
                + (alpha5*a(1,1,k1))*x1.at(k1) + (alpha5*a(1,2,k1))*x2.at(k1) + (alpha5*a(1,3,k1))*x3.at(k1)
                + (alpha5*b(1,k1));
        x2.at(k) = alpha1*x2.at(k-1) + alpha2*x2.at(k-2) + alpha3*x2.at(k-3) + alpha4*x2.at(k-4)
                + (alpha5*a(2,1,k1))*x1.at(k1) + (alpha5*a(2,2,k1))*x2.at(k1) + (alpha5*a(2,3,k1))*x3.at(k1)
                + (alpha5*b(2,k1));
        x3.at(k) = alpha1*x3.at(k-1) + alpha2*x3.at(k-2) + alpha3*x3.at(k-3) + alpha4*x3.at(k-4)
                + (alpha5*a(3,1,k1))*x1.at(k1) + (alpha5*a(3,2,k1))*x2.at(k1) + (alpha5*a(3,3,k1))*x3.at(k1)
                + (alpha5*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
}

void Example4::calculate3()
{
    h = 0.001;
    N = 1000;

    DoubleVector x01(N+1);
    DoubleVector x02(N+1);
    DoubleVector x03(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x01.at(i) = X1(i);
        x02.at(i) = X2(i);
        x03.at(i) = X3(i);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    IPrinter::printVector(14,10,x03);

    double alpha1 = +1.92;//+48.0/25.0;
    double alpha2 = -1.44;//-36.0/25.0;
    double alpha3 = +0.64;//+16.0/25.0;
    double alpha4 = -0.12;//-3.0/25.0;
    double alpha5 = +0.48*h;//-12.0/25.0;

    DoubleVector x1(N+1);
    DoubleVector x2(N+1);
    DoubleVector x3(N+1);

    x1.at(0) = X1(0); x1.at(1) = X1(1); x1.at(2) = X1(2); x1.at(3) = X1(3);
    x2.at(0) = X2(0); x2.at(1) = X2(1); x2.at(2) = X2(2); x2.at(3) = X2(3);
    x3.at(0) = X3(0); x3.at(1) = X3(1); x3.at(2) = X3(2); x3.at(3) = X3(3);

    DoubleCube p3(N+1, 3, 3);
    DoubleCube p2(N+1, 3, 3);
    DoubleCube p1(N+1, 3, 3);
    DoubleCube p0(N+1, 3, 3);
    DoubleMatrix q0(N+1, 3);

    for (unsigned int k=4; k<=N; k++)
    {
        if (k==4)
        {
            unsigned int k1 = k-1;

            p3.at(k,0,0) = 0.48*h*a(1,1,k1)+1.92; p3.at(k,0,1) = 0.48*h*a(1,2,k1);       p3.at(k,0,2) = 0.48*h*a(1,3,k1);
            p3.at(k,1,0) = 0.48*h*a(2,1,k1);      p3.at(k,1,1) = 0.48*h*a(2,2,k1)+1.92;  p3.at(k,1,2) = 0.48*h*a(2,3,k1);
            p3.at(k,2,0) = 0.48*h*a(3,1,k1);      p3.at(k,2,1) = 0.48*h*a(3,2,k1);       p3.at(k,2,2) = 0.48*h*a(3,3,k1)+1.92;

            p2.at(k,0,0) = -1.44; p2.at(k,0,1) = 0.0;   p2.at(k,0,2) = 0.0;
            p2.at(k,1,0) = 0.0;   p2.at(k,1,1) = -1.44; p2.at(k,1,2) = 0.0;
            p2.at(k,2,0) = 0.0;   p2.at(k,2,1) = 0.0;   p2.at(k,2,2) = -1.44;

            p1.at(k,0,0) = 0.64;  p1.at(k,0,1) = 0.0;   p1.at(k,0,2) = 0.0;
            p1.at(k,1,0) = 0.0;   p1.at(k,1,1) = 0.64;  p1.at(k,1,2) = 0.0;
            p1.at(k,2,0) = 0.0;   p1.at(k,2,1) = 0.0;   p1.at(k,2,2) = 0.64;

            p0.at(k,0,0) = -0.12; p0.at(k,0,1) = 0.0;   p0.at(k,0,2) = 0.0;
            p0.at(k,1,0) = 0.0;   p0.at(k,1,1) = -0.12; p0.at(k,1,2) = 0.0;
            p0.at(k,2,0) = 0.0;   p0.at(k,2,1) = 0.0;   p0.at(k,2,2) = -0.12;

            q0.at(k,0) = 0.48*h*b(1,k1);
            q0.at(k,1) = 0.48*h*b(2,k1);
            q0.at(k,2) = 0.48*h*b(3,k1);

            x1.at(k) = (p3.at(k,0,0)*x1.at(3)+p3.at(k,0,1)*x2.at(3)+p3.at(k,0,2)*x3.at(3))
                     + (p2.at(k,0,0)*x1.at(2)+p2.at(k,0,1)*x2.at(2)+p2.at(k,0,2)*x3.at(2))
                     + (p1.at(k,0,0)*x1.at(1)+p1.at(k,0,1)*x2.at(1)+p1.at(k,0,2)*x3.at(1))
                     + (p0.at(k,0,0)*x1.at(0)+p0.at(k,0,1)*x2.at(0)+p0.at(k,0,2)*x3.at(0)) + q0.at(k,0);
            x2.at(k) = (p3.at(k,1,0)*x1.at(3)+p3.at(k,1,1)*x2.at(3)+p3.at(k,1,2)*x3.at(3))
                     + (p2.at(k,1,0)*x1.at(2)+p2.at(k,1,1)*x2.at(2)+p2.at(k,1,2)*x3.at(2))
                     + (p1.at(k,1,0)*x1.at(1)+p1.at(k,1,1)*x2.at(1)+p1.at(k,1,2)*x3.at(1))
                     + (p0.at(k,1,0)*x1.at(0)+p0.at(k,1,1)*x2.at(0)+p0.at(k,1,2)*x3.at(0)) + q0.at(k,1);
            x3.at(k) = (p3.at(k,2,0)*x1.at(3)+p3.at(k,2,1)*x2.at(3)+p3.at(k,2,2)*x3.at(3))
                     + (p2.at(k,2,0)*x1.at(2)+p2.at(k,2,1)*x2.at(2)+p2.at(k,2,2)*x3.at(2))
                     + (p1.at(k,2,0)*x1.at(1)+p1.at(k,2,1)*x2.at(1)+p1.at(k,2,2)*x3.at(1))
                     + (p0.at(k,2,0)*x1.at(0)+p0.at(k,2,1)*x2.at(0)+p0.at(k,2,2)*x3.at(0)) + q0.at(k,2);
        }

        if (k==5)
        {
            unsigned int k1 = k-1;

            DoubleMatrix A51(3,3);
            A51.at(0,0) = 0.48*h*a(1,1,k1)+1.92; A51.at(0,1) = 0.48*h*a(1,2,k1);      A51.at(0,2) = 0.48*h*a(1,3,k1);
            A51.at(1,0) = 0.48*h*a(2,1,k1);      A51.at(1,1) = 0.48*h*a(2,2,k1)+1.92; A51.at(1,2) = 0.48*h*a(2,3,k1);
            A51.at(2,0) = 0.48*h*a(3,1,k1);      A51.at(2,1) = 0.48*h*a(3,2,k1);      A51.at(2,2) = 0.48*h*a(3,3,k1)+1.92;

            //p3
            p3.at(k,0,0) = A51.at(0,0)*p3.at(4,0,0) + A51.at(0,1)*p3.at(4,1,0) + A51.at(0,2)*p3.at(4,2,0) + (-1.44);
            p3.at(k,0,1) = A51.at(0,0)*p3.at(4,0,1) + A51.at(0,1)*p3.at(4,1,1) + A51.at(0,2)*p3.at(4,2,1) + (0.0);
            p3.at(k,0,2) = A51.at(0,0)*p3.at(4,0,2) + A51.at(0,1)*p3.at(4,1,2) + A51.at(0,2)*p3.at(4,2,2) + (0.0);
            p3.at(k,1,0) = A51.at(1,0)*p3.at(4,0,0) + A51.at(1,1)*p3.at(4,1,0) + A51.at(1,2)*p3.at(4,2,0) + (0.0);
            p3.at(k,1,1) = A51.at(1,0)*p3.at(4,0,1) + A51.at(1,1)*p3.at(4,1,1) + A51.at(1,2)*p3.at(4,2,1) + (-1.44);
            p3.at(k,1,2) = A51.at(1,0)*p3.at(4,0,2) + A51.at(1,1)*p3.at(4,1,2) + A51.at(1,2)*p3.at(4,2,2) + (0.0);
            p3.at(k,2,0) = A51.at(2,0)*p3.at(4,0,0) + A51.at(2,1)*p3.at(4,1,0) + A51.at(2,2)*p3.at(4,2,0) + (0.0);
            p3.at(k,2,1) = A51.at(2,0)*p3.at(4,0,1) + A51.at(2,1)*p3.at(4,1,1) + A51.at(2,2)*p3.at(4,2,1) + (0.0);
            p3.at(k,2,2) = A51.at(2,0)*p3.at(4,0,2) + A51.at(2,1)*p3.at(4,1,2) + A51.at(2,2)*p3.at(4,2,2) + (-1.44);

            //p2
            p2.at(k,0,0) = A51.at(0,0)*p2.at(4,0,0) + A51.at(0,1)*p2.at(4,1,0) + A51.at(0,2)*p2.at(4,2,0) + (0.64);
            p2.at(k,0,1) = A51.at(0,0)*p2.at(4,0,1) + A51.at(0,1)*p2.at(4,1,1) + A51.at(0,2)*p2.at(4,2,1) + (0.0);
            p2.at(k,0,2) = A51.at(0,0)*p2.at(4,0,2) + A51.at(0,1)*p2.at(4,1,2) + A51.at(0,2)*p2.at(4,2,2) + (0.0);
            p2.at(k,1,0) = A51.at(1,0)*p2.at(4,0,0) + A51.at(1,1)*p2.at(4,1,0) + A51.at(1,2)*p2.at(4,2,0) + (0.0);
            p2.at(k,1,1) = A51.at(1,0)*p2.at(4,0,1) + A51.at(1,1)*p2.at(4,1,1) + A51.at(1,2)*p2.at(4,2,1) + (0.64);
            p2.at(k,1,2) = A51.at(1,0)*p2.at(4,0,2) + A51.at(1,1)*p2.at(4,1,2) + A51.at(1,2)*p2.at(4,2,2) + (0.0);
            p2.at(k,2,0) = A51.at(2,0)*p2.at(4,0,0) + A51.at(2,1)*p2.at(4,1,0) + A51.at(2,2)*p2.at(4,2,0) + (0.0);
            p2.at(k,2,1) = A51.at(2,0)*p2.at(4,0,1) + A51.at(2,1)*p2.at(4,1,1) + A51.at(2,2)*p2.at(4,2,1) + (0.0);
            p2.at(k,2,2) = A51.at(2,0)*p2.at(4,0,2) + A51.at(2,1)*p2.at(4,1,2) + A51.at(2,2)*p2.at(4,2,2) + (0.64);

            //p1
            p1.at(k,0,0) = A51.at(0,0)*p1.at(4,0,0) + A51.at(0,1)*p1.at(4,1,0) + A51.at(0,2)*p1.at(4,2,0) + (-0.12);
            p1.at(k,0,1) = A51.at(0,0)*p1.at(4,0,1) + A51.at(0,1)*p1.at(4,1,1) + A51.at(0,2)*p1.at(4,2,1) + (0.0);
            p1.at(k,0,2) = A51.at(0,0)*p1.at(4,0,2) + A51.at(0,1)*p1.at(4,1,2) + A51.at(0,2)*p1.at(4,2,2) + (0.0);
            p1.at(k,1,0) = A51.at(1,0)*p1.at(4,0,0) + A51.at(1,1)*p1.at(4,1,0) + A51.at(1,2)*p1.at(4,2,0) + (0.0);
            p1.at(k,1,1) = A51.at(1,0)*p1.at(4,0,1) + A51.at(1,1)*p1.at(4,1,1) + A51.at(1,2)*p1.at(4,2,1) + (-0.12);
            p1.at(k,1,2) = A51.at(1,0)*p1.at(4,0,2) + A51.at(1,1)*p1.at(4,1,2) + A51.at(1,2)*p1.at(4,2,2) + (0.0);
            p1.at(k,2,0) = A51.at(2,0)*p1.at(4,0,0) + A51.at(2,1)*p1.at(4,1,0) + A51.at(2,2)*p1.at(4,2,0) + (0.0);
            p1.at(k,2,1) = A51.at(2,0)*p1.at(4,0,1) + A51.at(2,1)*p1.at(4,1,1) + A51.at(2,2)*p1.at(4,2,1) + (0.0);
            p1.at(k,2,2) = A51.at(2,0)*p1.at(4,0,2) + A51.at(2,1)*p1.at(4,1,2) + A51.at(2,2)*p1.at(4,2,2) + (-0.12);

            //p0
            p0.at(k,0,0) = A51.at(0,0)*p0.at(4,0,0) + A51.at(0,1)*p0.at(4,1,0) + A51.at(0,2)*p0.at(4,2,0);
            p0.at(k,0,1) = A51.at(0,0)*p0.at(4,0,1) + A51.at(0,1)*p0.at(4,1,1) + A51.at(0,2)*p0.at(4,2,1);
            p0.at(k,0,2) = A51.at(0,0)*p0.at(4,0,2) + A51.at(0,1)*p0.at(4,1,2) + A51.at(0,2)*p0.at(4,2,2);
            p0.at(k,1,0) = A51.at(1,0)*p0.at(4,0,0) + A51.at(1,1)*p0.at(4,1,0) + A51.at(1,2)*p0.at(4,2,0);
            p0.at(k,1,1) = A51.at(1,0)*p0.at(4,0,1) + A51.at(1,1)*p0.at(4,1,1) + A51.at(1,2)*p0.at(4,2,1);
            p0.at(k,1,2) = A51.at(1,0)*p0.at(4,0,2) + A51.at(1,1)*p0.at(4,1,2) + A51.at(1,2)*p0.at(4,2,2);
            p0.at(k,2,0) = A51.at(2,0)*p0.at(4,0,0) + A51.at(2,1)*p0.at(4,1,0) + A51.at(2,2)*p0.at(4,2,0);
            p0.at(k,2,1) = A51.at(2,0)*p0.at(4,0,1) + A51.at(2,1)*p0.at(4,1,1) + A51.at(2,2)*p0.at(4,2,1);
            p0.at(k,2,2) = A51.at(2,0)*p0.at(4,0,2) + A51.at(2,1)*p0.at(4,1,2) + A51.at(2,2)*p0.at(4,2,2);

            //q0
            q0.at(k,0) = A51.at(0,0)*q0.at(4,0) + A51.at(0,1)*q0.at(4,1) + A51.at(0,2)*q0.at(4,2) + 0.48*h*b(1,k1);
            q0.at(k,1) = A51.at(1,0)*q0.at(4,0) + A51.at(1,1)*q0.at(4,1) + A51.at(1,2)*q0.at(4,2) + 0.48*h*b(2,k1);
            q0.at(k,2) = A51.at(2,0)*q0.at(4,0) + A51.at(2,1)*q0.at(4,1) + A51.at(2,2)*q0.at(4,2) + 0.48*h*b(3,k1);

            x1.at(k) = (p3.at(k,0,0)*x1.at(3)+p3.at(k,0,1)*x2.at(3)+p3.at(k,0,2)*x3.at(3))
                     + (p2.at(k,0,0)*x1.at(2)+p2.at(k,0,1)*x2.at(2)+p2.at(k,0,2)*x3.at(2))
                     + (p1.at(k,0,0)*x1.at(1)+p1.at(k,0,1)*x2.at(1)+p1.at(k,0,2)*x3.at(1))
                     + (p0.at(k,0,0)*x1.at(0)+p0.at(k,0,1)*x2.at(0)+p0.at(k,0,2)*x3.at(0)) + q0.at(k,0);
            x2.at(k) = (p3.at(k,1,0)*x1.at(3)+p3.at(k,1,1)*x2.at(3)+p3.at(k,1,2)*x3.at(3))
                     + (p2.at(k,1,0)*x1.at(2)+p2.at(k,1,1)*x2.at(2)+p2.at(k,1,2)*x3.at(2))
                     + (p1.at(k,1,0)*x1.at(1)+p1.at(k,1,1)*x2.at(1)+p1.at(k,1,2)*x3.at(1))
                     + (p0.at(k,1,0)*x1.at(0)+p0.at(k,1,1)*x2.at(0)+p0.at(k,1,2)*x3.at(0)) + q0.at(k,1);
            x3.at(k) = (p3.at(k,2,0)*x1.at(3)+p3.at(k,2,1)*x2.at(3)+p3.at(k,2,2)*x3.at(3))
                     + (p2.at(k,2,0)*x1.at(2)+p2.at(k,2,1)*x2.at(2)+p2.at(k,2,2)*x3.at(2))
                     + (p1.at(k,2,0)*x1.at(1)+p1.at(k,2,1)*x2.at(1)+p1.at(k,2,2)*x3.at(1))
                     + (p0.at(k,2,0)*x1.at(0)+p0.at(k,2,1)*x2.at(0)+p0.at(k,2,2)*x3.at(0)) + q0.at(k,2);
        }

        if (k==6)
        {
            unsigned int k1 = k-1;

            DoubleMatrix A51(3,3);
            A51.at(0,0) = 0.48*h*a(1,1,k1)+1.92; A51.at(0,1) = 0.48*h*a(1,2,k1);      A51.at(0,2) = 0.48*h*a(1,3,k1);
            A51.at(1,0) = 0.48*h*a(2,1,k1);      A51.at(1,1) = 0.48*h*a(2,2,k1)+1.92; A51.at(1,2) = 0.48*h*a(2,3,k1);
            A51.at(2,0) = 0.48*h*a(3,1,k1);      A51.at(2,1) = 0.48*h*a(3,2,k1);      A51.at(2,2) = 0.48*h*a(3,3,k1)+1.92;

            //p3
            p3.at(k,0,0) = A51.at(0,0)*p3.at(4,0,0) + A51.at(0,1)*p3.at(4,1,0) + A51.at(0,2)*p3.at(4,2,0) + (-1.44);
            p3.at(k,0,1) = A51.at(0,0)*p3.at(4,0,1) + A51.at(0,1)*p3.at(4,1,1) + A51.at(0,2)*p3.at(4,2,1) + (0.0);
            p3.at(k,0,2) = A51.at(0,0)*p3.at(4,0,2) + A51.at(0,1)*p3.at(4,1,2) + A51.at(0,2)*p3.at(4,2,2) + (0.0);
            p3.at(k,1,0) = A51.at(1,0)*p3.at(4,0,0) + A51.at(1,1)*p3.at(4,1,0) + A51.at(1,2)*p3.at(4,2,0) + (0.0);
            p3.at(k,1,1) = A51.at(1,0)*p3.at(4,0,1) + A51.at(1,1)*p3.at(4,1,1) + A51.at(1,2)*p3.at(4,2,1) + (-1.44);
            p3.at(k,1,2) = A51.at(1,0)*p3.at(4,0,2) + A51.at(1,1)*p3.at(4,1,2) + A51.at(1,2)*p3.at(4,2,2) + (0.0);
            p3.at(k,2,0) = A51.at(2,0)*p3.at(4,0,0) + A51.at(2,1)*p3.at(4,1,0) + A51.at(2,2)*p3.at(4,2,0) + (0.0);
            p3.at(k,2,1) = A51.at(2,0)*p3.at(4,0,1) + A51.at(2,1)*p3.at(4,1,1) + A51.at(2,2)*p3.at(4,2,1) + (0.0);
            p3.at(k,2,2) = A51.at(2,0)*p3.at(4,0,2) + A51.at(2,1)*p3.at(4,1,2) + A51.at(2,2)*p3.at(4,2,2) + (-1.44);

            //p2
            p2.at(k,0,0) = A51.at(0,0)*p2.at(4,0,0) + A51.at(0,1)*p2.at(4,1,0) + A51.at(0,2)*p2.at(4,2,0) + (0.64);
            p2.at(k,0,1) = A51.at(0,0)*p2.at(4,0,1) + A51.at(0,1)*p2.at(4,1,1) + A51.at(0,2)*p2.at(4,2,1) + (0.0);
            p2.at(k,0,2) = A51.at(0,0)*p2.at(4,0,2) + A51.at(0,1)*p2.at(4,1,2) + A51.at(0,2)*p2.at(4,2,2) + (0.0);
            p2.at(k,1,0) = A51.at(1,0)*p2.at(4,0,0) + A51.at(1,1)*p2.at(4,1,0) + A51.at(1,2)*p2.at(4,2,0) + (0.0);
            p2.at(k,1,1) = A51.at(1,0)*p2.at(4,0,1) + A51.at(1,1)*p2.at(4,1,1) + A51.at(1,2)*p2.at(4,2,1) + (0.64);
            p2.at(k,1,2) = A51.at(1,0)*p2.at(4,0,2) + A51.at(1,1)*p2.at(4,1,2) + A51.at(1,2)*p2.at(4,2,2) + (0.0);
            p2.at(k,2,0) = A51.at(2,0)*p2.at(4,0,0) + A51.at(2,1)*p2.at(4,1,0) + A51.at(2,2)*p2.at(4,2,0) + (0.0);
            p2.at(k,2,1) = A51.at(2,0)*p2.at(4,0,1) + A51.at(2,1)*p2.at(4,1,1) + A51.at(2,2)*p2.at(4,2,1) + (0.0);
            p2.at(k,2,2) = A51.at(2,0)*p2.at(4,0,2) + A51.at(2,1)*p2.at(4,1,2) + A51.at(2,2)*p2.at(4,2,2) + (0.64);

            //p1
            p1.at(k,0,0) = A51.at(0,0)*p1.at(4,0,0) + A51.at(0,1)*p1.at(4,1,0) + A51.at(0,2)*p1.at(4,2,0) + (-0.12);
            p1.at(k,0,1) = A51.at(0,0)*p1.at(4,0,1) + A51.at(0,1)*p1.at(4,1,1) + A51.at(0,2)*p1.at(4,2,1) + (0.0);
            p1.at(k,0,2) = A51.at(0,0)*p1.at(4,0,2) + A51.at(0,1)*p1.at(4,1,2) + A51.at(0,2)*p1.at(4,2,2) + (0.0);
            p1.at(k,1,0) = A51.at(1,0)*p1.at(4,0,0) + A51.at(1,1)*p1.at(4,1,0) + A51.at(1,2)*p1.at(4,2,0) + (0.0);
            p1.at(k,1,1) = A51.at(1,0)*p1.at(4,0,1) + A51.at(1,1)*p1.at(4,1,1) + A51.at(1,2)*p1.at(4,2,1) + (-0.12);
            p1.at(k,1,2) = A51.at(1,0)*p1.at(4,0,2) + A51.at(1,1)*p1.at(4,1,2) + A51.at(1,2)*p1.at(4,2,2) + (0.0);
            p1.at(k,2,0) = A51.at(2,0)*p1.at(4,0,0) + A51.at(2,1)*p1.at(4,1,0) + A51.at(2,2)*p1.at(4,2,0) + (0.0);
            p1.at(k,2,1) = A51.at(2,0)*p1.at(4,0,1) + A51.at(2,1)*p1.at(4,1,1) + A51.at(2,2)*p1.at(4,2,1) + (0.0);
            p1.at(k,2,2) = A51.at(2,0)*p1.at(4,0,2) + A51.at(2,1)*p1.at(4,1,2) + A51.at(2,2)*p1.at(4,2,2) + (-0.12);

            //p0
            p0.at(k,0,0) = A51.at(0,0)*p0.at(4,0,0) + A51.at(0,1)*p0.at(4,1,0) + A51.at(0,2)*p0.at(4,2,0);
            p0.at(k,0,1) = A51.at(0,0)*p0.at(4,0,1) + A51.at(0,1)*p0.at(4,1,1) + A51.at(0,2)*p0.at(4,2,1);
            p0.at(k,0,2) = A51.at(0,0)*p0.at(4,0,2) + A51.at(0,1)*p0.at(4,1,2) + A51.at(0,2)*p0.at(4,2,2);
            p0.at(k,1,0) = A51.at(1,0)*p0.at(4,0,0) + A51.at(1,1)*p0.at(4,1,0) + A51.at(1,2)*p0.at(4,2,0);
            p0.at(k,1,1) = A51.at(1,0)*p0.at(4,0,1) + A51.at(1,1)*p0.at(4,1,1) + A51.at(1,2)*p0.at(4,2,1);
            p0.at(k,1,2) = A51.at(1,0)*p0.at(4,0,2) + A51.at(1,1)*p0.at(4,1,2) + A51.at(1,2)*p0.at(4,2,2);
            p0.at(k,2,0) = A51.at(2,0)*p0.at(4,0,0) + A51.at(2,1)*p0.at(4,1,0) + A51.at(2,2)*p0.at(4,2,0);
            p0.at(k,2,1) = A51.at(2,0)*p0.at(4,0,1) + A51.at(2,1)*p0.at(4,1,1) + A51.at(2,2)*p0.at(4,2,1);
            p0.at(k,2,2) = A51.at(2,0)*p0.at(4,0,2) + A51.at(2,1)*p0.at(4,1,2) + A51.at(2,2)*p0.at(4,2,2);

            //q0
            q0.at(k,0) = A51.at(0,0)*q0.at(4,0) + A51.at(0,1)*q0.at(4,1) + A51.at(0,2)*q0.at(4,2) + 0.48*h*b(1,k1);
            q0.at(k,1) = A51.at(1,0)*q0.at(4,0) + A51.at(1,1)*q0.at(4,1) + A51.at(1,2)*q0.at(4,2) + 0.48*h*b(2,k1);
            q0.at(k,2) = A51.at(2,0)*q0.at(4,0) + A51.at(2,1)*q0.at(4,1) + A51.at(2,2)*q0.at(4,2) + 0.48*h*b(3,k1);

            x1.at(k) = (p3.at(k,0,0)*x1.at(3)+p3.at(k,0,1)*x2.at(3)+p3.at(k,0,2)*x3.at(3))
                     + (p2.at(k,0,0)*x1.at(2)+p2.at(k,0,1)*x2.at(2)+p2.at(k,0,2)*x3.at(2))
                     + (p1.at(k,0,0)*x1.at(1)+p1.at(k,0,1)*x2.at(1)+p1.at(k,0,2)*x3.at(1))
                     + (p0.at(k,0,0)*x1.at(0)+p0.at(k,0,1)*x2.at(0)+p0.at(k,0,2)*x3.at(0)) + q0.at(k,0);
            x2.at(k) = (p3.at(k,1,0)*x1.at(3)+p3.at(k,1,1)*x2.at(3)+p3.at(k,1,2)*x3.at(3))
                     + (p2.at(k,1,0)*x1.at(2)+p2.at(k,1,1)*x2.at(2)+p2.at(k,1,2)*x3.at(2))
                     + (p1.at(k,1,0)*x1.at(1)+p1.at(k,1,1)*x2.at(1)+p1.at(k,1,2)*x3.at(1))
                     + (p0.at(k,1,0)*x1.at(0)+p0.at(k,1,1)*x2.at(0)+p0.at(k,1,2)*x3.at(0)) + q0.at(k,1);
            x3.at(k) = (p3.at(k,2,0)*x1.at(3)+p3.at(k,2,1)*x2.at(3)+p3.at(k,2,2)*x3.at(3))
                     + (p2.at(k,2,0)*x1.at(2)+p2.at(k,2,1)*x2.at(2)+p2.at(k,2,2)*x3.at(2))
                     + (p1.at(k,2,0)*x1.at(1)+p1.at(k,2,1)*x2.at(1)+p1.at(k,2,2)*x3.at(1))
                     + (p0.at(k,2,0)*x1.at(0)+p0.at(k,2,1)*x2.at(0)+p0.at(k,2,2)*x3.at(0)) + q0.at(k,2);
        }
    }

    printf("%.10f %.10f %.10f %.10f %.10f %.10f\n", x01.at(0), x01.at(1), x01.at(2), x01.at(3), x01.at(4), x01.at(5));
    printf("%.10f %.10f %.10f %.10f %.10f %.10f\n", x02.at(0), x02.at(1), x02.at(2), x02.at(3), x02.at(4), x02.at(5));
    printf("%.10f %.10f %.10f %.10f %.10f %.10f\n", x03.at(0), x03.at(1), x03.at(2), x03.at(3), x03.at(4), x03.at(5));
    printf("%.10f %.10f %.10f %.10f %.10f %.10f\n", x1.at(0), x1.at(1), x1.at(2), x1.at(3), x1.at(4), x1.at(5));
    printf("%.10f %.10f %.10f %.10f %.10f %.10f\n", x2.at(0), x2.at(1), x2.at(2), x2.at(3), x2.at(4), x2.at(5));
    printf("%.10f %.10f %.10f %.10f %.10f %.10f\n", x3.at(0), x3.at(1), x3.at(2), x3.at(3), x3.at(4), x3.at(5));
}

void Example4::calculate2()
{
    h = 0.00001;
    N = 100000;

    DoubleVector x01(N+1);
    DoubleVector x02(N+1);
    DoubleVector x03(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x01.at(i) = X1(i);
        x02.at(i) = X2(i);
        x03.at(i) = X3(i);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    IPrinter::printVector(14,10,x03);

    DoubleVector x1(N+1);
    DoubleVector x2(N+2);
    DoubleVector x3(N+2);

    x1.at(N-0) = X1(N-0); x2.at(N-0) = X2(N-0); x3.at(N-0) = X3(N-0);
    x1.at(N-1) = X1(N-1); x2.at(N-1) = X2(N-1); x3.at(N-1) = X3(N-1);
    x1.at(N-2) = X1(N-2); x2.at(N-2) = X2(N-2); x3.at(N-2) = X3(N-2);
    x1.at(N-3) = X1(N-3); x2.at(N-3) = X2(N-3); x3.at(N-3) = X3(N-3);

//    x1.at(0) = X1(0); x2.at(0) = X2(0); x3.at(0) = X3(0);
//    x1.at(1) = X1(1); x2.at(1) = X2(1); x3.at(1) = X3(1);
//    x1.at(2) = X1(2); x2.at(2) = X2(2); x3.at(2) = X3(2);
//    x1.at(3) = X1(3); x2.at(3) = X2(3); x3.at(3) = X3(3);

    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+1;
        double alpha1 = +1.92;//+48.0/25.0;
        double alpha2 = -1.44;//-36.0/25.0;
        double alpha3 = +0.64;//+16.0/25.0;
        double alpha4 = -0.12;//-3.0/25.0;
        double alpha5 = -0.48*h;//-12.0/25.0;

        x1.at(k) = alpha1*x1.at(k+1) + alpha2*x1.at(k+2) + alpha3*x1.at(k+3) + alpha4*x1.at(k+4)
                + (alpha5*a(1,1,k1))*x1.at(k1) + (alpha5*a(1,2,k1))*x2.at(k1) + (alpha5*a(1,3,k1))*x3.at(k1)
                + (alpha5*b(1,k1));
        x2.at(k) = alpha1*x2.at(k+1) + alpha2*x2.at(k+2) + alpha3*x2.at(k+3) + alpha4*x2.at(k+4)
                + (alpha5*a(2,1,k1))*x1.at(k1) + (alpha5*a(2,2,k1))*x2.at(k1) + (alpha5*a(2,3,k1))*x3.at(k1)
                + (alpha5*b(2,k1));
        x3.at(k) = alpha1*x3.at(k+1) + alpha2*x3.at(k+2) + alpha3*x3.at(k+3) + alpha4*x3.at(k+4)
                + (alpha5*a(3,1,k1))*x1.at(k1) + (alpha5*a(3,2,k1))*x2.at(k1) + (alpha5*a(3,3,k1))*x3.at(k1)
                + (alpha5*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+1;
        x1.at(k) = (48.0/25.0)*x1.at(k+1) + (-36.0/25.0)*x1.at(k+2) + (16.0/25.0)*x1.at(k+3) + (-3.0/25.0)*x1.at(k+4)
                + ((-12.0/25.0)*h*a(1,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(1,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(1,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(1,k1));
        x2.at(k) = (48.0/25.0)*x2.at(k+1) + (-36.0/25.0)*x2.at(k+2) + (16.0/25.0)*x2.at(k+3) + (-3.0/25.0)*x2.at(k+4)
                + ((-12.0/25.0)*h*a(2,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(2,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(2,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(2,k1));
        x3.at(k) = (48.0/25.0)*x3.at(k+1) + (-36.0/25.0)*x3.at(k+2) + (16.0/25.0)*x3.at(k+3) + (-3.0/25.0)*x3.at(k+4)
                + ((-12.0/25.0)*h*a(3,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(3,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(3,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+2;
        x1.at(k) = (48.0/25.0)*x1.at(k+1) + (-36.0/25.0)*x1.at(k+2) + (16.0/25.0)*x1.at(k+3) + (-3.0/25.0)*x1.at(k+4)
                + ((-12.0/25.0)*h*a(1,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(1,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(1,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(1,k1));
        x2.at(k) = (48.0/25.0)*x2.at(k+1) + (-36.0/25.0)*x2.at(k+2) + (16.0/25.0)*x2.at(k+3) + (-3.0/25.0)*x2.at(k+4)
                + ((-12.0/25.0)*h*a(2,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(2,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(2,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(2,k1));
        x3.at(k) = (48.0/25.0)*x3.at(k+1) + (-36.0/25.0)*x3.at(k+2) + (16.0/25.0)*x3.at(k+3) + (-3.0/25.0)*x3.at(k+4)
                + ((-12.0/25.0)*h*a(3,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(3,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(3,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+3;
        x1.at(k) = (48.0/25.0)*x1.at(k+1) + (-36.0/25.0)*x1.at(k+2) + (16.0/25.0)*x1.at(k+3) + (-3.0/25.0)*x1.at(k+4)
                + ((-12.0/25.0)*h*a(1,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(1,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(1,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(1,k1));
        x2.at(k) = (48.0/25.0)*x2.at(k+1) + (-36.0/25.0)*x2.at(k+2) + (16.0/25.0)*x2.at(k+3) + (-3.0/25.0)*x2.at(k+4)
                + ((-12.0/25.0)*h*a(2,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(2,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(2,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(2,k1));
        x3.at(k) = (48.0/25.0)*x3.at(k+1) + (-36.0/25.0)*x3.at(k+2) + (16.0/25.0)*x3.at(k+3) + (-3.0/25.0)*x3.at(k+4)
                + ((-12.0/25.0)*h*a(3,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(3,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(3,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
}

double Example4::X1(unsigned int k) const
{
    double t = k*h;
    return sin(2.0*t) + t*t;
}

double Example4::X2(unsigned int k) const
{
    double t = k*h;
    return 3.0*t;
}

double Example4::X3(unsigned int k) const
{
    double t = k*h;
    return cos(2.0*t) - sin(t);
}

double Example4::a(unsigned int i, unsigned int j, unsigned int k) const
{
    double t = k*h;

    if (i==1 && j==1) return +3.0;
    if (i==1 && j==2) return -t;
    if (i==1 && j==3) return +2.0;

    if (i==2 && j==1) return +1.0;
    if (i==2 && j==2) return +1.0;
    if (i==2 && j==3) return +3.0;

    if (i==3 && j==1) return -2.0;
    if (i==3 && j==2) return +t;
    if (i==3 && j==3) return +1.0;

    return 0.0;
}

double Example4::b(unsigned int i, unsigned int k) const
{
    double t = k*h;
    if (i==1) return 2.0*t + 2.0*sin(t) - 3.0*sin(2.0*t);
    if (i==2) return 3.0*sin(t) - sin(2.0*t) - 3.0*cos(2.0*t) - t*t - 3.0*t + 3.0;
    if (i==3) return sin(t) - cos(t) - cos(2.0*t) - t*t;
    return 0.0;
}
