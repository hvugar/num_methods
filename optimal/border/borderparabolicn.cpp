#include "borderparabolicn.h"
#include <math.h>

void BorderParabolicN::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderParabolicN bp;

    // Real solution
    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix ru(bp.M+1,bp.N+1);
        for (unsigned int i=0; i<=bp.M; i++)
        {
            for (unsigned int j=0; j<=bp.N; j++)
                ru.at(i,j) = bp.u(j,i);
        }
        IPrinter::printMatrix(14, 10, ru, 10, 10, NULL);
        ru.clear();
        IPrinter::printSeperatorLine();
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;

        DoubleMatrix u1;
        bp.calculateN(u1, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u1, 10, 10, NULL);
        u1.clear();
        IPrinter::printSeperatorLine();
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u2;
        bp.calculateN4L2RM(u2);
        IPrinter::printMatrix(14, 10, u2, 10, 10, NULL);
        u2.clear();
        IPrinter::printSeperatorLine();
    }

    //    {
    //        bp.hx = 0.01;
    //        bp.ht = 0.001;
    //        bp.N = 100;
    //        bp.M = 1000;
    //        DoubleMatrix u2;
    //        bp.calculateN4R2LM(u2);
    //        IPrinter::printMatrix(14, 10, u2, 10, 10, NULL);
    //        IPrinter::printSeperatorLine();
    //    }
}

double BorderParabolicN::initial(unsigned int i UNUSED_PARAM) const
{
    return u(i,0);
}

double BorderParabolicN::boundary(Boundary type UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t UNUSED_PARAM = j*ht;

#ifdef SAMPLE_1
    if (type == Left)  return 0.0;
    if (type == Right) return 2.0;
#endif
#ifdef SAMPLE_2
    if (type == Left)  return t;
    if (type == Right) return sin(20.0) + 20.0*cos(20.0) + t;
#endif
#ifdef SAMPLE_3
    if (type == Left)  return 12.0*t;
    if (type == Right) return 10.0*t*cos(10.0) + 2.0*t*exp(2.0*t);
#endif
    return 0.0;
}

double BorderParabolicN::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t UNUSED_PARAM = j*ht;
    double x UNUSED_PARAM = i*hx;
#ifdef SAMPLE_1
    return 2.0*t - 2.0*a*a;
#endif
#ifdef SAMPLE_2
    return x - a*a*(40.0*cos(20.0*x) - 400.0*x*sin(20.0*x));
#endif
#ifdef SAMPLE_3
    return (sin(10.0*x) + 2.0*x*exp(2.0*x*t) - a*a*(4.0*t*t*exp(2.0*x*t) - 100.0*t*sin(10.0*x)));
#endif
}

double BorderParabolicN::u(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
#ifdef SAMPLE_1
    return x*x + t*t;
#endif
#ifdef SAMPLE_2
    return x*sin(20.0*x) + t*x;
#endif
#ifdef SAMPLE_3
    return t*sin(10.0*x) + exp(2.0*x*t);
#endif
}

void BorderParabolicN::calculateN4L2RM(DoubleMatrix &u1)
{
    u1.resize(M+1, N+1);

    for (unsigned int i=0; i<=N; i++) u1.at(0,i) = initial(i);

    DoubleMatrix A(5,5,0.0);
    DoubleVector b(5,0.0);
    DoubleVector y(5,0.0);

    double beta1 = (a*a*ht)/(12.0*hx*hx);
    double beta2 = (a*a*ht)/(24.0*hx*hx);
    double gamma  = (a*a*ht)/hx;

    DoubleMatrix ems(N-3, 5);
    for (unsigned int m=1; m<=M; m++)
    {
        //0
        A.at(0,0) = -3.0*beta1 - 1.0;
        A.at(0,1) = -10.0*beta1;
        A.at(0,2) = +18.0*beta1;
        A.at(0,3) = -6.0*beta1;
        A.at(0,4) = +beta1;
        b.at(0)   = -u1.at(m-1,0) - ht*f(0,m) + gamma*boundary(Left,m);

        A.at(0,1) /= A.at(0,0);
        A.at(0,2) /= A.at(0,0);
        A.at(0,3) /= A.at(0,0);
        A.at(0,4) /= A.at(0,0);
        b.at(0)   /= A.at(0,0);
        A.at(0,0)  = 1.0;

        ems.at(0,0) = A.at(0,1);
        ems.at(0,1) = A.at(0,2);
        ems.at(0,2) = A.at(0,3);
        ems.at(0,3) = A.at(0,4);
        ems.at(0,4) = b.at(0);

        for (unsigned int n=0; n<=N-5; n++)
        {
            double g0 = +70.0*beta2 - 1.0;
            double g1 = -208.0*beta2;
            double g2 = +228.0*beta2;
            double g3 = -112.0*beta2;
            double g4 = +22.0*beta2;
            double fi = -u1.at(m-1,n) - ht*f(n,m);

            g1 /= -g0;
            g2 /= -g0;
            g3 /= -g0;
            g4 /= -g0;
            fi /= +g0;
            g0  = 1.0;

            double a00 = A.at(0,0);
            A.at(0,0) = A.at(0,1) + a00 * g1;
            A.at(0,1) = A.at(0,2) + a00 * g2;
            A.at(0,2) = A.at(0,3) + a00 * g3;
            A.at(0,3) = A.at(0,4) + a00 * g4;
            A.at(0,4) = 0.0;
            b.at(0)  = b.at(0) - a00*fi;

            A.at(0,1) /= A.at(0,0);
            A.at(0,2) /= A.at(0,0);
            A.at(0,3) /= A.at(0,0);
            A.at(0,4) /= A.at(0,0);
            b.at(0)   /= A.at(0,0);
            A.at(0,0)  = 1.0;

            ems.at(n+1,0) = A.at(0,1);
            ems.at(n+1,1) = A.at(0,2);
            ems.at(n+1,2) = A.at(0,3);
            ems.at(n+1,3) = A.at(0,4);
            ems.at(n+1,4) = b.at(0);
        }

        //N-3
        A.at(1,0) = +22.0*beta2;
        A.at(1,1) = -40.0*beta2 - 1.0;
        A.at(1,2) = +12.0*beta2;
        A.at(1,3) = +8.0*beta2;
        A.at(1,4) = -2.0*beta2;
        b.at(1)   = -u1.at(m-1,N-3) - ht*f(N-3,m);

        //N-2
        A.at(2,0) = -2.0*beta2;
        A.at(2,1) = +32.0*beta2;
        A.at(2,2) = -60.0*beta2 - 1.0;
        A.at(2,3) = +32.0*beta2;
        A.at(2,4) = -2.0*beta2;
        b.at(2)   = -u1.at(m-1,N-2) - ht*f(N-2,m);

        //N-1
        A.at(3,0) = -2.0*beta2;
        A.at(3,1) = +8.0*beta2;
        A.at(3,2) = +12.0*beta2;
        A.at(3,3) = -40.0*beta2 - 1.0;
        A.at(3,4) = +22.0*beta2;
        b.at(3)   = -u1.at(m-1,N-1) - ht*f(N-1,m);

        //N
        A.at(4,0) = -beta1;
        A.at(4,1) = +6*beta1;
        A.at(4,2) = -18.0*beta1;
        A.at(4,3) = +10.0*beta1;
        A.at(4,4) = +3.0*beta1 + 1.0;
        b.at(4)   = +u1.at(m-1,N) + ht*f(N,m) + gamma*boundary(Right,m);

        //IPrinter::print(A,A.rows(),A.cols(),18,10);

        GaussianElimination(A, b, y);
        //        IPrinter::print(y, y.size());

        u1.at(m, N-4) = y.at(0);
        u1.at(m, N-3) = y.at(1);
        u1.at(m, N-2) = y.at(2);
        u1.at(m, N-1) = y.at(3);
        u1.at(m, N-0) = y.at(4);
        //printf("%u, %14.10f %14.10f %14.10f %14.10f %14.10f\n", m, u(N-4,m), u(N-3,m), u(N-2,m), u(N-1,m), u(N,m));
        //printf("%u, %14.10f %14.10f %14.10f %14.10f %14.10f\n", m, u1.at(m,N-4), u1.at(m,N-3), u1.at(m,N-2), u1.at(m,N-1), u1.at(m,N));
        for (unsigned int n=N-5; n!=UINT32_MAX; n--)
        {
            u1.at(m,n) = -ems.at(n,0)*u1.at(m,n+1)
                    -ems.at(n,1)*u1.at(m,n+2)
                    -ems.at(n,2)*u1.at(m,n+3)
                    -ems.at(n,3)*u1.at(m,n+4)
                    +ems.at(n,4);
            //printf("%14.10f\n", u1.at(m,n));
        }

        //        IPrinter::printVector(u1.row(m));
        //        if (m>10) break;
    }
    IPrinter::printSeperatorLine();

    A.clear();
    b.clear();
    y.clear();
}

void BorderParabolicN::calculateN4R2LM(DoubleMatrix &u1)
{
    u1.resize(M+1, N+1);
    /* initial condition */
    for (unsigned int i=0; i<=N; i++) u1.at(0,i) = initial(i);
    /* border conditions */
    for (unsigned int j=1; j<=M; j++)
    {
        u1.at(j,0) = boundary(Left,j);
        u1.at(j,N) = boundary(Right, j);
    }

    //    double D41[5][5] =
    //    {
    //        {-25.0, +48.0, -36.0, +16.0, -3.0},
    //        {-3.0,  -10.0, +18.0, -6.0,  +1.0},
    //        {+1.0,  -8.0,  +0.0,  +8.0,  -1.0},
    //        {-1.0,  +6.0, -18.0,  +10.0, +3.0},
    //        {+3.0,  -16.0, +36.0, -48.0, +25.0}
    //    };

    //    double D42[5][5] =
    //    {
    //        {+70.0, -208.0, +228.0, -112.0, +22.0},
    //        {+22.0, -40.0,  +12.0,  +8.0,   -2.0},
    //        {-2.0,  +32.0,  -60.0,  +32.0,  -2.0},
    //        {-2.0,  +8.0,   +12.0,  -40.0,  +22.0},
    //        {+22.0, -112.0, +228.0, -208.0, +70.0}
    //    };

    double alpha = (a*a*ht)/(24.0*hx*hx);

    DoubleMatrix A1(4,4);
    DoubleVector b1(4);
    DoubleVector x(4);
    DoubleMatrix ems(N-4, 4);
    for (unsigned int m=1; m<=M; m++)
    {
        A1.at(0,0) = -40.0*alpha - 1.0;
        A1.at(0,1) = +12.0*alpha;
        A1.at(0,2) = +8.0*alpha;
        A1.at(0,3) = -2.0*alpha;
        b1.at(0)   = -u1.at(m-1,1) - (22.0*alpha)*u1.at(m,0) - ht*f(1,m);

        A1.at(1,0) = +32.0*alpha;
        A1.at(1,1) = -60.0*alpha - 1.0;
        A1.at(1,2) = +32.0*alpha;
        A1.at(1,3) = -2.0*alpha;
        b1.at(1)   = -u1.at(m-1,2) - (-2.0*alpha)*u1.at(m,0) - ht*f(2,m);

        A1[2][0] = +8.0*alpha;
        A1[2][1] = +12.0*alpha;
        A1[2][2] = -40.0*alpha - 1.0;
        A1[2][3] = +22.0*alpha;
        b1[2]    = -u1.at(m-1,3) - (-2.0*alpha)*u1.at(m,0) - ht*f(3,m);

        A1.at(3,0) = -2.0*alpha;
        A1.at(3,1) = +8.0*alpha;
        A1.at(3,2) = +12.0*alpha;
        A1.at(3,3) = -40.0*alpha - 1.0;
        b1.at(3)   = -u1.at(m-1,N-1) - (22.0*alpha)*u1.at(m,N) - ht*f(N-1,m);

        A1.at(3,0) /= A1.at(3,3);
        A1.at(3,1) /= A1.at(3,3);
        A1.at(3,2) /= A1.at(3,3);
        b1.at(3)   /= A1.at(3,3);
        A1.at(3,3) = 1.0;

        ems.at(N-5,0) = A1.at(3,0);
        ems.at(N-5,1) = A1.at(3,1);
        ems.at(N-5,2) = A1.at(3,2);
        ems.at(N-5,3) = b1.at(3);

        for (unsigned int n=N-1; n>=5; n--)
        {
            double g1 = +22.0*alpha;
            double g2 = -112.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -208.0*alpha;;
            double g5 = +70.0*alpha - 1.0;
            double g0  = -u1.at(m-1,n) - ht*f(n,m);

            g4 /= -g5;
            g3 /= -g5;
            g2 /= -g5;
            g1 /= -g5;
            g0 /= +g5;
            g5 = 1.0;

            A1.at(3,3) = A1.at(3,2) + g4;
            A1.at(3,2) = (A1.at(3,1) + g3)/A1.at(3,3);
            A1.at(3,1) = (A1.at(3,0) + g2)/A1.at(3,3);
            A1.at(3,0) = g1/A1.at(3,3);
            b1.at(3)   = (b1.at(3) - g0)/A1.at(3,3);
            A1.at(3,3) = 1.0;

            ems.at(n-5,0) = A1.at(3,0);
            ems.at(n-5,1) = A1.at(3,1);
            ems.at(n-5,2) = A1.at(3,2);
            ems.at(n-5,3) = b1.at(3);
        }

        GaussianElimination(A1, b1, x);

        u1.at(m, 1) = x.at(0);
        u1.at(m, 2) = x.at(1);
        u1.at(m, 3) = x.at(2);
        u1.at(m, 4) = x.at(3);
        for (unsigned int i=5; i<=N-1; i++)
        {
            u1.at(m,i) = -ems.at(i-4,2)*u1.at(m,i-1) - ems.at(i-4,1)*u1.at(m,i-2) - ems.at(i-4,0)*u1.at(m,i-3) + ems.at(i-4,3);
        }

        if (m==0)
        {
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u(1,m), u(2,m), u(3,m), u(4,m));
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, x[0], x[1], x[2], x[3]);
            //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u1.at(m,N-4), u1.at(m,N-3), u1.at(m,N-2), u1.at(m,N-1));
            IPrinter::printVector(14,10, u1.row(m));
            DoubleVector v(N+1); for (unsigned int n=0; n<=N; n++) v.at(n) = u(n,m);
            IPrinter::printVector(14,10, v);
            v.clear();
            IPrinter::printSeperatorLine();
        }
    }
    ems.clear();
    A1.clear();
    b1.clear();
    x.clear();
}
