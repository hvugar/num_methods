#include "borderparabolic.h"
#include <math.h>

#define NO_NORMALIZE

void BorderParabolic::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderParabolic bp;

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
        bp.hx = 0.001;
        bp.ht = 0.001;
        bp.N = 1000;
        bp.M = 1000;

        DoubleMatrix u1;
        bp.calculateU(u1, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u1, 10, 10, NULL);
        u1.clear();
        IPrinter::printSeperatorLine();
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.001;
        bp.N = 100;
        bp.M = 1000;
        DoubleMatrix u2;
        bp.calculateN4L2RM(u2);
        IPrinter::printMatrix(14, 10, u2, 10, 10, NULL);
        IPrinter::printSeperatorLine();
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.001;
        bp.N = 100;
        bp.M = 1000;
        DoubleMatrix u2;
        bp.calculateN4R2LM(u2);
        IPrinter::printMatrix(14, 10, u2, 10, 10, NULL);
        IPrinter::printSeperatorLine();
    }
}

double BorderParabolic::initial(unsigned int i UNUSED_PARAM) const
{
    return u(i,0);
}

double BorderParabolic::boundary(Boundary type UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    if (type == Left)  return u(0,j);
    if (type == Right) return u(N,j);
    return 0.0;
}

double BorderParabolic::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
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

double BorderParabolic::u(unsigned int i, unsigned int j) const
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

void BorderParabolic::calculateN4L2RM(DoubleMatrix &u1)
{
    u1.resize(M+1, N+1);
    /* initial condition */
    for (unsigned int i=0; i<=N; i++) u1.at(0,i) = initial(i);
    /* border conditions */
    for (unsigned int j=1; j<=M; j++)
    {
        u1.at(j,0) = boundary(Left, j);
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

    double D42[5][5] =
    {
        {+70.0, -208.0, +228.0, -112.0, +22.0},
        {+22.0, -40.0,  +12.0,  +8.0,   -2.0},
        {-2.0,  +32.0,  -60.0,  +32.0,  -2.0},
        {-2.0,  +8.0,   +12.0,  -40.0,  +22.0},
        {+22.0, -112.0, +228.0, -208.0, +70.0}
    };

    double alpha = (ht*a*a)/(24.0*hx*hx);

    DoubleMatrix A1(4,4);
    DoubleVector b1(4);
    DoubleVector x(4);

    DoubleMatrix ems(N-4, 4);
    for (unsigned int m=1; m<=M; m++)
    {
        A1.at(0,0) = (D42[1][1]*alpha - 1.0);
        A1.at(0,1) = (D42[1][2]*alpha)/A1.at(0,0);
        A1.at(0,2) = (D42[1][3]*alpha)/A1.at(0,0);
        A1.at(0,3) = (D42[1][4]*alpha)/A1.at(0,0);
        b1.at(0)   = (-u1.at(m-1,1) - (D42[1][0]*alpha)*u1.at(m,0) - ht*f(1,m))/A1.at(0,0);
        A1.at(0,0) = 1.0;

        ems.at(0,0) = A1.at(0,1);
        ems.at(0,1) = A1.at(0,2);
        ems.at(0,2) = A1.at(0,3);
        ems.at(0,3) = b1.at(0);

        // + * * * *
        for (unsigned int n=1; n<=N-5; n++)
        {
            double g1 = D42[0][0]*alpha-1.0;
            double g2 = D42[0][1]*alpha;
            double g3 = D42[0][2]*alpha;
            double g4 = D42[0][3]*alpha;
            double g5 = D42[0][4]*alpha;
            double g0 = u1.at(m-1,n) + ht*f(n,m);

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            g0 /= -g1;
            g1 = 1.0;

            A1.at(0,0) = (A1.at(0,1) + g2);
            A1.at(0,1) = (A1.at(0,2) + g3)/A1.at(0,0);
            A1.at(0,2) = (A1.at(0,3) + g4)/A1.at(0,0);
            A1.at(0,3) = g5/A1.at(0,0);
            b1.at(0)   = (b1.at(0) - g0)/A1.at(0,0);
            A1.at(0,0) = 1.0;

            ems.at(n,0) = A1.at(0,1);
            ems.at(n,1) = A1.at(0,2);
            ems.at(n,2) = A1.at(0,3);
            ems.at(n,3) = b1.at(0);
        }

        A1[1][0] = D42[1][0]*alpha;
        A1[1][1] = D42[1][1]*alpha - 1.0;
        A1[1][2] = D42[1][2]*alpha;
        A1[1][3] = D42[1][3]*alpha;
        b1[1]    = -u1.at(m-1,N-3) - (D42[1][4]*alpha)*u1.at(m,N) - ht*f(N-3,m);

        A1[2][0] = D42[2][0]*alpha;
        A1[2][1] = D42[2][1]*alpha;
        A1[2][2] = D42[2][2]*alpha - 1.0;
        A1[2][3] = D42[2][3]*alpha;
        b1[2]    = -u1.at(m-1,N-2) - (D42[2][4]*alpha)*u1.at(m,N) - ht*f(N-2,m);

        A1[3][0] = D42[3][0]*alpha;
        A1[3][1] = D42[3][1]*alpha;
        A1[3][2] = D42[3][2]*alpha;
        A1[3][3] = D42[3][3]*alpha - 1.0;
        b1[3]    = -u1.at(m-1,N-1) - (D42[3][4]*alpha)*u1.at(m,N) - ht*f(N-1,m);

        GaussianElimination(A1, b1, x);

        u1.at(m, N-1) = x.at(3);
        u1.at(m, N-2) = x.at(2);
        u1.at(m, N-3) = x.at(1);
        u1.at(m, N-4) = x.at(0);

        for (unsigned int i=N-5; i>=1; i--)
        {
            u1.at(m,i) = -ems.at(i-1,0)*u1.at(m,i+1) - ems.at(i-1,1)*u1.at(m,i+2) - ems.at(i-1,2)*u1.at(m,i+3) + ems.at(i-1,3);
        }

        if (m==0)
        {
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u(N-4,m), u(N-3,m), u(N-2,m), u(N-1,m));
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, x[0], x[1], x[2], x[3]);
            //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u1.at(m,N-4), u1.at(m,N-3), u1.at(m,N-2), u1.at(m,N-1));
            IPrinter::printVector(14,10, u1.row(m));
            DoubleVector v(N+1); for (unsigned int n=0; n<=N; n++) v.at(n) = u(n,m);
            IPrinter::printVector(14,10, v);
            v.clear();
            IPrinter::printSeperatorLine();
            break;
        }
    }
    ems.clear();
    A1.clear();
    b1.clear();
    x.clear();
}

void BorderParabolic::calculateN4R2LM(DoubleMatrix &u1)
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
