#include "borderparabolic.h"
#include <math.h>

#define NO_NORMALIZE

void BorderParabolic::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderParabolic bp;

    // Real solution
    {
        DoubleMatrix ur(bp.M+1,bp.N+1);
        for (unsigned int i=0; i<=bp.M; i++)
        {
            for (unsigned int j=0; j<=bp.N; j++)
                ur.at(i,j) = bp.u(j,i);
        }
        IPrinter::printMatrix(14, 10, ur, 10, 10, NULL);
        ur.clear();
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
        bp.ht = 0.00001;
        bp.N = 100;
        bp.M = 100000;

        DoubleMatrix u2;
        bp.calculateN4L2R(u2);
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

void BorderParabolic::calculateN4L2R(DoubleMatrix &u1)
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

    double alpha = (ht*a*a)/(24.0*hx*hx);
    for (unsigned int m=1; m<=M; m++)
    {
        double c1 = 40.0*alpha + 1.0;
        double c2 = (-12.0*alpha)/c1;
        double c3 = (-8.0*alpha)/c1;
        double c4 = (+2.0*alpha)/c1;
        double d  = (u1.at(m-1,1) + (22.0*alpha)*u1.at(m,0) + ht*f(1,m))/c1;
        c1 = 1.0;

        // + * * * *
        for (unsigned int n=1; n<=N-5; n++)
        {
            double q1 = +70.0*alpha-1.0;
            double q2 = -208.0*alpha;
            double q3 = +228.0*alpha;
            double q4 = -112.0*alpha;
            double q5 = +22.0*alpha;
            double e  = u1.at(m-1,n) + ht*f(n,m);

            q2 /= -q1;
            q3 /= -q1;
            q4 /= -q1;
            q5 /= -q1;
            e  /= -q1;
            q1 = 1.0;

            c1 = c2 + q2;
            c2 = (c3 + q3)/c1;
            c3 = (c4 + q4)/c1;
            c4 = q5/c1;
            d  = (d - e)/c1;
            c1 = 1.0;
        }

        /*
        // * + * * *
        for (unsigned int n=2; n<=N-4; n++)
        {
            double q1 = +22.0*alpha;
            double q2 = -40.0*alpha - 1.0;
            double q3 = +12.0*alpha;
            double q4 = +8.0*alpha;
            double q5 = -2.0*alpha;
            double e  = u1.at(m-1,n) + ht*f(n,m);

            q2 /= -q1;
            q3 /= -q1;
            q4 /= -q1;
            q5 /= -q1;
            e  /= -q1;
            q1 = 1.0;

            c1 = c2 + q2;
            c2 = (c3 + q3)/c1;
            c3 = (c4 + q4)/c1;
            c4 = q5/c1;
            d  = (d - e)/c1;
            c1 = 1.0;
        }
        */

        DoubleMatrix A1(4,4);
        DoubleVector b1(4);
        A1[0][0] = c1;
        A1[0][1] = c2;
        A1[0][2] = c3;
        A1[0][3] = c4;
        b1[0]    = d;

        A1[1][0] = -22.0*alpha;
        A1[1][1] = +40.0*alpha+1.0;
        A1[1][2] = -12.0*alpha;
        A1[1][3] = -8.0*alpha;
        b1[1]    = u1.at(m-1,N-3) - (2.0*alpha)*u1.at(m,N) + ht*f(N-3,m);

        A1[2][0] = +2.0*alpha;
        A1[2][1] = -32.0*alpha;
        A1[2][2] = +60.0*alpha+1.0;
        A1[2][3] = -32.0*alpha;
        b1[2]    = u1.at(m-1,N-2) - (2.0*alpha)*u1.at(m,N) + ht*f(N-2,m);

        A1[3][0] = +2.0*alpha;
        A1[3][1] = -8.0*alpha;
        A1[3][2] = -12.0*alpha;
        A1[3][3] = +40.0*alpha+1.0;
        b1[3]    = u1.at(m-1,N-1) + (22.0*alpha)*u1.at(m,N) + ht*f(N-1,m);

        //printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", c1, c2, c3, c4, d);

        DoubleVector x(4);
        GaussianElimination(A1, b1, x);
        A1.clear();
        b1.clear();

        u1.at(m, N-1) = x.at(3);
        u1.at(m, N-2) = x.at(2);
        u1.at(m, N-3) = x.at(1);
        u1.at(m, N-4) = x.at(0);

        {
            double d0 = -70.0*alpha+1.0;
            double d1 = -208.0*alpha;
            double d2 = +228.0*alpha;
            double d3 = -112.0*alpha;
            double d4 = +22.0*alpha;
            printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", d0, d1, d2, d3, d4, alpha);

            for (unsigned int i=N-5; i>=1; i--)
            {
                u1.at(m,i) = d1*u1.at(m,i+1) + d2*u1.at(m,i+2) + d3*u1.at(m,i+3) + d4*u1.at(m,i+4) + u1.at(m-1,i) + ht*f(i,m);
                //u1.at(m,i) = d1*u(i+1,m) + d2*u(i+2,m) + d3*u(i+3,m) + d4*u(i+4,m) + u(i,m-1) + ht*f(i,m);
                u1.at(m,i) /= d0;
            }
        }

        /*{
            double d0 = -22.0*alpha;
            double d1 = -40.0*alpha-1.0;
            double d2 = +12.0*alpha;
            double d3 = +8.0*alpha;
            double d4 = -2.0*alpha;
            for (unsigned int i=N-4; i>=2; i--)
            {
                u1.at(m,i-1) = d1*u1.at(m,i) + d2*u1.at(m,i+1) + d3*u1.at(m,i+2) + d4*u1.at(m,i+3) + u1.at(m-1,i) + ht*f(i,m);
                //u1.at(m,i-1) = d1*u(i,m)     + d2*u(i+1,m)       + d3*u(i+2,m)   + d4*u(i+3,m)     + u(i,m-1)     + ht*f(i,m);
                u1.at(m,i) /= d0;
            }
        }*/

        /*{
            double d0 = +2.0*alpha;
            double d1 = +32.0*alpha;
            double d2 = -60.0*alpha-1.0;
            double d3 = +32.0*alpha;
            double d4 = -2.0*alpha;
            for (unsigned int i=N-3; i>=3; i--)
            {
                u1.at(m,i-2) = d1*u1.at(m,i-1) + d2*u1.at(m,i) + d3*u1.at(m,i+1) + d4*u1.at(m,i+2) + u1.at(m-1,i) + ht*f(i,m);
                //u1.at(m,i) = d1*u1.at(m,i+1) + d2*u1.at(m,i+2) + d3*u1.at(m,i+3) + d4*u1.at(m,i+4) + u1.at(m-1,i+2) + ht*f(i+2,m);
                //u1.at(m,i-2) = d1*u(i-1,m)     + d2*u(i,m)     + d3*u(i+1,m)     + d4*u(i+2,m)     + u(i,m-1)     + ht*f(i,m);
                u1.at(m,i) /= d0;
            }
        }*/

        /*{
            double d0 = +2.0*alpha;
            double d1 = +8.0*alpha;
            double d2 = +12.0*alpha;
            double d3 = -40.0*alpha-1.0;
            double d4 = +22.0*alpha;
            for (unsigned int i=N-2; i>=4; i--)
            {
                //u1.at(m,i-3) = d1*u1.at(m,i-2) + d2*u1.at(m,i-1) + d3*u1.at(m,i) + d4*u1.at(m,i+1) + u1.at(m-1,i) + ht*f(i,m);
                //u1.at(m,i) = d1*u1.at(m,i+1) + d2*u1.at(m,i+2) + d3*u1.at(m,i+3) + d4*u1.at(m,i+4) + u1.at(m-1,i+2) + ht*f(i+2,m);
                u1.at(m,i-3) = d1*u(i-2,m)     + d2*u(i-1,m)     + d3*u(i,m)     + d4*u(i+1,m)     + u(i,m-1)     + ht*f(i,m);
                u1.at(m,i) /= d0;
            }
        }*/

        if (m==1)
        {
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u(N-4,m), u(N-3,m), u(N-2,m), u(N-1,m));
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, x[0], x[1], x[2], x[3]);
            //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u1.at(m,N-4), u1.at(m,N-3), u1.at(m,N-2), u1.at(m,N-1));
            IPrinter::printVector(14,10, u1.row(m));
            DoubleVector v(N+1); for (unsigned int n=0; n<=N; n++) v.at(n) = u(n,m);
            IPrinter::printVector(14,10, v);
            v.clear();
            IPrinter::printSeperatorLine();
        }
        if(m==1)
        {
            break;
        }

        A1.clear();
        b1.clear();
        x.clear();
    }
}

void BorderParabolic::calculateN4R2L(DoubleMatrix &u1)
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

    double alpha = (a*a*ht)/(24.0*hx*hx);

    for (unsigned int m=1; m<=M; m++)
    {
        DoubleMatrix A1(4,4);
        DoubleVector b1(4);

        A1[0][0] = +40.0*alpha+1.0;
        A1[0][1] = -12.0*alpha;
        A1[0][2] = -8.0*alpha;
        A1[0][3] = +2.0*alpha;
        b1[0]    = u1.at(m-1,1) + (22.0*alpha)*u1.at(m,0) + ht*f(1,m);

        A1[1][0] = -32.0*alpha;
        A1[1][1] = +60.0*alpha+1.0;
        A1[1][2] = -32.0*alpha;
        A1[1][3] = +2.0*alpha;
        b1[1]    = u1.at(m-1,2) - (2.0*alpha)*u1.at(m,0) + ht*f(2,m);

        A1[2][0] = -8.0*alpha;
        A1[2][1] = -12.0*alpha;
        A1[2][2] = +40.0*alpha+1.0;
        A1[2][3] = -22.0*alpha;
        b1[2]    = u1.at(m-1,3) - (2.0*alpha)*u1.at(m,0) + ht*f(3,m);

        double c1 = (40.0*alpha + 1.0);
        double c2 = (-12.0*alpha)/c1;
        double c3 = (-8.0*alpha)/c1;
        double c4 = +2.0*alpha/c1;
        double d  = (u1.at(m-1,N-1) + (22.0*alpha)*u1.at(m,N) + ht*f(N-1,m))/c1;
        c1 = 1.0;

        for (unsigned int n=N-2; n>=4; n--)
        {
            double q1 = -22.0*alpha;
            double q2 = -40.0*alpha - 1.0;
            double q3 = +12.0*alpha;
            double q4 = +8.0*alpha;;
            double q5 = -2.0*alpha;
            double e  = u1.at(m-1,n) + ht*f(n,m);

            q2 /= q1;
            q3 /= q1;
            q4 /= q1;
            q5 /= q1;
            e  /= q1;
            q1 = 1.0;

            c1 = c2 + q2;
            c2 = (c3 + q3)/c1;
            c3 = (c4 + q4)/c1;
            c4 = q5/c1;
            d  = (d - e)/c1;
            c1 = 1.0;
        }

        A1[3][0] = c4;
        A1[3][1] = c3;
        A1[3][2] = c2;
        A1[3][3] = c1;
        b1[3] = d;

        //printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", c1, c2, c3, c4, d);

        DoubleVector x(4);
        GaussianElimination(A1, b1, x);

        //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, x[0], x[1], x[2], x[3]);
        //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u(1,m), u(2,m), u(3,m), u(4,m));
        //return;

        u1.at(m, 1) = x.at(0);
        u1.at(m, 2) = x.at(1);
        u1.at(m, 3) = x.at(2);
        u1.at(m, 4) = x.at(3);

        {
            double d0 = -70.0*alpha+1.0;
            double d1 = -208.0*alpha;
            double d2 = +228.0*alpha;
            double d3 = -112.0*alpha;
            double d4 = +22.0*alpha;
            for (unsigned int i=5; i<=N-1; i++)
            {
                u1.at(m,i) = d1*u1.at(m,i-1) + d2*u1.at(m,i-2) + d3*u1.at(m,i-3) + d4*u1.at(m,i-4) + u1.at(m-1,i) + ht*f(i,m);
                //u1.at(m,i) = d1*u(i+1,m) + d2*u(i+2,m) + d3*u(i+3,m) + d4*u(i+4,m) + u(i,m-1) + ht*f(i,m);
                u1.at(m,i) /= d0;
            }
        }
        /*{
            double d0 = -22.0*alpha;
            double d1 = -40.0*alpha-1.0;
            double d2 = +12.0*alpha;
            double d3 = +8.0*alpha;
            double d4 = -2.0*alpha;
            for (unsigned int i=4; i<=N-2; i++)
            {
                u1.at(m,i+1) = d1*u1.at(m,i) + d2*u1.at(m,i-1) + d3*u1.at(m,i-2) + d4*u1.at(m,i-3) + u1.at(m-1,i) + ht*f(i,m);
                //u1.at(m,i-1) = d1*u(i,m)     + d2*u(i+1,m)       + d3*u(i+2,m)   + d4*u(i+3,m)     + u(i,m-1)     + ht*f(i,m);
                u1.at(m,i) /= d0;
            }
        }*/
        /*{
            double d0 = +2.0*alpha;
            double d1 = +32.0*alpha;
            double d2 = -60.0*alpha-1.0;
            double d3 = +32.0*alpha;
            double d4 = -2.0*alpha;
            for (unsigned int i=N-3; i>=3; i--)
            {
                u1.at(m,i-2) = d1*u1.at(m,i-1) + d2*u1.at(m,i) + d3*u1.at(m,i+1) + d4*u1.at(m,i+2) + u1.at(m-1,i) + ht*f(i,m);
                //u1.at(m,i) = d1*u1.at(m,i+1) + d2*u1.at(m,i+2) + d3*u1.at(m,i+3) + d4*u1.at(m,i+4) + u1.at(m-1,i+2) + ht*f(i+2,m);
                //u1.at(m,i-2) = d1*u(i-1,m)     + d2*u(i,m)     + d3*u(i+1,m)     + d4*u(i+2,m)     + u(i,m-1)     + ht*f(i,m);
                u1.at(m,i) /= d0;
            }
        }*/

        /*{
            double d0 = +2.0*alpha;
            double d1 = +8.0*alpha;
            double d2 = +12.0*alpha;
            double d3 = -40.0*alpha-1.0;
            double d4 = +22.0*alpha;
            for (unsigned int i=N-2; i>=4; i--)
            {
                //u1.at(m,i-3) = d1*u1.at(m,i-2) + d2*u1.at(m,i-1) + d3*u1.at(m,i) + d4*u1.at(m,i+1) + u1.at(m-1,i) + ht*f(i,m);
                //u1.at(m,i) = d1*u1.at(m,i+1) + d2*u1.at(m,i+2) + d3*u1.at(m,i+3) + d4*u1.at(m,i+4) + u1.at(m-1,i+2) + ht*f(i+2,m);
                u1.at(m,i-3) = d1*u(i-2,m)     + d2*u(i-1,m)     + d3*u(i,m)     + d4*u(i+1,m)     + u(i,m-1)     + ht*f(i,m);
                u1.at(m,i) /= d0;
            }
        }*/

        if (m==1)
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

        A1.clear();
        b1.clear();
        x.clear();
    }
}

void BorderParabolic::calculateN4(DoubleMatrix &u1)
{
    u1.resize(M+1, N+1);
    for (unsigned int i=0; i<=N; i++) u1.at(0,i) = initial(i);
    for (unsigned int j=1; j<=M; j++)
    {
        u1.at(j,0) = boundary(Left,j);
        u1.at(j,N) = boundary(Right, j);
    }

    double alpha = (ht*a*a)/(24.0*hx*hx);
    double norm = 0.0;

    DoubleMatrix A(N-1, N-1, 0.0);
    DoubleVector b(N-1, 0.0);

    for (unsigned int m=1; m<=M; m++)
    {
        //----------------------------------------------------------------------------------------------------------------
        A.at(0,0) = +40.0*alpha+1.0;
        A.at(0,1) = -12.0*alpha;
        A.at(0,2) = -8.0*alpha;
        A.at(0,3) = +2.0*alpha;
        b.at(0) = u1.at(m-1,1) + (22.0*alpha)*u1.at(m,0) + ht*f(1,m);
        norm = +40.0*alpha+1.0;

        A.at(0,0) /= norm;
        A.at(0,1) /= norm;
        A.at(0,2) /= norm;
        A.at(0,3) /= norm;
        b.at(0)   /= norm;
        //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, A.at(0,0), A.at(0,1), A.at(0,2), A.at(0,3), b.at(0));

        //----------------------------------------------------------------------------------------------------------------
        for (unsigned int n=2; n<=N-4; n++)
        {
            A.at(n-1,n-2) = -22.0*alpha;
            A.at(n-1,n-1) = +40.0*alpha+1.0;
            A.at(n-1,n-0) = -12.0*alpha;
            A.at(n-1,n+1) = -8.0*alpha;
            A.at(n-1,n+2) = +2.0*alpha;
            b.at(n-1) = u1.at(m-1,n) + ht*f(n,m);
            norm  = -22.0*alpha;

            A.at(n-1,n-2) /= norm;
            A.at(n-1,n-1) /= norm;
            A.at(n-1,n-0) /= norm;
            A.at(n-1,n+1) /= norm;
            A.at(n-1,n+2) /= norm;
            b.at(n-1)     /= norm;
        }

        //----------------------------------------------------------------------------------------------------------------
        A.at(N-4,N-5) = -22.0*alpha;
        A.at(N-4,N-4) = +40.0*alpha+1.0;
        A.at(N-4,N-3) = -12.0*alpha;
        A.at(N-4,N-2) = -8.0*alpha;
        b.at(N-4) = u1.at(m-1,N-3) - (2.0*alpha)*u1.at(m,N) + ht*f(N-3,m);
        norm  = -22.0*alpha;

        A.at(N-4,N-5) /= norm;
        A.at(N-4,N-4) /= norm;
        A.at(N-4,N-3) /= norm;
        A.at(N-4,N-2) /= norm;
        b.at(N-4)     /= norm;

        //----------------------------------------------------------------------------------------------------------------
        A.at(N-3,N-5) = +2.0*alpha;
        A.at(N-3,N-4) = -32.0*alpha;
        A.at(N-3,N-3) = +60.0*alpha+1.0;
        A.at(N-3,N-2) = -32.0*alpha;
        b.at(N-3) = u1.at(m-1,N-2) - (2.0*alpha)*u1.at(m,N) + ht*f(N-2,m);
        norm  = +2.0*alpha;

        A.at(N-3,N-5) /= norm;
        A.at(N-3,N-4) /= norm;
        A.at(N-3,N-3) /= norm;
        A.at(N-3,N-2) /= norm;
        b.at(N-3)     /= norm;

        //----------------------------------------------------------------------------------------------------------------
        A.at(N-2,N-5) = +2.0*alpha;
        A.at(N-2,N-4) = -8.0*alpha;
        A.at(N-2,N-3) = -12.0*alpha;
        A.at(N-2,N-2) = +40.0*alpha+1.0;
        b.at(N-2) = u1.at(m-1,N-1) + (22.0*alpha)*u1.at(m,N) + ht*f(N-1,m);
        norm  = +2.0*alpha;

        A.at(N-2,N-5) /= norm;
        A.at(N-2,N-4) /= norm;
        A.at(N-2,N-3) /= norm;
        A.at(N-2,N-2) /= norm;
        b.at(N-2)     /= norm;

        //printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", A.at(N-3,N-5), A.at(N-3,N-4), A.at(N-3,N-3), A.at(N-3,N-2), b.at(N-3));
        //printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", A.at(N-2,N-5), A.at(N-2,N-4), A.at(N-2,N-3), A.at(N-2,N-2), b.at(N-2));

        //FILE *file = fopen("A_Matrix.txt", "w");
        //IPrinter::print(A,A.rows(),A.cols(),14,10,file);
        //fclose(file);

        double z1 = A.at(0,0);
        double p1 = A.at(0,1);
        double q1 = A.at(0,2);
        double k1 = A.at(0,3);
        double r1 = b.at(0);

        //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, z1, p1, q1, k1, r1);
        //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, z2, p2, q2, k2, r2);
        //IPrinter::printSeperatorLine();

        for (unsigned int i=1; i<=N-5; i++)
        {
            //            double g1 = z1/A.at(i,i-1);
            //            z1 = -A.at(i,i-0) * g1 + p1;
            //            p1 = -A.at(i,i+1) * g1 + q1;
            //            q1 = -A.at(i,i+2) * g1 + k1;
            //            k1 = -A.at(i,i+3) * g1;
            //            r1 = -b.at(i) * g1 + r1;

            z1 = -A.at(i,i-0) + p1;
            p1 = -A.at(i,i+1) + q1;
            q1 = -A.at(i,i+2) + k1;
            k1 = -A.at(i,i+3);
            r1 = -b.at(i) + r1;

            //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f G1:%18.10f\n", i, z1, p1, q1, k1, r1, g1);

            double norm1 = z1;//sqrt(z1*z1 + p1*p1 + q1*q1 + k1*k1);
            z1 /= norm1;
            p1 /= norm1;
            q1 /= norm1;
            k1 /= norm1;
            r1 /= norm1;

            //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f N1: %14.10f\n", i, z1, p1, q1, k1, r1, norm1);
            //IPrinter::printSeperatorLine();
            // if (i>5) break;
        }

        //printf("%14.10f %14.10f\n", z1*u(N-4,m)+p1*u(N-3,m)+q1*u(N-2,m)+k1*u(N-1,m), r1);

        DoubleMatrix A1(4,4);
        DoubleVector b1(4);
        A1[0][0] = z1;            A1[0][1] = p1;            A1[0][2] = q1;            A1[0][3] = k1; b1[0] = r1;
        A1[1][0] = A.at(N-4,N-5); A1[1][1] = A.at(N-4,N-4); A1[1][2] = A.at(N-4,N-3); A1[1][3] = A.at(N-4,N-2); b1[1] = b.at(N-4);
        A1[2][0] = A.at(N-3,N-5); A1[2][1] = A.at(N-3,N-4); A1[2][2] = A.at(N-3,N-3); A1[2][3] = A.at(N-3,N-2); b1[2] = b.at(N-3);
        A1[3][0] = A.at(N-2,N-5); A1[3][1] = A.at(N-2,N-4); A1[3][2] = A.at(N-2,N-3); A1[3][3] = A.at(N-2,N-2); b1[3] = b.at(N-2);
        //printf("det A %14.10f\n", A1.determinant());

        //        DoubleVector x(N-1);
        //        GaussianElimination(A,b,x);
        //        for (unsigned int i=0; i<x.size(); i++)
        //        {
        //            u1.at(m,i+1) = x.at(i);
        //        }

        DoubleVector x(4);
        GaussianElimination(A1, b1, x);
        A1.clear();
        b1.clear();
        //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u(N-4,m), u(N-3,m), u(N-2,m), u(N-1,m));
        //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, x[0], x[1], x[2], x[3]);

        u1.at(m, N-1) = x.at(3);
        u1.at(m, N-2) = x.at(2);
        u1.at(m, N-3) = x.at(1);
        u1.at(m, N-4) = x.at(0);
        x.clear();
        //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u1.at(m,N-4), u1.at(m,N-3), u1.at(m,N-2), u1.at(m,N-1));
        //double k = ;

        double d0 = -70.0*alpha+1.0;
        double d1 = -208.0*alpha;
        double d2 = +228.0*alpha;
        double d3 = -112.0*alpha;
        double d4 = +22.0*alpha;

        //        double d0 = -22.0*alpha;
        //        double d1 = -40.0*alpha-1.0;
        //        double d2 = +12.0*alpha;
        //        double d3 = +8.0*alpha;
        //        double d4 = -2.0*alpha;

        for (unsigned int i=N-4; i>=1; i--)
        {
            u1.at(m,i) = d1*u1.at(m,i+1) + d2*u1.at(m,i+2) + d3*u1.at(m,i+3) + d4*u1.at(m,i+4) + u1.at(m-1,i) + ht*f(i,m);
            //u1.at(m,i) = d1*u(i+1,m) + d2*u(i+2,m) + d3*u(i+3,m) + d4*u(i+4,m) + u(i,m-1) + ht*f(i,m);
            //u1.at(m,i) = d1*u1.at(m,i+1) + d2*u1.at(m,i+2) + d3*u1.at(m,i+3) + d4*u1.at(m,i+4) + u1.at(m-1,i+1) + ht*f(i+1,m);
            //u1.at(m,i) = d1*u(i+1,m) + d2*u(i+2,m) + d3*u(i+3,m) + d4*u(i+4,m) + u(i+1,m-1) + ht*f(i+1,m);
            u1.at(m,i) /= d0;
        }
        //IPrinter::printVector(14,10,u1.row(m));
        //if (m>0) break;
    }
    //IPrinter::printSeperatorLine();
    A.clear();
    b.clear();
}

/*
void BorderParabolic::calculateN41(DoubleMatrix &u)
{
    u.resize(M+1, N+1);
    for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);
    for (unsigned int j=1; j<=M; j++)
    {
        u.at(j,0) = boundary(Left,j);
        u.at(j,N) = boundary(Right, j);
    }

    double alpha = (ht*a*a)/(24.0*hx*hx);

    for (unsigned int m=1; m<=M; m++)
    {
        DoubleMatrix A(N-1, N-1);
        DoubleVector b(N-1);

        A.at(0,0) = +40.0*alpha+1.0;
        A.at(0,1) = -12.0*alpha;
        A.at(0,2) = -8.0*alpha;
        A.at(0,3) = +2.0*alpha;
        b.at(0) = u.at(m-1,1) + (22.0*alpha)*u.at(m,0) + ht*f(1,m);

        A.at(1,0) = -32.0*alpha;
        A.at(1,1) = +60.0*alpha+1.0;
        A.at(1,2) = -32.0*alpha;
        A.at(1,3) = +2.0*alpha;
        b.at(1) = u.at(m-1,2) - (2.0*alpha)*u.at(m,0) + ht*f(2,m);

        for (unsigned int n=3; n<=N-3; n++)
        {
            A.at(n-1,n-3) = +2.0*alpha;
            A.at(n-1,n-2) = -32.0*alpha;
            A.at(n-1,n-1) = +60.0*alpha + 1.0;
            A.at(n-1,n-0) = -32.0*alpha;
            A.at(n-1,n+1) = +2.0*alpha;
            b.at(n-1) = u.at(m-1,n) + ht*f(n,m);
        }

        A.at(N-3,N-5) = +2.0*alpha;
        A.at(N-3,N-4) = -32.0*alpha;
        A.at(N-3,N-3) = +60.0*alpha+1.0;
        A.at(N-3,N-2) = -32.0*alpha;
        b.at(N-3) = u.at(m-1,N-2) - (2.0*alpha)*u.at(m,N) + ht*f(N-2,m);

        A.at(N-2,N-5) = +2.0*alpha;
        A.at(N-2,N-4) = -8.0*alpha;
        A.at(N-2,N-3) = -12.0*alpha;
        A.at(N-2,N-2) = +40.0*alpha+1.0;
        b.at(N-2) = u.at(m-1,N-1) + (22.0*alpha)*u.at(m,N) + ht*f(N-1,m);

        DoubleVector x(N-1);
        GaussianElimination(A,b,x);
        for (unsigned int i=0; i<x.size(); i++)
        {
            u.at(m,i+1) = x.at(i);
        }
        A.clear();
        b.clear();
        x.clear();
        //IPrinter::printVector(18,14,u.row(m));
        //IPrinter::printSeperatorLine();
    }
}

void BorderParabolic::calculateN42(DoubleMatrix &u)
{
    u.resize(M+1, N+1);
    for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);
    for (unsigned int j=1; j<=M; j++)
    {
        u.at(j,0) = boundary(Left,j);
        u.at(j,N) = boundary(Right, j);
    }

    double alpha = (ht*a*a)/(24.0*hx*hx);

    for (unsigned int m=1; m<=M; m++)
    {
        DoubleVector betta1(N+1);
        betta1[1] = +40.0*alpha+1.0;
        betta1[2] = -12.0*alpha;
        betta1[3] = -8.0*alpha;
        betta1[4] = +2.0*alpha;
        double eta1 = u.at(m-1,1) - (22.0*alpha)*u.at(m,0) + ht*f(1,m);

        for (unsigned int k=3; k<=N-3; k++)
        {
            double a1 = +16.0;
            double a2 = -(60.0*alpha+1.0)/(2.0*alpha);
            double a3 = +16.0;
            double a4 = -1.0;
            double a0 = (u.at(m-1,k) + ht*f(k,m))/(2.0*alpha);

            eta1 = eta1 - betta1[k-2]*a0;
            betta1[k-2+1] = betta1[k-2]*a1 + betta1[k-2+1];
            betta1[k-2+2] = betta1[k-2]*a2 + betta1[k-2+2];
            betta1[k-2+3] = betta1[k-2]*a3 + betta1[k-2+3];
            betta1[k-2+4] = betta1[k-2]*a4;
        }

        DoubleVector betta2(N+1);
        betta2[1] = -32.0*alpha;
        betta2[2] = +60.0*alpha+1.0;
        betta2[3] = -32.0*alpha;
        betta2[4] = +2.0*alpha;
        double eta2 = u.at(m-1,2) - (2.0*alpha)*u.at(m,0) + ht*f(2,m);

        for (unsigned int k=3; k<=N-3; k++)
        {
            double a1 = +16.0;
            double a2 = -(60.0*alpha+1.0)/(2.0*alpha);
            double a3 = +16.0;
            double a4 = -1.0;
            double a0 = (u.at(m-1,k) + ht*f(k,m))/(2.0*alpha);

            eta2 = eta2 - betta2[k-2]*a0;
            betta2[k-2+1] = betta2[k-2]*a1 + betta2[k-2+1];
            betta2[k-2+2] = betta2[k-2]*a2 + betta2[k-2+2];
            betta2[k-2+3] = betta2[k-2]*a3 + betta2[k-2+3];
            betta2[k-2+4] = betta2[k-2]*a4;
        }

        printf("%d\n", m);
        printf("%f %f %f %f\n", betta1[N-4], betta1[N-3], betta1[N-2], betta1[N-1]);
        printf("%f %f %f %f\n", betta2[N-4], betta2[N-3], betta2[N-2], betta2[N-1]);

        DoubleMatrix A(4,4);
        A.at(0,0) = betta1[N-4]; A.at(0,1) = betta1[N-3]; A.at(0,2) = betta1[N-2];     A.at(0,3) = betta1[N-1];
        A.at(1,0) = betta2[N-4]; A.at(1,1) = betta2[N-3]; A.at(1,2) = betta2[N-2];     A.at(1,3) = betta2[N-1];
        A.at(2,0) = +2.0*alpha;  A.at(2,1) = -32.0*alpha; A.at(2,2) = +60.0*alpha+1.0; A.at(2,3) = -32.0*alpha;
        A.at(3,0) = +2.0*alpha;  A.at(3,1) = -8.0*alpha;  A.at(3,2) = -12.0*alpha;     A.at(3,3) = +40.0*alpha+1.0;

        DoubleVector b(4);
        b.at(0) = eta1;
        b.at(1) = eta2;
        b.at(2) = u.at(m-1,N-2) - (2.0*alpha)*u.at(m,N) + ht*f(N-2,m);
        b.at(3) = u.at(m-1,N-1) + (22.0*alpha)*u.at(m,N) + ht*f(N-1,m);

        DoubleVector x(4);
        GaussianElimination(A,b,x);
        //printf("%f %f %f %f %f\n",x[0],x[1],x[2],x[3],A.determinant1());

        break;
    }
}

void BorderParabolic::calculateN43(DoubleMatrix &u1)
{
    u1.resize(M+1, N+1);
    for (unsigned int i=0; i<=N; i++) u1.at(0,i) = initial(i);
    for (unsigned int j=1; j<=M; j++)
    {
        u1.at(j,0) = boundary(Left,j);
        u1.at(j,N) = boundary(Right, j);
    }

    double alpha = (ht*a*a)/(24.0*hx*hx);
    printf("alpha: %14.10f\n", alpha);

    for (unsigned int m=1; m<=M; m++)
    {
        DoubleMatrix A(N-1, N-1, 0.0);
        DoubleVector b(N-1, 0.0);

        //        A.at(0,0) = +40.0*alpha+1.0;
        //        A.at(0,1) = -12.0*alpha;
        //        A.at(0,2) = -8.0*alpha;
        //        A.at(0,3) = +2.0*alpha;
        //        b.at(0) = u1.at(m-1,1) + (22.0*alpha)*u1.at(m,0) + ht*f(1,m);

        A.at(0,0) = 1.0;
        A.at(0,1) = -(12.0*alpha)/(+40.0*alpha+1.0);
        A.at(0,2) = -(8.0*alpha)/(+40.0*alpha+1.0);
        A.at(0,3) = +(2.0*alpha)/(+40.0*alpha+1.0);
        b.at(0) = (u1.at(m-1,1) + (22.0*alpha)*u1.at(m,0) + ht*f(1,m))/(+40.0*alpha+1.0);

        //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, A.at(0,0), A.at(0,1), A.at(0,2), A.at(0,3), b.at(0));

        //        printf("%14.10f %14.10f\n", A.at(0,0)*u(1,m)+A.at(0,1)*u(2,m)+A.at(0,2)*u(3,m)+A.at(0,3)*u(4,m), b.at(0));

        //        A.at(1,0) = -32.0*alpha;
        //        A.at(1,1) = +60.0*alpha+1.0;
        //        A.at(1,2) = -32.0*alpha;
        //        A.at(1,3) = +2.0*alpha;
        //        b.at(1) = u1.at(m-1,2) - (2.0*alpha)*u1.at(m,0) + ht*f(2,m);

        A.at(1,0) = 1.0;
        A.at(1,1) = +(60.0*alpha+1.0)/(-32.0*alpha);
        A.at(1,2) = -(32.0*alpha)/(-32.0*alpha);
        A.at(1,3) = +(2.0*alpha)/(-32.0*alpha);
        b.at(1) = (u1.at(m-1,2) - (2.0*alpha)*u1.at(m,0) + ht*f(2,m))/(-32.0*alpha);

        //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 1, A.at(1,0), A.at(1,1), A.at(1,2), A.at(1,3), b.at(1));
        //        printf("%14.10f %14.10f\n", A.at(1,0)*u(1,m)+A.at(1,1)*u(2,m)+A.at(1,2)*u(3,m)+A.at(1,3)*u(4,m), b.at(1));

        for (unsigned int n=3; n<=N-3; n++)
        {
            //            A.at(n-1,n-3) = +2.0*alpha;
            //            A.at(n-1,n-2) = -32.0*alpha;
            //            A.at(n-1,n-1) = +60.0*alpha + 1.0;
            //            A.at(n-1,n-0) = -32.0*alpha;
            //            A.at(n-1,n+1) = +2.0*alpha;
            //            b.at(n-1) = u1.at(m-1,n) + ht*f(n,m);

            A.at(n-1,n-3) = +1.0;
            A.at(n-1,n-2) = -16.0;
            A.at(n-1,n-1) = +(60.0*alpha + 1.0)/(2.0*alpha);
            A.at(n-1,n-0) = -16.0;
            A.at(n-1,n+1) = +1.0;
            b.at(n-1) = (u1.at(m-1,n) + ht*f(n,m))/(2.0*alpha);

            //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", n-1, A.at(n-1,n-3), A.at(n-1,n-2), A.at(n-1,n-1),
            //       A.at(n-1,n-0), A.at(n-1,n+1), b.at(n-1));


            //            printf("%14.10f %14.10f\n", A.at(n-1,n-3)*u(n-2,m)
            //                                       +A.at(n-1,n-2)*u(n-1,m)
            //                                       +A.at(n-1,n-1)*u(n+0,m)
            //                                       +A.at(n-1,n+0)*u(n+1,m)
            //                                       +A.at(n-1,n+1)*u(n+2,m), b.at(n-1));
        }

        //return;

        //        A.at(N-3,N-5) = +2.0*alpha;
        //        A.at(N-3,N-4) = -32.0*alpha;
        //        A.at(N-3,N-3) = +60.0*alpha+1.0;
        //        A.at(N-3,N-2) = -32.0*alpha;
        //        b.at(N-3) = u1.at(m-1,N-2) - (2.0*alpha)*u1.at(m,N) + ht*f(N-2,m);

        A.at(N-3,N-5) = 1.0;
        A.at(N-3,N-4) = -16.0;
        A.at(N-3,N-3) = (+60.0*alpha+1.0)/(2.0*alpha);
        A.at(N-3,N-2) = -16;
        b.at(N-3) = (u1.at(m-1,N-2) - (2.0*alpha)*u1.at(m,N) + ht*f(N-2,m))/(2.0*alpha);


        //        printf("%14.10f %14.10f\n", A.at(N-3,N-5)*u(N-4,m)
        //                                   +A.at(N-3,N-4)*u(N-3,m)
        //                                   +A.at(N-3,N-3)*u(N-2,m)
        //                                   +A.at(N-3,N-2)*u(N-1,m), b.at(N-3));

        //        A.at(N-2,N-5) = +2.0*alpha;
        //        A.at(N-2,N-4) = -8.0*alpha;
        //        A.at(N-2,N-3) = -12.0*alpha;
        //        A.at(N-2,N-2) = +40.0*alpha+1.0;
        //        b.at(N-2) = u1.at(m-1,N-1) + (22.0*alpha)*u1.at(m,N) + ht*f(N-1,m);

        A.at(N-2,N-5) = +1.0;
        A.at(N-2,N-4) = -4.0;
        A.at(N-2,N-3) = -6.0;
        A.at(N-2,N-2) = (+40.0*alpha+1.0)/(2.0*alpha);
        b.at(N-2) = (u1.at(m-1,N-1) + (22.0*alpha)*u1.at(m,N) + ht*f(N-1,m))/(2.0*alpha);

        //        printf("%14.10f %14.10f\n", A.at(N-2,N-5)*u(N-4,m)
        //                                   +A.at(N-2,N-4)*u(N-3,m)
        //                                   +A.at(N-2,N-3)*u(N-2,m)
        //                                   +A.at(N-2,N-2)*u(N-1,m), b.at(N-2));

        //        FILE *file = fopen("A_Matrix.txt", "w");
        //        IPrinter::print(A,A.rows(),A.cols(),14,10,file);
        //        fclose(file);

        double z1 = A.at(0,0);
        double p1 = A.at(0,1);
        double q1 = A.at(0,2);
        double k1 = A.at(0,3);
        double r1 = b.at(0);

        double z2 = A.at(1,0);
        double p2 = A.at(1,1);
        double q2 = A.at(1,2);
        double k2 = A.at(1,3);
        double r2 = b.at(1);

        //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, z1, p1, q1, k1, r1);
        //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, z2, p2, q2, k2, r2);
        //IPrinter::printSeperatorLine();

        for (unsigned int i=2; i<=N-4; i++)
        {
            //p1 = -p1/z1;
            //q1 = -q1/z1;
            //k1 = -k1/z1;
            //r1 = +r1/z1;
            //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", i, z1, p1, q1, k1, r1);

            z1 = -A.at(i,i-1) + p1;
            p1 = -A.at(i,i+0) + q1;
            q1 = -A.at(i,i+1) + k1;
            k1 = -A.at(i,i+2);
            r1 = -b.at(i) + r1;
            //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", i, z1, p1, q1, k1, r1);

            p1 /= z1;
            q1 /= z1;
            k1 /= z1;
            r1 /= z1;
            z1 /= z1;
            //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", i, z1, p1, q1, k1, r1);

            //            p2 = -p2/z2;
            //            q2 = -q2/z2;
            //            k2 = -k2/z2;
            //            r2 = +r2/z2;

            z2 = -A.at(i,i-1) + p2;
            p2 = -A.at(i,i+0) + q2;
            q2 = -A.at(i,i+1) + k2;
            k2 = -A.at(i,i+2);
            r2 = -b.at(i) + r2;
            //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", i, z2, p2, q2, k2, r2);

            p2 /= z2;
            q2 /= z2;
            k2 /= z2;
            r2 /= z2;
            z2 /= z2;
            //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", i, z2, p2, q2, k2, r2);

            //IPrinter::printSeperatorLine();
            //if (i>15)
            //            break;
        }

        printf("%14.10f %14.10f\n", z1*u(N-4,m)+p1*u(N-3,m)+q1*u(N-2,m)+k1*u(N-1,m), r1);

        DoubleMatrix A1(4,4);
        DoubleVector b1(4);
        A1[0][0] = -z1;            A1[0][1] = -p1;            A1[0][2] = -q1;            A1[0][3] = -k1; b1[0] = r1;
        A1[1][0] = -z2;            A1[1][1] = -p2;            A1[1][2] = -q2;            A1[1][3] = -k2; b1[1] = r2;
        A1[2][0] = A.at(N-3,N-5); A1[2][1] = A.at(N-3,N-4); A1[2][2] = A.at(N-3,N-3); A1[2][3] = A.at(N-3,N-2); b1[2] = b.at(N-3);
        A1[3][0] = A.at(N-2,N-5); A1[3][1] = A.at(N-2,N-4); A1[3][2] = A.at(N-2,N-3); A1[3][3] = A.at(N-2,N-2); b1[3] = b.at(N-2);
        DoubleVector x(4);
        GaussianElimination(A1, b1, x);
        printf("%18.10f %18.10f %18.10f %18.10f\n", x[0], x[1], x[2], x[3]);

        A.clear();
        b.clear();
        break;
    }
}

void BorderParabolic::calculateN44(DoubleMatrix &u1)
{
    u1.resize(M+1, N+1);
    for (unsigned int i=0; i<=N; i++) u1.at(0,i) = initial(i);
    for (unsigned int j=1; j<=M; j++)
    {
        u1.at(j,0) = boundary(Left,j);
        u1.at(j,N) = boundary(Right, j);
    }

    double alpha = (ht*a*a)/(24.0*hx*hx);

    for (unsigned int m=1; m<=M; m++)
    {
        DoubleMatrix A(N-1, N-1, 0.0);
        DoubleVector b(N-1, 0.0);

        A.at(0,0) = 1.0;
        A.at(0,1) = -(12.0*alpha)/(+40.0*alpha+1.0);
        A.at(0,2) = -(8.0*alpha)/(+40.0*alpha+1.0);
        A.at(0,3) = +(2.0*alpha)/(+40.0*alpha+1.0);
        b.at(0) = (u1.at(m-1,1) + (22.0*alpha)*u1.at(m,0) + ht*f(1,m))/(+40.0*alpha+1.0);

        A.at(1,0) = 1.0;
        A.at(1,1) = +(60.0*alpha+1.0)/(-32.0*alpha);
        A.at(1,2) = -(32.0*alpha)/(-32.0*alpha);
        A.at(1,3) = +(2.0*alpha)/(-32.0*alpha);
        b.at(1) = (u1.at(m-1,2) - (2.0*alpha)*u1.at(m,0) + ht*f(2,m))/(-32.0*alpha);

        for (unsigned int n=3; n<=N-3; n++)
        {
            A.at(n-1,n-3) = +1.0;
            A.at(n-1,n-2) = -16.0;
            A.at(n-1,n-1) = +(60.0*alpha + 1.0)/(2.0*alpha);
            A.at(n-1,n-0) = -16.0;
            A.at(n-1,n+1) = +1.0;
            b.at(n-1) = (u1.at(m-1,n) + ht*f(n,m))/(2.0*alpha);
        }

        A.at(N-3,N-5) = 1.0;
        A.at(N-3,N-4) = -16.0;
        A.at(N-3,N-3) = (+60.0*alpha+1.0)/(2.0*alpha);
        A.at(N-3,N-2) = -16;
        b.at(N-3) = (u1.at(m-1,N-2) - (2.0*alpha)*u1.at(m,N) + ht*f(N-2,m))/(2.0*alpha);

        A.at(N-2,N-5) = +1.0;
        A.at(N-2,N-4) = -4.0;
        A.at(N-2,N-3) = -6.0;
        A.at(N-2,N-2) = (+40.0*alpha+1.0)/(2.0*alpha);
        b.at(N-2) = (u1.at(m-1,N-1) + (22.0*alpha)*u1.at(m,N) + ht*f(N-1,m))/(2.0*alpha);

        FILE *file = fopen("A_Matrix.txt", "w");
        IPrinter::print(A,A.rows(),A.cols(),14,10,file);
        fclose(file);

        DoubleVector x(N-1);
        GaussianElimination(A,b,x);
        for (unsigned int i=0; i<x.size(); i++)
        {
            u1.at(m,i+1) = x.at(i);
        }
        A.clear();
        b.clear();
        x.clear();
    }
}

void BorderParabolic::calculateN45(DoubleMatrix &u1)
{
    u1.resize(M+1, N+1);
    for (unsigned int i=0; i<=N; i++) u1.at(0,i) = initial(i);
    for (unsigned int j=1; j<=M; j++)
    {
        u1.at(j,0) = boundary(Left,j);
        u1.at(j,N) = boundary(Right, j);
    }

    double alpha = (ht*a*a)/(24.0*hx*hx);
    //printf("alpha: %14.10f\n", alpha);
    double norm = 0.0;

    for (unsigned int m=1; m<=M; m++)
    {
        DoubleMatrix A(N-1, N-1, 0.0);
        DoubleVector b(N-1, 0.0);

        A.at(0,0) = +40.0*alpha+1.0;
        A.at(0,1) = -12.0*alpha;
        A.at(0,2) = -8.0*alpha;
        A.at(0,3) = +2.0*alpha;
        b.at(0) = u1.at(m-1,1) + (22.0*alpha)*u1.at(m,0) + ht*f(1,m);
        //norm = 1.0;
        //norm = sqrt(A.at(0,0)*A.at(0,0) + A.at(0,1)*A.at(0,1) + A.at(0,2)*A.at(0,2) + A.at(0,3)*A.at(0,3));
        norm = +40.0*alpha+1.0;

        A.at(0,0) /= norm;
        A.at(0,1) /= norm;
        A.at(0,2) /= norm;
        A.at(0,3) /= norm;
        b.at(0)   /= norm;
        printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, A.at(0,0), A.at(0,1), A.at(0,2), A.at(0,3), b.at(0));

        A.at(1,0) = -32.0*alpha;
        A.at(1,1) = +60.0*alpha+1.0;
        A.at(1,2) = -32.0*alpha;
        A.at(1,3) = +2.0*alpha;
        b.at(1) = u1.at(m-1,2) - (2.0*alpha)*u1.at(m,0) + ht*f(2,m);
        //norm = 1.0;
        //norm = sqrt(A.at(1,0)*A.at(1,0) + A.at(1,1)*A.at(1,1) + A.at(1,2)*A.at(1,2) + A.at(1,3)*A.at(1,3));
        norm = -32.0*alpha;

        A.at(1,0) /= norm;
        A.at(1,1) /= norm;
        A.at(1,2) /= norm;
        A.at(1,3) /= norm;
        b.at(1)   /= norm;
        printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 1, A.at(1,0), A.at(1,1), A.at(1,2), A.at(1,3), b.at(1));
        IPrinter::printSeperatorLine();

        for (unsigned int n=3; n<=N-3; n++)
        {
            A.at(n-1,n-3) = +2.0*alpha;
            A.at(n-1,n-2) = -32.0*alpha;
            A.at(n-1,n-1) = +60.0*alpha + 1.0;
            A.at(n-1,n-0) = -32.0*alpha;
            A.at(n-1,n+1) = +2.0*alpha;
            b.at(n-1) = u1.at(m-1,n) + ht*f(n,m);
            //norm = 1.0;
            //norm = sqrt(A.at(n-1,n-3)*A.at(n-1,n-3) + A.at(n-1,n-2)*A.at(n-1,n-2) + A.at(n-1,n-1)*A.at(n-1,n-1) + A.at(n-1,n-0)*A.at(n-1,n-0) + A.at(n-1,n+1)*A.at(n-1,n+1));
            norm  = +2.0*alpha;

            A.at(n-1,n-3) /= norm;
            A.at(n-1,n-2) /= norm;
            A.at(n-1,n-1) /= norm;
            A.at(n-1,n-0) /= norm;
            A.at(n-1,n+1) /= norm;
            b.at(n-1)     /= norm;
        }

        A.at(N-3,N-5) = +2.0*alpha;
        A.at(N-3,N-4) = -32.0*alpha;
        A.at(N-3,N-3) = +60.0*alpha+1.0;
        A.at(N-3,N-2) = -32.0*alpha;
        b.at(N-3) = u1.at(m-1,N-2) - (2.0*alpha)*u1.at(m,N) + ht*f(N-2,m);
        //norm = 1.0;
        //norm = sqrt(A.at(N-3,N-5)*A.at(N-3,N-5) + A.at(N-3,N-4)*A.at(N-3,N-4) + A.at(N-3,N-3)*A.at(N-3,N-3) + A.at(N-3,N-2)*A.at(N-3,N-2));
        norm  = +2.0*alpha;

        A.at(N-3,N-5) /= norm;
        A.at(N-3,N-4) /= norm;
        A.at(N-3,N-3) /= norm;
        A.at(N-3,N-2) /= norm;
        b.at(N-3)     /= norm;

        A.at(N-2,N-5) = +2.0*alpha;
        A.at(N-2,N-4) = -8.0*alpha;
        A.at(N-2,N-3) = -12.0*alpha;
        A.at(N-2,N-2) = +40.0*alpha+1.0;
        b.at(N-2) = u1.at(m-1,N-1) + (22.0*alpha)*u1.at(m,N) + ht*f(N-1,m);
        //norm = 1.0;
        //norm = sqrt(A.at(N-2,N-5)*A.at(N-2,N-5) + A.at(N-2,N-4)*A.at(N-2,N-4) + A.at(N-2,N-3)*A.at(N-2,N-3) + A.at(N-2,N-2)*A.at(N-2,N-2));
        norm  = +2.0*alpha;

        A.at(N-2,N-5) /= norm;
        A.at(N-2,N-4) /= norm;
        A.at(N-2,N-3) /= norm;
        A.at(N-2,N-2) /= norm;
        b.at(N-2)     /= norm;

        //printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", A.at(N-3,N-5), A.at(N-3,N-4), A.at(N-3,N-3), A.at(N-3,N-2), b.at(N-3));
        //printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", A.at(N-2,N-5), A.at(N-2,N-4), A.at(N-2,N-3), A.at(N-2,N-2), b.at(N-2));

        FILE *file = fopen("A_Matrix.txt", "w");
        IPrinter::print(A,A.rows(),A.cols(),14,10,file);
        fclose(file);

        double z1 = A.at(0,0);
        double p1 = A.at(0,1);
        double q1 = A.at(0,2);
        double k1 = A.at(0,3);
        double r1 = b.at(0);

        double z2 = A.at(1,0);
        double p2 = A.at(1,1);
        double q2 = A.at(1,2);
        double k2 = A.at(1,3);
        double r2 = b.at(1);

        //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, z1, p1, q1, k1, r1);
        //printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, z2, p2, q2, k2, r2);
        //IPrinter::printSeperatorLine();

        for (unsigned int i=2; i<=N-4; i++)
        {
            double g1 = z1/A.at(i,i-2);
            z1 = -A.at(i,i-1) * g1 + p1;
            p1 = -A.at(i,i+0) * g1 + q1;
            q1 = -A.at(i,i+1) * g1 + k1;
            k1 = -A.at(i,i+2) * g1;
            r1 = -b.at(i) * g1 + r1;

            double g2 = z2/A.at(i,i-2);
            z2 = -A.at(i,i-1) * g2 + p2;
            p2 = -A.at(i,i+0) * g2 + q2;
            q2 = -A.at(i,i+1) * g2 + k2;
            k2 = -A.at(i,i+2) * g2;
            r2 = -b.at(i) * g2 + r2;

            printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f G1:%18.10f\n", i, z1, p1, q1, k1, r1, g1);
            printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f G2:%18.10f\n", i, z2, p2, q2, k2, r2, g2);

            double norm1 = z1;//sqrt(z1*z1 + p1*p1 + q1*q1 + k1*k1);
            z1 /= norm1;
            p1 /= norm1;
            q1 /= norm1;
            k1 /= norm1;
            r1 /= norm1;
            double norm2 = z2;//sqrt(z2*z2 + p2*p2 + q2*q2 + k2*k2);
            z2 /= norm2;
            p2 /= norm2;
            q2 /= norm2;
            k2 /= norm2;
            r2 /= norm2;

            printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f N1: %14.10f\n", i, z1, p1, q1, k1, r1, norm1);
            printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f N2: %14.10f\n", i, z2, p2, q2, k2, r2, norm2);
            IPrinter::printSeperatorLine();
            // if (i>5) break;
        }

        //printf("%14.10f %14.10f\n", z1*u(N-4,m)+p1*u(N-3,m)+q1*u(N-2,m)+k1*u(N-1,m), r1);

        DoubleMatrix A1(4,4);
        DoubleVector b1(4);
        A1[0][0] = z1;            A1[0][1] = p1;            A1[0][2] = q1;            A1[0][3] = k1; b1[0] = r1;
        A1[1][0] = z2;            A1[1][1] = p2;            A1[1][2] = q2;            A1[1][3] = k2; b1[1] = r2;
        A1[2][0] = A.at(N-3,N-5); A1[2][1] = A.at(N-3,N-4); A1[2][2] = A.at(N-3,N-3); A1[2][3] = A.at(N-3,N-2); b1[2] = b.at(N-3);
        A1[3][0] = A.at(N-2,N-5); A1[3][1] = A.at(N-2,N-4); A1[3][2] = A.at(N-2,N-3); A1[3][3] = A.at(N-2,N-2); b1[3] = b.at(N-2);
        printf("det A %14.10f\n", A1.determinant());

        DoubleVector x(4);
        GaussianElimination(A1, b1, x);
        printf("%18.10f %18.10f %18.10f %18.10f\n", x[0], x[1], x[2], x[3]);

        for (unsigned int i=0; i<x.size(); i++)
        {
            u1.at(m,i+1) = x.at(i);
        }

        A.clear();
        b.clear();
        break;
    }
}
*/
