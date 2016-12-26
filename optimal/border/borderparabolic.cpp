#include "borderparabolic.h"
#include <math.h>

void BorderParabolic::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    DoubleMatrix u;
    BorderParabolic bp;
    bp.calculateU(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
    IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    IPrinter::printSeperatorLine();

    DoubleMatrix u1;
    bp.calculateN41(u1);
    IPrinter::printMatrix(14, 10, u1, 10, 10, NULL);
    IPrinter::printSeperatorLine();

//    DoubleMatrix u2;
//    bp.calculateN44(u2);
//    IPrinter::printMatrix(14, 10, u2, 10, 10, NULL);
//    IPrinter::printSeperatorLine();
}

double BorderParabolic::initial(unsigned int i UNUSED_PARAM) const
{
    double x = i*hx;
    return x*x;
    //return sin(x);
    //return x*x*x;
}

double BorderParabolic::boundary(Boundary type UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t = j*ht;
    if (type == Left)  return t*t;
    if (type == Right) return t*t + 1.0;
    //if (type == Left)  return exp(t)-1.0;
    //if (type == Right) return sin(1.0)+exp(t)-cos(2.0*t);
    //if (type == Left)  return t*t*t;
    //if (type == Right) return 1.0 + 2.0*t + t*t*t;
    return 0.0;
}

double BorderParabolic::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t UNUSED_PARAM = j*ht;
    double x UNUSED_PARAM = i*hx;
    return 2.0*t - 2.0*a*a;
    //return exp(t) - a*a*sin(x) + 2.0*x*sin(2.0*t*x) - 4.0*a*a*t*t*cos(2.0*x*t);
    //return 2.0*x + 3.0*t*t - 6.0*a*a*x;
}

double BorderParabolic::u(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    return x*x + t*t;
}

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

        printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, z1, p1, q1, k1, r1);
//        printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", 0, z2, p2, q2, k2, r2);
        IPrinter::printSeperatorLine();

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
            printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", i, z1, p1, q1, k1, r1);

            p1 /= z1;
            q1 /= z1;
            k1 /= z1;
            r1 /= z1;
            z1 /= z1;
            printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", i, z1, p1, q1, k1, r1);

//            p2 = -p2/z2;
//            q2 = -q2/z2;
//            k2 = -k2/z2;
//            r2 = +r2/z2;

            z2 = -A.at(i,i-1) + p2;
            p2 = -A.at(i,i+0) + q2;
            q2 = -A.at(i,i+1) + k2;
            k2 = -A.at(i,i+2);
            r2 = -b.at(i) + r2;
            printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", i, z2, p2, q2, k2, r2);

            p2 /= z2;
            q2 /= z2;
            k2 /= z2;
            r2 /= z2;
            z2 /= z2;
            printf("%4d %18.10f %18.10f %18.10f %18.10f %18.10f\n", i, z2, p2, q2, k2, r2);

            IPrinter::printSeperatorLine();
if (i>15)
            break;
        }

        printf("%14.10f %14.10f\n", z1*u(N-4,m)+p1*u(N-3,m)+q1*u(N-2,m)+k1*u(N-1,m), r1);

//        DoubleMatrix A1(4,4);
//        DoubleVector b1(4);
//        A1[0][0] = -z1;            A1[0][1] = -p1;            A1[0][2] = -q1;            A1[0][3] = -k1; b1[0] = r1;
//        A1[1][0] = -z2;            A1[1][1] = -p2;            A1[1][2] = -q2;            A1[1][3] = -k2; b1[1] = r2;
//        A1[2][0] = A.at(N-3,N-5); A1[2][1] = A.at(N-3,N-4); A1[2][2] = A.at(N-3,N-3); A1[2][3] = A.at(N-3,N-2); b1[2] = b.at(N-3);
//        A1[3][0] = A.at(N-2,N-5); A1[3][1] = A.at(N-2,N-4); A1[3][2] = A.at(N-2,N-3); A1[3][3] = A.at(N-2,N-2); b1[3] = b.at(N-2);
//        DoubleVector x(4);
//        GaussianElimination(A1, b1, x);
//        printf("%18.10f %18.10f %18.10f %18.10f\n", x[0], x[1], x[2], x[3]);

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

void BorderParabolic::calculateN6(DoubleMatrix &u UNUSED_PARAM)
{

}
