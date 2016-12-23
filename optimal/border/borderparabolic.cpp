#include "borderparabolic.h"

void BorderParabolic::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    DoubleMatrix u;
    BorderParabolic bp;
    bp.calculateU(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
    IPrinter::printMatrix(u, 10, 10, NULL);
    IPrinter::printSeperatorLine();

    DoubleMatrix u1;
    bp.calculateN42(u1);
    //IPrinter::printMatrix(u1, 10, 10, NULL);
}

double BorderParabolic::initial(unsigned int i UNUSED_PARAM) const
{
    double x = i*hx;
    return x*x;
}

double BorderParabolic::boundary(Boundary type UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t = j*ht;
    if (type == Left) return t*t;
    if (type == Right) return t*t + 1.0;
    return 0.0;
}

double BorderParabolic::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t = j*ht;
    return 2.0*t - 2.0*a*a;
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

        //        if (m<=2)
        {
            //            A.at(0,0) = +208.0*alpha;
            //            A.at(0,1) = -228.0*alpha;
            //            A.at(0,2) = +112.0*alpha;
            //            A.at(0,3) = -22.0*alpha;
            //            b.at(0) = u.at(m-1,0) + (70.0*alpha-1.0)*u.at(m,0) + ht*f(0,m);

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

            //            A.at(N,N-5) = -22.0*alpha;
            //            A.at(N,N-4) = +112.0*alpha;
            //            A.at(N,N-3) = -228.0*alpha;
            //            A.at(N,N-2) = +208.0*alpha;
            //            b.at(N-2) = u.at(m-1,N) + (70.0*alpha-1.0)*u.at(m,N) + ht*f(N,m);

            DoubleVector x(N-1);
            GaussianElimination(A,b,x);
            for (unsigned int i=0; i<x.size(); i++)
            {
                u.at(m,i+1) = x.at(i);
            }

            A.clear();
            b.clear();
            x.clear();
            //                        IPrinter::printVector(18,14,u.row(m));
            //                        IPrinter::printSeperatorLine();
        }
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

void BorderParabolic::calculateN6(DoubleMatrix &u)
{

}
