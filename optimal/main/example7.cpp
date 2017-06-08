#include "example7.h"

void Example7::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Example7 e;
    e.calculateK4_L_2_R();
//    e.calculateK4_R_2_L();
}

Example7::Example7()
{
    h = 0.01;
    N = 100;
    K = 4;
    w = 14;
    p = 10;
}

void Example7::calculateK4_L_2_R()
{
    DoubleVector rx(N+1);
    for (unsigned int k=0; k<=N; k++) rx[k] = f(k);
    IPrinter::printVector(w,p,rx);

    DoubleMatrix betta(N+1, N+1);
    DoubleVector eta(N+1);

    DoubleMatrix alpha(N+1, K+1);

    betta[0][0] = 1.0;
    betta[0][N/2] = 3.5;
    betta[0][N-4] = +4.5;
    betta[0][N-3] = +7.5;
    betta[0][N-2] = +1.5;
    betta[0][N-1] = +8.5;
    betta[0][N-0] = -2.5;
    eta[0] = 0.0; for (unsigned int k=0; k<=N; k++) eta[0] += betta[0][k]*rx[k];

    for (unsigned int k=0; k<=N-K; k++)
    {
        alpha[k][1] = +8.0;
        alpha[k][2] = +12.0*h*a(k+2);
        alpha[k][3] = -8.0;
        alpha[k][4] = +1.0;
        alpha[k][0] = +12.0*h*b(k+2);
    }

    IPrinter::printSeperatorLine();
    //IPrinter::print(alpha);

    for (unsigned int k=0; k<=N-K-1; k++)
    {
        betta[k+1][k+1] = betta[k][k+1] + alpha[k][1] /* betta[k][k]*/;
        betta[k+1][k+2] = betta[k][k+2] + alpha[k][2] /* betta[k][k]*/; betta[k+1][k+2] /= betta[k+1][k+1];
        betta[k+1][k+3] = betta[k][k+3] + alpha[k][3] /* betta[k][k]*/; betta[k+1][k+3] /= betta[k+1][k+1];
        betta[k+1][k+4] = betta[k][k+4] + alpha[k][4] /* betta[k][k]*/; betta[k+1][k+4] /= betta[k+1][k+1];
        for (unsigned int k2=k+5; k2<=N; k2++)
        {
            betta[k+1][k2] = betta[k][k2];
            betta[k+1][k2] /= betta[k+1][k+1];
        }
        eta[k+1]        = eta[k] - alpha[k][0] * betta[k][k];
        eta[k+1] /= betta[k+1][k+1];
        betta[k+1][k+1] = 1.0;
    }

    IPrinter::printSeperatorLine();
    //IPrinter::print(betta);
    IPrinter::printSeperatorLine();

    DoubleMatrix M(K+1,K+1);
    DoubleVector A(K+1);
    DoubleVector x(K+1);

    M[0][0] = betta[N-K][N-4];
    M[0][1] = betta[N-K][N-3];
    M[0][2] = betta[N-K][N-2];
    M[0][3] = betta[N-K][N-1];
    M[0][4] = betta[N-K][N-0];
    A[0] = eta[N-K];

//    M[0][0] = +1.0;
//    M[0][1] = -8.0;
//    M[0][2] = -12.0*h*a(N-2);
//    M[0][3] = +8.0;
//    M[0][4] = -1.0;
//    A[0] = 12.0*h*b(N-2);

    M[1][0] = -25.0 - 12.0*h*a(N-4);
    M[1][1] = +48.0;
    M[1][2] = -36.0;
    M[1][3] = +16.0;
    M[1][4] = -3.0;
    A[1] = 12.0*h*b(N-4);

    M[2][0] = -3.0;
    M[2][1] = -10.0 - 12.0*h*a(N-3);
    M[2][2] = +18.0;
    M[2][3] = -6.0;
    M[2][4] = +1.0;
    A[2] = 12.0*h*b(N-3);

    M[3][0] = -1.0;
    M[3][1] = +6.0;
    M[3][2] = -18.0;
    M[3][3] = +10.0 - 12.0*h*a(N-1);
    M[3][4] = +3.0;
    A[3] = 12.0*h*b(N-1);

    M[4][0] = +3.0;
    M[4][1] = -16.0;
    M[4][2] = +36.0;
    M[4][3] = -48.0;
    M[4][4] = +25.0 - 12.0*h*a(N-0);
    A[4] = 12.0*h*b(N-0);

    printf("det: %.10f\n", M.determinant());
    printf("%18.10f %18.10f\n", M[0][0]*rx[N-4] + M[0][1]*rx[N-3] + M[0][2]*rx[N-2] + M[0][3]*rx[N-1] + M[0][4]*rx[N-0], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*rx[N-4] + M[1][1]*rx[N-3] + M[1][2]*rx[N-2] + M[1][3]*rx[N-1] + M[1][4]*rx[N-0], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*rx[N-4] + M[2][1]*rx[N-3] + M[2][2]*rx[N-2] + M[2][3]*rx[N-1] + M[2][4]*rx[N-0], A[2]);
    printf("%18.10f %18.10f\n", M[3][0]*rx[N-4] + M[3][1]*rx[N-3] + M[3][2]*rx[N-2] + M[3][3]*rx[N-1] + M[3][4]*rx[N-0], A[3]);
    printf("%18.10f %18.10f\n", M[4][0]*rx[N-4] + M[4][1]*rx[N-3] + M[4][2]*rx[N-2] + M[4][3]*rx[N-1] + M[4][4]*rx[N-0], A[4]);

    IPrinter::printSeperatorLine();
//    IPrinter::print(M,5,5);
    IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    printf("%.10f %.10f %.10f %.10f %.10f\n", x[0], x[1], x[2], x[3], x[4]);

    printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2] + M[0][3]*x[3] + M[0][4]*x[4], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2] + M[1][3]*x[3] + M[1][4]*x[4], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2] + M[2][3]*x[3] + M[2][4]*x[4], A[2]);
    printf("%18.10f %18.10f\n", M[3][0]*x[0] + M[3][1]*x[1] + M[3][2]*x[2] + M[3][3]*x[3] + M[3][4]*x[4], A[3]);
    printf("%18.10f %18.10f\n", M[4][0]*x[0] + M[4][1]*x[1] + M[4][2]*x[2] + M[4][3]*x[3] + M[4][4]*x[4], A[4]);

    DoubleVector nx(N+1);
    nx[N-4] = rx[N-4];
    nx[N-3] = rx[N-3];
    nx[N-2] = rx[N-2];
    nx[N-1] = rx[N-1];
    nx[N-0] = rx[N-0];

    for (unsigned int k1=N-5; k1!=0; k1--)
    {
        for (unsigned int k2=N; k2!=0; k2--)
        {
            nx[k1] -= betta[k1][k2];
        }
        nx[k1] += eta[k1];
    }
    IPrinter::printVector(w,p,nx);
}

void Example7::calculateK4_R_2_L()
{
    DoubleVector rx(N+1);
    for (unsigned int k=0; k<=N; k++) rx[k] = f(k);
    IPrinter::printVector(w,p,rx);

    DoubleMatrix betta(N+1, N+1);
    DoubleVector eta(N+1);

    DoubleMatrix alpha(N+1, K+1);

    betta[0][0] = 1.0;
    betta[0][N] = -2.5;
    eta[0] = 0.0; for (unsigned int k=0; k<=N; k++) eta[0] += betta[0][k]*rx[k];

    for (unsigned int k=0; k<=N-K; k++)
    {
        alpha[k][1] = +8.0;
        alpha[k][2] = +12.0*h*a(k+2);
        alpha[k][3] = -8.0;
        alpha[k][4] = +1.0;
        alpha[k][0] = +12.0*h*b(k+2);
    }

    IPrinter::printSeperatorLine();
    //IPrinter::print(alpha);

    for (unsigned int k=0; k<=N-K-1; k++)
    {
        betta[k+1][k+1] = betta[k][k+1] + alpha[k][1] /* betta[k][k]*/;
        betta[k+1][k+2] = betta[k][k+2] + alpha[k][2] /* betta[k][k]*/; betta[k+1][k+2] /= betta[k+1][k+1];
        betta[k+1][k+3] = betta[k][k+3] + alpha[k][3] /* betta[k][k]*/; betta[k+1][k+3] /= betta[k+1][k+1];
        betta[k+1][k+4] = betta[k][k+4] + alpha[k][4] /* betta[k][k]*/; betta[k+1][k+4] /= betta[k+1][k+1];
        for (unsigned int k2=k+5; k2<=N; k2++)
        {
            betta[k+1][k2] = betta[k][k2];
            betta[k+1][k2] /= betta[k+1][k+1];
        }
        eta[k+1]        = eta[k] - alpha[k][0] * betta[k][k];
        eta[k+1] /= betta[k+1][k+1];
        betta[k+1][k+1] = 1.0;
    }

    IPrinter::printSeperatorLine();
    //IPrinter::print(betta);
    IPrinter::printSeperatorLine();

    DoubleMatrix M(K+1,K+1);
    DoubleVector A(K+1);
    DoubleVector x(K+1);

//    M[0][0] = betta[N-K][N-4];
//    M[0][1] = betta[N-K][N-3];
//    M[0][2] = betta[N-K][N-2];
//    M[0][3] = betta[N-K][N-1];
//    M[0][4] = betta[N-K][N-0];
//    A[0] = eta[N-K];

    M[0][0] = +1.0;
    M[0][1] = -8.0;
    M[0][2] = -12.0*h*a(2);
    M[0][3] = +8.0;
    M[0][4] = -1.0;
    A[0] = 12.0*h*b(2);

    M[1][0] = -25.0 - 12.0*h*a(0);
    M[1][1] = +48.0;
    M[1][2] = -36.0;
    M[1][3] = +16.0;
    M[1][4] = -3.0;
    A[1] = 12.0*h*b(0);

    M[2][0] = -3.0;
    M[2][1] = -10.0 - 12.0*h*a(1);
    M[2][2] = +18.0;
    M[2][3] = -6.0;
    M[2][4] = +1.0;
    A[2] = 12.0*h*b(1);

    M[3][0] = -1.0;
    M[3][1] = +6.0;
    M[3][2] = -18.0;
    M[3][3] = +10.0 - 12.0*h*a(3);
    M[3][4] = +3.0;
    A[3] = 12.0*h*b(3);

    M[4][0] = +3.0;
    M[4][1] = -16.0;
    M[4][2] = +36.0;
    M[4][3] = -48.0;
    M[4][4] = +25.0 - 12.0*h*a(4);
    A[4] = 12.0*h*b(4);

    printf("det: %.10f\n", M.determinant());
    printf("%18.10f %18.10f\n", M[0][0]*rx[0] + M[0][1]*rx[1] + M[0][2]*rx[2] + M[0][3]*rx[3] + M[0][4]*rx[4], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*rx[0] + M[1][1]*rx[1] + M[1][2]*rx[2] + M[1][3]*rx[3] + M[1][4]*rx[4], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*rx[0] + M[2][1]*rx[1] + M[2][2]*rx[2] + M[2][3]*rx[3] + M[2][4]*rx[4], A[2]);
    printf("%18.10f %18.10f\n", M[3][0]*rx[0] + M[3][1]*rx[1] + M[3][2]*rx[2] + M[3][3]*rx[3] + M[3][4]*rx[4], A[3]);
    printf("%18.10f %18.10f\n", M[4][0]*rx[0] + M[4][1]*rx[1] + M[4][2]*rx[2] + M[4][3]*rx[3] + M[4][4]*rx[4], A[4]);

    IPrinter::printSeperatorLine();
//    IPrinter::print(M,5,5);
    IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);
//    gaussianElimination(M.data(), A.data(), x.data(), x.size());

    printf("%.10f %.10f %.10f %.10f %.10f\n", x[0], x[1], x[2], x[3], x[4]);

    printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2] + M[0][3]*x[3] + M[0][4]*x[4], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2] + M[1][3]*x[3] + M[1][4]*x[4], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2] + M[2][3]*x[3] + M[2][4]*x[4], A[2]);
    printf("%18.10f %18.10f\n", M[3][0]*x[0] + M[3][1]*x[1] + M[3][2]*x[2] + M[3][3]*x[3] + M[3][4]*x[4], A[3]);
    printf("%18.10f %18.10f\n", M[4][0]*x[0] + M[4][1]*x[1] + M[4][2]*x[2] + M[4][3]*x[3] + M[4][4]*x[4], A[4]);

    DoubleVector nx(N+1);
    nx[N-4] = rx[N-4];
    nx[N-3] = rx[N-3];
    nx[N-2] = rx[N-2];
    nx[N-1] = rx[N-1];
    nx[N-0] = rx[N-0];

    for (unsigned int k1=N-5; k1!=0; k1--)
    {
        for (unsigned int k2=N; k2!=0; k2--)
        {
            nx[k1] -= betta[k1][k2];
        }
        nx[k1] += eta[k1];
    }
    IPrinter::printVector(w,p,nx);
}

double Example7::a(unsigned int k UNUSED_PARAM) const
{
    //double t = k*h;
    return 2.0;
}

double Example7::b(unsigned int k) const
{
    double t = k*h;
    return -2.0*sin(10.0*t) + 10.0*cos(10.0*t);
}

double Example7::f(unsigned int k) const
{
    double t = k*h;
    return sin(10.0*t);
}
