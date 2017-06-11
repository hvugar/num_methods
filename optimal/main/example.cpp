#include "example.h"

void Example::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Example e;
    e.calculate2R2LV2();
    e.calculate4R2LV2();
}

Example::Example()
{
    h = 0.01;
    N = 100;
    w = 14;
    p = 10;
}

void Example::calculate2R2LV1()
{
    DoubleVector rx(N+1);
    for (unsigned int k=0; k<=N; k++) rx[k] = f(k);
    IPrinter::printVector(w,p,rx);

    DoubleVector betta(N+1);
    betta[0] = 1.0;
    betta[N] = 3.5;
    double eta = 0.0;
    for (unsigned int i=0; i<=N; i++) eta += betta[i]*rx[i];

    DoubleMatrix alpha(N+1, 3);
    for (unsigned int k=0; k<=N-2; k++)
    {
        alpha[k][1] = -2.0*h*a(k+1);
        alpha[k][2] = +1.0;
        alpha[k][0] = -2.0*h*b(k+1);
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-2; k++)
    {
        betta[k+1] = betta[k+1] + alpha[k][1]*betta[k];
        betta[k+2] = betta[k+2] + alpha[k][2]*betta[k];
        eta        = eta        - alpha[k][0]*betta[k];

        //eta        /= betta[k+1];
        //betta[k+2] /= betta[k+1];
        //betta[k+1] = 1.0;
    }
    betta[N-2] = 0.0;

    DoubleMatrix M(3,3);
    DoubleVector A(3);
    DoubleVector x(3);

    M[0][0] = betta[N-2];
    M[0][1] = betta[N-1];
    M[0][2] = betta[N-0];
    A[0] = eta;

    M[1][0] = -3.0-2.0*h*a(N-2);
    M[1][1] = +4.0;
    M[1][2] = -1.0;
    A[1] = +2.0*h*b(N-2);

    M[2][0] = +1.0;
    M[2][1] = -4.0;
    M[2][2] = +3.0-2.0*h*a(N);
    A[2] = +2.0*h*b(N);

    GaussianElimination(M, A, x);

    printf("%.10f %.10f %.10f\n", x[0], x[1], x[2]);

    IPrinter::printSeperatorLine();
    IPrinter::print(M,3,3);
    IPrinter::printSeperatorLine();

    printf("det: %.10f\n", M.determinant());
    printf("%18.10f %18.10f\n", M[0][0]*rx[N-2] + M[0][1]*rx[N-1] + M[0][2]*rx[N-0], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*rx[N-2] + M[1][1]*rx[N-1] + M[1][2]*rx[N-0], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*rx[N-2] + M[2][1]*rx[N-1] + M[2][2]*rx[N-0], A[2]);
    IPrinter::printSeperatorLine();
    printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2], A[2]);
}

void Example::calculate2R2LV2()
{
    DoubleVector rx(N+1);
    for (unsigned int k=0; k<=N; k++) rx[k] = f(k);
    IPrinter::printVector(w,p,rx);

    FILE *file1 =fopen("data_rx.txt", "w");
    IPrinter::printVector(14,10,rx,NULL,rx.size(),0,0,file1);
    fclose(file1);

    DoubleVector betta(N+1);
    betta[0] = 1.0;
    betta[N] = 3.5;
    double eta = betta[0]*rx[0]+betta[N]*rx[N];
    //double eta = 0.0;
    //for (unsigned int i=0; i<=N; i++) eta += betta[i]*rx[i];

    DoubleMatrix alpha(N+1, 3);
    for (unsigned int k=0; k<=N-2; k++)
    {
        double m = +3.0 + 2.0*h*a(k);
        alpha[k][1] = +4.0/m;
        alpha[k][2] = -1.0/m;
        alpha[k][0] = -2.0*h*b(k)/m;
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-2; k++)
    {
        betta[k+1] = betta[k+1] + alpha[k][1]*betta[k];
        betta[k+2] = betta[k+2] + alpha[k][2]*betta[k];
        eta        = eta        - alpha[k][0]*betta[k];

        //eta        /= betta[k+1];
        //betta[k+2] /= betta[k+1];
        //betta[k+1] = 1.0;
    }
    betta[N-2] = 0.0;

    DoubleMatrix M(3,3);
    DoubleVector A(3);
    DoubleVector x(3);

    M[0][0] = betta[N-2];
    M[0][1] = betta[N-1];
    M[0][2] = betta[N-0];
    A[0] = eta;

    M[1][0] = +1.0;
    M[1][1] = +2.0*h*a(N-1);
    M[1][2] = -1.0;
    A[1] = -2.0*h*b(N-1);

    M[2][0] = +1.0;
    M[2][1] = -4.0;
    M[2][2] = +3.0-2.0*h*a(N);
    A[2] = +2.0*h*b(N);

    printf("det: %.10f\n", M.determinant());
    printf("%18.10f %18.10f\n", M[0][0]*rx[N-2] + M[0][1]*rx[N-1] + M[0][2]*rx[N-0], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*rx[N-2] + M[1][1]*rx[N-1] + M[1][2]*rx[N-0], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*rx[N-2] + M[2][1]*rx[N-1] + M[2][2]*rx[N-0], A[2]);

    IPrinter::printSeperatorLine();
    IPrinter::print(M,3,3);
    IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    printf("%.10f %.10f %.10f\n", x[0], x[1], x[2]);

    IPrinter::printSeperatorLine();
    printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2], A[2]);

    DoubleVector nx(N+1);
    nx[N-0] = x[2];
    nx[N-1] = x[1];
    nx[N-2] = x[0];

    IPrinter::printSeperatorLine();

    printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
    printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
    printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);

    betta[N-1] -= alpha[N-2][1]*betta[N-2];
    betta[N-0] -= alpha[N-2][2]*betta[N-2];
    eta        += alpha[N-2][0]*betta[N-2];

    for (unsigned int k=N-2; k>0; k--)
    {
        betta[k+0] -= alpha[k-1][1]*betta[k-1];
        betta[k+1] -= alpha[k-1][2]*betta[k-1];
        eta        += alpha[k-1][0]*betta[k-1];

        nx[k-1] = eta;
        for (unsigned int i=k; i<=N; i++)
        {
            nx[k-1] -= betta[i]*nx[i];
        }
        nx[k-1] /= betta[k-1];
        //printf("%d %.10f %.10f\n", k-1, rx[k-1], nx[k-1]);
    }

    IPrinter::printVector(w,p,nx);

    FILE *file =fopen("data_k2.txt", "w");
    IPrinter::printVector(14,10,nx,NULL,nx.size(),0,0,file);
    fclose(file);
}

void Example::calculate4R2LV1()
{
    DoubleVector rx(N+1);
    for (unsigned int k=0; k<=N; k++) rx[k] = f(k);
    IPrinter::printVector(w,p,rx);

    DoubleVector betta(N+1);

    DoubleMatrix alpha(N+1, 5);

    betta[0] = +1.0;
    betta[N] = +1.5;
    double eta = betta[0]*rx[0]+betta[N]*rx[N];
    //double eta = 0.0;
    //for (unsigned int k=0; k<=N; k++) eta += betta[0]*rx[k];

    for (unsigned int k=0; k<=N-4; k++)
    {
        alpha[k][1] = +8.0;
        alpha[k][2] = +12.0*h*a(k+2);
        alpha[k][3] = -8.0;
        alpha[k][4] = +1.0;
        alpha[k][0] = +12.0*h*b(k+2);
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-4; k++)
    {
        betta[k+1] = betta[k+1] + alpha[k][1]*betta[k];
        betta[k+2] = betta[k+2] + alpha[k][2]*betta[k];
        betta[k+3] = betta[k+3] + alpha[k][3]*betta[k];
        betta[k+4] = betta[k+4] + alpha[k][4]*betta[k];
        eta        = eta        - alpha[k][0]*betta[k];
    }
    betta[N-4] = 0.0;

    IPrinter::printSeperatorLine();
//    IPrinter::print(betta);
    IPrinter::printSeperatorLine();

    DoubleMatrix M(5,5);
    DoubleVector A(5);
    DoubleVector x(5);

    M[0][0] = betta[N-4];
    M[0][1] = betta[N-3];
    M[0][2] = betta[N-2];
    M[0][3] = betta[N-1];
    M[0][4] = betta[N-0];
    A[0]    = eta;

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
    IPrinter::print(M,5,5);
    IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    printf("%.10f %.10f %.10f %.10f %.10f\n", x[0], x[1], x[2], x[3], x[4]);

//    printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2] + M[0][3]*x[3] + M[0][4]*x[4], A[0]);
//    printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2] + M[1][3]*x[3] + M[1][4]*x[4], A[1]);
//    printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2] + M[2][3]*x[3] + M[2][4]*x[4], A[2]);
//    printf("%18.10f %18.10f\n", M[3][0]*x[0] + M[3][1]*x[1] + M[3][2]*x[2] + M[3][3]*x[3] + M[3][4]*x[4], A[3]);
//    printf("%18.10f %18.10f\n", M[4][0]*x[0] + M[4][1]*x[1] + M[4][2]*x[2] + M[4][3]*x[3] + M[4][4]*x[4], A[4]);

//    DoubleVector nx(N+1);
//    nx[N-4] = rx[N-4];
//    nx[N-3] = rx[N-3];
//    nx[N-2] = rx[N-2];
//    nx[N-1] = rx[N-1];
//    nx[N-0] = rx[N-0];

//    for (unsigned int k1=N-5; k1!=0; k1--)
//    {
//        for (unsigned int k2=N; k2!=0; k2--)
//        {
//            nx[k1] -= betta[k1][k2];
//        }
//        nx[k1] += eta[k1];
//    }
//    IPrinter::printVector(w,p,nx);
}

void Example::calculate4R2LV2()
{
    DoubleVector rx(N+1);
    for (unsigned int k=0; k<=N; k++) rx[k] = f(k);
    IPrinter::printVector(w,p,rx);

//    DoubleMatrix betta(N+1, N+1, 0.0);
    DoubleVector betta(N+1, 0.0);

    DoubleMatrix alpha(N+1, 5);

//    betta[0][0] = +1.0;
//    betta[0][N] = +1.5;
//    double eta = betta[0][0]*rx[0]+betta[0][N]*rx[N];
    //for (unsigned int k=0; k<=N; k++) eta += betta[0]*rx[k];
    betta[0] = +1.0;
    betta[N] = +1.5;
    double eta = betta[0]*rx[0]+betta[N]*rx[N];

    for (unsigned int k=0; k<=N-4; k++)
    {
        double m = -25.0 - 12.0*h*a(k);
        alpha[k][1] = -48.0/m;
        alpha[k][2] = +36.0/m;
        alpha[k][3] = -16.0/m;
        alpha[k][4] = +3.0/m;
        alpha[k][0] = +12.0*h*b(k)/m;
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-4; k++)
    {
//        betta[k+1][k+1] = betta[k+1][k+1] + alpha[k][1]*betta[k][k];
//        betta[k+1][k+2] = betta[k+1][k+2] + alpha[k][2]*betta[k][k];
//        betta[k+1][k+3] = betta[k+1][k+3] + alpha[k][3]*betta[k][k];
//        betta[k+1][k+4] = betta[k+1][k+4] + alpha[k][4]*betta[k][k];
//        eta             = eta           - alpha[k][0]*betta[k][k];

        betta[k+1] = betta[k+1] + alpha[k][1]*betta[k];
        betta[k+2] = betta[k+2] + alpha[k][2]*betta[k];
        betta[k+3] = betta[k+3] + alpha[k][3]*betta[k];
        betta[k+4] = betta[k+4] + alpha[k][4]*betta[k];
        eta        = eta        - alpha[k][0]*betta[k];
    }
//    betta[N-3][N-4] = 0.0;
//    betta[N-4] = 0.0;

    IPrinter::printSeperatorLine();
//    IPrinter::print(betta);
    IPrinter::printSeperatorLine();

    DoubleMatrix M(5,5);
    DoubleVector A(5);
    DoubleVector x(5);

//    M[0][0] = betta[N-3][N-4];
//    M[0][1] = betta[N-3][N-3];
//    M[0][2] = betta[N-3][N-2];
//    M[0][3] = betta[N-3][N-1];
//    M[0][4] = betta[N-3][N-0];
//    A[0]    = eta;
    M[0][0] = 0.0;//betta[N-4];
    M[0][1] = betta[N-3];
    M[0][2] = betta[N-2];
    M[0][3] = betta[N-1];
    M[0][4] = betta[N-0];
    A[0]    = eta;

    M[1][0] = -3.0;
    M[1][1] = -10.0 - 12.0*h*a(N-3);
    M[1][2] = +18.0;
    M[1][3] = -6.0;
    M[1][4] = +1.0;
    A[1] = 12.0*h*b(N-3);

    M[2][0] = +1.0;
    M[2][1] = -8.0;
    M[2][2] = +0.0 -12.0*h*a(N-2);
    M[2][3] = +8.0;
    M[2][4] = -1.0;
    A[2] = 12.0*h*b(N-2);

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
    IPrinter::print(M,5,5);
    IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    printf("%.10f %.10f %.10f %.10f %.10f\n", x[0], x[1], x[2], x[3], x[4]);

    printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2] + M[0][3]*x[3] + M[0][4]*x[4], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2] + M[1][3]*x[3] + M[1][4]*x[4], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2] + M[2][3]*x[3] + M[2][4]*x[4], A[2]);
    printf("%18.10f %18.10f\n", M[3][0]*x[0] + M[3][1]*x[1] + M[3][2]*x[2] + M[3][3]*x[3] + M[3][4]*x[4], A[3]);
    printf("%18.10f %18.10f\n", M[4][0]*x[0] + M[4][1]*x[1] + M[4][2]*x[2] + M[4][3]*x[3] + M[4][4]*x[4], A[4]);

    DoubleVector nx(N+1);
    nx[N-0] = x[4];
    nx[N-1] = x[3];
    nx[N-2] = x[2];
    nx[N-3] = x[1];
    nx[N-4] = x[0];

    IPrinter::printSeperatorLine();

    printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
    printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
    printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);
    printf("%d %.10f %.10f\n", N-3, rx[N-3], nx[N-3]);
    printf("%d %.10f %.10f\n", N-4, rx[N-4], nx[N-4]);

    betta[N-3] -= alpha[N-4][1]*betta[N-4];
    betta[N-2] -= alpha[N-4][2]*betta[N-4];
    betta[N-1] -= alpha[N-4][3]*betta[N-4];
    betta[N-0] -= alpha[N-4][4]*betta[N-4];
    eta        += alpha[N-4][0]*betta[N-4];

    for (unsigned int k=N-4; k>0; k--)
    {
        betta[k+0] -= alpha[k-1][1]*betta[k-1];
        betta[k+1] -= alpha[k-1][2]*betta[k-1];
        betta[k+2] -= alpha[k-1][3]*betta[k-1];
        betta[k+3] -= alpha[k-1][4]*betta[k-1];
        eta        += alpha[k-1][0]*betta[k-1];

        nx[k-1] = eta;
        for (unsigned int i=k; i<=N; i++)
        {
            nx[k-1] -= betta[i]*nx[i];
        }
        nx[k-1] /= betta[k-1];
        //printf("%d %.10f %.10f\n", k-1, rx[k-1], nx[k-1]);
    }

    IPrinter::printVector(w,p,nx);

    FILE *file =fopen("data_k4.txt", "w");
    IPrinter::printVector(14,10,nx,NULL,nx.size(),0,0,file);
    fclose(file);
}

void Example::calculate6R2LV2()
{
    DoubleVector rx(N+1);
    for (unsigned int k=0; k<=N; k++) rx[k] = f(k);
    IPrinter::printVector(w,p,rx);

//    DoubleMatrix betta(N+1, N+1, 0.0);
    DoubleVector betta(N+1, 0.0);

    DoubleMatrix alpha(N+1, 5);

//    betta[0][0] = +1.0;
//    betta[0][N] = +1.5;
//    double eta = betta[0][0]*rx[0]+betta[0][N]*rx[N];
    //for (unsigned int k=0; k<=N; k++) eta += betta[0]*rx[k];
    betta[0] = +1.0;
    betta[N] = +1.5;
    double eta = betta[0]*rx[0]+betta[N]*rx[N];

    for (unsigned int k=0; k<=N-4; k++)
    {
        double m = -25.0 - 12.0*h*a(k);
        alpha[k][1] = -48.0/m;
        alpha[k][2] = +36.0/m;
        alpha[k][3] = -16.0/m;
        alpha[k][4] = +3.0/m;
        alpha[k][0] = +12.0*h*b(k)/m;
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-4; k++)
    {
//        betta[k+1][k+1] = betta[k+1][k+1] + alpha[k][1]*betta[k][k];
//        betta[k+1][k+2] = betta[k+1][k+2] + alpha[k][2]*betta[k][k];
//        betta[k+1][k+3] = betta[k+1][k+3] + alpha[k][3]*betta[k][k];
//        betta[k+1][k+4] = betta[k+1][k+4] + alpha[k][4]*betta[k][k];
//        eta             = eta           - alpha[k][0]*betta[k][k];

        betta[k+1] = betta[k+1] + alpha[k][1]*betta[k];
        betta[k+2] = betta[k+2] + alpha[k][2]*betta[k];
        betta[k+3] = betta[k+3] + alpha[k][3]*betta[k];
        betta[k+4] = betta[k+4] + alpha[k][4]*betta[k];
        eta        = eta        - alpha[k][0]*betta[k];
    }
//    betta[N-3][N-4] = 0.0;
//    betta[N-4] = 0.0;

    IPrinter::printSeperatorLine();
//    IPrinter::print(betta);
    IPrinter::printSeperatorLine();

    DoubleMatrix M(5,5);
    DoubleVector A(5);
    DoubleVector x(5);

//    M[0][0] = betta[N-3][N-4];
//    M[0][1] = betta[N-3][N-3];
//    M[0][2] = betta[N-3][N-2];
//    M[0][3] = betta[N-3][N-1];
//    M[0][4] = betta[N-3][N-0];
//    A[0]    = eta;
    M[0][0] = 0.0;//betta[N-4];
    M[0][1] = betta[N-3];
    M[0][2] = betta[N-2];
    M[0][3] = betta[N-1];
    M[0][4] = betta[N-0];
    A[0]    = eta;

    M[1][0] = -3.0;
    M[1][1] = -10.0 - 12.0*h*a(N-3);
    M[1][2] = +18.0;
    M[1][3] = -6.0;
    M[1][4] = +1.0;
    A[1] = 12.0*h*b(N-3);

    M[2][0] = +1.0;
    M[2][1] = -8.0;
    M[2][2] = +0.0 -12.0*h*a(N-2);
    M[2][3] = +8.0;
    M[2][4] = -1.0;
    A[2] = 12.0*h*b(N-2);

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
    IPrinter::print(M,5,5);
    IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    printf("%.10f %.10f %.10f %.10f %.10f\n", x[0], x[1], x[2], x[3], x[4]);

    printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2] + M[0][3]*x[3] + M[0][4]*x[4], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2] + M[1][3]*x[3] + M[1][4]*x[4], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2] + M[2][3]*x[3] + M[2][4]*x[4], A[2]);
    printf("%18.10f %18.10f\n", M[3][0]*x[0] + M[3][1]*x[1] + M[3][2]*x[2] + M[3][3]*x[3] + M[3][4]*x[4], A[3]);
    printf("%18.10f %18.10f\n", M[4][0]*x[0] + M[4][1]*x[1] + M[4][2]*x[2] + M[4][3]*x[3] + M[4][4]*x[4], A[4]);

    DoubleVector nx(N+1);
    nx[N-0] = x[4];
    nx[N-1] = x[3];
    nx[N-2] = x[2];
    nx[N-3] = x[1];
    nx[N-4] = x[0];

    IPrinter::printSeperatorLine();

    printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
    printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
    printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);
    printf("%d %.10f %.10f\n", N-3, rx[N-3], nx[N-3]);
    printf("%d %.10f %.10f\n", N-4, rx[N-4], nx[N-4]);

    betta[N-3] -= alpha[N-4][1]*betta[N-4];
    betta[N-2] -= alpha[N-4][2]*betta[N-4];
    betta[N-1] -= alpha[N-4][3]*betta[N-4];
    betta[N-0] -= alpha[N-4][4]*betta[N-4];
    eta        += alpha[N-4][0]*betta[N-4];

    for (unsigned int k=N-4; k>0; k--)
    {
        betta[k+0] -= alpha[k-1][1]*betta[k-1];
        betta[k+1] -= alpha[k-1][2]*betta[k-1];
        betta[k+2] -= alpha[k-1][3]*betta[k-1];
        betta[k+3] -= alpha[k-1][4]*betta[k-1];
        eta        += alpha[k-1][0]*betta[k-1];

        nx[k-1] = eta;
        for (unsigned int i=k; i<=N; i++)
        {
            nx[k-1] -= betta[i]*nx[i];
        }
        nx[k-1] /= betta[k-1];
        //printf("%d %.10f %.10f\n", k-1, rx[k-1], nx[k-1]);
    }

    IPrinter::printVector(w,p,nx);

    FILE *file =fopen("data_k4.txt", "w");
    IPrinter::printVector(14,10,nx,NULL,nx.size(),0,0,file);
    fclose(file);
}

double Example::a(unsigned int k UNUSED_PARAM ) const
{
    double t = k*h;
    return t;
//    return 2.0;
}

double Example::b(unsigned int k) const
{
    double t = k*h;
    return -t*sin(10.0*t) + 10.0*cos(10.0*t);
    //return -2.0*sin(10.0*t) + 10.0*cos(10.0*t);
}

double Example::f(unsigned int k) const
{
    double t = k*h;
    return sin(10.0*t);
}
