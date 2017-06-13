#include "singledifequ.h"

void SingleDifEquation::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    SingleDifEquation e;

    DoubleVector rx;
    e.calculateRX(rx);

    e.calculate2R2LV1(rx);
//    e.calculate2R2LV3(rx);

//    e.calculate4R2LV1(rx);
//    e.calculate6R2LV1(rx);
}

SingleDifEquation::SingleDifEquation()
{
    h = 0.01;
    N = 100;
    w = 14;
    p = 10;
}

void SingleDifEquation::calculateRX(DoubleVector &rx)
{
    rx.clear();
    rx.resize(N+1);
    for (unsigned int k=0; k<=N; k++) rx[k] = f(k);
    IPrinter::printVector(w,p,rx);

    FILE *file1 =fopen("data_rx.txt", "w");
    IPrinter::printVector(14,10,rx,"rx0",rx.size(),0,0,file1);
    fclose(file1);
}

void SingleDifEquation::calculate2R2LV1(const DoubleVector &rx)
{
    DoubleVector betta(N+1);
    DoubleMatrix alpha(N+1, 3);

    betta[0] = +1.0;
    betta[N/2] = +2.5;
    betta[N] = +1.5;
    double eta = betta[0]*rx[0]+betta[N/2]*rx[N/2]+betta[N]*rx[N];

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
    }

    DoubleMatrix M(3,3);
    DoubleVector A(3);
    DoubleVector x(3);

    M[0][0] = 0.0;
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

    //printf("det: %.10f\n", M.determinant());
    //printf("%18.10f %18.10f\n", M[0][0]*rx[N-2] + M[0][1]*rx[N-1] + M[0][2]*rx[N-0], A[0]);
    //printf("%18.10f %18.10f\n", M[1][0]*rx[N-2] + M[1][1]*rx[N-1] + M[1][2]*rx[N-0], A[1]);
    //printf("%18.10f %18.10f\n", M[2][0]*rx[N-2] + M[2][1]*rx[N-1] + M[2][2]*rx[N-0], A[2]);

    //IPrinter::printSeperatorLine();
    //IPrinter::print(M,3,3);
    //IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    //printf("%.10f %.10f %.10f\n", x[0], x[1], x[2]);

    //IPrinter::printSeperatorLine();
    //printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2], A[0]);
    //printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2], A[1]);
    //printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2], A[2]);

    DoubleVector nx(N+1);
    nx[N-0] = x[2];
    nx[N-1] = x[1];
    nx[N-2] = x[0];

    //IPrinter::printSeperatorLine();
    //printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
    //printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
    //printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);

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
    }

    IPrinter::printVector(w,p,nx);

    FILE *file =fopen("data_rx.txt", "a");
    IPrinter::printVector(14,10,nx,"nv2",nx.size(),0,0,file);
    fclose(file);
}

void SingleDifEquation::calculate2R2LV3(const DoubleVector &rx)
{
    DoubleVector betta(N+1);
    DoubleMatrix alpha(N+1, 3);

    betta[0] = +1.0;
    betta[N/2] = +2.5;
    betta[N] = +1.5;
    double eta = betta[0]*rx[0]+betta[N/2]*rx[N/2]+betta[N]*rx[N];

    for (unsigned int k=0; k<=N-2; k++)
    {
        alpha[k][1] = +4.0;
        alpha[k][2] = -3.0 + 2.0*h*a(k+2);
        alpha[k][0] = +2.0*h*b(k+2);
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-2; k++)
    {
        betta[k+1] = betta[k+1] + alpha[k][1]*betta[k];
        betta[k+2] = betta[k+2] + alpha[k][2]*betta[k];
        eta        = eta        - alpha[k][0]*betta[k];
    }

    DoubleMatrix M(3,3);
    DoubleVector A(3);
    DoubleVector x(3);

    M[0][0] = -3.0 - 2.0*h*a(N-2);
    M[0][1] = +4.0;
    M[0][2] = -1.0;
    A[0] = 2.0*h*b(N-2);

    M[1][0] = +1.0;
    M[1][1] = +2.0*h*a(N-1);
    M[1][2] = -1.0;
    A[1] = -2.0*h*b(N-1);

    M[2][0] = 0.0;
    M[2][1] = betta[N-1];
    M[2][2] = betta[N-0];
    A[2] = eta;

    printf("det: %.10f\n", M.determinant());
    printf("%18.10f %18.10f\n", M[0][0]*rx[N-2] + M[0][1]*rx[N-1] + M[0][2]*rx[N-0], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*rx[N-2] + M[1][1]*rx[N-1] + M[1][2]*rx[N-0], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*rx[N-2] + M[2][1]*rx[N-1] + M[2][2]*rx[N-0], A[2]);

    //IPrinter::printSeperatorLine();
    //IPrinter::print(M,3,3);
    //IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    printf("%.10f %.10f %.10f\n", x[0], x[1], x[2]);
    return;

    IPrinter::printSeperatorLine();
    printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2], A[0]);
    printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2], A[1]);
    printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2], A[2]);

    DoubleVector nx(N+1);
    nx[N-0] = x[2];
    nx[N-1] = x[1];
    nx[N-2] = x[0];

    //IPrinter::printSeperatorLine();
    //printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
    //printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
    //printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);

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
    }

    IPrinter::printVector(w,p,nx);

    FILE *file =fopen("data_rx.txt", "a");
    IPrinter::printVector(14,10,nx,"nv2",nx.size(),0,0,file);
    fclose(file);
}

void SingleDifEquation::calculate4R2LV1(const DoubleVector &rx)
{
    DoubleVector betta(N+1, 0.0);
    DoubleMatrix alpha(N+1, 5);

    betta[0] = +1.0;
    betta[N/2] = +2.5;
    betta[N] = +1.5;
    double eta = betta[0]*rx[0]+betta[N/2]*rx[N/2]+betta[N]*rx[N];

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
        betta[k+1] = betta[k+1] + alpha[k][1]*betta[k];
        betta[k+2] = betta[k+2] + alpha[k][2]*betta[k];
        betta[k+3] = betta[k+3] + alpha[k][3]*betta[k];
        betta[k+4] = betta[k+4] + alpha[k][4]*betta[k];
        eta        = eta        - alpha[k][0]*betta[k];
    }

    DoubleMatrix M(5,5);
    DoubleVector A(5);
    DoubleVector x(5);

    M[0][0] = 0.0;
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

    //printf("det: %.10f\n", M.determinant());
    //printf("%18.10f %18.10f\n", M[0][0]*rx[N-4] + M[0][1]*rx[N-3] + M[0][2]*rx[N-2] + M[0][3]*rx[N-1] + M[0][4]*rx[N-0], A[0]);
    //printf("%18.10f %18.10f\n", M[1][0]*rx[N-4] + M[1][1]*rx[N-3] + M[1][2]*rx[N-2] + M[1][3]*rx[N-1] + M[1][4]*rx[N-0], A[1]);
    //printf("%18.10f %18.10f\n", M[2][0]*rx[N-4] + M[2][1]*rx[N-3] + M[2][2]*rx[N-2] + M[2][3]*rx[N-1] + M[2][4]*rx[N-0], A[2]);
    //printf("%18.10f %18.10f\n", M[3][0]*rx[N-4] + M[3][1]*rx[N-3] + M[3][2]*rx[N-2] + M[3][3]*rx[N-1] + M[3][4]*rx[N-0], A[3]);
    //printf("%18.10f %18.10f\n", M[4][0]*rx[N-4] + M[4][1]*rx[N-3] + M[4][2]*rx[N-2] + M[4][3]*rx[N-1] + M[4][4]*rx[N-0], A[4]);

    //IPrinter::printSeperatorLine();
    //IPrinter::print(M,5,5);
    //IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    //printf("%.10f %.10f %.10f %.10f %.10f\n", x[0], x[1], x[2], x[3], x[4]);

    //IPrinter::printSeperatorLine();
    //printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2] + M[0][3]*x[3] + M[0][4]*x[4], A[0]);
    //printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2] + M[1][3]*x[3] + M[1][4]*x[4], A[1]);
    //printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2] + M[2][3]*x[3] + M[2][4]*x[4], A[2]);
    //printf("%18.10f %18.10f\n", M[3][0]*x[0] + M[3][1]*x[1] + M[3][2]*x[2] + M[3][3]*x[3] + M[3][4]*x[4], A[3]);
    //printf("%18.10f %18.10f\n", M[4][0]*x[0] + M[4][1]*x[1] + M[4][2]*x[2] + M[4][3]*x[3] + M[4][4]*x[4], A[4]);

    DoubleVector nx(N+1);
    nx[N-0] = x[4];
    nx[N-1] = x[3];
    nx[N-2] = x[2];
    nx[N-3] = x[1];
    nx[N-4] = x[0];

    //IPrinter::printSeperatorLine();
    //printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
    //printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
    //printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);
    //printf("%d %.10f %.10f\n", N-3, rx[N-3], nx[N-3]);
    //printf("%d %.10f %.10f\n", N-4, rx[N-4], nx[N-4]);

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
    }

    IPrinter::printVector(w,p,nx);

    FILE *file =fopen("data_rx.txt", "a");
    IPrinter::printVector(14,10,nx,"nv4",nx.size(),0,0,file);
    fclose(file);
}

void SingleDifEquation::calculate6R2LV1(const DoubleVector &rx)
{
    DoubleVector betta(N+1, 0.0);
    DoubleMatrix alpha(N+1, 7);

    betta[0] = +1.0;
    betta[N/2] = +2.5;
    betta[N] = +1.5;
    double eta = betta[0]*rx[0]+betta[N/2]*rx[N/2]+betta[N]*rx[N];

    for (unsigned int k=0; k<=N-6; k++)
    {
        double m = -147.0 - 60.0*h*a(k);
        alpha[k][1] = -360.0/m;
        alpha[k][2] = +450.0/m;
        alpha[k][3] = -400.0/m;
        alpha[k][4] = +225.0/m;
        alpha[k][5] = -72.0/m;
        alpha[k][6] = +10.0/m;
        alpha[k][0] = +60.0*h*b(k)/m;
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-6; k++)
    {
        betta[k+1] = betta[k+1] + alpha[k][1]*betta[k];
        betta[k+2] = betta[k+2] + alpha[k][2]*betta[k];
        betta[k+3] = betta[k+3] + alpha[k][3]*betta[k];
        betta[k+4] = betta[k+4] + alpha[k][4]*betta[k];
        betta[k+5] = betta[k+5] + alpha[k][5]*betta[k];
        betta[k+6] = betta[k+6] + alpha[k][6]*betta[k];
        eta        = eta        - alpha[k][0]*betta[k];
    }

    DoubleMatrix M(7,7);
    DoubleVector A(7);
    DoubleVector x(7);

    M[0][0] = 0.0;
    M[0][1] = betta[N-5];
    M[0][2] = betta[N-4];
    M[0][3] = betta[N-3];
    M[0][4] = betta[N-2];
    M[0][5] = betta[N-1];
    M[0][6] = betta[N-0];
    A[0]    = eta;

    M[1][0] = -10.0;
    M[1][1] = -77.0 - 60.0*h*a(N-5);
    M[1][2] = +150.0;
    M[1][3] = -100.0;
    M[1][4] = +50.0;
    M[1][5] = -15.0;
    M[1][6] = +2.0;
    A[1] = 60.0*h*b(N-5);

    M[2][0] = +2.0;
    M[2][1] = -24.0;
    M[2][2] = -35.0 - 60.0*h*a(N-4);
    M[2][3] = +80.0;
    M[2][4] = -30.0;
    M[2][5] = +8.0;
    M[2][6] = -1.0;
    A[2] = 60.0*h*b(N-4);

    M[3][0] = -1.0;
    M[3][1] = +9.0;
    M[3][2] = -45.0;
    M[3][3] = -60.0*h*a(N-3);
    M[3][4] = +45.0;
    M[3][5] = -9.0;
    M[3][6] = +1.0;
    A[3] = 60.0*h*b(N-3);

    M[4][0] = +1.0;
    M[4][1] = -8.0;
    M[4][2] = +30.0;
    M[4][3] = -80.0;
    M[4][4] = +35.0 - 60.0*h*a(N-2);
    M[4][5] = +24.0;
    M[4][6] = -2.0;
    A[4] = 60.0*h*b(N-2);

    M[5][0] = -2.0;
    M[5][1] = +15.0;
    M[5][2] = -50.0;
    M[5][3] = +100.0;
    M[5][4] = -150.0;
    M[5][5] = +77.0 - 60.0*h*a(N-1);
    M[5][6] = +10.0;
    A[5] = 60.0*h*b(N-1);

    M[6][0] = +10.0;
    M[6][1] = -72.0;
    M[6][2] = +225.0;
    M[6][3] = -400.0;
    M[6][4] = +450.0;
    M[6][5] = -360.0;
    M[6][6] = +147.0 - 60.0*h*a(N-0);
    A[6] = 60.0*h*b(N-0);

    //printf("det: %.10f\n", M.determinant());
    //printf("%18.10f %18.10f\n", M[0][0]*rx[N-6] + M[0][1]*rx[N-5] + M[0][2]*rx[N-4] + M[0][3]*rx[N-3] + M[0][4]*rx[N-2] + M[0][5]*rx[N-1] + M[0][6]*rx[N-0], A[0]);
    //printf("%18.10f %18.10f\n", M[1][0]*rx[N-6] + M[1][1]*rx[N-5] + M[1][2]*rx[N-4] + M[1][3]*rx[N-3] + M[1][4]*rx[N-2] + M[1][5]*rx[N-1] + M[1][6]*rx[N-0], A[1]);
    //printf("%18.10f %18.10f\n", M[2][0]*rx[N-6] + M[2][1]*rx[N-5] + M[2][2]*rx[N-4] + M[2][3]*rx[N-3] + M[2][4]*rx[N-2] + M[2][5]*rx[N-1] + M[2][6]*rx[N-0], A[2]);
    //printf("%18.10f %18.10f\n", M[3][0]*rx[N-6] + M[3][1]*rx[N-5] + M[3][2]*rx[N-4] + M[3][3]*rx[N-3] + M[3][4]*rx[N-2] + M[3][5]*rx[N-1] + M[3][6]*rx[N-0], A[3]);
    //printf("%18.10f %18.10f\n", M[4][0]*rx[N-6] + M[4][1]*rx[N-5] + M[4][2]*rx[N-4] + M[4][3]*rx[N-3] + M[4][4]*rx[N-2] + M[4][5]*rx[N-1] + M[4][6]*rx[N-0], A[4]);
    //printf("%18.10f %18.10f\n", M[5][0]*rx[N-6] + M[5][1]*rx[N-5] + M[5][2]*rx[N-4] + M[5][3]*rx[N-3] + M[5][4]*rx[N-2] + M[5][5]*rx[N-1] + M[5][6]*rx[N-0], A[5]);
    //printf("%18.10f %18.10f\n", M[6][0]*rx[N-6] + M[6][1]*rx[N-5] + M[6][2]*rx[N-4] + M[6][3]*rx[N-3] + M[6][4]*rx[N-2] + M[6][5]*rx[N-1] + M[6][6]*rx[N-0], A[6]);

    //IPrinter::printSeperatorLine();
    //IPrinter::print(M,7,7);
    //IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    //printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]);

    //IPrinter::printSeperatorLine();
    //printf("%18.10f %18.10f\n", M[0][0]*x[N-6] + M[0][1]*x[N-5] + M[0][2]*x[N-4] + M[0][3]*x[N-3] + M[0][4]*x[N-2] + M[0][5]*x[N-1] + M[0][6]*x[N-0], A[0]);
    //printf("%18.10f %18.10f\n", M[1][0]*x[N-6] + M[1][1]*x[N-5] + M[1][2]*x[N-4] + M[1][3]*x[N-3] + M[1][4]*x[N-2] + M[1][5]*x[N-1] + M[1][6]*x[N-0], A[1]);
    //printf("%18.10f %18.10f\n", M[2][0]*x[N-6] + M[2][1]*x[N-5] + M[2][2]*x[N-4] + M[2][3]*x[N-3] + M[2][4]*x[N-2] + M[2][5]*x[N-1] + M[2][6]*x[N-0], A[2]);
    //printf("%18.10f %18.10f\n", M[3][0]*x[N-6] + M[3][1]*x[N-5] + M[3][2]*x[N-4] + M[3][3]*x[N-3] + M[3][4]*x[N-2] + M[3][5]*x[N-1] + M[3][6]*x[N-0], A[3]);
    //printf("%18.10f %18.10f\n", M[4][0]*x[N-6] + M[4][1]*x[N-5] + M[4][2]*x[N-4] + M[4][3]*x[N-3] + M[4][4]*x[N-2] + M[4][5]*x[N-1] + M[4][6]*x[N-0], A[4]);
    //printf("%18.10f %18.10f\n", M[5][0]*x[N-6] + M[5][1]*x[N-5] + M[5][2]*x[N-4] + M[5][3]*x[N-3] + M[5][4]*x[N-2] + M[5][5]*x[N-1] + M[5][6]*x[N-0], A[5]);
    //printf("%18.10f %18.10f\n", M[6][0]*x[N-6] + M[6][1]*x[N-5] + M[6][2]*x[N-4] + M[6][3]*x[N-3] + M[6][4]*x[N-2] + M[6][5]*x[N-1] + M[6][6]*x[N-0], A[6]);

    DoubleVector nx(N+1);
    nx[N-0] = x[6];
    nx[N-1] = x[5];
    nx[N-2] = x[4];
    nx[N-3] = x[3];
    nx[N-4] = x[2];
    nx[N-5] = x[1];
    nx[N-6] = x[0];

    //IPrinter::printSeperatorLine();
    //printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
    //printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
    //printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);
    //printf("%d %.10f %.10f\n", N-3, rx[N-3], nx[N-3]);
    //printf("%d %.10f %.10f\n", N-4, rx[N-4], nx[N-4]);
    //printf("%d %.10f %.10f\n", N-5, rx[N-5], nx[N-5]);
    //printf("%d %.10f %.10f\n", N-6, rx[N-6], nx[N-6]);

    betta[N-5] -= alpha[N-6][1]*betta[N-6];
    betta[N-4] -= alpha[N-6][2]*betta[N-6];
    betta[N-3] -= alpha[N-6][3]*betta[N-6];
    betta[N-2] -= alpha[N-6][4]*betta[N-6];
    betta[N-1] -= alpha[N-6][5]*betta[N-6];
    betta[N-0] -= alpha[N-6][6]*betta[N-6];
    eta        += alpha[N-6][0]*betta[N-6];

    for (unsigned int k=N-6; k>0; k--)
    {
        betta[k+0] -= alpha[k-1][1]*betta[k-1];
        betta[k+1] -= alpha[k-1][2]*betta[k-1];
        betta[k+2] -= alpha[k-1][3]*betta[k-1];
        betta[k+3] -= alpha[k-1][4]*betta[k-1];
        betta[k+4] -= alpha[k-1][5]*betta[k-1];
        betta[k+5] -= alpha[k-1][6]*betta[k-1];
        eta        += alpha[k-1][0]*betta[k-1];

        nx[k-1] = eta;
        for (unsigned int i=k; i<=N; i++)
        {
            nx[k-1] -= betta[i]*nx[i];
        }
        nx[k-1] /= betta[k-1];
    }

    IPrinter::printVector(w,p,nx);

    FILE *file =fopen("data_rx.txt", "a");
    IPrinter::printVector(14,10,nx,"nv6",nx.size(),0,0,file);
    fclose(file);
}

void SingleDifEquation::calculate2R2LV2()
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

void SingleDifEquation::calculate4R2LV2()
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

double SingleDifEquation::a(unsigned int k UNUSED_PARAM ) const
{
    double t = k*h;

#ifdef SAMPLE_1
    return 2.0;
#endif
#ifdef SAMPLE_2
    return t;
#endif
}

double SingleDifEquation::b(unsigned int k) const
{
    double t = k*h;
#ifdef SAMPLE_1
    return -2.0*sin(10.0*t) + 10.0*cos(10.0*t);
#endif
#ifdef SAMPLE_2
    return (1.0-t*t)*sin(10.0*t*t) + 20.0*t*t*cos(10.0*t*t);
#endif
}

double SingleDifEquation::f(unsigned int k) const
{
    double t = k*h;
#ifdef SAMPLE_1
    return sin(10.0*t);
#endif
#ifdef SAMPLE_2
    return t*sin(10.0*t*t);
#endif
}
