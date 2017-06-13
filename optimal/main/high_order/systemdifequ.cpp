#include "systemdifequ.h"

void SystemDifEquation::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    SystemDifEquation e;

    DoubleMatrix rx;
    e.calculateRX(rx);
    IPrinter::printSeperatorLine();

    e.calculate2R2LV1(rx);
    e.calculate4R2LV1(rx);
//    e.calculate6R2LV1(rx);
}

SystemDifEquation::SystemDifEquation()
{
    h = 0.01;
    N = 100;
    w = 14;
    p = 10;
}

void SystemDifEquation::calculate2R2LV1(const DoubleMatrix &rx)
{
    DoubleMatrix* betta = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) betta[i].resize(3,3,0.0);

    betta[0].randomData();
    betta[N].randomData();

//    betta[0][0][0] = +1.0; betta[0][0][1] = +8.0; betta[0][0][2] = -1.0;
//    betta[0][1][0] = +5.0; betta[0][1][1] = +3.0; betta[0][1][2] = +7.0;
//    betta[0][2][0] = -2.0; betta[0][2][1] = +2.0; betta[0][2][2] = -3.0;

//    betta[N][0][0] = +7.5; betta[N][0][1] = -1.5; betta[N][0][2] = +1.5;
//    betta[N][1][0] = -2.5; betta[N][1][1] = +2.5; betta[N][1][2] = -5.5;
//    betta[N][2][0] = +4.5; betta[N][2][1] = +5.5; betta[N][2][2] = -3.5;

    DoubleVector eta = DoubleVector(betta[0]*rx.col(0)) + DoubleVector(betta[N]*rx.col(N));
    IPrinter::printSeperatorLine();
    //IPrinter::print(eta, eta.size());

    DoubleMatrix* alpha1 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha1[i].resize(3,3,0.0);
    DoubleMatrix* alpha2 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha2[i].resize(3,3,0.0);
    DoubleMatrix* alpha0 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha0[i].resize(3,1,0.0);

    DoubleMatrix AI(3,3);
    for (unsigned int k=0; k<=N-2; k++)
    {
        AI[0][0] = 2.0*h*a(1,1,k) + 3.0; AI[0][1] = 2.0*h*a(1,2,k) + 0.0; AI[0][2] = 2.0*h*a(1,3,k) + 0.0;
        AI[1][0] = 2.0*h*a(2,1,k) + 0.0; AI[1][1] = 2.0*h*a(2,2,k) + 3.0; AI[1][2] = 2.0*h*a(2,3,k) + 0.0;
        AI[2][0] = 2.0*h*a(3,1,k) + 0.0; AI[2][1] = 2.0*h*a(3,2,k) + 0.0; AI[2][2] = 2.0*h*a(3,3,k) + 3.0;
        AI.inverse();

        alpha1[k][0][0] = +4.0; alpha1[k][0][1] = +0.0; alpha1[k][0][2] = +0.0;
        alpha1[k][1][0] = +0.0; alpha1[k][1][1] = +4.0; alpha1[k][1][2] = +0.0;
        alpha1[k][2][0] = +0.0; alpha1[k][2][1] = +0.0; alpha1[k][2][2] = +4.0;

        alpha2[k][0][0] = -1.0; alpha2[k][0][1] = +0.0; alpha2[k][0][2] = +0.0;
        alpha2[k][1][0] = +0.0; alpha2[k][1][1] = -1.0; alpha2[k][1][2] = +0.0;
        alpha2[k][2][0] = +0.0; alpha2[k][2][1] = +0.0; alpha2[k][2][2] = -1.0;

        alpha0[k][0][0] = -2.0*h*b(1,k);
        alpha0[k][1][0] = -2.0*h*b(2,k);
        alpha0[k][2][0] = -2.0*h*b(3,k);

        alpha1[k] = AI*alpha1[k];
        alpha2[k] = AI*alpha2[k];
        alpha0[k] = AI*alpha0[k];
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-2; k++)
    {
        betta[k+1] = betta[k+1] + betta[k]*alpha1[k];
        betta[k+2] = betta[k+2] + betta[k]*alpha2[k];
        eta        = eta        - betta[k]*alpha0[k];
    }

    DoubleMatrix M(9,9);
    DoubleVector A(9);
    DoubleVector x(9);

    M[0][0] = 0.0; M[0][1] = 0.0; M[0][2] = 0.0; M[0][3] = betta[N-1][0][0]; M[0][4] = betta[N-1][0][1]; M[0][5] = betta[N-1][0][2]; M[0][6] = betta[N-0][0][0]; M[0][7] = betta[N-0][0][1]; M[0][8] = betta[N-0][0][2]; A[0] = eta[0];
    M[1][0] = 0.0; M[1][1] = 0.0; M[1][2] = 0.0; M[1][3] = betta[N-1][1][0]; M[1][4] = betta[N-1][1][1]; M[1][5] = betta[N-1][1][2]; M[1][6] = betta[N-0][1][0]; M[1][7] = betta[N-0][1][1]; M[1][8] = betta[N-0][1][2]; A[1] = eta[1];
    M[2][0] = 0.0; M[2][1] = 0.0; M[2][2] = 0.0; M[2][3] = betta[N-1][2][0]; M[2][4] = betta[N-1][2][1]; M[2][5] = betta[N-1][2][2]; M[2][6] = betta[N-0][2][0]; M[2][7] = betta[N-0][2][1]; M[2][8] = betta[N-0][2][2]; A[2] = eta[2];

    M[3][0] = +1.0; M[3][1] = +0.0; M[3][2] = +0.0; M[3][3] = +2.0*h*a(1,1,N-1); M[3][4] = +2.0*h*a(1,2,N-1); M[3][5] = +2.0*h*a(1,3,N-1); M[3][6] = -1.0; M[3][7] = +0.0; M[3][8] = +0.0; A[3] = -2.0*h*b(1,N-1);
    M[4][0] = +0.0; M[4][1] = +1.0; M[4][2] = +0.0; M[4][3] = +2.0*h*a(2,1,N-1); M[4][4] = +2.0*h*a(2,2,N-1); M[4][5] = +2.0*h*a(2,3,N-1); M[4][6] = +0.0; M[4][7] = -1.0; M[4][8] = +0.0; A[4] = -2.0*h*b(2,N-1);
    M[5][0] = +0.0; M[5][1] = +0.0; M[5][2] = +1.0; M[5][3] = +2.0*h*a(3,1,N-1); M[5][4] = +2.0*h*a(3,2,N-1); M[5][5] = +2.0*h*a(3,3,N-1); M[5][6] = +0.0; M[5][7] = +0.0; M[5][8] = -1.0; A[5] = -2.0*h*b(3,N-1);

    M[6][0] = +1.0; M[6][1] = +0.0; M[6][2] = +0.0; M[6][3] = -4.0; M[6][4] = +0.0; M[6][5] = +0.0; M[6][6] = +3.0-2.0*h*a(1,1,N); M[6][7] = +0.0-2.0*h*a(1,2,N); M[6][8] = +0.0-2.0*h*a(3,1,N); A[6] = +2.0*h*b(1,N);
    M[7][0] = +0.0; M[7][1] = +1.0; M[7][2] = +0.0; M[7][3] = +0.0; M[7][4] = -4.0; M[7][5] = +0.0; M[7][6] = +0.0-2.0*h*a(2,1,N); M[7][7] = +3.0-2.0*h*a(2,2,N); M[7][8] = +0.0-2.0*h*a(3,2,N); A[7] = +2.0*h*b(2,N);
    M[8][0] = +0.0; M[8][1] = +0.0; M[8][2] = +1.0; M[8][3] = +0.0; M[8][4] = +0.0; M[8][5] = -4.0; M[8][6] = +0.0-2.0*h*a(3,1,N); M[8][7] = +0.0-2.0*h*a(3,2,N); M[8][8] = +3.0-2.0*h*a(3,3,N); A[8] = +2.0*h*b(3,N);

    //printf("det: %.10f\n", M.determinant());
    //printf("%18.10f %18.10f\n", M[0][0]*rx[0][N-2]+M[0][1]*rx[1][N-2]+M[0][2]*rx[2][N-2]+M[0][3]*rx[0][N-1]+M[0][4]*rx[1][N-1]+M[0][5]*rx[2][N-1]+M[0][6]*rx[0][N-0]+M[0][7]*rx[1][N-0]+M[0][8]*rx[2][N-0], A[0]);
    //printf("%18.10f %18.10f\n", M[1][0]*rx[0][N-2]+M[1][1]*rx[1][N-2]+M[1][2]*rx[2][N-2]+M[1][3]*rx[0][N-1]+M[1][4]*rx[1][N-1]+M[1][5]*rx[2][N-1]+M[1][6]*rx[0][N-0]+M[1][7]*rx[1][N-0]+M[1][8]*rx[2][N-0], A[1]);
    //printf("%18.10f %18.10f\n", M[2][0]*rx[0][N-2]+M[2][1]*rx[1][N-2]+M[2][2]*rx[2][N-2]+M[2][3]*rx[0][N-1]+M[2][4]*rx[1][N-1]+M[2][5]*rx[2][N-1]+M[2][6]*rx[0][N-0]+M[2][7]*rx[1][N-0]+M[2][8]*rx[2][N-0], A[2]);
    //printf("%18.10f %18.10f\n", M[3][0]*rx[0][N-2]+M[3][1]*rx[1][N-2]+M[3][2]*rx[2][N-2]+M[3][3]*rx[0][N-1]+M[3][4]*rx[1][N-1]+M[3][5]*rx[2][N-1]+M[3][6]*rx[0][N-0]+M[3][7]*rx[1][N-0]+M[3][8]*rx[2][N-0], A[3]);
    //printf("%18.10f %18.10f\n", M[4][0]*rx[0][N-2]+M[4][1]*rx[1][N-2]+M[4][2]*rx[2][N-2]+M[4][3]*rx[0][N-1]+M[4][4]*rx[1][N-1]+M[4][5]*rx[2][N-1]+M[4][6]*rx[0][N-0]+M[4][7]*rx[1][N-0]+M[4][8]*rx[2][N-0], A[4]);
    //printf("%18.10f %18.10f\n", M[5][0]*rx[0][N-2]+M[5][1]*rx[1][N-2]+M[5][2]*rx[2][N-2]+M[5][3]*rx[0][N-1]+M[5][4]*rx[1][N-1]+M[5][5]*rx[2][N-1]+M[5][6]*rx[0][N-0]+M[5][7]*rx[1][N-0]+M[5][8]*rx[2][N-0], A[5]);
    //printf("%18.10f %18.10f\n", M[6][0]*rx[0][N-2]+M[6][1]*rx[1][N-2]+M[6][2]*rx[2][N-2]+M[6][3]*rx[0][N-1]+M[6][4]*rx[1][N-1]+M[6][5]*rx[2][N-1]+M[6][6]*rx[0][N-0]+M[6][7]*rx[1][N-0]+M[6][8]*rx[2][N-0], A[6]);
    //printf("%18.10f %18.10f\n", M[7][0]*rx[0][N-2]+M[7][1]*rx[1][N-2]+M[7][2]*rx[2][N-2]+M[7][3]*rx[0][N-1]+M[7][4]*rx[1][N-1]+M[7][5]*rx[2][N-1]+M[7][6]*rx[0][N-0]+M[7][7]*rx[1][N-0]+M[7][8]*rx[2][N-0], A[7]);
    //printf("%18.10f %18.10f\n", M[8][0]*rx[0][N-2]+M[8][1]*rx[1][N-2]+M[8][2]*rx[2][N-2]+M[8][3]*rx[0][N-1]+M[8][4]*rx[1][N-1]+M[8][5]*rx[2][N-1]+M[8][6]*rx[0][N-0]+M[8][7]*rx[1][N-0]+M[8][8]*rx[2][N-0], A[8]);

    //IPrinter::printSeperatorLine();
    //IPrinter::print(M,9,9);
    //IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

    //printf("%.10f %.10f %.10f\n", x[0], x[1], x[2]);
    //printf("%.10f %.10f %.10f\n", x[3], x[4], x[5]);
    //printf("%.10f %.10f %.10f\n", x[6], x[7], x[8]);

    //IPrinter::printSeperatorLine();
    //printf("%18.10f %18.10f\n", M[0][0]*x[0]+M[0][1]*x[1]+M[0][2]*x[2]+M[0][3]*x[3]+M[0][4]*x[4]+M[0][5]*x[5]+M[0][6]*x[6]+M[0][7]*x[7]+M[0][8]*x[8], A[0]);
    //printf("%18.10f %18.10f\n", M[1][0]*x[0]+M[1][1]*x[1]+M[1][2]*x[2]+M[1][3]*x[3]+M[1][4]*x[4]+M[1][5]*x[5]+M[1][6]*x[6]+M[1][7]*x[7]+M[1][8]*x[8], A[1]);
    //printf("%18.10f %18.10f\n", M[2][0]*x[0]+M[2][1]*x[1]+M[2][2]*x[2]+M[2][3]*x[3]+M[2][4]*x[4]+M[2][5]*x[5]+M[2][6]*x[6]+M[2][7]*x[7]+M[2][8]*x[8], A[2]);
    //printf("%18.10f %18.10f\n", M[3][0]*x[0]+M[3][1]*x[1]+M[3][2]*x[2]+M[3][3]*x[3]+M[3][4]*x[4]+M[3][5]*x[5]+M[3][6]*x[6]+M[3][7]*x[7]+M[3][8]*x[8], A[3]);
    //printf("%18.10f %18.10f\n", M[4][0]*x[0]+M[4][1]*x[1]+M[4][2]*x[2]+M[4][3]*x[3]+M[4][4]*x[4]+M[4][5]*x[5]+M[4][6]*x[6]+M[4][7]*x[7]+M[4][8]*x[8], A[4]);
    //printf("%18.10f %18.10f\n", M[5][0]*x[0]+M[5][1]*x[1]+M[5][2]*x[2]+M[5][3]*x[3]+M[5][4]*x[4]+M[5][5]*x[5]+M[5][6]*x[6]+M[5][7]*x[7]+M[5][8]*x[8], A[5]);
    //printf("%18.10f %18.10f\n", M[6][0]*x[0]+M[6][1]*x[1]+M[6][2]*x[2]+M[6][3]*x[3]+M[6][4]*x[4]+M[6][5]*x[5]+M[6][6]*x[6]+M[6][7]*x[7]+M[6][8]*x[8], A[6]);
    //printf("%18.10f %18.10f\n", M[7][0]*x[0]+M[7][1]*x[1]+M[7][2]*x[2]+M[7][3]*x[3]+M[7][4]*x[4]+M[7][5]*x[5]+M[7][6]*x[6]+M[7][7]*x[7]+M[7][8]*x[8], A[7]);
    //printf("%18.10f %18.10f\n", M[8][0]*x[0]+M[8][1]*x[1]+M[8][2]*x[2]+M[8][3]*x[3]+M[8][4]*x[4]+M[8][5]*x[5]+M[8][6]*x[6]+M[8][7]*x[7]+M[8][8]*x[8], A[8]);

    DoubleMatrix nx(3, N+1, 0.0);
    nx[2][N-0] = x[8]; nx[2][N-1] = x[5]; nx[2][N-2] = x[2];
    nx[1][N-0] = x[7]; nx[1][N-1] = x[4]; nx[1][N-2] = x[1];
    nx[0][N-0] = x[6]; nx[0][N-1] = x[3]; nx[0][N-2] = x[0];

//    //IPrinter::printSeperatorLine();
//    //printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
//    //printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
//    //printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);

    betta[N-1] = betta[N-1] - betta[N-2]*alpha1[N-2];
    betta[N-0] = betta[N-0] - betta[N-2]*alpha2[N-2];
    eta        = eta + DoubleVector(betta[N-2]*alpha0[N-2]);

    for (unsigned int k=N-2; k>0; k--)
    {
        betta[k+0] = betta[k+0] - betta[k-1]*alpha1[k-1];
        betta[k+1] = betta[k+1] - betta[k-1]*alpha2[k-1];
        eta        = eta + DoubleVector(betta[k-1]*alpha0[k-1]);

        nx[0][k-1] = eta[0];
        nx[1][k-1] = eta[1];
        nx[2][k-1] = eta[2];
        for (unsigned int i=k; i<=N; i++)
        {
            DoubleVector col = nx.col(k-1);
            DoubleVector aaa = DoubleVector(betta[i]*nx.col(i));
            col[0] -= aaa[0];
            col[1] -= aaa[1];
            col[2] -= aaa[2];

            //DoubleVector col = nx.col(k-1)-DoubleVector(betta[i]*nx.col(i));
            nx.setColumn(k-1, col);
        }
        DoubleMatrix bbb = betta[k-1];
        bbb.inverse();
        nx.setColumn(k-1, DoubleVector(bbb*nx.col(k-1)));
    }

    IPrinter::printVector(w,p,nx.row(0));
    IPrinter::printVector(w,p,nx.row(1));
    IPrinter::printVector(w,p,nx.row(2));

//    FILE *file =fopen("data_rx.txt", "a");
//    IPrinter::printVector(14,10,nx,"nv2",nx.size(),0,0,file);
//    fclose(file);
}

void SystemDifEquation::calculate4R2LV1(const DoubleMatrix &rx)
{
    DoubleMatrix* betta = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) betta[i].resize(3,3,0.0);

    betta[0].randomData();
    betta[N].randomData();

//    betta[0][0][0] = +1.0; betta[0][0][1] = +8.0; betta[0][0][2] = -1.0;
//    betta[0][1][0] = +5.0; betta[0][1][1] = +3.0; betta[0][1][2] = +7.0;
//    betta[0][2][0] = -2.0; betta[0][2][1] = +2.0; betta[0][2][2] = -3.0;

//    betta[N][0][0] = +7.5; betta[N][0][1] = -1.5; betta[N][0][2] = +1.5;
//    betta[N][1][0] = -2.5; betta[N][1][1] = +2.5; betta[N][1][2] = -5.5;
//    betta[N][2][0] = +4.5; betta[N][2][1] = +5.5; betta[N][2][2] = -3.5;

    DoubleVector eta = DoubleVector(betta[0]*rx.col(0)) + DoubleVector(betta[N]*rx.col(N));
    IPrinter::printSeperatorLine();
    //IPrinter::print(eta, eta.size());

    DoubleMatrix* alpha1 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha1[i].resize(3,3,0.0);
    DoubleMatrix* alpha2 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha2[i].resize(3,3,0.0);
    DoubleMatrix* alpha3 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha3[i].resize(3,3,0.0);
    DoubleMatrix* alpha4 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha4[i].resize(3,3,0.0);
    DoubleMatrix* alpha0 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha0[i].resize(3,1,0.0);

    DoubleMatrix AI(3,3);
    for (unsigned int k=0; k<=N-4; k++)
    {
        AI[0][0] = -12.0*h*a(1,1,k) - 25.0; AI[0][1] = -12.0*h*a(1,2,k) + 0.0;  AI[0][2] = -12.0*h*a(1,3,k) + 0.0;
        AI[1][0] = -12.0*h*a(2,1,k) + 0.0;  AI[1][1] = -12.0*h*a(2,2,k) - 25.0; AI[1][2] = -12.0*h*a(2,3,k) + 0.0;
        AI[2][0] = -12.0*h*a(3,1,k) + 0.0;  AI[2][1] = -12.0*h*a(3,2,k) + 0.0;  AI[2][2] = -12.0*h*a(3,3,k) - 25.0;
        AI.inverse();

        alpha1[k][0][0] = -48.0; alpha1[k][0][1] =  +0.0; alpha1[k][0][2] =  +0.0;
        alpha1[k][1][0] =  +0.0; alpha1[k][1][1] = -48.0; alpha1[k][1][2] =  +0.0;
        alpha1[k][2][0] =  +0.0; alpha1[k][2][1] =  +0.0; alpha1[k][2][2] = -48.0;

        alpha2[k][0][0] = +36.0; alpha2[k][0][1] = +0.0;  alpha2[k][0][2] = +0.0;
        alpha2[k][1][0] =  +0.0; alpha2[k][1][1] = +36.0; alpha2[k][1][2] = +0.0;
        alpha2[k][2][0] =  +0.0; alpha2[k][2][1] = +0.0;  alpha2[k][2][2] = +36.0;

        alpha3[k][0][0] = -16.0; alpha3[k][0][1] = +0.0;  alpha3[k][0][2] = +0.0;
        alpha3[k][1][0] =  +0.0; alpha3[k][1][1] = -16.0; alpha3[k][1][2] = +0.0;
        alpha3[k][2][0] =  +0.0; alpha3[k][2][1] = +0.0;  alpha3[k][2][2] = -16.0;

        alpha4[k][0][0] =  +3.0; alpha4[k][0][1] = +0.0; alpha4[k][0][2] = +0.0;
        alpha4[k][1][0] =  +0.0; alpha4[k][1][1] = +3.0; alpha4[k][1][2] = +0.0;
        alpha4[k][2][0] =  +0.0; alpha4[k][2][1] = +0.0; alpha4[k][2][2] = +3.0;

        alpha0[k][0][0] = +12.0*h*b(1,k);
        alpha0[k][1][0] = +12.0*h*b(2,k);
        alpha0[k][2][0] = +12.0*h*b(3,k);

        alpha1[k] = AI*alpha1[k];
        alpha2[k] = AI*alpha2[k];
        alpha3[k] = AI*alpha3[k];
        alpha4[k] = AI*alpha4[k];
        alpha0[k] = AI*alpha0[k];
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-4; k++)
    {
        betta[k+1] = betta[k+1] + betta[k]*alpha1[k];
        betta[k+2] = betta[k+2] + betta[k]*alpha2[k];
        betta[k+3] = betta[k+3] + betta[k]*alpha3[k];
        betta[k+4] = betta[k+4] + betta[k]*alpha4[k];
        eta        = eta        - betta[k]*alpha0[k];
    }

    DoubleMatrix M(15,15);
    DoubleVector A(15);
    DoubleVector x(15);

    M[0][0]  = 0.0;              M[0][1] = 0.0;               M[0][2] = 0.0;
    M[0][3]  = betta[N-3][0][0]; M[0][4]  = betta[N-3][0][1]; M[0][5]  = betta[N-3][0][2];
    M[0][6]  = betta[N-2][0][0]; M[0][7]  = betta[N-2][0][1]; M[0][8]  = betta[N-2][0][2];
    M[0][9]  = betta[N-1][0][0]; M[0][10] = betta[N-1][0][1]; M[0][11] = betta[N-1][0][2];
    M[0][12] = betta[N-0][0][0]; M[0][13] = betta[N-0][0][1]; M[0][14] = betta[N-0][0][2]; A[0] = eta[0];

    M[1][0]  = 0.0;              M[1][1]  = 0.0;              M[1][2]  = 0.0;
    M[1][3]  = betta[N-3][1][0]; M[1][4]  = betta[N-3][1][1]; M[1][5]  = betta[N-3][1][2];
    M[1][6]  = betta[N-2][1][0]; M[1][7]  = betta[N-2][1][1]; M[1][8]  = betta[N-2][1][2];
    M[1][9]  = betta[N-1][1][0]; M[1][10] = betta[N-1][1][1]; M[1][11] = betta[N-1][1][2];
    M[1][12] = betta[N-0][1][0]; M[1][13] = betta[N-0][1][1]; M[1][14] = betta[N-0][1][2]; A[1] = eta[1];

    M[2][0]  = 0.0;              M[2][1]  = 0.0;              M[2][2]  = 0.0;
    M[2][3]  = betta[N-3][2][0]; M[2][4]  = betta[N-3][2][1]; M[2][5]  = betta[N-3][2][2];
    M[2][6]  = betta[N-2][2][0]; M[2][7]  = betta[N-2][2][1]; M[2][8]  = betta[N-2][2][2];
    M[2][9]  = betta[N-1][2][0]; M[2][10] = betta[N-1][2][1]; M[2][11] = betta[N-1][2][2];
    M[2][12] = betta[N-0][2][0]; M[2][13] = betta[N-0][2][1]; M[2][14] = betta[N-0][2][2]; A[2] = eta[2];

    M[3][0] = -3.0; M[3][1] = +0.0; M[3][2] = +0.0; M[3][3] = -12.0*h*a(1,1,N-3)-10.0; M[3][4] = -12.0*h*a(1,2,N-3);      M[3][5] = -12.0*h*a(1,3,N-3);      M[3][6] = +18.0; M[3][7] =  +0.0; M[3][8] =  +0.0; M[3][9] = -6.0; M[3][10] = +0.0; M[3][11] = +0.0; M[3][12] = +1.0; M[3][13] = +0.0; M[3][14] = +0.0; A[3] = +12.0*h*b(1,N-3);
    M[4][0] = +0.0; M[4][1] = -3.0; M[4][2] = +0.0; M[4][3] = -12.0*h*a(2,1,N-3);      M[4][4] = -12.0*h*a(2,2,N-3)-10.0; M[4][5] = -12.0*h*a(2,3,N-3);      M[4][6] =  +0.0; M[4][7] = +18.0; M[4][8] =  +0.0; M[4][9] = +0.0; M[4][10] = -6.0; M[4][11] = +0.0; M[4][12] = +0.0; M[4][13] = +1.0; M[4][14] = +0.0; A[4] = +12.0*h*b(2,N-3);
    M[5][0] = +0.0; M[5][1] = +0.0; M[5][2] = -3.0; M[5][3] = -12.0*h*a(3,1,N-3);      M[5][4] = -12.0*h*a(3,2,N-3);      M[5][5] = -12.0*h*a(3,3,N-3)-10.0; M[5][6] =  +0.0; M[5][7] =  +0.0; M[5][8] = +18.0; M[5][9] = +0.0; M[5][10] = +0.0; M[5][11] = -6.0; M[5][12] = +0.0; M[5][13] = +0.0; M[5][14] = +1.0; A[5] = +12.0*h*b(3,N-3);

    M[6][0] = +1.0; M[6][1] = +0.0; M[6][2] = +0.0; M[6][3] = -8.0; M[6][4] = +0.0; M[6][5] = +0.0; M[6][6] = -12.0*h*a(1,1,N-2); M[6][7] = -12.0*h*a(1,2,N-2); M[6][8] = -12.0*h*a(1,3,N-2); M[6][9] = +8.0; M[6][10] = +0.0; M[6][11] = +0.0; M[6][12] = -1.0; M[6][13] = +0.0; M[6][14] = +0.0; A[6] = +12.0*h*b(1,N-2);
    M[7][0] = +0.0; M[7][1] = +1.0; M[7][2] = +0.0; M[7][3] = +0.0; M[7][4] = -8.0; M[7][5] = +0.0; M[7][6] = -12.0*h*a(2,1,N-2); M[7][7] = -12.0*h*a(2,2,N-2); M[7][8] = -12.0*h*a(2,3,N-2); M[7][9] = +0.0; M[7][10] = +8.0; M[7][11] = +0.0; M[7][12] = +0.0; M[7][13] = -1.0; M[7][14] = +0.0; A[7] = +12.0*h*b(2,N-2);
    M[8][0] = +0.0; M[8][1] = +0.0; M[8][2] = +1.0; M[8][3] = +0.0; M[8][4] = +0.0; M[8][5] = -8.0; M[8][6] = -12.0*h*a(3,1,N-2); M[8][7] = -12.0*h*a(3,2,N-2); M[8][8] = -12.0*h*a(3,3,N-2); M[8][9] = +0.0; M[8][10] = +0.0; M[8][11] = +8.0; M[8][12] = +0.0; M[8][13] = +0.0; M[8][14] = -1.0; A[8] = +12.0*h*b(3,N-2);

    M[9][0]  = -1.0; M[9][1]  = +0.0; M[9][2]  = +0.0; M[9][3]  = +6.0; M[9][4]  = +0.0; M[9][5]  = +0.0; M[9][6]  = -18.0; M[9][7]  = +0.0;  M[9][8]  = +0.0;  M[9][9]  = -12.0*h*a(1,1,N-1)+10.0;  M[9][10]  = -12.0*h*a(1,2,N-1);      M[9][11]  = -12.0*h*a(1,3,N-1);      M[9][12]  = +3.0; M[9][13]  = +0.0; M[9][14]  = +0.0; A[9]  = +12.0*h*b(1,N-1);
    M[10][0] = +0.0; M[10][1] = -1.0; M[10][2] = +0.0; M[10][3] = +0.0; M[10][4] = +6.0; M[10][5] = +0.0; M[10][6] = +0.0;  M[10][7] = -18.0; M[10][8] = +0.0;  M[10][9] = -12.0*h*a(2,1,N-1);       M[10][10] = -12.0*h*a(2,2,N-1)+10.0; M[10][11] = -12.0*h*a(2,3,N-1);      M[10][12] = +0.0; M[10][13] = +3.0; M[10][14] = +0.0; A[10] = +12.0*h*b(2,N-1);
    M[11][0] = +0.0; M[11][1] = +0.0; M[11][2] = -1.0; M[11][3] = +0.0; M[11][4] = +0.0; M[11][5] = +6.0; M[11][6] = +0.0;  M[11][7] = +0.0;  M[11][8] = -18.0; M[11][9] = -12.0*h*a(3,1,N-1);       M[11][10] = -12.0*h*a(3,2,N-1);      M[11][11] = -12.0*h*a(3,3,N-1)+10.0; M[11][12] = +0.0; M[11][13] = +0.0; M[11][14] = +3.0; A[11] = +12.0*h*b(3,N-1);

    M[12][0] = +3.0; M[12][1] = +0.0; M[12][2] = +0.0; M[12][3] = -16.0; M[12][4] = +0.0;  M[12][5] = +0.0;  M[12][6] = +36.0; M[12][7] = +0.0;  M[12][8] = +0.0;  M[12][9] = -48.0; M[12][10] = +0.0;  M[12][11] = +0.0;  M[12][12] = -12.0*h*a(1,1,N-0)+25.0; M[12][13] = -12.0*h*a(1,2,N-0);      M[12][14] = -12.0*h*a(1,3,N-0);      A[12] = +12.0*h*b(1,N-0);
    M[13][0] = +0.0; M[13][1] = +3.0; M[13][2] = +0.0; M[13][3] = +0.0;  M[13][4] = -16.0; M[13][5] = +0.0;  M[13][6] = +0.0;  M[13][7] = +36.0; M[13][8] = +0.0;  M[13][9] = +0.0;  M[13][10] = -48.0; M[13][11] = +0.0;  M[13][12] = -12.0*h*a(2,1,N-0);      M[13][13] = -12.0*h*a(2,2,N-0)+25.0; M[13][14] = -12.0*h*a(2,3,N-0);      A[13] = +12.0*h*b(2,N-0);
    M[14][0] = +0.0; M[14][1] = +0.0; M[14][2] = +3.0; M[14][3] = +0.0;  M[14][4] = +0.0;  M[14][5] = -16.0; M[14][6] = +0.0;  M[14][7] = +0.0;  M[14][8] = +36.0; M[14][9] = +0.0;  M[14][10] = +0.0;  M[14][11] = -48.0; M[14][12] = -12.0*h*a(3,1,N-0);      M[14][13] = -12.0*h*a(3,2,N-0);      M[14][14] = -12.0*h*a(3,3,N-0)+25.0; A[14] = +12.0*h*b(3,N-0);

//    printf("det: %.10f\n", M.determinant());

//    printf("%18.10f %18.10f\n", M[0][0]*rx[0][N-4]+ M[0][1]*rx[1][N-4]+ M[0][2]*rx[2][N-4]+
//                                M[0][3]*rx[0][N-3]+ M[0][4]*rx[1][N-3]+ M[0][5]*rx[2][N-3]+
//                                M[0][6]*rx[0][N-2]+ M[0][7]*rx[1][N-2]+ M[0][8]*rx[2][N-2]+
//                                M[0][9]*rx[0][N-1]+ M[0][10]*rx[2][N-1]+M[0][11]*rx[2][N-1]+
//                                M[0][12]*rx[0][N-0]+M[0][13]*rx[2][N-0]+M[0][14]*rx[2][N-0], A[0]);

//    printf("%18.10f %18.10f\n", M[1][0]*rx[0][N-4]+M[1][1]*rx[1][N-4]+M[1][2]*rx[2][N-4]+
//                                M[1][3]*rx[0][N-3]+M[1][4]*rx[1][N-3]+M[1][5]*rx[2][N-3]+
//                                M[1][6]*rx[0][N-2]+M[1][7]*rx[1][N-2]+M[1][8]*rx[2][N-2]+
//                                M[1][9]*rx[0][N-1]+M[1][10]*rx[1][N-1]+M[1][11]*rx[2][N-1]+
//                                M[1][12]*rx[0][N-0]+M[1][13]*rx[1][N-0]+M[1][14]*rx[2][N-0], A[1]);

//    printf("%18.10f %18.10f\n", M[2][0]*rx[0][N-4]+M[2][1]*rx[1][N-4]+M[2][2]*rx[2][N-4]+
//                                M[2][3]*rx[0][N-3]+M[2][4]*rx[1][N-3]+M[2][5]*rx[2][N-3]+
//                                M[2][6]*rx[0][N-2]+M[2][7]*rx[1][N-2]+M[2][8]*rx[2][N-2]+
//                                M[2][9]*rx[0][N-1]+M[2][10]*rx[1][N-1]+M[2][11]*rx[2][N-1]+
//                                M[2][12]*rx[0][N-0]+M[2][13]*rx[1][N-0]+M[2][14]*rx[2][N-0], A[2]);

//    printf("%18.10f %18.10f\n", M[3][0]*rx[0][N-4] +M[3][1]*rx[1][N-4]+ M[3][2]*rx[2][N-4]+
//                                M[3][3]*rx[0][N-3] +M[3][4]*rx[1][N-3]+ M[3][5]*rx[2][N-3]+
//                                M[3][6]*rx[0][N-2] +M[3][7]*rx[1][N-2]+ M[3][8]*rx[2][N-2]+
//                                M[3][9]*rx[0][N-1] +M[3][10]*rx[1][N-1]+M[3][11]*rx[2][N-1]+
//                                M[3][12]*rx[0][N-0]+M[3][13]*rx[1][N-0]+M[3][14]*rx[2][N-0], A[3]);

//    printf("%18.10f %18.10f\n", M[4][0]*rx[0][N-4]+ M[4][1]*rx[1][N-4]+ M[4][2]*rx[2][N-4]+
//                                M[4][3]*rx[0][N-3]+ M[4][4]*rx[1][N-3]+ M[4][5]*rx[2][N-3]+
//                                M[4][6]*rx[0][N-2]+ M[4][7]*rx[1][N-2]+ M[4][8]*rx[2][N-2]+
//                                M[4][9]*rx[0][N-1]+ M[4][10]*rx[1][N-1]+M[4][11]*rx[2][N-1]+
//                                M[4][12]*rx[0][N-0]+M[4][13]*rx[1][N-0]+M[4][14]*rx[2][N-0], A[4]);

//    printf("%18.10f %18.10f\n", M[5][0]*rx[0][N-4]+ M[5][1]*rx[1][N-4]+ M[5][2]*rx[2][N-4]+
//                                M[5][3]*rx[0][N-3]+ M[5][4]*rx[1][N-3]+ M[5][5]*rx[2][N-3]+
//                                M[5][6]*rx[0][N-2]+ M[5][7]*rx[1][N-2]+ M[5][8]*rx[2][N-2]+
//                                M[5][9]*rx[0][N-1]+ M[5][10]*rx[1][N-1]+M[5][11]*rx[2][N-1]+
//                                M[5][12]*rx[0][N-0]+M[5][13]*rx[1][N-0]+M[5][14]*rx[2][N-0], A[5]);

//    printf("%18.10f %18.10f\n", M[6][0]*rx[0][N-4]+ M[6][1]*rx[1][N-4]+ M[6][2]*rx[2][N-4]+
//                                M[6][3]*rx[0][N-3]+ M[6][4]*rx[1][N-3]+ M[6][5]*rx[2][N-3]+
//                                M[6][6]*rx[0][N-2]+ M[6][7]*rx[1][N-2]+ M[6][8]*rx[2][N-2]+
//                                M[6][9]*rx[0][N-1]+ M[6][10]*rx[1][N-1]+M[6][11]*rx[2][N-1]+
//                                M[6][12]*rx[0][N-0]+M[6][13]*rx[1][N-0]+M[6][14]*rx[2][N-0], A[6]);

//    printf("%18.10f %18.10f\n", M[7][0]*rx[0][N-4]+M[7][1]*rx[1][N-4]+M[7][2]*rx[2][N-4]+
//                                M[7][3]*rx[0][N-3]+M[7][4]*rx[1][N-3]+M[7][5]*rx[2][N-3]+
//                                M[7][6]*rx[0][N-2]+M[7][7]*rx[1][N-2]+M[7][8]*rx[2][N-2]+
//                                M[7][9]*rx[0][N-1]+M[7][10]*rx[1][N-1]+M[7][11]*rx[2][N-1]+
//                                M[7][12]*rx[0][N-0]+M[7][13]*rx[1][N-0]+M[7][14]*rx[2][N-0], A[7]);

//    printf("%18.10f %18.10f\n", M[8][0]*rx[0][N-4]+M[8][1]*rx[1][N-4]+M[8][2]*rx[2][N-4]+
//                                M[8][3]*rx[0][N-3]+M[8][4]*rx[1][N-3]+M[8][5]*rx[2][N-3]+
//                                M[8][6]*rx[0][N-2]+M[8][7]*rx[1][N-2]+M[8][8]*rx[2][N-2]+
//                                M[8][9]*rx[0][N-1]+M[8][10]*rx[1][N-1]+M[8][11]*rx[2][N-1]+
//                                M[8][12]*rx[0][N-0]+M[8][13]*rx[1][N-0]+M[8][14]*rx[2][N-0], A[8]);

//    printf("%18.10f %18.10f\n", M[9][0]*rx[0][N-4]+ M[9][1]*rx[1][N-4]+ M[9][2]*rx[2][N-4]+
//                                M[9][3]*rx[0][N-3]+ M[9][4]*rx[1][N-3]+ M[9][5]*rx[2][N-3]+
//                                M[9][6]*rx[0][N-2]+ M[9][7]*rx[1][N-2]+ M[9][8]*rx[2][N-2]+
//                                M[9][9]*rx[0][N-1]+ M[9][10]*rx[1][N-1]+M[9][11]*rx[2][N-1]+
//                                M[9][12]*rx[0][N-0]+M[9][13]*rx[1][N-0]+M[9][14]*rx[2][N-0], A[9]);

//    printf("%18.10f %18.10f\n", M[10][0]*rx[0][N-4]+M[10][1]*rx[1][N-4]+M[10][2]*rx[2][N-4]+
//                                M[10][3]*rx[0][N-3]+M[10][4]*rx[1][N-3]+M[10][5]*rx[2][N-3]+
//                                M[10][6]*rx[0][N-2]+M[10][7]*rx[1][N-2]+M[10][8]*rx[2][N-2]+
//                                M[10][9]*rx[0][N-1]+M[10][10]*rx[1][N-1]+M[10][11]*rx[2][N-1]+
//                                M[10][12]*rx[0][N-0]+M[10][13]*rx[1][N-0]+M[10][14]*rx[2][N-0], A[10]);

//    printf("%18.10f %18.10f\n", M[11][0]*rx[0][N-4]+M[11][1]*rx[1][N-4]+M[11][2]*rx[2][N-4]+
//                                M[11][3]*rx[0][N-3]+M[11][4]*rx[1][N-3]+M[11][5]*rx[2][N-3]+
//                                M[11][6]*rx[0][N-2]+M[11][7]*rx[1][N-2]+M[11][8]*rx[2][N-2]+
//                                M[11][9]*rx[0][N-1]+M[11][10]*rx[1][N-1]+M[11][11]*rx[2][N-1]+
//                                M[11][12]*rx[0][N-0]+M[11][13]*rx[1][N-0]+M[11][14]*rx[2][N-0], A[11]);

//    printf("%18.10f %18.10f\n", M[12][0]*rx[0][N-4]+M[12][1]*rx[1][N-4]+M[12][2]*rx[2][N-4]+
//                                M[12][3]*rx[0][N-3]+M[12][4]*rx[1][N-3]+M[12][5]*rx[2][N-3]+
//                                M[12][6]*rx[0][N-2]+M[12][7]*rx[1][N-2]+M[12][8]*rx[2][N-2]+
//                                M[12][9]*rx[0][N-1]+M[12][10]*rx[1][N-1]+M[12][11]*rx[2][N-1]+
//                                M[12][12]*rx[0][N-0]+M[12][13]*rx[1][N-0]+M[12][14]*rx[2][N-0], A[12]);

//    printf("%18.10f %18.10f\n", M[13][0]*rx[0][N-4]+M[13][1]*rx[1][N-4]+M[13][2]*rx[2][N-4]+
//                                M[13][3]*rx[0][N-3]+M[13][4]*rx[1][N-3]+M[13][5]*rx[2][N-3]+
//                                M[13][6]*rx[0][N-2]+M[13][7]*rx[1][N-2]+M[13][8]*rx[2][N-2]+
//                                M[13][9]*rx[0][N-1]+M[13][10]*rx[1][N-1]+M[13][11]*rx[2][N-1]+
//                                M[13][12]*rx[0][N-0]+M[13][13]*rx[1][N-0]+M[13][14]*rx[2][N-0], A[13]);

//    printf("%18.10f %18.10f\n", M[14][0]*rx[0][N-4]+M[14][1]*rx[1][N-4]+M[14][2]*rx[2][N-4]+
//                                M[14][3]*rx[0][N-3]+M[14][4]*rx[1][N-3]+M[14][5]*rx[2][N-3]+
//                                M[14][6]*rx[0][N-2]+M[14][7]*rx[1][N-2]+M[14][8]*rx[2][N-2]+
//                                M[14][9]*rx[0][N-1]+M[14][10]*rx[1][N-1]+M[14][11]*rx[2][N-1]+
//                                M[14][12]*rx[0][N-0]+M[14][13]*rx[1][N-0]+M[14][14]*rx[2][N-0], A[14]);

//    IPrinter::printSeperatorLine();
//    IPrinter::print(M,15,15);
//    IPrinter::printSeperatorLine();

    GaussianElimination(M, A, x);

//    printf("%.10f %.10f %.10f\n", x[0],  x[1],  x[2]);
//    printf("%.10f %.10f %.10f\n", x[3],  x[4],  x[5]);
//    printf("%.10f %.10f %.10f\n", x[6],  x[7],  x[8]);
//    printf("%.10f %.10f %.10f\n", x[9],  x[10], x[11]);
//    printf("%.10f %.10f %.10f\n", x[12], x[13], x[14]);

//    IPrinter::printSeperatorLine();
//    printf("%18.10f %18.10f\n", M[0][0]*x[0]+M[0][1]*x[1]+M[0][2]*x[2]+M[0][3]*x[3]+M[0][4]*x[4]+M[0][5]*x[5]+M[0][6]*x[6]+M[0][7]*x[7]+M[0][8]*x[8], A[0]);
//    printf("%18.10f %18.10f\n", M[1][0]*x[0]+M[1][1]*x[1]+M[1][2]*x[2]+M[1][3]*x[3]+M[1][4]*x[4]+M[1][5]*x[5]+M[1][6]*x[6]+M[1][7]*x[7]+M[1][8]*x[8], A[1]);
//    printf("%18.10f %18.10f\n", M[2][0]*x[0]+M[2][1]*x[1]+M[2][2]*x[2]+M[2][3]*x[3]+M[2][4]*x[4]+M[2][5]*x[5]+M[2][6]*x[6]+M[2][7]*x[7]+M[2][8]*x[8], A[2]);
//    printf("%18.10f %18.10f\n", M[3][0]*x[0]+M[3][1]*x[1]+M[3][2]*x[2]+M[3][3]*x[3]+M[3][4]*x[4]+M[3][5]*x[5]+M[3][6]*x[6]+M[3][7]*x[7]+M[3][8]*x[8], A[3]);
//    printf("%18.10f %18.10f\n", M[4][0]*x[0]+M[4][1]*x[1]+M[4][2]*x[2]+M[4][3]*x[3]+M[4][4]*x[4]+M[4][5]*x[5]+M[4][6]*x[6]+M[4][7]*x[7]+M[4][8]*x[8], A[4]);
//    printf("%18.10f %18.10f\n", M[5][0]*x[0]+M[5][1]*x[1]+M[5][2]*x[2]+M[5][3]*x[3]+M[5][4]*x[4]+M[5][5]*x[5]+M[5][6]*x[6]+M[5][7]*x[7]+M[5][8]*x[8], A[5]);
//    printf("%18.10f %18.10f\n", M[6][0]*x[0]+M[6][1]*x[1]+M[6][2]*x[2]+M[6][3]*x[3]+M[6][4]*x[4]+M[6][5]*x[5]+M[6][6]*x[6]+M[6][7]*x[7]+M[6][8]*x[8], A[6]);
//    printf("%18.10f %18.10f\n", M[7][0]*x[0]+M[7][1]*x[1]+M[7][2]*x[2]+M[7][3]*x[3]+M[7][4]*x[4]+M[7][5]*x[5]+M[7][6]*x[6]+M[7][7]*x[7]+M[7][8]*x[8], A[7]);
//    printf("%18.10f %18.10f\n", M[8][0]*x[0]+M[8][1]*x[1]+M[8][2]*x[2]+M[8][3]*x[3]+M[8][4]*x[4]+M[8][5]*x[5]+M[8][6]*x[6]+M[8][7]*x[7]+M[8][8]*x[8], A[8]);
//    printf("%18.10f %18.10f\n", M[9][0]*x[0]+M[9][1]*x[1]+M[9][2]*x[2]+M[9][3]*x[3]+M[9][4]*x[4]+M[9][5]*x[5]+M[9][6]*x[6]+M[8][7]*x[7]+M[8][8]*x[8], A[8]);
//    printf("%18.10f %18.10f\n", M[8][0]*x[0]+M[8][1]*x[1]+M[8][2]*x[2]+M[8][3]*x[3]+M[8][4]*x[4]+M[8][5]*x[5]+M[8][6]*x[6]+M[8][7]*x[7]+M[8][8]*x[8], A[8]);
//    printf("%18.10f %18.10f\n", M[8][0]*x[0]+M[8][1]*x[1]+M[8][2]*x[2]+M[8][3]*x[3]+M[8][4]*x[4]+M[8][5]*x[5]+M[8][6]*x[6]+M[8][7]*x[7]+M[8][8]*x[8], A[8]);
//    printf("%18.10f %18.10f\n", M[8][0]*x[0]+M[8][1]*x[1]+M[8][2]*x[2]+M[8][3]*x[3]+M[8][4]*x[4]+M[8][5]*x[5]+M[8][6]*x[6]+M[8][7]*x[7]+M[8][8]*x[8], A[8]);
//    printf("%18.10f %18.10f\n", M[8][0]*x[0]+M[8][1]*x[1]+M[8][2]*x[2]+M[8][3]*x[3]+M[8][4]*x[4]+M[8][5]*x[5]+M[8][6]*x[6]+M[8][7]*x[7]+M[8][8]*x[8], A[8]);
//    printf("%18.10f %18.10f\n", M[8][0]*x[0]+M[8][1]*x[1]+M[8][2]*x[2]+M[8][3]*x[3]+M[8][4]*x[4]+M[8][5]*x[5]+M[8][6]*x[6]+M[8][7]*x[7]+M[8][8]*x[8], A[8]);

    DoubleMatrix nx(3, N+1, 0.0);
    nx[0][N-4] = x[0]; nx[0][N-3] = x[3]; nx[0][N-2] = x[6]; nx[0][N-1] = x[9];  nx[0][N-0] = x[12];
    nx[1][N-4] = x[1]; nx[1][N-3] = x[4]; nx[1][N-2] = x[7]; nx[1][N-1] = x[10]; nx[1][N-0] = x[13];
    nx[2][N-4] = x[2]; nx[2][N-3] = x[5]; nx[2][N-2] = x[8]; nx[2][N-1] = x[11]; nx[2][N-0] = x[14];

//    //IPrinter::printSeperatorLine();
//    //printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
//    //printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
//    //printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);

    betta[N-3] = betta[N-3] - betta[N-4]*alpha1[N-4];
    betta[N-2] = betta[N-2] - betta[N-4]*alpha2[N-4];
    betta[N-1] = betta[N-1] - betta[N-4]*alpha3[N-4];
    betta[N-0] = betta[N-0] - betta[N-4]*alpha4[N-4];
    eta        = eta + DoubleVector(betta[N-4]*alpha0[N-4]);

    for (unsigned int k=N-4; k>0; k--)
    {
        betta[k+0] = betta[k+0] - betta[k-1]*alpha1[k-1];
        betta[k+1] = betta[k+1] - betta[k-1]*alpha2[k-1];
        betta[k+2] = betta[k+2] - betta[k-1]*alpha3[k-1];
        betta[k+3] = betta[k+3] - betta[k-1]*alpha4[k-1];
        eta        = eta + DoubleVector(betta[k-1]*alpha0[k-1]);

        nx[0][k-1] = eta[0];
        nx[1][k-1] = eta[1];
        nx[2][k-1] = eta[2];
        for (unsigned int i=k; i<=N; i++)
        {
            DoubleVector col = nx.col(k-1);
            DoubleVector aaa = DoubleVector(betta[i]*nx.col(i));
            col[0] -= aaa[0];
            col[1] -= aaa[1];
            col[2] -= aaa[2];

            //DoubleVector col = nx.col(k-1)-DoubleVector(betta[i]*nx.col(i));
            nx.setColumn(k-1, col);
        }
        DoubleMatrix bbb = betta[k-1];
        bbb.inverse();
        nx.setColumn(k-1, DoubleVector(bbb*nx.col(k-1)));
    }

    IPrinter::printVector(w,p,nx.row(0));
    IPrinter::printVector(w,p,nx.row(1));
    IPrinter::printVector(w,p,nx.row(2));

//    FILE *file =fopen("data_rx.txt", "a");
//    IPrinter::printVector(14,10,nx,"nv2",nx.size(),0,0,file);
//    fclose(file);
}

void SystemDifEquation::calculateRX(DoubleMatrix &rx)
{
    rx.clear();
    rx.resize(3, N+1, 0.0);

    for (unsigned int k=0; k<=N; k++)
    {
        rx[0][k] = f(1,k);
        rx[1][k] = f(2,k);
        rx[2][k] = f(3,k);
    }
    IPrinter::printVector(w,p,rx.row(0));
    IPrinter::printVector(w,p,rx.row(1));
    IPrinter::printVector(w,p,rx.row(2));

//    //    FILE *file1 =fopen("data_rx.txt", "w");
//    //    IPrinter::printVector(14,10,rx,"rx0",rx.size(),0,0,file1);
//    //    fclose(file1);
}

double SystemDifEquation::a(unsigned int i, unsigned int j, unsigned int k) const
{
    if (i==1 && j==1) return +2.0;
    if (i==1 && j==2) return +3.0;
    if (i==1 && j==3) return -1.0;

    if (i==2 && j==1) return +4.0;
    if (i==2 && j==2) return +6.0;
    if (i==2 && j==3) return -2.0;

    if (i==3 && j==1) return -1.0;
    if (i==3 && j==2) return +1.0;
    if (i==3 && j==3) return -1.0;

    return NAN;
}

double SystemDifEquation::b(unsigned int i, unsigned int k) const
{
    double t = k*h;

#ifdef SAMPLE_1
    if (i==1) return -(+2.0*sin(20.0*t*t) + 3.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+40.0*t*cos(20.0*t*t));
    if (i==2) return -(+4.0*sin(20.0*t*t) + 6.0*(cos(10.0*t) - sin(20.0*t)) - 2.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (-10.0*sin(10.0*t) - 20.0*cos(20.0*t));
    if (i==3) return -(-1.0*sin(20.0*t*t) + 1.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+3.0*t*t - 16.0*cos(8.0*t)*sin(8.0*t));
#endif
#ifdef SAMPLE_2
    if (i==1) return -(+2.0*sin(t*t) + 3.0*(cos(t) - sin(t)) - 1.0*(t*t*t - sin(t)*sin(t))) + (+2.0*t*cos(t*t));
    if (i==2) return -(+4.0*sin(t*t) + 6.0*(cos(t) - sin(t)) - 2.0*(t*t*t - sin(t)*sin(t))) + (-sin(t) - cos(t));
    if (i==3) return -(-1.0*sin(t*t) + 1.0*(cos(t) - sin(t)) - 1.0*(t*t*t - sin(t)*sin(t))) + (+3.0*t*t - 2.0*cos(t)*sin(t));
#endif
    return NAN;
}

double SystemDifEquation::f(unsigned int i, unsigned int k) const
{
    double t = k*h;
#ifdef SAMPLE_1
    if (i==1) return sin(20.0*t*t);
    if (i==2) return cos(10.0*t) - sin(20.0*t);
    if (i==3) return t*t*t - sin(8.0*t)*sin(8.0*t);
#endif
#ifdef SAMPLE_2
    if (i==1) return sin(t*t);
    if (i==2) return cos(t) - sin(t);
    if (i==3) return t*t*t - sin(t)*sin(t);
#endif
    return NAN;
}
