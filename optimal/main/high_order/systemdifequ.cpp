#include "systemdifequ.h"

double norm_1(const DoubleVector &v1, const DoubleVector &v2)
{
    double mod = 0.0;
    for (unsigned int i=0; i<v1.size(); i++)
        mod += (v1[i]-v2[i])*(v1[i]-v2[i]);
    return sqrt(mod);
}

double Determinant(const DoubleMatrix &m)
{
    return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
          -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
          +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][1]);
}

void inverse(const DoubleMatrix &src, DoubleMatrix &dst)
{
    dst.clear();
    dst.resize(src.rows(), src.cols());
    double det = Determinant(src);

    for (unsigned int r=0; r<src.rows(); r++)
    {
        for (unsigned int c=0; c<src.cols(); c++)
        {
            DoubleMatrix minor = src.minor(r,c);
            dst[r][c] = minor.determinant();
            if ((r+c)%2==1) dst[r][c] *= -1.0;
        }
    }

    dst.transpose();

    for (unsigned int r=0; r<src.rows(); r++)
    {
        for (unsigned int c=0; c<src.cols(); c++)
        {
            dst[r][c] *= 1.0/det;
        }
    }
}

void SystemDifEquation::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    SystemDifEquation e;

    DoubleMatrix rx;
    e.calculateRX(rx);
    IPrinter::printSeperatorLine();

    e.calculate2R2LV1(rx);
    e.calculate4R2LV1(rx);
    e.calculate6R2LV1(rx);
}

SystemDifEquation::SystemDifEquation()
{
    h = 0.01;
    N = 100;
    w = 16;
    p = 12;
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

//        DoubleMatrix ccc = AI;
//        AI.inverse();
//        DoubleMatrix ddd = ccc*AI;
//        IPrinter::print(ddd,ddd.rows(),ddd.cols());

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

    printf("det: %.10f\n", M.determinant());
    //IPrinter::printSeperatorLine();
    //IPrinter::print(M,9,9);
    //IPrinter::printSeperatorLine();
    //for (unsigned int i=0; i<9; i++)
    //{
    //    printf("%4d %18.10f %18.10f\n", i,
    //           M[i][0]*rx[0][N-2]+ M[i][1]*rx[1][N-2]+ M[i][2]*rx[2][N-2]+
    //           M[i][3]*rx[0][N-1]+ M[i][4]*rx[1][N-1]+ M[i][5]*rx[2][N-1]+
    //           M[i][6]*rx[0][N-0]+ M[i][7]*rx[1][N-0]+ M[i][8]*rx[2][N-0], A[i]);
    //}

    GaussianElimination(M, A, x);

    DoubleMatrix nx(3, N+1, 0.0);
    nx[0][N-2] = x[0]; nx[0][N-1] = x[3]; nx[0][N-0] = x[6];
    nx[1][N-2] = x[1]; nx[1][N-1] = x[4]; nx[1][N-0] = x[7];
    nx[2][N-2] = x[2]; nx[2][N-1] = x[5]; nx[2][N-0] = x[8];

//    IPrinter::printSeperatorLine();
//    printf("%18.12f %18.12f %18.12f\n", nx[0][N-2], nx[0][N-1], nx[0][N-0]);
//    printf("%18.12f %18.12f %18.12f\n", rx[0][N-2], rx[0][N-1], rx[0][N-0]);
//    printf("%18.12f %18.12f %18.12f\n", nx[1][N-2], nx[1][N-1], nx[1][N-0]);
//    printf("%18.12f %18.12f %18.12f\n", rx[1][N-2], rx[1][N-1], rx[1][N-0]);
//    printf("%18.12f %18.12f %18.12f\n", nx[2][N-2], nx[2][N-1], nx[2][N-0]);
//    printf("%18.12f %18.12f %18.12f\n", rx[2][N-2], rx[2][N-1], rx[2][N-0]);
//    IPrinter::printSeperatorLine();

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
            nx.setColumn(k-1, col);
        }
        DoubleMatrix bbb = betta[k-1];
        bbb.inverse();
        nx.setColumn(k-1, DoubleVector(bbb*nx.col(k-1)));
    }

    IPrinter::printVector(w,p,nx.row(0));
    IPrinter::printVector(w,p,nx.row(1));
    IPrinter::printVector(w,p,nx.row(2));

    printf("mod1: %.12f\n", norm_1(rx.row(0), nx.row(0)));
    printf("mod2: %.12f\n", norm_1(rx.row(1), nx.row(1)));
    printf("mod3: %.12f\n", norm_1(rx.row(2), nx.row(2)));

    FILE *file =fopen("data_rx.txt", "a");
    IPrinter::printVector(14,10,nx.row(0),"nv21",N+1,0,0,file);
    IPrinter::printVector(14,10,nx.row(1),"nv22",N+1,0,0,file);
    IPrinter::printVector(14,10,nx.row(2),"nv22",N+1,0,0,file);
    fclose(file);
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

    printf("det: %.10f\n", M.determinant());
    //IPrinter::printSeperatorLine();
    //IPrinter::print(M,15,15,10,6);
    //IPrinter::printSeperatorLine();
    //for (unsigned int i=0; i<15; i++)
    //{
    //    printf("%4d %18.10f %18.10f\n", i,
    //           M[i][0]*rx[0][N-4]+ M[i][1]*rx[1][N-4]+ M[i][2]*rx[2][N-4]+
    //           M[i][3]*rx[0][N-3]+ M[i][4]*rx[1][N-3]+ M[i][5]*rx[2][N-3]+
    //           M[i][6]*rx[0][N-2]+ M[i][7]*rx[1][N-2]+ M[i][8]*rx[2][N-2]+
    //           M[i][9]*rx[0][N-1]+ M[i][10]*rx[1][N-1]+M[i][11]*rx[2][N-1]+
    //           M[i][12]*rx[0][N-0]+M[i][13]*rx[1][N-0]+M[i][14]*rx[2][N-0], A[i]);
    //}

    GaussianElimination(M, A, x);

    DoubleMatrix nx(3, N+1, 0.0);
    nx[0][N-4] = x[0]; nx[0][N-3] = x[3]; nx[0][N-2] = x[6]; nx[0][N-1] = x[9];  nx[0][N-0] = x[12];
    nx[1][N-4] = x[1]; nx[1][N-3] = x[4]; nx[1][N-2] = x[7]; nx[1][N-1] = x[10]; nx[1][N-0] = x[13];
    nx[2][N-4] = x[2]; nx[2][N-3] = x[5]; nx[2][N-2] = x[8]; nx[2][N-1] = x[11]; nx[2][N-0] = x[14];

    IPrinter::printSeperatorLine();
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f\n", nx[0][N-4], nx[0][N-3], nx[0][N-2], nx[0][N-1], nx[0][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f\n", rx[0][N-4], rx[0][N-3], rx[0][N-2], rx[0][N-1], rx[0][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f\n", nx[1][N-4], nx[1][N-3], nx[1][N-2], nx[1][N-1], nx[1][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f\n", rx[1][N-4], rx[1][N-3], rx[1][N-2], rx[1][N-1], rx[1][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f\n", nx[2][N-4], nx[2][N-3], nx[2][N-2], nx[2][N-1], nx[2][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f\n", rx[2][N-4], rx[2][N-3], rx[2][N-2], rx[2][N-1], rx[2][N-0]);
    IPrinter::printSeperatorLine();

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
            nx.setColumn(k-1, col);
        }
        DoubleMatrix bbb = betta[k-1];
        bbb.inverse();
        nx.setColumn(k-1, DoubleVector(bbb*nx.col(k-1)));
    }

    IPrinter::printVector(w,p,nx.row(0));
    IPrinter::printVector(w,p,nx.row(1));
    IPrinter::printVector(w,p,nx.row(2));

    printf("mod1: %.12f\n", norm_1(rx.row(0), nx.row(0)));
    printf("mod2: %.12f\n", norm_1(rx.row(1), nx.row(1)));
    printf("mod3: %.12f\n", norm_1(rx.row(2), nx.row(2)));

    FILE *file =fopen("data_rx.txt", "a");
    IPrinter::printVector(14,10,nx.row(0),"nv41",N+1,0,0,file);
    IPrinter::printVector(14,10,nx.row(1),"nv42",N+1,0,0,file);
    IPrinter::printVector(14,10,nx.row(2),"nv42",N+1,0,0,file);
    fclose(file);

}

void SystemDifEquation::calculate6R2LV1(const DoubleMatrix &rx)
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

    DoubleMatrix* alpha5 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha5[i].resize(3,3,0.0);

    DoubleMatrix* alpha6 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha6[i].resize(3,3,0.0);

    DoubleMatrix* alpha0 = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) alpha0[i].resize(3,1,0.0);

    DoubleMatrix AI(3,3);
    for (unsigned int k=0; k<=N-6; k++)
    {
        AI[0][0] = -60.0*h*a(1,1,k) - 147.0; AI[0][1] = -60.0*h*a(1,2,k) + 0.0;   AI[0][2] = -60.0*h*a(1,3,k) + 0.0;
        AI[1][0] = -60.0*h*a(2,1,k) + 0.0;   AI[1][1] = -60.0*h*a(2,2,k) - 147.0; AI[1][2] = -60.0*h*a(2,3,k) + 0.0;
        AI[2][0] = -60.0*h*a(3,1,k) + 0.0;   AI[2][1] = -60.0*h*a(3,2,k) + 0.0;   AI[2][2] = -60.0*h*a(3,3,k) - 147.0;
        AI.inverse();

        alpha1[k][0][0] = -360.0; alpha1[k][0][1] =   +0.0; alpha1[k][0][2] =  +0.0;
        alpha1[k][1][0] =   +0.0; alpha1[k][1][1] = -360.0; alpha1[k][1][2] =  +0.0;
        alpha1[k][2][0] =   +0.0; alpha1[k][2][1] =   +0.0; alpha1[k][2][2] = -360.0;

        alpha2[k][0][0] = +450.0; alpha2[k][0][1] = +0.0;   alpha2[k][0][2] = +0.0;
        alpha2[k][1][0] =   +0.0; alpha2[k][1][1] = +450.0; alpha2[k][1][2] = +0.0;
        alpha2[k][2][0] =   +0.0; alpha2[k][2][1] =  +0.0;  alpha2[k][2][2] = +450.0;

        alpha3[k][0][0] = -400.0; alpha3[k][0][1] =   +0.0; alpha3[k][0][2] = +0.0;
        alpha3[k][1][0] =   +0.0; alpha3[k][1][1] = -400.0; alpha3[k][1][2] = +0.0;
        alpha3[k][2][0] =   +0.0; alpha3[k][2][1] =   +0.0; alpha3[k][2][2] = -400.0;

        alpha4[k][0][0] = +225.0; alpha4[k][0][1] =   +0.0; alpha4[k][0][2] = +0.0;
        alpha4[k][1][0] =   +0.0; alpha4[k][1][1] = +225.0; alpha4[k][1][2] = +0.0;
        alpha4[k][2][0] =   +0.0; alpha4[k][2][1] =   +0.0; alpha4[k][2][2] = +225.0;

        alpha5[k][0][0] = -72.0;  alpha5[k][0][1] =  +0.0;  alpha5[k][0][2] =  +0.0;
        alpha5[k][1][0] =  +0.0;  alpha5[k][1][1] = -72.0;  alpha5[k][1][2] =  +0.0;
        alpha5[k][2][0] =  +0.0;  alpha5[k][2][1] =  +0.0;  alpha5[k][2][2] = -72.0;

        alpha6[k][0][0] =  +10.0; alpha6[k][0][1] =  +0.0;  alpha6[k][0][2] = +0.0;
        alpha6[k][1][0] =   +0.0; alpha6[k][1][1] = +10.0;  alpha6[k][1][2] = +0.0;
        alpha6[k][2][0] =   +0.0; alpha6[k][2][1] =  +0.0;  alpha6[k][2][2] = +10.0;

        alpha0[k][0][0] = +60.0*h*b(1,k);
        alpha0[k][1][0] = +60.0*h*b(2,k);
        alpha0[k][2][0] = +60.0*h*b(3,k);

        alpha1[k] = AI*alpha1[k];
        alpha2[k] = AI*alpha2[k];
        alpha3[k] = AI*alpha3[k];
        alpha4[k] = AI*alpha4[k];
        alpha5[k] = AI*alpha5[k];
        alpha6[k] = AI*alpha6[k];
        alpha0[k] = AI*alpha0[k];
    }

    IPrinter::printSeperatorLine();

    for (unsigned int k=0; k<=N-6; k++)
    {
        betta[k+1] = betta[k+1] + betta[k]*alpha1[k];
        betta[k+2] = betta[k+2] + betta[k]*alpha2[k];
        betta[k+3] = betta[k+3] + betta[k]*alpha3[k];
        betta[k+4] = betta[k+4] + betta[k]*alpha4[k];
        betta[k+5] = betta[k+5] + betta[k]*alpha5[k];
        betta[k+6] = betta[k+6] + betta[k]*alpha6[k];
        eta        = eta        - betta[k]*alpha0[k];
    }

    DoubleMatrix M(21,21);
    DoubleVector A(21);
    DoubleVector x(21);

    M[0][0]  = 0.0;              M[0][1] = 0.0;               M[0][2] = 0.0;
    M[0][3]  = betta[N-5][0][0]; M[0][4]  = betta[N-5][0][1]; M[0][5]  = betta[N-5][0][2];
    M[0][6]  = betta[N-4][0][0]; M[0][7]  = betta[N-4][0][1]; M[0][8]  = betta[N-4][0][2];
    M[0][9]  = betta[N-3][0][0]; M[0][10] = betta[N-3][0][1]; M[0][11] = betta[N-3][0][2];
    M[0][12] = betta[N-2][0][0]; M[0][13] = betta[N-2][0][1]; M[0][14] = betta[N-2][0][2];
    M[0][15] = betta[N-1][0][0]; M[0][16] = betta[N-1][0][1]; M[0][17] = betta[N-1][0][2];
    M[0][18] = betta[N-0][0][0]; M[0][19] = betta[N-0][0][1]; M[0][20] = betta[N-0][0][2]; A[0] = eta[0];

    M[1][0]  = 0.0;              M[1][1]  = 0.0;              M[1][2]  = 0.0;
    M[1][3]  = betta[N-5][1][0]; M[1][4]  = betta[N-5][1][1]; M[1][5]  = betta[N-5][1][2];
    M[1][6]  = betta[N-4][1][0]; M[1][7]  = betta[N-4][1][1]; M[1][8]  = betta[N-4][1][2];
    M[1][9]  = betta[N-3][1][0]; M[1][10] = betta[N-3][1][1]; M[1][11] = betta[N-3][1][2];
    M[1][12] = betta[N-2][1][0]; M[1][13] = betta[N-2][1][1]; M[1][14] = betta[N-2][1][2];
    M[1][15] = betta[N-1][1][0]; M[1][16] = betta[N-1][1][1]; M[1][17] = betta[N-1][1][2];
    M[1][18] = betta[N-0][1][0]; M[1][19] = betta[N-0][1][1]; M[1][20] = betta[N-0][1][2]; A[1] = eta[1];

    M[2][0]  = 0.0;              M[2][1]  = 0.0;              M[2][2]  = 0.0;
    M[2][3]  = betta[N-5][2][0]; M[2][4]  = betta[N-5][2][1]; M[2][5]  = betta[N-5][2][2];
    M[2][6]  = betta[N-4][2][0]; M[2][7]  = betta[N-4][2][1]; M[2][8]  = betta[N-4][2][2];
    M[2][9]  = betta[N-3][2][0]; M[2][10] = betta[N-3][2][1]; M[2][11] = betta[N-3][2][2];
    M[2][12] = betta[N-2][2][0]; M[2][13] = betta[N-2][2][1]; M[2][14] = betta[N-2][2][2];
    M[2][15] = betta[N-1][2][0]; M[2][16] = betta[N-1][2][1]; M[2][17] = betta[N-1][2][2];
    M[2][18] = betta[N-0][2][0]; M[2][19] = betta[N-0][2][1]; M[2][20] = betta[N-0][2][2]; A[2] = eta[2];

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    M[3][0]  = -10.0;                   M[3][1]  = +0.0;                    M[3][2]  = +0.0;
    M[3][3]  = -60.0*h*a(1,1,N-5)-77.0; M[3][4]  = -60.0*h*a(1,2,N-5);      M[3][5]  = -60.0*h*a(1,3,N-5);
    M[3][6]  = +150.0;                  M[3][7]  = +0.0;                    M[3][8]  = +0.0;
    M[3][9]  = -100.0;                  M[3][10] = +0.0;                    M[3][11] = +0.0;
    M[3][12] = +50.0;                   M[3][13] = +0.0;                    M[3][14] = +0.0;
    M[3][15] = -15.0;                   M[3][16] = +0.0;                    M[3][17] = +0.0;
    M[3][18] = +2.0;                    M[3][19] = +0.0;                    M[3][20] = +0.0; A[3] = +60.0*h*b(1,N-5);

    M[4][0]  = +0.0;                    M[4][1]  = -10.0;                   M[4][2]  = +0.0;
    M[4][3]  = -60.0*h*a(2,1,N-5);      M[4][4]  = -60.0*h*a(2,2,N-5)-77.0; M[4][5]  = -60.0*h*a(2,3,N-5);
    M[4][6]  = +0.0;                    M[4][7]  = +150.0;                  M[4][8]  = +0.0;
    M[4][9]  = +0.0;                    M[4][10] = -100.0;                  M[4][11] = +0.0;
    M[4][12] = +0.0;                    M[4][13] = +50.0;                   M[4][14] = +0.0;
    M[4][15] = +0.0;                    M[4][16] = -15.0;                   M[4][17] = +0.0;
    M[4][18] = +0.0;                    M[4][19] = +2.0;                    M[4][20] = +0.0; A[4] = +60.0*h*b(2,N-5);

    M[5][0]  = +0.0;                    M[5][1]  = +0.0;                    M[5][2]  = -10.0;
    M[5][3]  = -60.0*h*a(3,1,N-5);      M[5][4]  = -60.0*h*a(3,2,N-5);      M[5][5]  = -60.0*h*a(3,3,N-5)-77.0;
    M[5][6]  = +0.0;                    M[5][7]  = +0.0;                    M[5][8]  = +150.0;
    M[5][9]  = +0.0;                    M[5][10] = +0.0;                    M[5][11] = -100.0;
    M[5][12] = +0.0;                    M[5][13] = +0.0;                    M[5][14] = +50.0;
    M[5][15] = +0.0;                    M[5][16] = +0.0;                    M[5][17] = -15.0;
    M[5][18] = +0.0;                    M[5][19] = +0.0;                    M[5][20] = +2.0; A[5] = +60.0*h*b(3,N-5);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    M[6][0]  = +2.0;                    M[6][1]  = +0.0;                    M[6][2]  = +0.0;
    M[6][3]  = -24.0;                   M[6][4]  = +0.0;                    M[6][5]  = +0.0;
    M[6][6]  = -60.0*h*a(1,1,N-4)-35.0; M[6][7]  = -60.0*h*a(1,2,N-4);      M[6][8]  = -60.0*h*a(1,3,N-4);
    M[6][9]  = +80.0;                   M[6][10] = +0.0;                    M[6][11] = +0.0;
    M[6][12] = -30.0;                   M[6][13] = +0.0;                    M[6][14] = +0.0;
    M[6][15] = +8.0;                    M[6][16] = +0.0;                    M[6][17] = +0.0;
    M[6][18] = -1.0;                    M[6][19] = +0.0;                    M[6][20] = +0.0; A[6] = +60.0*h*b(1,N-4);

    M[7][0]  = +0.0;                    M[7][1]  = +2.0;                    M[7][2]  = +0.0;
    M[7][3]  = +0.0;                    M[7][4]  = -24.0;                   M[7][5]  = +0.0;
    M[7][6]  = -60.0*h*a(2,1,N-4);      M[7][7]  = -60.0*h*a(2,2,N-4)-35.0; M[7][8]  = -60.0*h*a(2,3,N-4);
    M[7][9]  = +0.0;                    M[7][10] = +80.0;                   M[7][11] = +0.0;
    M[7][12] = +0.0;                    M[7][13] = -30.0;                   M[7][14] = +0.0;
    M[7][15] = +0.0;                    M[7][16] = +8.0;                    M[7][17] = +0.0;
    M[7][18] = +0.0;                    M[7][19] = -1.0;                    M[7][20] = +0.0; A[7] = +60.0*h*b(2,N-4);

    M[8][0]  = +0.0;                    M[8][1]  = +0.0;                    M[8][2]  = +2.0;
    M[8][3]  = +0.0;                    M[8][4]  = +0.0;                    M[8][5]  = -24.0;
    M[8][6]  = -60.0*h*a(3,1,N-4);      M[8][7]  = -60.0*h*a(3,2,N-4);      M[8][8]  = -60.0*h*a(3,3,N-4)-35.0;
    M[8][9]  = +0.0;                    M[8][10] = +0.0;                    M[8][11] = +80.0;
    M[8][12] = +0.0;                    M[8][13] = +0.0;                    M[8][14] = -30.0;
    M[8][15] = +0.0;                    M[8][16] = +0.0;                    M[8][17] = +8.0;
    M[8][18] = +0.0;                    M[8][19] = +0.0;                    M[8][20] = -1.0; A[8] = +60.0*h*b(3,N-4);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    M[9][0]  = -1.0;                    M[9][1]  = +0.0;                    M[9][2]  = +0.0;
    M[9][3]  = +9.0;                    M[9][4]  = +0.0;                    M[9][5]  = +0.0;
    M[9][6]  = -45.0;                   M[9][7]  = +0.0;                    M[9][8]  = +0.0;
    M[9][9]  = -60.0*h*a(1,1,N-3);      M[9][10] = -60.0*h*a(1,2,N-3);      M[9][11] = -60.0*h*a(1,3,N-3);
    M[9][12] = +45.0;                   M[9][13] = +0.0;                    M[9][14] = +0.0;
    M[9][15] = -9.0;                    M[9][16] = +0.0;                    M[9][17] = +0.0;
    M[9][18] = +1.0;                    M[9][19] = +0.0;                    M[9][20] = +0.0; A[9] = +60.0*h*b(1,N-3);

    M[10][0]  = +0.0;                   M[10][1]  = -1.0;                   M[10][2]  = +0.0;
    M[10][3]  = +0.0;                   M[10][4]  = +9.0;                   M[10][5]  = +0.0;
    M[10][6]  = +0.0;                   M[10][7]  = -45.0;                  M[10][8]  = +0.0;
    M[10][9]  = -60.0*h*a(2,1,N-3);     M[10][10] = -60.0*h*a(2,2,N-3);     M[10][11] = -60.0*h*a(2,3,N-3);
    M[10][12] = +0.0;                   M[10][13] = +45.0;                  M[10][14] = +0.0;
    M[10][15] = +0.0;                   M[10][16] = -9.0;                   M[10][17] = +0.0;
    M[10][18] = +0.0;                   M[10][19] = +1.0;                   M[10][20] = +0.0; A[10] = +60.0*h*b(2,N-3);

    M[11][0]  = +0.0;                   M[11][1]  = +0.0;                   M[11][2]  = -1.0;
    M[11][3]  = +0.0;                   M[11][4]  = +0.0;                   M[11][5]  = +9.0;
    M[11][6]  = +0.0;                   M[11][7]  = +0.0;                   M[11][8]  = -45.0;
    M[11][9]  = -60.0*h*a(3,1,N-3);     M[11][10] = -60.0*h*a(3,2,N-3);     M[11][11] = -60.0*h*a(3,3,N-3);
    M[11][12] = +0.0;                   M[11][13] = +0.0;                   M[11][14] = +45.0;
    M[11][15] = +0.0;                   M[11][16] = +0.0;                   M[11][17] = -9.0;
    M[11][18] = +0.0;                   M[11][19] = +0.0;                   M[11][20] = +1.0; A[11] = +60.0*h*b(3,N-3);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    M[12][0]  = +1.0;                    M[12][1]  = +0.0;                    M[12][2]  = +0.0;
    M[12][3]  = -8.0;                    M[12][4]  = +0.0;                    M[12][5]  = +0.0;
    M[12][6]  = +30.0;                   M[12][7]  = +0.0;                    M[12][8]  = +0.0;
    M[12][9]  = -80.0;                   M[12][10] = +0.0;                    M[12][11] = +0.0;
    M[12][12] = -60.0*h*a(1,1,N-2)+35.0; M[12][13] = -60.0*h*a(1,2,N-2);      M[12][14] = -60.0*h*a(1,3,N-2);
    M[12][15] = +24.0;                   M[12][16] = +0.0;                    M[12][17] = +0.0;
    M[12][18] = -2.0;                    M[12][19] = +0.0;                    M[12][20] = +0.0; A[12] = +60.0*h*b(1,N-2);

    M[13][0]  = +0.0;                    M[13][1]  = +1.0;                    M[13][2]  = +0.0;
    M[13][3]  = +0.0;                    M[13][4]  = -8.0;                    M[13][5]  = +0.0;
    M[13][6]  = +0.0;                    M[13][7]  = +30.0;                   M[13][8]  = +0.0;
    M[13][9]  = +0.0;                    M[13][10] = -80.0;                   M[13][11] = +0.0;
    M[13][12] = -60.0*h*a(2,1,N-2);      M[13][13] = -60.0*h*a(2,2,N-2)+35.0; M[13][14] = -60.0*h*a(2,3,N-2);
    M[13][15] = +0.0;                    M[13][16] = +24.0;                   M[13][17] = +0.0;
    M[13][18] = +0.0;                    M[13][19] = -2.0;                    M[13][20] = +0.0; A[13] = +60.0*h*b(2,N-2);

    M[14][0]  = +0.0;                    M[14][1]  = +0.0;                    M[14][2]  = +1.0;
    M[14][3]  = +0.0;                    M[14][4]  = +0.0;                    M[14][5]  = -8.0;
    M[14][6]  = +0.0;                    M[14][7]  = +0.0;                    M[14][8]  = +30.0;
    M[14][9]  = +0.0;                    M[14][10] = +0.0;                    M[14][11] = -80.0;
    M[14][12] = -60.0*h*a(3,1,N-2);      M[14][13] = -60.0*h*a(3,2,N-2);      M[14][14] = -60.0*h*a(3,3,N-2)+35.0;
    M[14][15] = +0.0;                    M[14][16] = +0.0;                    M[14][17] = +24.0;
    M[14][18] = +0.0;                    M[14][19] = +0.0;                    M[14][20] = -2.0; A[14] = +60.0*h*b(3,N-2);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    M[15][0]  = -2.0;                    M[15][1]  = +0.0;                    M[15][2]  = +0.0;
    M[15][3]  = +15.0;                   M[15][4]  = +0.0;                    M[15][5]  = +0.0;
    M[15][6]  = -50.0;                   M[15][7]  = +0.0;                    M[15][8]  = +0.0;
    M[15][9]  = +100.0;                  M[15][10] = +0.0;                    M[15][11] = +0.0;
    M[15][12] = -150.0;                  M[15][13] = +0.0;                    M[15][14] = +0.0;
    M[15][15] = -60.0*h*a(1,1,N-1)+77.0; M[15][16] = -60.0*h*a(1,2,N-1);      M[15][17] = -60.0*h*a(1,3,N-1);
    M[15][18] = +10.0;                   M[15][19] = +0.0;                    M[15][20] = +0.0; A[15] = +60.0*h*b(1,N-1);

    M[16][0]  = +0.0;                    M[16][1]  = -2.0;                    M[16][2]  = +0.0;
    M[16][3]  = +0.0;                    M[16][4]  = +15.0;                   M[16][5]  = +0.0;
    M[16][6]  = +0.0;                    M[16][7]  = -50.0;                   M[16][8]  = +0.0;
    M[16][9]  = +0.0;                    M[16][10] = +100.0;                  M[16][11] = +0.0;
    M[16][12] = +0.0;                    M[16][13] = -150.0;                  M[16][14] = +0.0;
    M[16][15] = -60.0*h*a(2,1,N-1);      M[16][16] = -60.0*h*a(2,2,N-1)+77.0; M[16][17] = -60.0*h*a(2,3,N-1);
    M[16][18] = +0.0;                    M[16][19] = +10.0;                   M[16][20] = +0.0; A[16] = +60.0*h*b(2,N-1);

    M[17][0]  = +0.0;                    M[17][1]  = +0.0;                    M[17][2]  = -2.0;
    M[17][3]  = +0.0;                    M[17][4]  = +0.0;                    M[17][5]  = +15.0;
    M[17][6]  = +0.0;                    M[17][7]  = +0.0;                    M[17][8]  = -50.0;
    M[17][9]  = +0.0;                    M[17][10] = +0.0;                    M[17][11] = +100.0;
    M[17][12] = +0.0;                    M[17][13] = +0.0;                    M[17][14] = -150.0;
    M[17][15] = -60.0*h*a(3,1,N-1);      M[17][16] = -60.0*h*a(3,2,N-1);      M[17][17] = -60.0*h*a(3,3,N-1)+77.0;
    M[17][18] = +0.0;                    M[17][19] = +0.0;                    M[17][20] = +10.0; A[17] = +60.0*h*b(3,N-1);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    M[18][0]  = +10.0;                    M[18][1]  = +0.0;                   M[18][2]  = +0.0;
    M[18][3]  = -72.0;                    M[18][4]  = +0.0;                   M[18][5]  = +0.0;
    M[18][6]  = +225.0;                   M[18][7]  = +0.0;                   M[18][8]  = +0.0;
    M[18][9]  = -400.0;                   M[18][10] = +0.0;                   M[18][11] = +0.0;
    M[18][12] = +450.0;                   M[18][13] = +0.0;                   M[18][14] = +0.0;
    M[18][15] = -360.0;                   M[18][16] = +0.0;                   M[18][17] = +0.0;
    M[18][18] = -60.0*h*a(1,1,N-0)+147.0; M[18][19] = -60.0*h*a(1,2,N-0);     M[18][20] = -60.0*h*a(1,3,N-0); A[18] = +60.0*h*b(1,N-0);

    M[19][0]  = +0.0;                    M[19][1]  = +10.0;                    M[19][2]  = +0.0;
    M[19][3]  = +0.0;                    M[19][4]  = -72.0;                    M[19][5]  = +0.0;
    M[19][6]  = +0.0;                    M[19][7]  = +225.0;                   M[19][8]  = +0.0;
    M[19][9]  = +0.0;                    M[19][10] = -400.0;                   M[19][11] = +0.0;
    M[19][12] = +0.0;                    M[19][13] = +450.0;                   M[19][14] = +0.0;
    M[19][15] = +0.0;                    M[19][16] = -360.0;                   M[19][17] = +0.0;
    M[19][18] = -60.0*h*a(2,1,N-0);      M[19][19] = -60.0*h*a(2,2,N-0)+147.0; M[19][20] = -60.0*h*a(2,3,N-0); A[19] = +60.0*h*b(2,N-0);

    M[20][0]  = +0.0;                    M[20][1]  = +0.0;                    M[20][2]  = +10.0;
    M[20][3]  = +0.0;                    M[20][4]  = +0.0;                    M[20][5]  = -72.0;
    M[20][6]  = +0.0;                    M[20][7]  = +0.0;                    M[20][8]  = +225.0;
    M[20][9]  = +0.0;                    M[20][10] = +0.0;                    M[20][11] = -400.0;
    M[20][12] = +0.0;                    M[20][13] = +0.0;                    M[20][14] = +450.0;
    M[20][15] = +0.0;                    M[20][16] = +0.0;                    M[20][17] = -360.0;
    M[20][18] = -60.0*h*a(3,1,N-0);      M[20][19] = -60.0*h*a(3,2,N-0);      M[20][20] = -60.0*h*a(3,3,N-0)+147; A[20] = +60.0*h*b(3,N-0);


    printf("det: %.10f\n", M.determinant());
    //IPrinter::printSeperatorLine();
    //IPrinter::print(M,21,21,7,2);
    //IPrinter::printSeperatorLine();
    //for (unsigned int i=0; i<21; i++)
    //{
    //    printf("%4d %18.10f %18.10f\n", i,
    //           M[i][0]*rx[0][N-6]+ M[i][1]*rx[1][N-6]+ M[i][2]*rx[2][N-6]+
    //           M[i][3]*rx[0][N-5]+ M[i][4]*rx[1][N-5]+ M[i][5]*rx[2][N-5]+
    //           M[i][6]*rx[0][N-4]+ M[i][7]*rx[1][N-4]+ M[i][8]*rx[2][N-4]+
    //           M[i][9]*rx[0][N-3]+ M[i][10]*rx[1][N-3]+M[i][11]*rx[2][N-3]+
    //           M[i][12]*rx[0][N-2]+M[i][13]*rx[1][N-2]+M[i][14]*rx[2][N-2]+
    //           M[i][15]*rx[0][N-1]+M[i][16]*rx[1][N-1]+M[i][17]*rx[2][N-1]+
    //           M[i][18]*rx[0][N-0]+M[i][19]*rx[1][N-0]+M[i][20]*rx[2][N-0], A[i]);
    //}

    GaussianElimination(M, A, x);

    DoubleMatrix nx(3, N+1, 0.0);
    nx[0][N-6] = x[0]; nx[0][N-5] = x[3]; nx[0][N-4] = x[6]; nx[0][N-3] = x[9];  nx[0][N-2] = x[12]; nx[0][N-1] = x[15]; nx[0][N-0] = x[18];
    nx[1][N-6] = x[1]; nx[1][N-5] = x[4]; nx[1][N-4] = x[7]; nx[1][N-3] = x[10]; nx[1][N-2] = x[13]; nx[1][N-1] = x[16]; nx[1][N-0] = x[19];
    nx[2][N-6] = x[2]; nx[2][N-5] = x[5]; nx[2][N-4] = x[8]; nx[2][N-3] = x[11]; nx[2][N-2] = x[14]; nx[2][N-1] = x[17]; nx[2][N-0] = x[20];

    IPrinter::printSeperatorLine();
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f %18.12f %18.12f\n", nx[0][N-6], nx[0][N-5], nx[0][N-4], nx[0][N-3], nx[0][N-2], nx[0][N-1], nx[0][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f %18.12f %18.12f\n", rx[0][N-6], rx[0][N-5], rx[0][N-4], rx[0][N-3], rx[0][N-2], rx[0][N-1], rx[0][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f %18.12f %18.12f\n", nx[1][N-6], nx[1][N-5], nx[1][N-4], nx[1][N-3], nx[1][N-2], nx[1][N-1], nx[1][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f %18.12f %18.12f\n", rx[1][N-6], rx[1][N-5], rx[1][N-4], rx[1][N-3], rx[1][N-2], rx[1][N-1], rx[1][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f %18.12f %18.12f\n", nx[2][N-6], nx[2][N-5], nx[2][N-4], nx[2][N-3], nx[2][N-2], nx[2][N-1], nx[2][N-0]);
    printf("%18.12f %18.12f %18.12f %18.12f %18.12f %18.12f %18.12f\n", rx[2][N-6], rx[2][N-5], rx[2][N-4], rx[2][N-3], rx[2][N-2], rx[2][N-1], rx[2][N-0]);
    IPrinter::printSeperatorLine();

    betta[N-5] = betta[N-5] - betta[N-6]*alpha1[N-6];
    betta[N-4] = betta[N-4] - betta[N-6]*alpha2[N-6];
    betta[N-3] = betta[N-3] - betta[N-6]*alpha3[N-6];
    betta[N-2] = betta[N-2] - betta[N-6]*alpha4[N-6];
    betta[N-1] = betta[N-1] - betta[N-6]*alpha5[N-6];
    betta[N-0] = betta[N-0] - betta[N-6]*alpha6[N-6];
    eta        = eta + DoubleVector(betta[N-6]*alpha0[N-6]);

    for (unsigned int k=N-6; k>0; k--)
    {
        betta[k+0] = betta[k+0] - betta[k-1]*alpha1[k-1];
        betta[k+1] = betta[k+1] - betta[k-1]*alpha2[k-1];
        betta[k+2] = betta[k+2] - betta[k-1]*alpha3[k-1];
        betta[k+3] = betta[k+3] - betta[k-1]*alpha4[k-1];
        betta[k+4] = betta[k+4] - betta[k-1]*alpha5[k-1];
        betta[k+5] = betta[k+5] - betta[k-1]*alpha6[k-1];
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
            nx.setColumn(k-1, col);
        }
        DoubleMatrix bbb = betta[k-1];
        bbb.inverse();
        nx.setColumn(k-1, DoubleVector(bbb*nx.col(k-1)));
    }

    IPrinter::printVector(w,p,nx.row(0));
    IPrinter::printVector(w,p,nx.row(1));
    IPrinter::printVector(w,p,nx.row(2));

    printf("mod1: %.12f\n", norm_1(rx.row(0), nx.row(0)));
    printf("mod2: %.12f\n", norm_1(rx.row(1), nx.row(1)));
    printf("mod3: %.12f\n", norm_1(rx.row(2), nx.row(2)));

    FILE *file =fopen("data_rx.txt", "a");
    IPrinter::printVector(14,10,nx.row(0),"nv61",N+1,0,0,file);
    IPrinter::printVector(14,10,nx.row(1),"nv62",N+1,0,0,file);
    IPrinter::printVector(14,10,nx.row(2),"nv62",N+1,0,0,file);
    fclose(file);
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

    FILE *file =fopen("data_rx.txt", "w");
    IPrinter::printVector(14,10,rx.row(0),"rx01",N+1,0,0,file);
    IPrinter::printVector(14,10,rx.row(1),"rx02",N+1,0,0,file);
    IPrinter::printVector(14,10,rx.row(2),"rx02",N+1,0,0,file);
    fclose(file);
}

double SystemDifEquation::a(unsigned int i, unsigned int j, unsigned int k UNUSED_PARAM) const
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
#ifdef SAMPLE_3
    if (i==1) return -(+2.0*t*t + 3.0*t - 1.0*t*t*t) + (+2.0*t);
    if (i==2) return -(+4.0*t*t + 6.0*t - 2.0*t*t*t) + (+1.0);
    if (i==3) return -(-1.0*t*t + 1.0*t - 1.0*t*t*t) + (+3.0*t*t);
#endif
#ifdef SAMPLE_4
    if (i==1) return -(+2.0*1.0*t + 3.0*2.0*t - 1.0*3.0*t) + (+1.0);
    if (i==2) return -(+4.0*1.0*t + 6.0*2.0*t - 2.0*3.0*t) + (+2.0);
    if (i==3) return -(-1.0*1.0*t + 1.0*2.0*t - 1.0*3.0*t) + (+3.0);
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
#ifdef SAMPLE_3
    if (i==1) return t*t;
    if (i==2) return t;
    if (i==3) return t*t*t;
#endif
#ifdef SAMPLE_4
    if (i==1) return 1.0*t;
    if (i==2) return 2.0*t;
    if (i==3) return 3.0*t;
#endif
    return NAN;
}
