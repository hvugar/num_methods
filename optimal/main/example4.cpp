#include "example4.h"

//double my_rand1()
//{
//    return ((rand() % 1000) + 1) / 1000.0;
//}

FILE *file;

void Example4::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    file = fopen("data.txt", "w");

    Example4 e;
    e.h = 0.01;
    e.N = 100;
    e.F = e.N/10;
    e.n = 3;
    e.K = 4;
    e.w = 14;
    e.p = 8;
    e.L = 4;

    const unsigned int s[][4] =
    {
        {0, 1*e.F, 8*e.F, 10*e.F},
        {0, 2*e.F, 9*e.F, 10*e.F},
        {0, 4*e.F, 5*e.F, 10*e.F},
        {0, 3*e.F, 6*e.F, 10*e.F}
    };

    IPrinter::printSeperatorLine(NULL,'-', stdout);
    printf("s0 %8u %8u %8u %8u\n", s[0][0], s[0][1], s[0][2], s[0][3]);
    printf("s1 %8u %8u %8u %8u\n", s[1][0], s[1][1], s[1][2], s[1][3]);
    printf("s2 %8u %8u %8u %8u\n", s[2][0], s[2][1], s[2][2], s[2][3]);
    printf("s3 %8u %8u %8u %8u\n", s[3][0], s[3][1], s[3][2], s[3][3]);

    //--------------------------------------------------------------------------
    fputs("Real process solution:\n",stdout);
    IPrinter::printSeperatorLine(NULL,'-', stdout);
    DoubleMatrix rx;
    e.calculateRX(rx);
    IPrinter::printVector(e.w,e.p,rx.row(0),"x1: ");
    IPrinter::printVector(e.w,e.p,rx.row(1),"x2: ");
    IPrinter::printVector(e.w,e.p,rx.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);

    //--------------------------------------------------------------------------
    fputs("Initial first 3 points:\n",stdout);
    DoubleMatrix nx;
    e.calculateNX(rx,nx);
    IPrinter::printVector(e.w,e.p,nx.row(0),"x1: ");
    IPrinter::printVector(e.w,e.p,nx.row(1),"x2: ");
    IPrinter::printVector(e.w,e.p,nx.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);

    //--------------------------------------------------------------------------

    fputs("Method #1",file);
    IPrinter::printSeperatorLine(NULL,'-', file);
    e.calculateM1(s, rx, nx);
    IPrinter::printSeperatorLine(NULL,'-', file);

    //puts("Method #2");
    //IPrinter::printSeperatorLine();
    //e.calculateM2(s, rx);

    fclose(file);
}

Example4::Example4()
{
}

void Example4::init(std::vector<DoubleVector> &rx)
{
    DoubleVector tx01;
    DoubleVector tx02;
    DoubleVector tx03;

    rx.resize(N+1);

    for (unsigned int i=0; i<=N; i++)
    {
        DoubleVector &px = rx.at(i);
        px.resize(n);

        px.at(0) = fx1(i);
        px.at(1) = fx2(i);
        px.at(2) = fx3(i);

        tx01 << px.at(0);
        tx02 << px.at(1);
        tx03 << px.at(2);
    }

    IPrinter::printVector(w,p,tx01,NULL,tx01.size(),0,0,file); tx01.clear();
    IPrinter::printVector(w,p,tx02,NULL,tx02.size(),0,0,file); tx02.clear();
    IPrinter::printVector(w,p,tx03,NULL,tx03.size(),0,0,file); tx03.clear();
}

void Example4::initAMatrices(std::vector<DoubleMatrix> &A)
{
    A.clear();
    A.resize(5);
    A[0].resize(n,1);
    A[1].resize(n,n);
    A[2].resize(n,n);
    A[3].resize(n,n);
    A[4].resize(n,n);
}

void Example4::updateAMatrices(std::vector<DoubleMatrix> &A, unsigned int k)
{
    A[1].at(0,0) = 0.48*h*a(1,1,k-1)+1.92; A[1].at(0,1) = 0.48*h*a(1,2,k-1);      A[1].at(0,2) = 0.48*h*a(1,3,k-1);
    A[1].at(1,0) = 0.48*h*a(2,1,k-1);      A[1].at(1,1) = 0.48*h*a(2,2,k-1)+1.92; A[1].at(1,2) = 0.48*h*a(2,3,k-1);
    A[1].at(2,0) = 0.48*h*a(3,1,k-1);      A[1].at(2,1) = 0.48*h*a(3,2,k-1);      A[1].at(2,2) = 0.48*h*a(3,3,k-1)+1.92;

    A[2].at(0,0) = -1.44; A[2].at(0,1) = +0.00; A[2].at(0,2) = +0.00;
    A[2].at(1,0) = +0.00; A[2].at(1,1) = -1.44; A[2].at(1,2) = +0.00;
    A[2].at(2,0) = +0.00; A[2].at(2,1) = +0.00; A[2].at(2,2) = -1.44;

    A[3].at(0,0) = +0.64; A[3].at(0,1) = +0.00; A[3].at(0,2) = +0.00;
    A[3].at(1,0) = +0.00; A[3].at(1,1) = +0.64; A[3].at(1,2) = +0.00;
    A[3].at(2,0) = +0.00; A[3].at(2,1) = +0.00; A[3].at(2,2) = +0.64;

    A[4].at(0,0) = -0.12; A[4].at(0,1) = +0.00; A[4].at(0,2) = +0.00;
    A[4].at(1,0) = +0.00; A[4].at(1,1) = -0.12; A[4].at(1,2) = +0.00;
    A[4].at(2,0) = +0.00; A[4].at(2,1) = +0.00; A[4].at(2,2) = -0.12;

    A[0].at(0,0) = 0.48*h*b(1,k-1);
    A[0].at(1,0) = 0.48*h*b(2,k-1);
    A[0].at(2,0) = 0.48*h*b(3,k-1);

//    double m1 = 25.0 - 12.0*h*a(1,1,k);
//    double m2 = 25.0 - 12.0*h*a(2,2,k);
//    double m3 = 25.0 - 12.0*h*a(3,3,k);

//    A[1].at(0,0) = (48.0/m1);              A[1].at(0,1) = (12.0*h*a(1,2,k-1))/m1; A[1].at(0,2) = (12.0*h*a(1,3,k-1))/m1;
//    A[1].at(1,0) = (12.0*h*a(2,1,k-1))/m2; A[1].at(1,1) = (48.0/m2);              A[1].at(1,2) = (12.0*h*a(2,3,k-1))/m2;
//    A[1].at(2,0) = (12.0*h*a(3,1,k-1))/m3; A[1].at(2,1) = (12.0*h*a(3,2,k-1))/m3; A[1].at(2,2) = (48.0/m3);

//    A[2].at(0,0) = -36.0/m1; A[2].at(0,1) = +0.00;    A[2].at(0,2) = +0.00;
//    A[2].at(1,0) = +0.00;    A[2].at(1,1) = -36.0/m2; A[2].at(1,2) = +0.00;
//    A[2].at(2,0) = +0.00;    A[2].at(2,1) = +0.00;    A[2].at(2,2) = -36.0/m3;

//    A[3].at(0,0) = +16.0/m1; A[3].at(0,1) = +0.00;    A[3].at(0,2) = +0.00;
//    A[3].at(1,0) = +0.00;    A[3].at(1,1) = +16.0/m2; A[3].at(1,2) = +0.00;
//    A[3].at(2,0) = +0.00;    A[3].at(2,1) = +0.00;    A[3].at(2,2) = +16.0/m3;

//    A[4].at(0,0) = -3.0/m1; A[4].at(0,1) = +0.00;   A[4].at(0,2) = +0.00;
//    A[4].at(1,0) = +0.00;   A[4].at(1,1) = -3.0/m2; A[4].at(1,2) = +0.00;
//    A[4].at(2,0) = +0.00;   A[4].at(2,1) = +0.00;   A[4].at(2,2) = -3.0/m3;

//    A[0].at(0,0) = (12.0*h*b(1,k-1))/m1;
//    A[0].at(1,0) = (12.0*h*b(2,k-1))/m2;
//    A[0].at(2,0) = (12.0*h*b(3,k-1))/m3;
}

void Example4::clearAMatrices(std::vector<DoubleMatrix> &A)
{
    A[1].clear();
    A[2].clear();
    A[3].clear();
    A[4].clear();
    A[0].clear();
    A.clear();
}

void Example4::calculateNX(const std::vector<DoubleVector> &rx, DoubleVector &x1, DoubleVector &x2, DoubleVector &x3, std::vector<DoubleVector> &nx)
{
    nx.clear();
    nx.resize(N+1);

    nx.at(0) = rx.at(0);
    nx.at(1) = rx.at(1);
    nx.at(2) = rx.at(2);
    nx.at(3) = rx.at(3);

    for (unsigned int i=0; i<4; i++)
    {
        x1 << nx.at(i).at(0);
        x2 << nx.at(i).at(1);
        x3 << nx.at(i).at(2);
    }

    std::vector<DoubleMatrix> A;
    initAMatrices(A);
    for (unsigned int k=K; k<=N; k++)
    {
        updateAMatrices(A,k);
        nx.at(k) = A[1]*nx.at(k-1) + A[2]*nx.at(k-2) + A[3]*nx.at(k-3) + A[4]*nx.at(k-4) + A[0];

        x1 << nx.at(k).at(0);
        x2 << nx.at(k).at(1);
        x3 << nx.at(k).at(2);
    }
    clearAMatrices(A);
}

void Example4::calculateM1(const unsigned int s[][4], const DoubleMatrix &rx UNUSED_PARAM, const DoubleMatrix &nx)
{
    DoubleMatrix M(K*n, K*n, 0.0);
    DoubleVector B(K*n,0.0);

    calculateM1BE(0,s[0],L,nx,M,B);
    calculateM1BE(1,s[1],L,nx,M,B);
    calculateM1BE(2,s[2],L,nx,M,B);
    calculateM1BE(3,s[3],L,nx,M,B);

    printf("det: %14.10f\n",M.determinant1());

    DoubleVector x(n*K);
    GaussianElimination(M,B,x);

    IPrinter::printSeperatorLine(NULL,'-',stdout);
    IPrinter::print(x,x.size(),w,p,stdout);
    IPrinter::printSeperatorLine(NULL,'-',stdout);

    DoubleMatrix rx1(n,K);
    rx1.at(0,0) = x[0]; rx1.at(0,1) = x[3]; rx1.at(0,2) = x[6]; rx1.at(0,3) = x[9];
    rx1.at(1,0) = x[1]; rx1.at(1,1) = x[4]; rx1.at(1,2) = x[7]; rx1.at(1,3) = x[10];
    rx1.at(2,0) = x[2]; rx1.at(2,1) = x[5]; rx1.at(2,2) = x[8]; rx1.at(2,3) = x[11];

    DoubleMatrix nx1;
    calculateNX(rx1,nx1);
    IPrinter::printVector(w,p,nx1.row(0),"x1: ");
    IPrinter::printVector(w,p,nx1.row(1),"x2: ");
    IPrinter::printVector(w,p,nx1.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
}

void Example4::calculateM2(const unsigned int s[][4], const std::vector<DoubleVector> &rx UNUSED_PARAM)
{
    /* real solution vectors */
    std::vector<DoubleMatrix> P3(N+1);
    std::vector<DoubleMatrix> P2(N+1);
    std::vector<DoubleMatrix> P1(N+1);
    std::vector<DoubleMatrix> P0(N+1);
    std::vector<DoubleVector> Q(N+1);
    calculatePQ(P3,P2,P1,P0,Q);

    /* numerical solution vectors */
    std::vector<DoubleVector> nx(N+1);
    calculateNS(nx,P3,P2,P1,P0,Q);

    /* find x0, x1, x2, x3 */

    DoubleMatrix M(K*n, K*n,0.0);
    DoubleVector B(K*n);

    calculateM2BE(0,s[0],L,nx,M,B,P3,P2,P1,P0,Q);
    calculateM2BE(1,s[1],L,nx,M,B,P3,P2,P1,P0,Q);
    calculateM2BE(2,s[2],L,nx,M,B,P3,P2,P1,P0,Q);
    calculateM2BE(3,s[3],L,nx,M,B,P3,P2,P1,P0,Q);

    DoubleVector x(n*K);
    GaussianElimination(M, B, x);

    puts("---------------------------------------------------");
    printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]);
    IPrinter::print(x,10,12,8);
    printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", fx1(3), fx2(3), fx3(3), fx1(2), fx2(2), fx3(2), fx1(1), fx2(1), fx3(1), fx1(0), fx2(0), fx3(0));
    puts("---------------------------------------------------");

    {
        std::vector<DoubleVector> nx1(N+1);
        nx1.at(0).resize(n);
        nx1.at(1).resize(n);
        nx1.at(2).resize(n);
        nx1.at(3).resize(n);

        nx1.at(3).at(0) = x[0]; nx1.at(2).at(0) = x[3]; nx1.at(1).at(0) = x[6]; nx1.at(0).at(0) = x[9];
        nx1.at(3).at(1) = x[1]; nx1.at(2).at(1) = x[4]; nx1.at(1).at(1) = x[7]; nx1.at(0).at(1) = x[10];
        nx1.at(3).at(2) = x[2]; nx1.at(2).at(2) = x[5]; nx1.at(1).at(2) = x[8]; nx1.at(0).at(2) = x[11];

        for (unsigned int k=K; k<=N; k++)
        {
            nx1.at(k) = P3[k]*nx1.at(3) + P2[k]*nx1.at(2) + P1[k]*nx1.at(1) + P0[k]*nx1.at(0) + Q[k];
        }

        DoubleVector x1;
        DoubleVector x2;
        DoubleVector x3;
        for (unsigned int i=0; i<=N; i++)
        {
            DoubleVector &px = nx1.at(i);
            x1 << px.at(0);
            x2 << px.at(1);
            x3 << px.at(2);
        }
        IPrinter::printVector(w,p,x1); x1.clear();
        IPrinter::printVector(w,p,x2); x2.clear();
        IPrinter::printVector(w,p,x3); x3.clear();
        puts("--------------------------------------------------------------------------------------");
    }

    for (unsigned int i=0; i<=N; i++)
    {
        P3[i].clear();
        P2[i].clear();
        P1[i].clear();
        P0[i].clear();
        Q[i].clear();
    }
    P3.clear();
    P2.clear();
    P1.clear();
    P0.clear();
    Q.clear();
}

void Example4::calculateM1BE(unsigned int c, const unsigned int s[], unsigned int L, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B)
{
    std::vector<DoubleMatrix> betta(N+1);
    DoubleVector eta(n,0.0);

    std::vector<DoubleMatrix> GAMMA(L);
    for (unsigned int i=0; i<L; i++)
    {
        GAMMA[i].resize(n,n,1.0);
//        GAMMA[i].randomData();

        for (unsigned int i1=0; i1<n; i1++)
        {
            for (unsigned int i2=0; i2<n; i2++)
            {
                GAMMA[i].at(i1,i2) = (double)(rand() % RAND_MAX) * 300.0;
            }
        }
        IPrinter::print(GAMMA[i],n,n,14,6);
        puts("---");

        eta = eta + DoubleVector(GAMMA[i]*nx.col(s[i]));
    }

    printf("%14.10f %14.10f %14.10f\n", eta.at(0), eta.at(1), eta.at(2));

    std::vector<DoubleMatrix> A;
    initAMatrices(A);
    for (unsigned int k=N; k>=K; k--)
    {
        updateAMatrices(A,k);

        if (k==N)
        {
            betta[N] = GAMMA[L-1];
            betta[N-1].resize(n,n,0.0);
            betta[N-2].resize(n,n,0.0);
            betta[N-3].resize(n,n,0.0);
        }

        betta[k-1] = betta[k]*A[1] + betta[k-1];
        betta[k-2] = betta[k]*A[2] + betta[k-2];
        betta[k-3] = betta[k]*A[3] + betta[k-3];
        betta[k-4] = betta[k]*A[4];
        eta        = eta - betta[k]*A[0];

        for (unsigned int i=0; i<L; i++)
        {
            if (k==(s[i]+K))
            {
                betta[k-4] = betta[k-4] + GAMMA[i];
            }
        }
    }
    clearAMatrices(A);

    for (unsigned int i=0; i<n; i++)
    {
        for (unsigned int j=0; j<n; j++)
        {
            M[c*n+i][0*n+j] = betta[0][i][j];
            M[c*n+i][1*n+j] = betta[1][i][j];
            M[c*n+i][2*n+j] = betta[2][i][j];
            M[c*n+i][3*n+j] = betta[3][i][j];
        }
        B.at(c*n+i) = eta.at(i);
    }

    eta.clear();
    for (unsigned int i=0; i<L; i++) GAMMA[i].clear();
    GAMMA.clear();
    for (unsigned int i=0; i<=N; i++) betta[i].clear();
    betta.clear();
}

void Example4::calculateM2BE(unsigned int c, const unsigned int s[], unsigned int L, const std::vector<DoubleVector> &nx, DoubleMatrix &M, DoubleVector &B,
                             const std::vector<DoubleMatrix> &P3, const std::vector<DoubleMatrix> &P2, const std::vector<DoubleMatrix> &P1, const std::vector<DoubleMatrix> &P0,
                             const std::vector<DoubleVector> &Q)
{
    std::vector<DoubleMatrix> GAMMA(L);
    for (unsigned int i=0; i<L; i++)
    {
        GAMMA[i].resize(n,n,0.0);
        GAMMA[i].randomData();
    }

    DoubleMatrix U3(n,n,0.0);
    DoubleMatrix U2(n,n,0.0);
    DoubleMatrix U1(n,n,0.0);
    DoubleMatrix U0(n,n,0.0);
    DoubleMatrix V0(n,1,0.0);

    for (unsigned int i=0; i<L; i++)
    {
        if (s[i] == 3)
            U3 = U3 + GAMMA[i];
        else if (s[i] == 2)
            U2 = U2 + GAMMA[i];
        else if (s[i] == 1)
            U1 = U1 + GAMMA[i];
        else if (s[i] == 0)
            U0 = U0 + GAMMA[i];
        else
        {
            U3 = U3 + GAMMA[i]*P3[s[i]];
            U2 = U2 + GAMMA[i]*P2[s[i]];
            U1 = U1 + GAMMA[i]*P1[s[i]];
            U0 = U0 + GAMMA[i]*P0[s[i]];
            V0 = V0 + GAMMA[i]*Q[s[i]];
        }
    }

    DoubleVector B1(n,0.0);
    for (unsigned int i=0; i<L; i++)
    {
        B1 = B1 + DoubleVector(GAMMA[i]*nx.at(s[i]));
    }
    B1 = B1 - V0;

    for (unsigned int i=0; i<n; i++)
    {
        for (unsigned int j=0; j<n; j++)
        {
            M[c*n+i][0*n+j] = U3[i][j];
            M[c*n+i][1*n+j] = U2[i][j];
            M[c*n+i][2*n+j] = U1[i][j];
            M[c*n+i][3*n+j] = U0[i][j];
        }
        B.at(c*n+i) = B1.at(i);
    }
}

void Example4::calculateRX(DoubleMatrix &rx)
{
    if (rx.empty())
    {
        rx.resize(n,N+1,0.0);
        for (unsigned int i=0; i<=N; i++)
        {
            rx.at(0,i) = fx1(i);
            rx.at(1,i) = fx2(i);
            rx.at(2,i) = fx3(i);
        }
    }
}

void Example4::calculateNX(const DoubleMatrix &rx, DoubleMatrix &nx)
{
    if (nx.empty())
    {
        nx.resize(n,N+1,0.0);

        nx.setColumn(0, rx.col(0));
        nx.setColumn(1, rx.col(1));
        nx.setColumn(2, rx.col(2));
        nx.setColumn(3, rx.col(3));

        std::vector<DoubleMatrix> A;
        initAMatrices(A);
        for (unsigned int k=K; k<=N; k++)
        {
            updateAMatrices(A,k);
            DoubleVector ck = A[1]*nx.col(k-1) + A[2]*nx.col(k-2) + A[3]*nx.col(k-3) + A[4]*nx.col(k-4) + A[0];
            nx.setColumn(k,ck);
        }
        clearAMatrices(A);
    }
}


void Example4::calculateNS(std::vector<DoubleVector> &nx, const std::vector<DoubleMatrix> &P3, const std::vector<DoubleMatrix> &P2,
                           const std::vector<DoubleMatrix> &P1, const std::vector<DoubleMatrix> &P0, const std::vector<DoubleVector> &Q)
{
    nx.at(0).resize(n);
    nx.at(1).resize(n);
    nx.at(2).resize(n);
    nx.at(3).resize(n);

    nx.at(0).at(0) = fx1(0); nx.at(1).at(0) = fx1(1); nx.at(2).at(0) = fx1(2); nx.at(3).at(0) = fx1(3);
    nx.at(0).at(1) = fx2(0); nx.at(1).at(1) = fx2(1); nx.at(2).at(1) = fx2(2); nx.at(3).at(1) = fx2(3);
    nx.at(0).at(2) = fx3(0); nx.at(1).at(2) = fx3(1); nx.at(2).at(2) = fx3(2); nx.at(3).at(2) = fx3(3);
    for (unsigned int k=K; k<=N; k++)
    {
        nx.at(k) = P3[k]*nx.at(3) + P2[k]*nx.at(2) + P1[k]*nx.at(1) + P0[k]*nx.at(0) + Q[k];
    }

    DoubleVector x1;
    DoubleVector x2;
    DoubleVector x3;
    for (unsigned int i=0; i<=N; i++)
    {
        DoubleVector &px = nx.at(i);
        x1 << px.at(0);
        x2 << px.at(1);
        x3 << px.at(2);
    }
    IPrinter::printVector(w,p,x1); x1.clear();
    IPrinter::printVector(w,p,x2); x2.clear();
    IPrinter::printVector(w,p,x3); x3.clear();
    puts("--------------------------------------------------------------------------------------");
}

void Example4::calculatePQ(std::vector<DoubleMatrix> &P3, std::vector<DoubleMatrix> &P2, std::vector<DoubleMatrix> &P1, std::vector<DoubleMatrix> &P0, std::vector<DoubleVector> &Q)
{
    /* calculating P,Q matrices */
    std::vector<DoubleMatrix> A;
    initAMatrices(A);
    for (unsigned int k=K; k<=N; k++)
    {
        updateAMatrices(A,k);

        if (k==K)
        {
            P3[k] = A[1];
            P2[k] = A[2];
            P1[k] = A[3];
            P0[k] = A[4];
            Q[k]  = A[0];
        }
        else if (k==K+1)
        {
            P3[k] = A[1]*P3[k-1] + A[2];
            P2[k] = A[1]*P2[k-1] + A[3];
            P1[k] = A[1]*P1[k-1] + A[4];
            P0[k] = A[1]*P0[k-1];
            Q[k]  = A[1]*Q[k-1] + A[0];
        }
        else if (k==K+2)
        {
            P3[k] = A[1]*P3[k-1] + A[2]*P3[k-2] + A[3];
            P2[k] = A[1]*P2[k-1] + A[2]*P2[k-2] + A[4];
            P1[k] = A[1]*P1[k-1] + A[2]*P1[k-2];
            P0[k] = A[1]*P0[k-1] + A[2]*P0[k-2];
            Q[k]  = A[1]*DoubleMatrix(Q[k-1]) + A[2]*DoubleMatrix(Q[k-2]) + A[0];
        }
        else if (k==K+3)
        {
            P3[k] = A[1]*P3[k-1] + A[2]*P3[k-2] + A[3]*P3[k-3] + A[4];
            P2[k] = A[1]*P2[k-1] + A[2]*P2[k-2] + A[3]*P2[k-3];
            P1[k] = A[1]*P1[k-1] + A[2]*P1[k-2] + A[3]*P1[k-3];
            P0[k] = A[1]*P0[k-1] + A[2]*P0[k-2] + A[3]*P0[k-3];
            Q[k]  = A[1]*DoubleMatrix(Q[k-1]) + A[2]*DoubleMatrix(Q[k-2]) + A[3]*DoubleMatrix(Q[k-3]) + A[0];
        }
        if (k>=2*K)
        {
            P3[k] = A[1]*P3[k-1] + A[2]*P3[k-2] + A[3]*P3[k-3] + A[4]*P3[k-4];
            P2[k] = A[1]*P2[k-1] + A[2]*P2[k-2] + A[3]*P2[k-3] + A[4]*P2[k-4];
            P1[k] = A[1]*P1[k-1] + A[2]*P1[k-2] + A[3]*P1[k-3] + A[4]*P1[k-4];
            P0[k] = A[1]*P0[k-1] + A[2]*P0[k-2] + A[3]*P0[k-3] + A[4]*P0[k-4];
            Q[k]  = A[1]*DoubleMatrix(Q[k-1]) + A[2]*DoubleMatrix(Q[k-2]) + A[3]*DoubleMatrix(Q[k-3]) + A[4]*DoubleMatrix(Q[k-4]) + A[0];
        }
    }
    A[4].clear();
    A[3].clear();
    A[2].clear();
    A[1].clear();
    A[0].clear();
    A.clear();
    /* calculating P,Q matrices */
}

double Example4::fx1(unsigned int k) const
{
    double t = k*h;
#ifdef SAMPLE_1
    return sin(2.0*t) + t*t;
#endif
#ifdef SAMPLE_2
    return t*t+t;
#endif
}

double Example4::fx2(unsigned int k) const
{
    double t = k*h;
#ifdef SAMPLE_1
    return 3.0*t;
#endif
#ifdef SAMPLE_2
    return 2.0*t;
#endif
}

double Example4::fx3(unsigned int k) const
{
    double t = k*h;
#ifdef SAMPLE_1
    return cos(2.0*t) - sin(t);
#endif
#ifdef SAMPLE_2
    return 3.0*t*t;
#endif
}

double Example4::a(unsigned int i, unsigned int j, unsigned int k) const
{
    double t = k*h;
    C_UNUSED(t);

#ifdef SAMPLE_1
    if (i==1 && j==1) return +3.0;
    if (i==1 && j==2) return -t;
    if (i==1 && j==3) return +2.0;

    if (i==2 && j==1) return +1.0;
    if (i==2 && j==2) return +1.0;
    if (i==2 && j==3) return +3.0;

    if (i==3 && j==1) return -2.0;
    if (i==3 && j==2) return +t;
    if (i==3 && j==3) return +1.0;
#endif

#ifdef SAMPLE_2
    if (i==1 && j==1) return -3.0;
    if (i==1 && j==2) return +1.0;
    if (i==1 && j==3) return +1.0;

    if (i==2 && j==1) return +3.0;
    if (i==2 && j==2) return -1.5;
    if (i==2 && j==3) return -1.0;

    if (i==3 && j==1) return +6.0;
    if (i==3 && j==2) return +3.0;
    if (i==3 && j==3) return -2.0;
#endif

    return 0.0;
}

double Example4::b(unsigned int i, unsigned int k) const
{
    double t = k*h;
    C_UNUSED(t);

#ifdef SAMPLE_1
    if (i==1) return 2.0*t + 2.0*sin(t) - 3.0*sin(2.0*t);
    if (i==2) return 3.0*sin(t) - sin(2.0*t) - 3.0*cos(2.0*t) - t*t - 3.0*t + 3.0;
    if (i==3) return sin(t) - cos(t) - cos(2.0*t) - t*t;
#endif
#ifdef SAMPLE_2
    if (i==1) return 3.0*t + 1.0;
    if (i==2) return 2.0;
    if (i==3) return -6.0*t;
#endif
    return 0.0;
}
