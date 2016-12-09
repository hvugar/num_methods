#include "example4.h"

double my_rand1()
{
    return ((rand() % 1000) + 1) / 1000.0;
}

void Example4::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Example4 e;
    //e.calculate1();
    puts("Method #1");
    e.calculateM1();
    puts("Method #2");
    e.calculateM2();
}

Example4::Example4()
{
    h = 0.01;
    N = 100;
    //F = N/100;
    n = 3;
    K = 4;
    w = 12;
    p = 8;
}

void Example4::calculateM1()
{
    /* real solution vectors */
    std::vector<DoubleVector> rx(N+1);
    calculateRS(rx);

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

    unsigned int L=4;

    puts("Condition #1");
    unsigned int s0[] = {0,10,80,100};
    calculateM1BE(0,s0,L,nx,M,B,P3,P2,P1,P0,Q);

    puts("Condition #2");
    unsigned int s1[] = {0,20,90,100};
    calculateM1BE(1,s1,L,nx,M,B,P3,P2,P1,P0,Q);

    puts("Condition #3");
    unsigned int s2[] = {0,40,50,100};
    calculateM1BE(2,s2,L,nx,M,B,P3,P2,P1,P0,Q);

    puts("Condition #3");
    unsigned int s3[] = {0,15,75,100};
    calculateM1BE(3,s3,L,nx,M,B,P3,P2,P1,P0,Q);

    DoubleVector x(n*K);
    GaussianElimination(M, B, x);

    puts("---------------------------------------------------");
    printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]);
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

void Example4::calculateM2()
{
    /* real solution vectors */
    std::vector<DoubleVector> rx(N+1);
    calculateRS(rx);

    std::vector<DoubleVector> nx(N+1);
    {
        std::vector<DoubleMatrix> A(5);
        A[0].resize(n,1);
        A[1].resize(n,n);
        A[2].resize(n,n);
        A[3].resize(n,n);
        A[4].resize(n,n);

        DoubleVector tx01;
        DoubleVector tx02;
        DoubleVector tx03;

        nx.at(0).resize(n);
        nx.at(1).resize(n);
        nx.at(2).resize(n);
        nx.at(3).resize(n);

        nx.at(0).at(0) = fx1(0); nx.at(1).at(0) = fx1(1); nx.at(2).at(0) = fx1(2); nx.at(3).at(0) = fx1(3);
        nx.at(0).at(1) = fx2(0); nx.at(1).at(1) = fx2(1); nx.at(2).at(1) = fx2(2); nx.at(3).at(1) = fx2(3);
        nx.at(0).at(2) = fx3(0); nx.at(1).at(2) = fx3(1); nx.at(2).at(2) = fx3(2); nx.at(3).at(2) = fx3(3);

        tx01 << nx.at(0).at(0); tx01 << nx.at(1).at(0); tx01 << nx.at(2).at(0); tx01 << nx.at(3).at(0);
        tx02 << nx.at(0).at(1); tx02 << nx.at(1).at(1); tx02 << nx.at(2).at(1); tx02 << nx.at(3).at(1);
        tx03 << nx.at(0).at(2); tx03 << nx.at(1).at(2); tx03 << nx.at(2).at(2); tx03 << nx.at(3).at(2);

        for (unsigned int k=K; k<=N; k++)
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

            nx.at(k) = A[1]*nx.at(k-1) + A[2]*nx.at(k-2) + A[3]*nx.at(k-3) + A[4]*nx.at(k-4) + A[0];

            tx01 << nx.at(k).at(0);
            tx02 << nx.at(k).at(1);
            tx03 << nx.at(k).at(2);
        }

        IPrinter::printVector(w,p,tx01); tx01.clear();
        IPrinter::printVector(w,p,tx02); tx02.clear();
        IPrinter::printVector(w,p,tx03); tx03.clear();
        puts("---------------------------------------------------------------------------------");

        A[1].clear();
        A[2].clear();
        A[3].clear();
        A[4].clear();
        A[0].clear();
        A.clear();
    }

    DoubleMatrix M(K*n, K*n, 0.0);
    DoubleVector B(K*n,0.0);

    unsigned int L=4;

    puts("Condition #1");
    unsigned int s0[] = {0,10,80,100};
    calculateM2BE(0,s0,L,nx,M,B);

    puts("Condition #2");
    unsigned int s1[] = {0,20,90,100};
    calculateM2BE(1,s1,L,nx,M,B);

    puts("Condition #3");
    unsigned int s2[] = {0,40,50,100};
    calculateM2BE(2,s2,L,nx,M,B);

    puts("Condition #3");
    unsigned int s3[] = {0,15,75,100};
    calculateM2BE(3,s3,L,nx,M,B);

    DoubleVector x(n*K);
    GaussianElimination(M,B,x);

    puts("---------------------------------------------------");
    printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]);
    printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", fx1(3), fx2(3), fx3(3), fx1(2), fx2(2), fx3(2), fx1(1), fx2(1), fx3(1), fx1(0), fx2(0), fx3(0));
    puts("---------------------------------------------------");

    {
        std::vector<DoubleVector> nx1(N+1);
        std::vector<DoubleMatrix> A(5);
        A[0].resize(n,1);
        A[1].resize(n,n);
        A[2].resize(n,n);
        A[3].resize(n,n);
        A[4].resize(n,n);

        DoubleVector tx01;
        DoubleVector tx02;
        DoubleVector tx03;

        nx1.at(0).resize(n);
        nx1.at(1).resize(n);
        nx1.at(2).resize(n);
        nx1.at(3).resize(n);

        nx1.at(3).at(0) = x[0]; nx1.at(2).at(0) = x[3]; nx1.at(1).at(0) = x[6]; nx1.at(0).at(0) = x[9];
        nx1.at(3).at(1) = x[1]; nx1.at(2).at(1) = x[4]; nx1.at(1).at(1) = x[7]; nx1.at(0).at(1) = x[10];
        nx1.at(3).at(2) = x[2]; nx1.at(2).at(2) = x[5]; nx1.at(1).at(2) = x[8]; nx1.at(0).at(2) = x[11];

        tx01 << nx1.at(0).at(0); tx01 << nx1.at(1).at(0); tx01 << nx1.at(2).at(0); tx01 << nx1.at(3).at(0);
        tx02 << nx1.at(0).at(1); tx02 << nx1.at(1).at(1); tx02 << nx1.at(2).at(1); tx02 << nx1.at(3).at(1);
        tx03 << nx1.at(0).at(2); tx03 << nx1.at(1).at(2); tx03 << nx1.at(2).at(2); tx03 << nx1.at(3).at(2);

        for (unsigned int k=K; k<=N; k++)
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

            nx1.at(k) = A[1]*nx1.at(k-1) + A[2]*nx1.at(k-2) + A[3]*nx1.at(k-3) + A[4]*nx1.at(k-4) + A[0];

            tx01 << nx1.at(k).at(0);
            tx02 << nx1.at(k).at(1);
            tx03 << nx1.at(k).at(2);
        }

        IPrinter::printVector(w,p,tx01); tx01.clear();
        IPrinter::printVector(w,p,tx02); tx02.clear();
        IPrinter::printVector(w,p,tx03); tx03.clear();
        puts("---------------------------------------------------------------------------------");

        A[1].clear();
        A[2].clear();
        A[3].clear();
        A[4].clear();
        A[0].clear();
        A.clear();
    }
}

/**
 * @brief Example4::calculateM2BE
 * @param c Condition number
 * @param s Indexes
 * @param L Count of indexes
 */
void Example4::calculateM2BE(unsigned int c, unsigned int s[], unsigned int L, const std::vector<DoubleVector> &rx, DoubleMatrix &M, DoubleVector &B)
{
    {
        std::vector<DoubleMatrix> betta(N+1);
        DoubleVector eta(n,0.0);

        std::vector<DoubleMatrix> GAMMA(L);
        for (unsigned int i=0; i<L; i++)
        {
            GAMMA[i].resize(n,n,0.0);
            GAMMA[i].randomData();

            eta = eta + DoubleVector(GAMMA[i]*rx.at(s[i]));
        }

        std::vector<DoubleMatrix> A(5);
        A[0].resize(n,1);
        A[1].resize(n,n);
        A[2].resize(n,n);
        A[3].resize(n,n);
        A[4].resize(n,n);

        for (unsigned int k=N; k>=K; k--)
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

            if (k==N)
            {
                betta[N] = GAMMA[3];
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

        for (unsigned int i=0; i<n; i++)
        {
            for (unsigned int j=0; j<n; j++)
            {
                M[c*n+i][0*n+j] = betta[3][i][j];
                M[c*n+i][1*n+j] = betta[2][i][j];
                M[c*n+i][2*n+j] = betta[1][i][j];
                M[c*n+i][3*n+j] = betta[0][i][j];
            }
            B.at(c*n+i) = eta.at(i);
        }

        A[4].clear();
        A[3].clear();
        A[2].clear();
        A[1].clear();
        A[0].clear();
        A.clear();

        for (unsigned int i=0; i<L; i++)
        {
            GAMMA[i].clear();
        }
        GAMMA.clear();
        eta.clear();

        for (unsigned int i=0; i<=N; i++) betta.clear();
        betta.clear();
    }
}

void Example4::calculateRS(std::vector<DoubleVector> &rx)
{
    DoubleVector tx01;
    DoubleVector tx02;
    DoubleVector tx03;
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

    IPrinter::printVector(w,p,tx01); tx01.clear();
    IPrinter::printVector(w,p,tx02); tx02.clear();
    IPrinter::printVector(w,p,tx03); tx03.clear();
    puts("---------------------------------------------------------------------------------");
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
    std::vector<DoubleMatrix> A(5);
    A[0].resize(n,1);
    A[1].resize(n,n);
    A[2].resize(n,n);
    A[3].resize(n,n);
    A[4].resize(n,n);
    for (unsigned int k=K; k<=N; k++)
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

void Example4::calculateM1BE(unsigned int c, unsigned int s[], unsigned int L, const std::vector<DoubleVector> &nx, DoubleMatrix &M, DoubleVector &B,
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
