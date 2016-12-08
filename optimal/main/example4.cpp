#include "example4.h"

double my_rand1()
{
    return ((rand() % 1000) + 1) / 1000.0;
}

void Example4::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Example4 e;
    //e.calculate1();
    e.calculateM1();
}

Example4::Example4()
{
    h = 0.001;
    N = 1000;
    n = 3;
    K = 4;
    w = 18;
    p = 14;
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

void Example4::calculateXK(DoubleVector &x, const std::vector<DoubleVector> &rx,
                           const std::vector<DoubleMatrix> &P3, const std::vector<DoubleMatrix> &P2,
                           const std::vector<DoubleMatrix> &P1, const std::vector<DoubleMatrix> &P0,
                           const std::vector<DoubleVector> &Q UNUSED_PARAM)
{
    DoubleMatrix M(K*n, K*n,0.0);
    DoubleVector b(K*n);

    //    DoubleMatrix G000(n,n,0.0); G000.randomData();
    //    DoubleMatrix G001(n,n,0.0); G001.randomData();
    //    DoubleMatrix G002(n,n,0.0); G002.randomData();
    //    DoubleMatrix G003(n,n,0.0); G003.randomData();
    //    DoubleMatrix G100(n,n,0.0); G100.randomData();
    //    DoubleMatrix G200(n,n,0.0); G200.randomData();
    //    DoubleMatrix G250(n,n,0.0); G250.randomData();
    //    DoubleMatrix G300(n,n,0.0); G300.randomData();
    //    DoubleMatrix G450(n,n,0.0); G450.randomData();
    //    DoubleMatrix G800(n,n,0.0); G800.randomData();
    //    DoubleMatrix G900(n,n,0.0); G900.randomData();
    //    DoubleMatrix G00N(n,n,0.0); G00N.randomData();

    //M.randomData();
    {
        DoubleMatrix G000(n,n,0.0); G000.randomData();
        DoubleMatrix G001(n,n,0.0); G001.randomData();
        DoubleMatrix G002(n,n,0.0); G002.randomData();
        DoubleMatrix G003(n,n,0.0); G003.randomData();
        DoubleMatrix G100(n,n,0.0); G100.randomData();
        DoubleMatrix G900(n,n,0.0); G900.randomData();
        DoubleMatrix G00N(n,n,0.0); G00N.randomData();

        DoubleMatrix U3 = G100*P3[100] + G900*P3[900] + G00N*P3[N];// + G003;
        DoubleMatrix U2 = G100*P2[100] + G900*P2[900] + G00N*P2[N] + G002;
        DoubleMatrix U1 = G100*P1[100] + G900*P1[900] + G00N*P1[N] + G001;
        DoubleMatrix U0 = G100*P0[100] + G900*P0[900] + G00N*P0[N] + G000;
        DoubleMatrix V0 = G100*Q[100]  + G900*Q[900]  + G00N*Q[N];

        DoubleVector B1 = G000*rx.at(0) + G001*rx.at(1) + G002*rx.at(2) + G100*rx.at(100) + G900*rx.at(900) + G00N*rx.at(N) - V0;
        DoubleVector B2 = U3*rx.at(3) + U2*rx.at(2) + U1*rx.at(1) + U0*rx.at(0);
        printf("*** %14.10f %14.10f %14.10f\n", B1.at(0), B1.at(1), B1.at(2));
        printf("*** %14.10f %14.10f %14.10f\n", B2.at(0), B2.at(1), B2.at(2));

        for (unsigned int i=0; i<n; i++)
        {
            for (unsigned int j=0; j<n; j++)
            {
                M[j][i+0*n] = U3[j][i];
                M[j][i+1*n] = U2[j][i];
                M[j][i+2*n] = U1[j][i];
                M[j][i+3*n] = U0[j][i];
            }
            b[i] = B2[i];
        }
    }

    {
        DoubleMatrix G000(n,n,0.0); G000.randomData();
        DoubleMatrix G001(n,n,0.0); G001.randomData();
        DoubleMatrix G200(n,n,0.0); G200.randomData();
        DoubleMatrix G300(n,n,0.0); G300.randomData();
        DoubleMatrix G00N(n,n,0.0); G00N.randomData();

        DoubleMatrix U3 = G200*P3[200] + G300*P3[300] + G00N*P3[N];
        DoubleMatrix U2 = G200*P2[200] + G300*P2[300] + G00N*P2[N];
        DoubleMatrix U1 = G200*P1[200] + G300*P1[300] + G00N*P1[N] + G001;
        DoubleMatrix U0 = G200*P0[200] + G300*P0[300] + G00N*P0[N] + G000;
        DoubleMatrix V0 = G200*Q[200]  + G300*Q[300]  + G00N*Q[N];

        DoubleVector B1 = G000*rx.at(0) + G001*rx.at(1) + G200*rx.at(200) + G300*rx.at(300) + G00N*rx.at(N) - V0;
        DoubleVector B2 = U3*rx.at(3) + U2*rx.at(2) + U1*rx.at(1) + U0*rx.at(0);
        printf("*** %14.10f %14.10f %14.10f\n", B1.at(0), B1.at(1), B1.at(2));
        printf("*** %14.10f %14.10f %14.10f\n", B2.at(0), B2.at(1), B2.at(2));

        for (unsigned int i=0; i<n; i++)
        {
            for (unsigned int j=0; j<n; j++)
            {
                M[n+j][i+0*n] = U3[j][i];
                M[n+j][i+1*n] = U2[j][i];
                M[n+j][i+2*n] = U1[j][i];
                M[n+j][i+3*n] = U0[j][i];
            }
            b[n+i] = B2[i];
        }
    }

    {
        DoubleMatrix G000(n,n,0.0); G000.randomData();
        DoubleMatrix G800(n,n,0.0); G800.randomData();
        DoubleMatrix G900(n,n,0.0); G900.randomData();
        DoubleMatrix G00N(n,n,0.0); G00N.randomData();

        DoubleMatrix U3 = G800*P3[800] + G900*P3[900] + G00N*P3[N];
        DoubleMatrix U2 = G800*P2[800] + G900*P2[900] + G00N*P2[N];
        DoubleMatrix U1 = G800*P1[800] + G900*P1[900] + G00N*P1[N];
        DoubleMatrix U0 = G800*P0[800] + G900*P0[900] + G00N*P0[N] + G000;
        DoubleMatrix V0 = G800*Q[800]  + G900*Q[900]  + G00N*Q[N];

        DoubleVector B1 = G000*rx.at(0) + G800*rx.at(800) + G900*rx.at(900) + G00N*rx.at(N) - V0;
        DoubleVector B2 = U3*rx.at(3) + U2*rx.at(2) + U1*rx.at(1) + U0*rx.at(0);
        printf("*** %14.10f %14.10f %14.10f\n", B1.at(0), B1.at(1), B1.at(2));
        printf("*** %14.10f %14.10f %14.10f\n", B2.at(0), B2.at(1), B2.at(2));

        for (unsigned int i=0; i<n; i++)
        {
            for (unsigned int j=0; j<n; j++)
            {
                M[2*n+j][i+0*n] = U3[j][i];
                M[2*n+j][i+1*n] = U2[j][i];
                M[2*n+j][i+2*n] = U1[j][i];
                M[2*n+j][i+3*n] = U0[j][i];
            }
            b[2*n+i] = B2[i];
        }
    }

    {
        DoubleMatrix G000(n,n,0.0); G000.randomData();
        DoubleMatrix G001(n,n,0.0); G001.randomData();
        DoubleMatrix G250(n,n,0.0); G250.randomData();
        DoubleMatrix G450(n,n,0.0); G450.randomData();
        DoubleMatrix G00N(n,n,0.0); G00N.randomData();

        DoubleMatrix U3 = G250*P3[250] + G450*P3[450] + G00N*P3[N];
        DoubleMatrix U2 = G250*P2[250] + G450*P2[450] + G00N*P2[N];
        DoubleMatrix U1 = G250*P1[250] + G450*P1[450] + G00N*P1[N];// + G001;
        DoubleMatrix U0 = G250*P0[250] + G450*P0[450] + G00N*P0[N];// + G000;
        DoubleMatrix V0 = G250*Q[250]  + G450*Q[450]  + G00N*Q[N];

        DoubleVector B1 = G250*rx.at(250) + G450*rx.at(450) + G00N*rx.at(N) - V0;
        DoubleVector B2 = U3*rx.at(3) + U2*rx.at(2) + U1*rx.at(1) + U0*rx.at(0);
        printf("*** %14.10f %14.10f %14.10f\n", B1.at(0), B1.at(1), B1.at(2));
        printf("*** %14.10f %14.10f %14.10f\n", B2.at(0), B2.at(1), B2.at(2));

        for (unsigned int i=0; i<n; i++)
        {
            for (unsigned int j=0; j<n; j++)
            {
                M[3*n+j][i+0*n] = U3[j][i];
                M[3*n+j][i+1*n] = U2[j][i];
                M[3*n+j][i+2*n] = U1[j][i];
                M[3*n+j][i+3*n] = U0[j][i];
            }
            b[3*n+i] = B2[i];
        }
    }

    //    puts("+++");
    //    IPrinter::print(M,M.rows(),M.cols(),w,p);
    //    puts("+++");

    GaussianElimination(M, b, x);

    //    system("pause");
    //    printf("Determinat: %14.10f\n", M.determinant());
    //    system("pause");

    //    puts("**********************************************************");
    //    for (unsigned int i=0; i<12; i++)
    //    {
    //        for (unsigned int j=0; j<12; j++)
    //        {
    //            DoubleVector v(12);
    //            for (unsigned int k=0; k<12; k++) v.at(k) = M.at(i,k)/M.at(j,k);
    //            IPrinter::printVector(18,14,v);
    //            printf("%d %d**********************************************************\n", i,j);
    //            system("pause");
    //        }
    //    }

    //    DoubleMatrix M1 = M;
    //    DoubleVector b1 = b;
    //    //M1.randomData();
    //    puts("..............................................................");
    //    IPrinter::print(M1,12,12,w,p);
    //    for (unsigned int k=0; k<11; k++)
    //    {
    //        for (unsigned int j=k+1; j<12; j++)
    //        {
    //            printf("%d %d\n", k,j);
    //            double c = -M1.at(j,k)/M1.at(k,k);
    //            for (unsigned int i=k; i<12; i++) M1.at(j, i) = M1.at(j, i) + M1.at(k,i)*c;
    //            b1[j] = b1[j] + b1[k] *c;
    //            puts("..............................................................");
    //            system("pause");
    //            IPrinter::print(M1,12,12,18,14);
    //        }
    //    }
    //    puts("..............................................................");
    //    printf("%18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f\n", b1[0], b1[1], b1[2], b1[3], b1[4], b1[5], b1[6], b1[7], b1[8], b1[9], b1[10], b1[11]);

    //    xx[11] = (b1[11]) / M1.at(11,11);
    //    xx[10] = (b1[10] - M1.at(10,11)*xx[11]) / M1.at(10,10);
    //    xx[9]  = (b1[9] - M1.at(9,11)*xx[11] - M1.at(9,10)*xx[10]) / M1.at(9,9);
    //    xx[8]  = (b1[8] - /*M1.at(8,11)*xx[11] - M1.at(8,10)*xx[10] -*/ M1.at(8,9)*xx[9]) / M1.at(8,8);
    //    xx[7]  = (b1[7] - M1.at(7,11)*xx[11] - M1.at(7,10)*xx[10] - M1.at(7,9)*xx[9] - M1.at(7,8)*xx[8]) / M1.at(7,7);


    //    for (unsigned int i=11; i!=UINT32_MAX; i--)
    //    {
    //        for (unsigned int j=11; j>i; j--) b1[i] -= (M1.at(i,j) * xx[j]);
    //        xx[i] = b1[i] / M1.at(i,i);
    //        printf("%d %14.10f\n", i, xx[i]);
    //    }

    //    puts("+++");
    //    for (unsigned int i=0; i<M.rows(); i++)
    //    {
    //        for (unsigned int j=0; j<M.cols(); j++)
    //        {
    //            printf("%14.10f ",M[i][j]);
    //        }
    //        printf("%14.10f\n", b[i]);
    //    }
    //    puts("+++");
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
    DoubleVector x(n*K);
    calculateXK(x,rx,P3,P2,P1,P0,Q);

    puts("---------------------------------------------------");
    printf("%18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]);
    printf("%18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f\n", fx1(3), fx2(3), fx3(3), fx1(2), fx2(2), fx3(2), fx1(1), fx2(1), fx3(1), fx1(0), fx2(0), fx3(0));
    puts("---------------------------------------------------");

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

    std::vector<DoubleMatrix> P3(N+1);
    std::vector<DoubleMatrix> P2(N+1);
    std::vector<DoubleMatrix> P1(N+1);
    std::vector<DoubleMatrix> P0(N+1);
    std::vector<DoubleVector> Q(N+1);

    std::vector<DoubleMatrix> A(5);
    A[0].resize(n,1);
    A[1].resize(n,n);
    A[2].resize(n,n);
    A[3].resize(n,n);
    A[4].resize(n,n);

    {
        DoubleMatrix G000(n,n,0.0); G000.randomData();
        DoubleMatrix G001(n,n,0.0); G001.randomData();
        DoubleMatrix G002(n,n,0.0); G002.randomData();
        DoubleMatrix G003(n,n,0.0); G003.randomData();
        DoubleMatrix G100(n,n,0.0); G100.randomData();
        DoubleMatrix G900(n,n,0.0); G900.randomData();
        DoubleMatrix G00N(n,n,0.0); G00N.randomData();

        for (unsigned int k=N; k!=UINT32_MAX; k--)
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
                P3[k] = G00N;
                P2[k].resize(n,n,0.0);
                P1[k].resize(n,n,0.0);
                P0[k].resize(n,n,0.0);
            }
        }
    }
}

void Example4::calculate1()
{
    DoubleVector x01(N+1);
    DoubleVector x02(N+1);
    DoubleVector x03(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x01.at(i) = fx1(i);
        x02.at(i) = fx2(i);
        x03.at(i) = fx3(i);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    IPrinter::printVector(14,10,x03);

    DoubleVector x1(N+1);
    DoubleVector x2(N+2);
    DoubleVector x3(N+2);

    x1.at(0) = fx1(0); x2.at(0) = fx2(0); x3.at(0) = fx3(0);
    x1.at(1) = fx1(1); x2.at(1) = fx2(1); x3.at(1) = fx3(1);
    x1.at(2) = fx1(2); x2.at(2) = fx2(2); x3.at(2) = fx3(2);
    x1.at(3) = fx1(3); x2.at(3) = fx2(3); x3.at(3) = fx3(3);

    puts("---");
    for (unsigned int k=4; k<=N; k++)
    {
        unsigned int k1 = k-1;
        double alpha1 = +1.92;//+48.0/25.0;
        double alpha2 = -1.44;//-36.0/25.0;
        double alpha3 = +0.64;//+16.0/25.0;
        double alpha4 = -0.12;//-3.0/25.0;
        double alpha5 = +0.48*h;//12.0/25.0;

        x1.at(k) = alpha1*x1.at(k-1) + alpha2*x1.at(k-2) + alpha3*x1.at(k-3) + alpha4*x1.at(k-4)
                + (alpha5*a(1,1,k1))*x1.at(k1) + (alpha5*a(1,2,k1))*x2.at(k1) + (alpha5*a(1,3,k1))*x3.at(k1)
                + (alpha5*b(1,k1));
        x2.at(k) = alpha1*x2.at(k-1) + alpha2*x2.at(k-2) + alpha3*x2.at(k-3) + alpha4*x2.at(k-4)
                + (alpha5*a(2,1,k1))*x1.at(k1) + (alpha5*a(2,2,k1))*x2.at(k1) + (alpha5*a(2,3,k1))*x3.at(k1)
                + (alpha5*b(2,k1));
        x3.at(k) = alpha1*x3.at(k-1) + alpha2*x3.at(k-2) + alpha3*x3.at(k-3) + alpha4*x3.at(k-4)
                + (alpha5*a(3,1,k1))*x1.at(k1) + (alpha5*a(3,2,k1))*x2.at(k1) + (alpha5*a(3,3,k1))*x3.at(k1)
                + (alpha5*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
}

void Example4::calculate2()
{
    h = 0.00001;
    N = 100000;

    DoubleVector x01(N+1);
    DoubleVector x02(N+1);
    DoubleVector x03(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x01.at(i) = fx1(i);
        x02.at(i) = fx2(i);
        x03.at(i) = fx3(i);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    IPrinter::printVector(14,10,x03);

    DoubleVector x1(N+1);
    DoubleVector x2(N+2);
    DoubleVector x3(N+2);

    x1.at(N-0) = fx1(N-0); x2.at(N-0) = fx2(N-0); x3.at(N-0) = fx3(N-0);
    x1.at(N-1) = fx1(N-1); x2.at(N-1) = fx2(N-1); x3.at(N-1) = fx3(N-1);
    x1.at(N-2) = fx1(N-2); x2.at(N-2) = fx2(N-2); x3.at(N-2) = fx3(N-2);
    x1.at(N-3) = fx1(N-3); x2.at(N-3) = fx2(N-3); x3.at(N-3) = fx3(N-3);

    //    x1.at(0) = X1(0); x2.at(0) = X2(0); x3.at(0) = X3(0);
    //    x1.at(1) = X1(1); x2.at(1) = X2(1); x3.at(1) = X3(1);
    //    x1.at(2) = X1(2); x2.at(2) = X2(2); x3.at(2) = X3(2);
    //    x1.at(3) = X1(3); x2.at(3) = X2(3); x3.at(3) = X3(3);

    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+1;
        double alpha1 = +1.92;//+48.0/25.0;
        double alpha2 = -1.44;//-36.0/25.0;
        double alpha3 = +0.64;//+16.0/25.0;
        double alpha4 = -0.12;//-3.0/25.0;
        double alpha5 = -0.48*h;//-12.0/25.0;

        x1.at(k) = alpha1*x1.at(k+1) + alpha2*x1.at(k+2) + alpha3*x1.at(k+3) + alpha4*x1.at(k+4)
                + (alpha5*a(1,1,k1))*x1.at(k1) + (alpha5*a(1,2,k1))*x2.at(k1) + (alpha5*a(1,3,k1))*x3.at(k1)
                + (alpha5*b(1,k1));
        x2.at(k) = alpha1*x2.at(k+1) + alpha2*x2.at(k+2) + alpha3*x2.at(k+3) + alpha4*x2.at(k+4)
                + (alpha5*a(2,1,k1))*x1.at(k1) + (alpha5*a(2,2,k1))*x2.at(k1) + (alpha5*a(2,3,k1))*x3.at(k1)
                + (alpha5*b(2,k1));
        x3.at(k) = alpha1*x3.at(k+1) + alpha2*x3.at(k+2) + alpha3*x3.at(k+3) + alpha4*x3.at(k+4)
                + (alpha5*a(3,1,k1))*x1.at(k1) + (alpha5*a(3,2,k1))*x2.at(k1) + (alpha5*a(3,3,k1))*x3.at(k1)
                + (alpha5*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+1;
        x1.at(k) = (48.0/25.0)*x1.at(k+1) + (-36.0/25.0)*x1.at(k+2) + (16.0/25.0)*x1.at(k+3) + (-3.0/25.0)*x1.at(k+4)
                + ((-12.0/25.0)*h*a(1,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(1,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(1,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(1,k1));
        x2.at(k) = (48.0/25.0)*x2.at(k+1) + (-36.0/25.0)*x2.at(k+2) + (16.0/25.0)*x2.at(k+3) + (-3.0/25.0)*x2.at(k+4)
                + ((-12.0/25.0)*h*a(2,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(2,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(2,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(2,k1));
        x3.at(k) = (48.0/25.0)*x3.at(k+1) + (-36.0/25.0)*x3.at(k+2) + (16.0/25.0)*x3.at(k+3) + (-3.0/25.0)*x3.at(k+4)
                + ((-12.0/25.0)*h*a(3,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(3,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(3,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+2;
        x1.at(k) = (48.0/25.0)*x1.at(k+1) + (-36.0/25.0)*x1.at(k+2) + (16.0/25.0)*x1.at(k+3) + (-3.0/25.0)*x1.at(k+4)
                + ((-12.0/25.0)*h*a(1,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(1,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(1,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(1,k1));
        x2.at(k) = (48.0/25.0)*x2.at(k+1) + (-36.0/25.0)*x2.at(k+2) + (16.0/25.0)*x2.at(k+3) + (-3.0/25.0)*x2.at(k+4)
                + ((-12.0/25.0)*h*a(2,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(2,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(2,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(2,k1));
        x3.at(k) = (48.0/25.0)*x3.at(k+1) + (-36.0/25.0)*x3.at(k+2) + (16.0/25.0)*x3.at(k+3) + (-3.0/25.0)*x3.at(k+4)
                + ((-12.0/25.0)*h*a(3,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(3,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(3,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+3;
        x1.at(k) = (48.0/25.0)*x1.at(k+1) + (-36.0/25.0)*x1.at(k+2) + (16.0/25.0)*x1.at(k+3) + (-3.0/25.0)*x1.at(k+4)
                + ((-12.0/25.0)*h*a(1,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(1,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(1,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(1,k1));
        x2.at(k) = (48.0/25.0)*x2.at(k+1) + (-36.0/25.0)*x2.at(k+2) + (16.0/25.0)*x2.at(k+3) + (-3.0/25.0)*x2.at(k+4)
                + ((-12.0/25.0)*h*a(2,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(2,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(2,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(2,k1));
        x3.at(k) = (48.0/25.0)*x3.at(k+1) + (-36.0/25.0)*x3.at(k+2) + (16.0/25.0)*x3.at(k+3) + (-3.0/25.0)*x3.at(k+4)
                + ((-12.0/25.0)*h*a(3,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(3,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*a(3,3,k1))*x3.at(k1)
                + ((-12.0/25.0)*h*b(3,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");
}

void Example4::calculate3()
{
    h = 0.001;
    N = 1000;

    DoubleVector x01(N+1);
    DoubleVector x02(N+1);
    DoubleVector x03(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x01.at(i) = fx1(i);
        x02.at(i) = fx2(i);
        x03.at(i) = fx3(i);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    IPrinter::printVector(14,10,x03);
    puts("---");

    DoubleVector x1(N+1);
    DoubleVector x2(N+1);
    DoubleVector x3(N+1);

    x1.at(0) = fx1(0); x1.at(1) = fx1(1); x1.at(2) = fx1(2); x1.at(3) = fx1(3);
    x2.at(0) = fx2(0); x2.at(1) = fx2(1); x2.at(2) = fx2(2); x2.at(3) = fx2(3);
    x3.at(0) = fx3(0); x3.at(1) = fx3(1); x3.at(2) = fx3(2); x3.at(3) = fx3(3);

    DoubleCube p3(N+1, 3, 3);
    DoubleCube p2(N+1, 3, 3);
    DoubleCube p1(N+1, 3, 3);
    DoubleCube p0(N+1, 3, 3);
    DoubleMatrix q0(N+1, 3);

    for (unsigned int k=4; k<=N; k++)
    {
        if (k==4)
        {
            unsigned int k1 = k-1;

            DoubleMatrix A41(3,3);
            A41.at(0,0) = 0.48*h*a(1,1,k1)+1.92; A41.at(0,1) = 0.48*h*a(1,2,k1);      A41.at(0,2) = 0.48*h*a(1,3,k1);
            A41.at(1,0) = 0.48*h*a(2,1,k1);      A41.at(1,1) = 0.48*h*a(2,2,k1)+1.92; A41.at(1,2) = 0.48*h*a(2,3,k1);
            A41.at(2,0) = 0.48*h*a(3,1,k1);      A41.at(2,1) = 0.48*h*a(3,2,k1);      A41.at(2,2) = 0.48*h*a(3,3,k1)+1.92;

            DoubleMatrix A42(3,3);
            A42.at(0,0) = -1.44; A42.at(0,1) = +0.00; A42.at(0,2) = +0.00;
            A42.at(1,0) = +0.00; A42.at(1,1) = -1.44; A42.at(1,2) = +0.00;
            A42.at(2,0) = +0.00; A42.at(2,1) = +0.00; A42.at(2,2) = -1.44;

            DoubleMatrix A43(3,3);
            A43.at(0,0) = +0.64; A43.at(0,1) = +0.00; A43.at(0,2) = +0.00;
            A43.at(1,0) = +0.00; A43.at(1,1) = +0.64; A43.at(1,2) = +0.00;
            A43.at(2,0) = +0.00; A43.at(2,1) = +0.00; A43.at(2,2) = +0.64;

            DoubleMatrix A44(3,3);
            A44.at(0,0) = -0.12; A44.at(0,1) = +0.00; A44.at(0,2) = +0.00;
            A44.at(1,0) = +0.00; A44.at(1,1) = -0.12; A44.at(1,2) = +0.00;
            A44.at(2,0) = +0.00; A44.at(2,1) = +0.00; A44.at(2,2) = -0.12;

            //            p3.at(k,0,0) = A41(0,0); p3.at(k,0,1) = A41(0,1);  p3.at(k,0,2) = A41(0,2);
            //            p3.at(k,1,0) = A41(1,0); p3.at(k,1,1) = A41(1,1);  p3.at(k,1,2) = A41(1,2);
            //            p3.at(k,2,0) = A41(2,0); p3.at(k,2,1) = A41(2,1);  p3.at(k,2,2) = A41(2,2);

            for (unsigned int i=0; i<3; i++)
            {
                for (unsigned int j=0; j<3; j++)
                {
                    p3.at(k,i,j) = A41(i,j);
                    p2.at(k,i,j) = A42(i,j);
                    p1.at(k,i,j) = A43(i,j);
                    p0.at(k,i,j) = A44(i,j);
                }
            }

            //            p3.at(k,0,0) = 0.48*h*a(1,1,k1)+1.92; p3.at(k,0,1) = 0.48*h*a(1,2,k1);       p3.at(k,0,2) = 0.48*h*a(1,3,k1);
            //            p3.at(k,1,0) = 0.48*h*a(2,1,k1);      p3.at(k,1,1) = 0.48*h*a(2,2,k1)+1.92;  p3.at(k,1,2) = 0.48*h*a(2,3,k1);
            //            p3.at(k,2,0) = 0.48*h*a(3,1,k1);      p3.at(k,2,1) = 0.48*h*a(3,2,k1);       p3.at(k,2,2) = 0.48*h*a(3,3,k1)+1.92;

            //            p2.at(k,0,0) = -1.44; p2.at(k,0,1) = 0.0;   p2.at(k,0,2) = 0.0;
            //            p2.at(k,1,0) = 0.0;   p2.at(k,1,1) = -1.44; p2.at(k,1,2) = 0.0;
            //            p2.at(k,2,0) = 0.0;   p2.at(k,2,1) = 0.0;   p2.at(k,2,2) = -1.44;

            //            p1.at(k,0,0) = 0.64;  p1.at(k,0,1) = 0.0;   p1.at(k,0,2) = 0.0;
            //            p1.at(k,1,0) = 0.0;   p1.at(k,1,1) = 0.64;  p1.at(k,1,2) = 0.0;
            //            p1.at(k,2,0) = 0.0;   p1.at(k,2,1) = 0.0;   p1.at(k,2,2) = 0.64;

            //            p0.at(k,0,0) = -0.12; p0.at(k,0,1) = 0.0;   p0.at(k,0,2) = 0.0;
            //            p0.at(k,1,0) = 0.0;   p0.at(k,1,1) = -0.12; p0.at(k,1,2) = 0.0;
            //            p0.at(k,2,0) = 0.0;   p0.at(k,2,1) = 0.0;   p0.at(k,2,2) = -0.12;

            q0.at(k,0) = 0.48*h*b(1,k1);
            q0.at(k,1) = 0.48*h*b(2,k1);
            q0.at(k,2) = 0.48*h*b(3,k1);
        }

        if (k==5)
        {
            unsigned int k1 = k-1;

            DoubleMatrix A51(3,3);
            A51.at(0,0) = 0.48*h*a(1,1,k1)+1.92; A51.at(0,1) = 0.48*h*a(1,2,k1);      A51.at(0,2) = 0.48*h*a(1,3,k1);
            A51.at(1,0) = 0.48*h*a(2,1,k1);      A51.at(1,1) = 0.48*h*a(2,2,k1)+1.92; A51.at(1,2) = 0.48*h*a(2,3,k1);
            A51.at(2,0) = 0.48*h*a(3,1,k1);      A51.at(2,1) = 0.48*h*a(3,2,k1);      A51.at(2,2) = 0.48*h*a(3,3,k1)+1.92;

            DoubleMatrix A52(3,3);
            A52.at(0,0) = -1.44; A52.at(0,1) = +0.00; A52.at(0,2) = +0.00;
            A52.at(1,0) = +0.00; A52.at(1,1) = -1.44; A52.at(1,2) = +0.00;
            A52.at(2,0) = +0.00; A52.at(2,1) = +0.00; A52.at(2,2) = -1.44;

            DoubleMatrix A53(3,3);
            A53.at(0,0) = +0.64; A53.at(0,1) = +0.00; A53.at(0,2) = +0.00;
            A53.at(1,0) = +0.00; A53.at(1,1) = +0.64; A53.at(1,2) = +0.00;
            A53.at(2,0) = +0.00; A53.at(2,1) = +0.00; A53.at(2,2) = +0.64;

            DoubleMatrix A54(3,3);
            A54.at(0,0) = -0.12; A54.at(0,1) = +0.00; A54.at(0,2) = +0.00;
            A54.at(1,0) = +0.00; A54.at(1,1) = -0.12; A54.at(1,2) = +0.00;
            A54.at(2,0) = +0.00; A54.at(2,1) = +0.00; A54.at(2,2) = -0.12;

            DoubleMatrix P43(3,3);
            DoubleMatrix P42(3,3);
            DoubleMatrix P41(3,3);
            for (unsigned int i=0; i<3; i++)
            {
                for (unsigned int j=0; j<3; j++)
                {
                    P43.at(i,j) = p3.at(4,i,j);
                    P42.at(i,j) = p2.at(4,i,j);
                    P41.at(i,j) = p1.at(4,i,j);
                }
            }

            //            DoubleMatrix P53 = A51,P43 + A52;
            //            DoubleMatrix P52 = A51,P42 + A53;
            //            DoubleMatrix P51 = A51,P41 + A54;

            //p3
            p3.at(k,0,0) = A51.at(0,0)*p3.at(4,0,0) + A51.at(0,1)*p3.at(4,1,0) + A51.at(0,2)*p3.at(4,2,0) + (-1.44);
            p3.at(k,0,1) = A51.at(0,0)*p3.at(4,0,1) + A51.at(0,1)*p3.at(4,1,1) + A51.at(0,2)*p3.at(4,2,1) + (0.0);
            p3.at(k,0,2) = A51.at(0,0)*p3.at(4,0,2) + A51.at(0,1)*p3.at(4,1,2) + A51.at(0,2)*p3.at(4,2,2) + (0.0);
            p3.at(k,1,0) = A51.at(1,0)*p3.at(4,0,0) + A51.at(1,1)*p3.at(4,1,0) + A51.at(1,2)*p3.at(4,2,0) + (0.0);
            p3.at(k,1,1) = A51.at(1,0)*p3.at(4,0,1) + A51.at(1,1)*p3.at(4,1,1) + A51.at(1,2)*p3.at(4,2,1) + (-1.44);
            p3.at(k,1,2) = A51.at(1,0)*p3.at(4,0,2) + A51.at(1,1)*p3.at(4,1,2) + A51.at(1,2)*p3.at(4,2,2) + (0.0);
            p3.at(k,2,0) = A51.at(2,0)*p3.at(4,0,0) + A51.at(2,1)*p3.at(4,1,0) + A51.at(2,2)*p3.at(4,2,0) + (0.0);
            p3.at(k,2,1) = A51.at(2,0)*p3.at(4,0,1) + A51.at(2,1)*p3.at(4,1,1) + A51.at(2,2)*p3.at(4,2,1) + (0.0);
            p3.at(k,2,2) = A51.at(2,0)*p3.at(4,0,2) + A51.at(2,1)*p3.at(4,1,2) + A51.at(2,2)*p3.at(4,2,2) + (-1.44);

            //p2
            p2.at(k,0,0) = A51.at(0,0)*p2.at(4,0,0) + A51.at(0,1)*p2.at(4,1,0) + A51.at(0,2)*p2.at(4,2,0) + (0.64);
            p2.at(k,0,1) = A51.at(0,0)*p2.at(4,0,1) + A51.at(0,1)*p2.at(4,1,1) + A51.at(0,2)*p2.at(4,2,1) + (0.0);
            p2.at(k,0,2) = A51.at(0,0)*p2.at(4,0,2) + A51.at(0,1)*p2.at(4,1,2) + A51.at(0,2)*p2.at(4,2,2) + (0.0);
            p2.at(k,1,0) = A51.at(1,0)*p2.at(4,0,0) + A51.at(1,1)*p2.at(4,1,0) + A51.at(1,2)*p2.at(4,2,0) + (0.0);
            p2.at(k,1,1) = A51.at(1,0)*p2.at(4,0,1) + A51.at(1,1)*p2.at(4,1,1) + A51.at(1,2)*p2.at(4,2,1) + (0.64);
            p2.at(k,1,2) = A51.at(1,0)*p2.at(4,0,2) + A51.at(1,1)*p2.at(4,1,2) + A51.at(1,2)*p2.at(4,2,2) + (0.0);
            p2.at(k,2,0) = A51.at(2,0)*p2.at(4,0,0) + A51.at(2,1)*p2.at(4,1,0) + A51.at(2,2)*p2.at(4,2,0) + (0.0);
            p2.at(k,2,1) = A51.at(2,0)*p2.at(4,0,1) + A51.at(2,1)*p2.at(4,1,1) + A51.at(2,2)*p2.at(4,2,1) + (0.0);
            p2.at(k,2,2) = A51.at(2,0)*p2.at(4,0,2) + A51.at(2,1)*p2.at(4,1,2) + A51.at(2,2)*p2.at(4,2,2) + (0.64);

            //p1
            p1.at(k,0,0) = A51.at(0,0)*p1.at(4,0,0) + A51.at(0,1)*p1.at(4,1,0) + A51.at(0,2)*p1.at(4,2,0) + (-0.12);
            p1.at(k,0,1) = A51.at(0,0)*p1.at(4,0,1) + A51.at(0,1)*p1.at(4,1,1) + A51.at(0,2)*p1.at(4,2,1) + (0.0);
            p1.at(k,0,2) = A51.at(0,0)*p1.at(4,0,2) + A51.at(0,1)*p1.at(4,1,2) + A51.at(0,2)*p1.at(4,2,2) + (0.0);
            p1.at(k,1,0) = A51.at(1,0)*p1.at(4,0,0) + A51.at(1,1)*p1.at(4,1,0) + A51.at(1,2)*p1.at(4,2,0) + (0.0);
            p1.at(k,1,1) = A51.at(1,0)*p1.at(4,0,1) + A51.at(1,1)*p1.at(4,1,1) + A51.at(1,2)*p1.at(4,2,1) + (-0.12);
            p1.at(k,1,2) = A51.at(1,0)*p1.at(4,0,2) + A51.at(1,1)*p1.at(4,1,2) + A51.at(1,2)*p1.at(4,2,2) + (0.0);
            p1.at(k,2,0) = A51.at(2,0)*p1.at(4,0,0) + A51.at(2,1)*p1.at(4,1,0) + A51.at(2,2)*p1.at(4,2,0) + (0.0);
            p1.at(k,2,1) = A51.at(2,0)*p1.at(4,0,1) + A51.at(2,1)*p1.at(4,1,1) + A51.at(2,2)*p1.at(4,2,1) + (0.0);
            p1.at(k,2,2) = A51.at(2,0)*p1.at(4,0,2) + A51.at(2,1)*p1.at(4,1,2) + A51.at(2,2)*p1.at(4,2,2) + (-0.12);

            //p0
            p0.at(k,0,0) = A51.at(0,0)*p0.at(4,0,0) + A51.at(0,1)*p0.at(4,1,0) + A51.at(0,2)*p0.at(4,2,0);
            p0.at(k,0,1) = A51.at(0,0)*p0.at(4,0,1) + A51.at(0,1)*p0.at(4,1,1) + A51.at(0,2)*p0.at(4,2,1);
            p0.at(k,0,2) = A51.at(0,0)*p0.at(4,0,2) + A51.at(0,1)*p0.at(4,1,2) + A51.at(0,2)*p0.at(4,2,2);
            p0.at(k,1,0) = A51.at(1,0)*p0.at(4,0,0) + A51.at(1,1)*p0.at(4,1,0) + A51.at(1,2)*p0.at(4,2,0);
            p0.at(k,1,1) = A51.at(1,0)*p0.at(4,0,1) + A51.at(1,1)*p0.at(4,1,1) + A51.at(1,2)*p0.at(4,2,1);
            p0.at(k,1,2) = A51.at(1,0)*p0.at(4,0,2) + A51.at(1,1)*p0.at(4,1,2) + A51.at(1,2)*p0.at(4,2,2);
            p0.at(k,2,0) = A51.at(2,0)*p0.at(4,0,0) + A51.at(2,1)*p0.at(4,1,0) + A51.at(2,2)*p0.at(4,2,0);
            p0.at(k,2,1) = A51.at(2,0)*p0.at(4,0,1) + A51.at(2,1)*p0.at(4,1,1) + A51.at(2,2)*p0.at(4,2,1);
            p0.at(k,2,2) = A51.at(2,0)*p0.at(4,0,2) + A51.at(2,1)*p0.at(4,1,2) + A51.at(2,2)*p0.at(4,2,2);

            //q0
            q0.at(k,0) = A51.at(0,0)*q0.at(4,0) + A51.at(0,1)*q0.at(4,1) + A51.at(0,2)*q0.at(4,2) + 0.48*h*b(1,k1);
            q0.at(k,1) = A51.at(1,0)*q0.at(4,0) + A51.at(1,1)*q0.at(4,1) + A51.at(1,2)*q0.at(4,2) + 0.48*h*b(2,k1);
            q0.at(k,2) = A51.at(2,0)*q0.at(4,0) + A51.at(2,1)*q0.at(4,1) + A51.at(2,2)*q0.at(4,2) + 0.48*h*b(3,k1);

            //            printf("%14.10f %14.10f %14.10f\n", P53.at(0,0),P53.at(0,1),P53.at(0,2));
            //            printf("%14.10f %14.10f %14.10f\n", p3.at(k,0,0),p3.at(k,0,1),p3.at(k,0,2));
            //            printf("%14.10f %14.10f %14.10f\n", P52.at(0,0),P52.at(0,1),P52.at(0,2));
            //            printf("%14.10f %14.10f %14.10f\n", p2.at(k,0,0),p2.at(k,0,1),p2.at(k,0,2));
            //            printf("%14.10f %14.10f %14.10f\n", P51.at(0,0),P51.at(0,1),P51.at(0,2));
            //            printf("%14.10f %14.10f %14.10f\n", p1.at(k,0,0),p1.at(k,0,1),p1.at(k,0,2));
        }

        if (k==6)
        {
            unsigned int k1 = k-1;

            DoubleMatrix A61(3,3);
            A61.at(0,0) = 0.48*h*a(1,1,k1)+1.92; A61.at(0,1) = 0.48*h*a(1,2,k1);      A61.at(0,2) = 0.48*h*a(1,3,k1);
            A61.at(1,0) = 0.48*h*a(2,1,k1);      A61.at(1,1) = 0.48*h*a(2,2,k1)+1.92; A61.at(1,2) = 0.48*h*a(2,3,k1);
            A61.at(2,0) = 0.48*h*a(3,1,k1);      A61.at(2,1) = 0.48*h*a(3,2,k1);      A61.at(2,2) = 0.48*h*a(3,3,k1)+1.92;

            DoubleMatrix A62(3,3);
            A62.at(0,0) = -1.44; A62.at(0,1) = 0.0;   A62.at(0,2) = 0.0;
            A62.at(1,0) = 0.0;   A62.at(1,1) = -1.44; A62.at(1,2) = 0.0;
            A62.at(2,0) = 0.0;   A62.at(2,1) = 0.0;   A62.at(2,2) = -1.44;

            //p3
            p3.at(k,0,0) = A61.at(0,0)*p3.at(5,0,0) + A61.at(0,1)*p3.at(5,1,0) + A61.at(0,2)*p3.at(5,2,0) + A62.at(0,0)*p3.at(4,0,0) + A62.at(0,1)*p3.at(4,1,0) + A62.at(0,2)*p3.at(4,2,0) + (0.64);
            p3.at(k,0,1) = A61.at(0,0)*p3.at(5,0,1) + A61.at(0,1)*p3.at(5,1,1) + A61.at(0,2)*p3.at(5,2,1) + A62.at(0,0)*p3.at(4,0,1) + A62.at(0,1)*p3.at(4,1,1) + A62.at(0,2)*p3.at(4,2,1) + (0.00);
            p3.at(k,0,2) = A61.at(0,0)*p3.at(5,0,2) + A61.at(0,1)*p3.at(5,1,2) + A61.at(0,2)*p3.at(5,2,2) + A62.at(0,0)*p3.at(4,0,2) + A62.at(0,1)*p3.at(4,1,2) + A62.at(0,2)*p3.at(4,2,2) + (0.00);

            p3.at(k,1,0) = A61.at(1,0)*p3.at(5,0,0) + A61.at(1,1)*p3.at(5,1,0) + A61.at(1,2)*p3.at(5,2,0) + A62.at(1,0)*p3.at(4,0,0) + A62.at(1,1)*p3.at(4,1,0) + A62.at(1,2)*p3.at(4,2,0) + (0.00);
            p3.at(k,1,1) = A61.at(1,0)*p3.at(5,0,1) + A61.at(1,1)*p3.at(5,1,1) + A61.at(1,2)*p3.at(5,2,1) + A62.at(1,0)*p3.at(4,0,1) + A62.at(1,1)*p3.at(4,1,1) + A62.at(1,2)*p3.at(4,2,1) + (0.64);
            p3.at(k,1,2) = A61.at(1,0)*p3.at(5,0,2) + A61.at(1,1)*p3.at(5,1,2) + A61.at(1,2)*p3.at(5,2,2) + A62.at(1,0)*p3.at(4,0,2) + A62.at(1,1)*p3.at(4,1,2) + A62.at(1,2)*p3.at(4,2,2) + (0.00);

            p3.at(k,2,0) = A61.at(2,0)*p3.at(5,0,0) + A61.at(2,1)*p3.at(5,1,0) + A61.at(2,2)*p3.at(5,2,0) + A62.at(2,0)*p3.at(4,0,0) + A62.at(2,1)*p3.at(4,1,0) + A62.at(2,2)*p3.at(4,2,0) + (0.00);
            p3.at(k,2,1) = A61.at(2,0)*p3.at(5,0,1) + A61.at(2,1)*p3.at(5,1,1) + A61.at(2,2)*p3.at(5,2,1) + A62.at(2,0)*p3.at(4,0,1) + A62.at(2,1)*p3.at(4,1,1) + A62.at(2,2)*p3.at(4,2,1) + (0.00);
            p3.at(k,2,2) = A61.at(2,0)*p3.at(5,0,2) + A61.at(2,1)*p3.at(5,1,2) + A61.at(2,2)*p3.at(5,2,2) + A62.at(2,0)*p3.at(4,0,2) + A62.at(2,1)*p3.at(4,1,2) + A62.at(2,2)*p3.at(4,2,2) + (0.64);

            //p2
            p2.at(k,0,0) = A61.at(0,0)*p2.at(5,0,0) + A61.at(0,1)*p2.at(5,1,0) + A61.at(0,2)*p2.at(5,2,0) + A62.at(0,0)*p2.at(4,0,0) + A62.at(0,1)*p2.at(4,1,0) + A62.at(0,2)*p2.at(4,2,0) + (-0.12);
            p2.at(k,0,1) = A61.at(0,0)*p2.at(5,0,1) + A61.at(0,1)*p2.at(5,1,1) + A61.at(0,2)*p2.at(5,2,1) + A62.at(0,0)*p2.at(4,0,1) + A62.at(0,1)*p2.at(4,1,1) + A62.at(0,2)*p2.at(4,2,1) + (+0.00);
            p2.at(k,0,2) = A61.at(0,0)*p2.at(5,0,2) + A61.at(0,1)*p2.at(5,1,2) + A61.at(0,2)*p2.at(5,2,2) + A62.at(0,0)*p2.at(4,0,2) + A62.at(0,1)*p2.at(4,1,2) + A62.at(0,2)*p2.at(4,2,2) + (+0.00);

            p2.at(k,1,0) = A61.at(1,0)*p2.at(5,0,0) + A61.at(1,1)*p2.at(5,1,0) + A61.at(1,2)*p2.at(5,2,0) + A62.at(1,0)*p2.at(4,0,0) + A62.at(1,1)*p2.at(4,1,0) + A62.at(1,2)*p2.at(4,2,0) + (+0.00);
            p2.at(k,1,1) = A61.at(1,0)*p2.at(5,0,1) + A61.at(1,1)*p2.at(5,1,1) + A61.at(1,2)*p2.at(5,2,1) + A62.at(1,0)*p2.at(4,0,1) + A62.at(1,1)*p2.at(4,1,1) + A62.at(1,2)*p2.at(4,2,1) + (-0.12);
            p2.at(k,1,2) = A61.at(1,0)*p2.at(5,0,2) + A61.at(1,1)*p2.at(5,1,2) + A61.at(1,2)*p2.at(5,2,2) + A62.at(1,0)*p2.at(4,0,2) + A62.at(1,1)*p2.at(4,1,2) + A62.at(1,2)*p2.at(4,2,2) + (+0.00);

            p2.at(k,2,0) = A61.at(2,0)*p2.at(5,0,0) + A61.at(2,1)*p2.at(5,1,0) + A61.at(2,2)*p2.at(5,2,0) + A62.at(2,0)*p2.at(4,0,0) + A62.at(2,1)*p2.at(4,1,0) + A62.at(2,2)*p2.at(4,2,0) + (+0.00);
            p2.at(k,2,1) = A61.at(2,0)*p2.at(5,0,1) + A61.at(2,1)*p2.at(5,1,1) + A61.at(2,2)*p2.at(5,2,1) + A62.at(2,0)*p2.at(4,0,1) + A62.at(2,1)*p2.at(4,1,1) + A62.at(2,2)*p2.at(4,2,1) + (+0.00);
            p2.at(k,2,2) = A61.at(2,0)*p2.at(5,0,2) + A61.at(2,1)*p2.at(5,1,2) + A61.at(2,2)*p2.at(5,2,2) + A62.at(2,0)*p2.at(4,0,2) + A62.at(2,1)*p2.at(4,1,2) + A62.at(2,2)*p2.at(4,2,2) + (-0.12);

            //p1
            p1.at(k,0,0) = A61.at(0,0)*p1.at(5,0,0) + A61.at(0,1)*p1.at(5,1,0) + A61.at(0,2)*p1.at(5,2,0) + A62.at(0,0)*p1.at(4,0,0) + A62.at(0,1)*p1.at(4,1,0) + A62.at(0,2)*p1.at(4,2,0) + (+0.00);
            p1.at(k,0,1) = A61.at(0,0)*p1.at(5,0,1) + A61.at(0,1)*p1.at(5,1,1) + A61.at(0,2)*p1.at(5,2,1) + A62.at(0,0)*p1.at(4,0,1) + A62.at(0,1)*p1.at(4,1,1) + A62.at(0,2)*p1.at(4,2,1) + (+0.00);
            p1.at(k,0,2) = A61.at(0,0)*p1.at(5,0,2) + A61.at(0,1)*p1.at(5,1,2) + A61.at(0,2)*p1.at(5,2,2) + A62.at(0,0)*p1.at(4,0,2) + A62.at(0,1)*p1.at(4,1,2) + A62.at(0,2)*p1.at(4,2,2) + (+0.00);

            p1.at(k,1,0) = A61.at(1,0)*p1.at(5,0,0) + A61.at(1,1)*p1.at(5,1,0) + A61.at(1,2)*p1.at(5,2,0) + A62.at(1,0)*p1.at(4,0,0) + A62.at(1,1)*p1.at(4,1,0) + A62.at(1,2)*p1.at(4,2,0) + (+0.00);
            p1.at(k,1,1) = A61.at(1,0)*p1.at(5,0,1) + A61.at(1,1)*p1.at(5,1,1) + A61.at(1,2)*p1.at(5,2,1) + A62.at(1,0)*p1.at(4,0,1) + A62.at(1,1)*p1.at(4,1,1) + A62.at(1,2)*p1.at(4,2,1) + (+0.00);
            p1.at(k,1,2) = A61.at(1,0)*p1.at(5,0,2) + A61.at(1,1)*p1.at(5,1,2) + A61.at(1,2)*p1.at(5,2,2) + A62.at(1,0)*p1.at(4,0,2) + A62.at(1,1)*p1.at(4,1,2) + A62.at(1,2)*p1.at(4,2,2) + (+0.00);

            p1.at(k,2,0) = A61.at(2,0)*p1.at(5,0,0) + A61.at(2,1)*p1.at(5,1,0) + A61.at(2,2)*p1.at(5,2,0) + A62.at(2,0)*p1.at(4,0,0) + A62.at(2,1)*p1.at(4,1,0) + A62.at(2,2)*p1.at(4,2,0) + (+0.00);
            p1.at(k,2,1) = A61.at(2,0)*p1.at(5,0,1) + A61.at(2,1)*p1.at(5,1,1) + A61.at(2,2)*p1.at(5,2,1) + A62.at(2,0)*p1.at(4,0,1) + A62.at(2,1)*p1.at(4,1,1) + A62.at(2,2)*p1.at(4,2,1) + (+0.00);
            p1.at(k,2,2) = A61.at(2,0)*p1.at(5,0,2) + A61.at(2,1)*p1.at(5,1,2) + A61.at(2,2)*p1.at(5,2,2) + A62.at(2,0)*p1.at(4,0,2) + A62.at(2,1)*p1.at(4,1,2) + A62.at(2,2)*p1.at(4,2,2) + (+0.00);

            //p0
            p0.at(k,0,0) = A61.at(0,0)*p0.at(5,0,0) + A61.at(0,1)*p0.at(5,1,0) + A61.at(0,2)*p0.at(5,2,0) + A62.at(0,0)*p0.at(4,0,0) + A62.at(0,1)*p0.at(4,1,0) + A62.at(0,2)*p0.at(4,2,0) + (+0.00);
            p0.at(k,0,1) = A61.at(0,0)*p0.at(5,0,1) + A61.at(0,1)*p0.at(5,1,1) + A61.at(0,2)*p0.at(5,2,1) + A62.at(0,0)*p0.at(4,0,1) + A62.at(0,1)*p0.at(4,1,1) + A62.at(0,2)*p0.at(4,2,1) + (+0.00);
            p0.at(k,0,2) = A61.at(0,0)*p0.at(5,0,2) + A61.at(0,1)*p0.at(5,1,2) + A61.at(0,2)*p0.at(5,2,2) + A62.at(0,0)*p0.at(4,0,2) + A62.at(0,1)*p0.at(4,1,2) + A62.at(0,2)*p0.at(4,2,2) + (+0.00);

            p0.at(k,1,0) = A61.at(1,0)*p0.at(5,0,0) + A61.at(1,1)*p0.at(5,1,0) + A61.at(1,2)*p0.at(5,2,0) + A62.at(1,0)*p0.at(4,0,0) + A62.at(1,1)*p0.at(4,1,0) + A62.at(1,2)*p0.at(4,2,0) + (+0.00);
            p0.at(k,1,1) = A61.at(1,0)*p0.at(5,0,1) + A61.at(1,1)*p0.at(5,1,1) + A61.at(1,2)*p0.at(5,2,1) + A62.at(1,0)*p0.at(4,0,1) + A62.at(1,1)*p0.at(4,1,1) + A62.at(1,2)*p0.at(4,2,1) + (+0.00);
            p0.at(k,1,2) = A61.at(1,0)*p0.at(5,0,2) + A61.at(1,1)*p0.at(5,1,2) + A61.at(1,2)*p0.at(5,2,2) + A62.at(1,0)*p0.at(4,0,2) + A62.at(1,1)*p0.at(4,1,2) + A62.at(1,2)*p0.at(4,2,2) + (+0.00);

            p0.at(k,2,0) = A61.at(2,0)*p0.at(5,0,0) + A61.at(2,1)*p0.at(5,1,0) + A61.at(2,2)*p0.at(5,2,0) + A62.at(2,0)*p0.at(4,0,0) + A62.at(2,1)*p0.at(4,1,0) + A62.at(2,2)*p0.at(4,2,0) + (+0.00);
            p0.at(k,2,1) = A61.at(2,0)*p0.at(5,0,1) + A61.at(2,1)*p0.at(5,1,1) + A61.at(2,2)*p0.at(5,2,1) + A62.at(2,0)*p0.at(4,0,1) + A62.at(2,1)*p0.at(4,1,1) + A62.at(2,2)*p0.at(4,2,1) + (+0.00);
            p0.at(k,2,2) = A61.at(2,0)*p0.at(5,0,2) + A61.at(2,1)*p0.at(5,1,2) + A61.at(2,2)*p0.at(5,2,2) + A62.at(2,0)*p0.at(4,0,2) + A62.at(2,1)*p0.at(4,1,2) + A62.at(2,2)*p0.at(4,2,2) + (+0.00);

            //q0
            q0.at(k,0) = A61.at(0,0)*q0.at(5,0) + A61.at(0,1)*q0.at(5,1) + A61.at(0,2)*q0.at(5,2) + A62.at(0,0)*q0.at(4,0) + A62.at(0,1)*q0.at(4,1) + A62.at(0,2)*q0.at(4,2) + 0.48*h*b(1,k1);
            q0.at(k,1) = A61.at(1,0)*q0.at(5,0) + A61.at(1,1)*q0.at(5,1) + A61.at(1,2)*q0.at(5,2) + A62.at(1,0)*q0.at(4,0) + A62.at(1,1)*q0.at(4,1) + A62.at(1,2)*q0.at(4,2) + 0.48*h*b(2,k1);
            q0.at(k,2) = A61.at(2,0)*q0.at(5,0) + A61.at(2,1)*q0.at(5,1) + A61.at(2,2)*q0.at(5,2) + A62.at(2,0)*q0.at(4,0) + A62.at(2,1)*q0.at(4,1) + A62.at(2,2)*q0.at(4,2) + 0.48*h*b(3,k1);
        }

        if (k==7)
        {
            unsigned int k1 = k-1;

            DoubleMatrix A71(3,3);
            A71.at(0,0) = 0.48*h*a(1,1,k1)+1.92; A71.at(0,1) = 0.48*h*a(1,2,k1);      A71.at(0,2) = 0.48*h*a(1,3,k1);
            A71.at(1,0) = 0.48*h*a(2,1,k1);      A71.at(1,1) = 0.48*h*a(2,2,k1)+1.92; A71.at(1,2) = 0.48*h*a(2,3,k1);
            A71.at(2,0) = 0.48*h*a(3,1,k1);      A71.at(2,1) = 0.48*h*a(3,2,k1);      A71.at(2,2) = 0.48*h*a(3,3,k1)+1.92;

            DoubleMatrix A72(3,3);
            A72.at(0,0) = -1.44; A72.at(0,1) = 0.0;   A72.at(0,2) = 0.0;
            A72.at(1,0) = 0.0;   A72.at(1,1) = -1.44; A72.at(1,2) = 0.0;
            A72.at(2,0) = 0.0;   A72.at(2,1) = 0.0;   A72.at(2,2) = -1.44;

            DoubleMatrix A73(3,3);
            A73.at(0,0) = 0.64; A73.at(0,1) = 0.0;  A73.at(0,2) = 0.0;
            A73.at(1,0) = 0.0;  A73.at(1,1) = 0.64; A73.at(1,2) = 0.0;
            A73.at(2,0) = 0.0;  A73.at(2,1) = 0.0;  A73.at(2,2) = 0.64;

            //p3
            p3.at(k,0,0) = A71.at(0,0)*p3.at(6,0,0) + A71.at(0,1)*p3.at(6,1,0) + A71.at(0,2)*p3.at(6,2,0) + A72.at(0,0)*p3.at(5,0,0) + A72.at(0,1)*p3.at(5,1,0) + A72.at(0,2)*p3.at(5,2,0) + A73.at(0,0)*p3.at(4,0,0) + A73.at(0,1)*p3.at(4,1,0) + A73.at(0,2)*p3.at(4,2,0) + (-0.12);
            p3.at(k,0,1) = A71.at(0,0)*p3.at(6,0,1) + A71.at(0,1)*p3.at(6,1,1) + A71.at(0,2)*p3.at(6,2,1) + A72.at(0,0)*p3.at(5,0,1) + A72.at(0,1)*p3.at(5,1,1) + A72.at(0,2)*p3.at(5,2,1) + A73.at(0,0)*p3.at(4,0,1) + A73.at(0,1)*p3.at(4,1,1) + A73.at(0,2)*p3.at(4,2,1) + (+0.00);
            p3.at(k,0,2) = A71.at(0,0)*p3.at(6,0,2) + A71.at(0,1)*p3.at(6,1,2) + A71.at(0,2)*p3.at(6,2,2) + A72.at(0,0)*p3.at(5,0,2) + A72.at(0,1)*p3.at(5,1,2) + A72.at(0,2)*p3.at(5,2,2) + A73.at(0,0)*p3.at(4,0,2) + A73.at(0,1)*p3.at(4,1,2) + A73.at(0,2)*p3.at(4,2,2) + (+0.00);

            p3.at(k,1,0) = A71.at(1,0)*p3.at(6,0,0) + A71.at(1,1)*p3.at(6,1,0) + A71.at(1,2)*p3.at(6,2,0) + A72.at(1,0)*p3.at(5,0,0) + A72.at(1,1)*p3.at(5,1,0) + A72.at(1,2)*p3.at(5,2,0) + A73.at(1,0)*p3.at(4,0,0) + A73.at(1,1)*p3.at(4,1,0) + A73.at(1,2)*p3.at(4,2,0) + (+0.00);
            p3.at(k,1,1) = A71.at(1,0)*p3.at(6,0,1) + A71.at(1,1)*p3.at(6,1,1) + A71.at(1,2)*p3.at(6,2,1) + A72.at(1,0)*p3.at(5,0,1) + A72.at(1,1)*p3.at(5,1,1) + A72.at(1,2)*p3.at(5,2,1) + A73.at(1,0)*p3.at(4,0,1) + A73.at(1,1)*p3.at(4,1,1) + A73.at(1,2)*p3.at(4,2,1) + (-0.12);
            p3.at(k,1,2) = A71.at(1,0)*p3.at(6,0,2) + A71.at(1,1)*p3.at(6,1,2) + A71.at(1,2)*p3.at(6,2,2) + A72.at(1,0)*p3.at(5,0,2) + A72.at(1,1)*p3.at(5,1,2) + A72.at(1,2)*p3.at(5,2,2) + A73.at(1,0)*p3.at(4,0,2) + A73.at(1,1)*p3.at(4,1,2) + A73.at(1,2)*p3.at(4,2,2) + (+0.00);

            p3.at(k,2,0) = A71.at(2,0)*p3.at(6,0,0) + A71.at(2,1)*p3.at(6,1,0) + A71.at(2,2)*p3.at(6,2,0) + A72.at(2,0)*p3.at(5,0,0) + A72.at(2,1)*p3.at(5,1,0) + A72.at(2,2)*p3.at(5,2,0) + A73.at(2,0)*p3.at(4,0,0) + A73.at(2,1)*p3.at(4,1,0) + A73.at(2,2)*p3.at(4,2,0) + (+0.00);
            p3.at(k,2,1) = A71.at(2,0)*p3.at(6,0,1) + A71.at(2,1)*p3.at(6,1,1) + A71.at(2,2)*p3.at(6,2,1) + A72.at(2,0)*p3.at(5,0,1) + A72.at(2,1)*p3.at(5,1,1) + A72.at(2,2)*p3.at(5,2,1) + A73.at(2,0)*p3.at(4,0,1) + A73.at(2,1)*p3.at(4,1,1) + A73.at(2,2)*p3.at(4,2,1) + (+0.00);
            p3.at(k,2,2) = A71.at(2,0)*p3.at(6,0,2) + A71.at(2,1)*p3.at(6,1,2) + A71.at(2,2)*p3.at(6,2,2) + A72.at(2,0)*p3.at(5,0,2) + A72.at(2,1)*p3.at(5,1,2) + A72.at(2,2)*p3.at(5,2,2) + A73.at(2,0)*p3.at(4,0,2) + A73.at(2,1)*p3.at(4,1,2) + A73.at(2,2)*p3.at(4,2,2) + (-0.12);

            //p2
            p2.at(k,0,0) = A71.at(0,0)*p2.at(6,0,0) + A71.at(0,1)*p2.at(6,1,0) + A71.at(0,2)*p2.at(6,2,0) + A72.at(0,0)*p2.at(5,0,0) + A72.at(0,1)*p2.at(5,1,0) + A72.at(0,2)*p2.at(5,2,0) + A73.at(0,0)*p2.at(4,0,0) + A73.at(0,1)*p2.at(4,1,0) + A73.at(0,2)*p2.at(4,2,0) + (+0.00);
            p2.at(k,0,1) = A71.at(0,0)*p2.at(6,0,1) + A71.at(0,1)*p2.at(6,1,1) + A71.at(0,2)*p2.at(6,2,1) + A72.at(0,0)*p2.at(5,0,1) + A72.at(0,1)*p2.at(5,1,1) + A72.at(0,2)*p2.at(5,2,1) + A73.at(0,0)*p2.at(4,0,1) + A73.at(0,1)*p2.at(4,1,1) + A73.at(0,2)*p2.at(4,2,1) + (+0.00);
            p2.at(k,0,2) = A71.at(0,0)*p2.at(6,0,2) + A71.at(0,1)*p2.at(6,1,2) + A71.at(0,2)*p2.at(6,2,2) + A72.at(0,0)*p2.at(5,0,2) + A72.at(0,1)*p2.at(5,1,2) + A72.at(0,2)*p2.at(5,2,2) + A73.at(0,0)*p2.at(4,0,2) + A73.at(0,1)*p2.at(4,1,2) + A73.at(0,2)*p2.at(4,2,2) + (+0.00);

            p2.at(k,1,0) = A71.at(1,0)*p2.at(6,0,0) + A71.at(1,1)*p2.at(6,1,0) + A71.at(1,2)*p2.at(6,2,0) + A72.at(1,0)*p2.at(5,0,0) + A72.at(1,1)*p2.at(5,1,0) + A72.at(1,2)*p2.at(5,2,0) + A73.at(1,0)*p2.at(4,0,0) + A73.at(1,1)*p2.at(4,1,0) + A73.at(1,2)*p2.at(4,2,0) + (+0.00);
            p2.at(k,1,1) = A71.at(1,0)*p2.at(6,0,1) + A71.at(1,1)*p2.at(6,1,1) + A71.at(1,2)*p2.at(6,2,1) + A72.at(1,0)*p2.at(5,0,1) + A72.at(1,1)*p2.at(5,1,1) + A72.at(1,2)*p2.at(5,2,1) + A73.at(1,0)*p2.at(4,0,1) + A73.at(1,1)*p2.at(4,1,1) + A73.at(1,2)*p2.at(4,2,1) + (+0.00);
            p2.at(k,1,2) = A71.at(1,0)*p2.at(6,0,2) + A71.at(1,1)*p2.at(6,1,2) + A71.at(1,2)*p2.at(6,2,2) + A72.at(1,0)*p2.at(5,0,2) + A72.at(1,1)*p2.at(5,1,2) + A72.at(1,2)*p2.at(5,2,2) + A73.at(1,0)*p2.at(4,0,2) + A73.at(1,1)*p2.at(4,1,2) + A73.at(1,2)*p2.at(4,2,2) + (+0.00);

            p2.at(k,2,0) = A71.at(2,0)*p2.at(6,0,0) + A71.at(2,1)*p2.at(6,1,0) + A71.at(2,2)*p2.at(6,2,0) + A72.at(2,0)*p2.at(5,0,0) + A72.at(2,1)*p2.at(5,1,0) + A72.at(2,2)*p2.at(5,2,0) + A73.at(2,0)*p2.at(4,0,0) + A73.at(2,1)*p2.at(4,1,0) + A73.at(2,2)*p2.at(4,2,0) + (+0.00);
            p2.at(k,2,1) = A71.at(2,0)*p2.at(6,0,1) + A71.at(2,1)*p2.at(6,1,1) + A71.at(2,2)*p2.at(6,2,1) + A72.at(2,0)*p2.at(5,0,1) + A72.at(2,1)*p2.at(5,1,1) + A72.at(2,2)*p2.at(5,2,1) + A73.at(2,0)*p2.at(4,0,1) + A73.at(2,1)*p2.at(4,1,1) + A73.at(2,2)*p2.at(4,2,1) + (+0.00);
            p2.at(k,2,2) = A71.at(2,0)*p2.at(6,0,2) + A71.at(2,1)*p2.at(6,1,2) + A71.at(2,2)*p2.at(6,2,2) + A72.at(2,0)*p2.at(5,0,2) + A72.at(2,1)*p2.at(5,1,2) + A72.at(2,2)*p2.at(5,2,2) + A73.at(2,0)*p2.at(4,0,2) + A73.at(2,1)*p2.at(4,1,2) + A73.at(2,2)*p2.at(4,2,2) + (+0.00);

            //p1
            p1.at(k,0,0) = A71.at(0,0)*p1.at(6,0,0) + A71.at(0,1)*p1.at(6,1,0) + A71.at(0,2)*p1.at(6,2,0) + A72.at(0,0)*p1.at(5,0,0) + A72.at(0,1)*p1.at(5,1,0) + A72.at(0,2)*p1.at(5,2,0) + A73.at(0,0)*p1.at(4,0,0) + A73.at(0,1)*p1.at(4,1,0) + A73.at(0,2)*p1.at(4,2,0) + (+0.00);
            p1.at(k,0,1) = A71.at(0,0)*p1.at(6,0,1) + A71.at(0,1)*p1.at(6,1,1) + A71.at(0,2)*p1.at(6,2,1) + A72.at(0,0)*p1.at(5,0,1) + A72.at(0,1)*p1.at(5,1,1) + A72.at(0,2)*p1.at(5,2,1) + A73.at(0,0)*p1.at(4,0,1) + A73.at(0,1)*p1.at(4,1,1) + A73.at(0,2)*p1.at(4,2,1) + (+0.00);
            p1.at(k,0,2) = A71.at(0,0)*p1.at(6,0,2) + A71.at(0,1)*p1.at(6,1,2) + A71.at(0,2)*p1.at(6,2,2) + A72.at(0,0)*p1.at(5,0,2) + A72.at(0,1)*p1.at(5,1,2) + A72.at(0,2)*p1.at(5,2,2) + A73.at(0,0)*p1.at(4,0,2) + A73.at(0,1)*p1.at(4,1,2) + A73.at(0,2)*p1.at(4,2,2) + (+0.00);

            p1.at(k,1,0) = A71.at(1,0)*p1.at(6,0,0) + A71.at(1,1)*p1.at(6,1,0) + A71.at(1,2)*p1.at(6,2,0) + A72.at(1,0)*p1.at(5,0,0) + A72.at(1,1)*p1.at(5,1,0) + A72.at(1,2)*p1.at(5,2,0) + A73.at(1,0)*p1.at(4,0,0) + A73.at(1,1)*p1.at(4,1,0) + A73.at(1,2)*p1.at(4,2,0) + (+0.00);
            p1.at(k,1,1) = A71.at(1,0)*p1.at(6,0,1) + A71.at(1,1)*p1.at(6,1,1) + A71.at(1,2)*p1.at(6,2,1) + A72.at(1,0)*p1.at(5,0,1) + A72.at(1,1)*p1.at(5,1,1) + A72.at(1,2)*p1.at(5,2,1) + A73.at(1,0)*p1.at(4,0,1) + A73.at(1,1)*p1.at(4,1,1) + A73.at(1,2)*p1.at(4,2,1) + (+0.00);
            p1.at(k,1,2) = A71.at(1,0)*p1.at(6,0,2) + A71.at(1,1)*p1.at(6,1,2) + A71.at(1,2)*p1.at(6,2,2) + A72.at(1,0)*p1.at(5,0,2) + A72.at(1,1)*p1.at(5,1,2) + A72.at(1,2)*p1.at(5,2,2) + A73.at(1,0)*p1.at(4,0,2) + A73.at(1,1)*p1.at(4,1,2) + A73.at(1,2)*p1.at(4,2,2) + (+0.00);

            p1.at(k,2,0) = A71.at(2,0)*p1.at(6,0,0) + A71.at(2,1)*p1.at(6,1,0) + A71.at(2,2)*p1.at(6,2,0) + A72.at(2,0)*p1.at(5,0,0) + A72.at(2,1)*p1.at(5,1,0) + A72.at(2,2)*p1.at(5,2,0) + A73.at(2,0)*p1.at(4,0,0) + A73.at(2,1)*p1.at(4,1,0) + A73.at(2,2)*p1.at(4,2,0) + (+0.00);
            p1.at(k,2,1) = A71.at(2,0)*p1.at(6,0,1) + A71.at(2,1)*p1.at(6,1,1) + A71.at(2,2)*p1.at(6,2,1) + A72.at(2,0)*p1.at(5,0,1) + A72.at(2,1)*p1.at(5,1,1) + A72.at(2,2)*p1.at(5,2,1) + A73.at(2,0)*p1.at(4,0,1) + A73.at(2,1)*p1.at(4,1,1) + A73.at(2,2)*p1.at(4,2,1) + (+0.00);
            p1.at(k,2,2) = A71.at(2,0)*p1.at(6,0,2) + A71.at(2,1)*p1.at(6,1,2) + A71.at(2,2)*p1.at(6,2,2) + A72.at(2,0)*p1.at(5,0,2) + A72.at(2,1)*p1.at(5,1,2) + A72.at(2,2)*p1.at(5,2,2) + A73.at(2,0)*p1.at(4,0,2) + A73.at(2,1)*p1.at(4,1,2) + A73.at(2,2)*p1.at(4,2,2) + (+0.00);

            //p0
            p0.at(k,0,0) = A71.at(0,0)*p0.at(6,0,0) + A71.at(0,1)*p0.at(6,1,0) + A71.at(0,2)*p0.at(6,2,0) + A72.at(0,0)*p0.at(5,0,0) + A72.at(0,1)*p0.at(5,1,0) + A72.at(0,2)*p0.at(5,2,0) + A73.at(0,0)*p0.at(4,0,0) + A73.at(0,1)*p0.at(4,1,0) + A73.at(0,2)*p0.at(4,2,0) + (+0.00);
            p0.at(k,0,1) = A71.at(0,0)*p0.at(6,0,1) + A71.at(0,1)*p0.at(6,1,1) + A71.at(0,2)*p0.at(6,2,1) + A72.at(0,0)*p0.at(5,0,1) + A72.at(0,1)*p0.at(5,1,1) + A72.at(0,2)*p0.at(5,2,1) + A73.at(0,0)*p0.at(4,0,1) + A73.at(0,1)*p0.at(4,1,1) + A73.at(0,2)*p0.at(4,2,1) + (+0.00);
            p0.at(k,0,2) = A71.at(0,0)*p0.at(6,0,2) + A71.at(0,1)*p0.at(6,1,2) + A71.at(0,2)*p0.at(6,2,2) + A72.at(0,0)*p0.at(5,0,2) + A72.at(0,1)*p0.at(5,1,2) + A72.at(0,2)*p0.at(5,2,2) + A73.at(0,0)*p0.at(4,0,2) + A73.at(0,1)*p0.at(4,1,2) + A73.at(0,2)*p0.at(4,2,2) + (+0.00);

            p0.at(k,1,0) = A71.at(1,0)*p0.at(6,0,0) + A71.at(1,1)*p0.at(6,1,0) + A71.at(1,2)*p0.at(6,2,0) + A72.at(1,0)*p0.at(5,0,0) + A72.at(1,1)*p0.at(5,1,0) + A72.at(1,2)*p0.at(5,2,0) + A73.at(1,0)*p0.at(4,0,0) + A73.at(1,1)*p0.at(4,1,0) + A73.at(1,2)*p0.at(4,2,0) + (+0.00);
            p0.at(k,1,1) = A71.at(1,0)*p0.at(6,0,1) + A71.at(1,1)*p0.at(6,1,1) + A71.at(1,2)*p0.at(6,2,1) + A72.at(1,0)*p0.at(5,0,1) + A72.at(1,1)*p0.at(5,1,1) + A72.at(1,2)*p0.at(5,2,1) + A73.at(1,0)*p0.at(4,0,1) + A73.at(1,1)*p0.at(4,1,1) + A73.at(1,2)*p0.at(4,2,1) + (+0.00);
            p0.at(k,1,2) = A71.at(1,0)*p0.at(6,0,2) + A71.at(1,1)*p0.at(6,1,2) + A71.at(1,2)*p0.at(6,2,2) + A72.at(1,0)*p0.at(5,0,2) + A72.at(1,1)*p0.at(5,1,2) + A72.at(1,2)*p0.at(5,2,2) + A73.at(1,0)*p0.at(4,0,2) + A73.at(1,1)*p0.at(4,1,2) + A73.at(1,2)*p0.at(4,2,2) + (+0.00);

            p0.at(k,2,0) = A71.at(2,0)*p0.at(6,0,0) + A71.at(2,1)*p0.at(6,1,0) + A71.at(2,2)*p0.at(6,2,0) + A72.at(2,0)*p0.at(5,0,0) + A72.at(2,1)*p0.at(5,1,0) + A72.at(2,2)*p0.at(5,2,0) + A73.at(2,0)*p0.at(4,0,0) + A73.at(2,1)*p0.at(4,1,0) + A73.at(2,2)*p0.at(4,2,0) + (+0.00);
            p0.at(k,2,1) = A71.at(2,0)*p0.at(6,0,1) + A71.at(2,1)*p0.at(6,1,1) + A71.at(2,2)*p0.at(6,2,1) + A72.at(2,0)*p0.at(5,0,1) + A72.at(2,1)*p0.at(5,1,1) + A72.at(2,2)*p0.at(5,2,1) + A73.at(2,0)*p0.at(4,0,1) + A73.at(2,1)*p0.at(4,1,1) + A73.at(2,2)*p0.at(4,2,1) + (+0.00);
            p0.at(k,2,2) = A71.at(2,0)*p0.at(6,0,2) + A71.at(2,1)*p0.at(6,1,2) + A71.at(2,2)*p0.at(6,2,2) + A72.at(2,0)*p0.at(5,0,2) + A72.at(2,1)*p0.at(5,1,2) + A72.at(2,2)*p0.at(5,2,2) + A73.at(2,0)*p0.at(4,0,2) + A73.at(2,1)*p0.at(4,1,2) + A73.at(2,2)*p0.at(4,2,2) + (+0.00);

            //q0
            q0.at(k,0) = A71.at(0,0)*q0.at(6,0) + A71.at(0,1)*q0.at(6,1) + A71.at(0,2)*q0.at(6,2) + A72.at(0,0)*q0.at(5,0) + A72.at(0,1)*q0.at(5,1) + A72.at(0,2)*q0.at(5,2) + A73.at(0,0)*q0.at(4,0) + A73.at(0,1)*q0.at(4,1) + A73.at(0,2)*q0.at(4,2) + 0.48*h*b(1,k1);
            q0.at(k,1) = A71.at(1,0)*q0.at(6,0) + A71.at(1,1)*q0.at(6,1) + A71.at(1,2)*q0.at(6,2) + A72.at(1,0)*q0.at(5,0) + A72.at(1,1)*q0.at(5,1) + A72.at(1,2)*q0.at(5,2) + A73.at(1,0)*q0.at(4,0) + A73.at(1,1)*q0.at(4,1) + A73.at(1,2)*q0.at(4,2) + 0.48*h*b(2,k1);
            q0.at(k,2) = A71.at(2,0)*q0.at(6,0) + A71.at(2,1)*q0.at(6,1) + A71.at(2,2)*q0.at(6,2) + A72.at(2,0)*q0.at(5,0) + A72.at(2,1)*q0.at(5,1) + A72.at(2,2)*q0.at(5,2) + A73.at(2,0)*q0.at(4,0) + A73.at(2,1)*q0.at(4,1) + A73.at(2,2)*q0.at(4,2) + 0.48*h*b(3,k1);
        }

        if (k>=8)
        {
            unsigned int k1 = k-1;

            DoubleMatrix A71(3,3);
            A71.at(0,0) = 0.48*h*a(1,1,k1)+1.92; A71.at(0,1) = 0.48*h*a(1,2,k1);      A71.at(0,2) = 0.48*h*a(1,3,k1);
            A71.at(1,0) = 0.48*h*a(2,1,k1);      A71.at(1,1) = 0.48*h*a(2,2,k1)+1.92; A71.at(1,2) = 0.48*h*a(2,3,k1);
            A71.at(2,0) = 0.48*h*a(3,1,k1);      A71.at(2,1) = 0.48*h*a(3,2,k1);      A71.at(2,2) = 0.48*h*a(3,3,k1)+1.92;

            DoubleMatrix A72(3,3);
            A72.at(0,0) = -1.44; A72.at(0,1) = 0.0;   A72.at(0,2) = 0.0;
            A72.at(1,0) = 0.0;   A72.at(1,1) = -1.44; A72.at(1,2) = 0.0;
            A72.at(2,0) = 0.0;   A72.at(2,1) = 0.0;   A72.at(2,2) = -1.44;

            DoubleMatrix A73(3,3);
            A73.at(0,0) = 0.64; A73.at(0,1) = 0.0;  A73.at(0,2) = 0.0;
            A73.at(1,0) = 0.0;  A73.at(1,1) = 0.64; A73.at(1,2) = 0.0;
            A73.at(2,0) = 0.0;  A73.at(2,1) = 0.0;  A73.at(2,2) = 0.64;

            DoubleMatrix A74(3,3);
            A74.at(0,0) = -0.12; A74.at(0,1) = +0.00; A74.at(0,2) = +0.00;
            A74.at(1,0) = +0.00; A74.at(1,1) = -0.12; A74.at(1,2) = +0.00;
            A74.at(2,0) = +0.00; A74.at(2,1) = +0.00; A74.at(2,2) = -0.12;

            //p3
            p3.at(k,0,0) = A71.at(0,0)*p3.at(k-1,0,0) + A71.at(0,1)*p3.at(k-1,1,0) + A71.at(0,2)*p3.at(k-1,2,0) + A72.at(0,0)*p3.at(k-2,0,0) + A72.at(0,1)*p3.at(k-2,1,0) + A72.at(0,2)*p3.at(k-2,2,0) + A73.at(0,0)*p3.at(k-3,0,0) + A73.at(0,1)*p3.at(k-3,1,0) + A73.at(0,2)*p3.at(k-3,2,0) + A74.at(0,0)*p3.at(k-4,0,0) + A74.at(0,1)*p3.at(k-4,1,0) + A74.at(0,2)*p3.at(k-4,2,0);
            p3.at(k,0,1) = A71.at(0,0)*p3.at(k-1,0,1) + A71.at(0,1)*p3.at(k-1,1,1) + A71.at(0,2)*p3.at(k-1,2,1) + A72.at(0,0)*p3.at(k-2,0,1) + A72.at(0,1)*p3.at(k-2,1,1) + A72.at(0,2)*p3.at(k-2,2,1) + A73.at(0,0)*p3.at(k-3,0,1) + A73.at(0,1)*p3.at(k-3,1,1) + A73.at(0,2)*p3.at(k-3,2,1) + A74.at(0,0)*p3.at(k-4,0,1) + A74.at(0,1)*p3.at(k-4,1,1) + A74.at(0,2)*p3.at(k-4,2,1);
            p3.at(k,0,2) = A71.at(0,0)*p3.at(k-1,0,2) + A71.at(0,1)*p3.at(k-1,1,2) + A71.at(0,2)*p3.at(k-1,2,2) + A72.at(0,0)*p3.at(k-2,0,2) + A72.at(0,1)*p3.at(k-2,1,2) + A72.at(0,2)*p3.at(k-2,2,2) + A73.at(0,0)*p3.at(k-3,0,2) + A73.at(0,1)*p3.at(k-3,1,2) + A73.at(0,2)*p3.at(k-3,2,2) + A74.at(0,0)*p3.at(k-4,0,2) + A74.at(0,1)*p3.at(k-4,1,2) + A74.at(0,2)*p3.at(k-4,2,2);

            p3.at(k,1,0) = A71.at(1,0)*p3.at(k-1,0,0) + A71.at(1,1)*p3.at(k-1,1,0) + A71.at(1,2)*p3.at(k-1,2,0) + A72.at(1,0)*p3.at(k-2,0,0) + A72.at(1,1)*p3.at(k-2,1,0) + A72.at(1,2)*p3.at(k-2,2,0) + A73.at(1,0)*p3.at(k-3,0,0) + A73.at(1,1)*p3.at(k-3,1,0) + A73.at(1,2)*p3.at(k-3,2,0) + A74.at(1,0)*p3.at(k-4,0,0) + A74.at(1,1)*p3.at(k-4,1,0) + A74.at(1,2)*p3.at(k-4,2,0);
            p3.at(k,1,1) = A71.at(1,0)*p3.at(k-1,0,1) + A71.at(1,1)*p3.at(k-1,1,1) + A71.at(1,2)*p3.at(k-1,2,1) + A72.at(1,0)*p3.at(k-2,0,1) + A72.at(1,1)*p3.at(k-2,1,1) + A72.at(1,2)*p3.at(k-2,2,1) + A73.at(1,0)*p3.at(k-3,0,1) + A73.at(1,1)*p3.at(k-3,1,1) + A73.at(1,2)*p3.at(k-3,2,1) + A74.at(1,0)*p3.at(k-4,0,1) + A74.at(1,1)*p3.at(k-4,1,1) + A74.at(1,2)*p3.at(k-4,2,1);
            p3.at(k,1,2) = A71.at(1,0)*p3.at(k-1,0,2) + A71.at(1,1)*p3.at(k-1,1,2) + A71.at(1,2)*p3.at(k-1,2,2) + A72.at(1,0)*p3.at(k-2,0,2) + A72.at(1,1)*p3.at(k-2,1,2) + A72.at(1,2)*p3.at(k-2,2,2) + A73.at(1,0)*p3.at(k-3,0,2) + A73.at(1,1)*p3.at(k-3,1,2) + A73.at(1,2)*p3.at(k-3,2,2) + A74.at(1,0)*p3.at(k-4,0,2) + A74.at(1,1)*p3.at(k-4,1,2) + A74.at(1,2)*p3.at(k-4,2,2);

            p3.at(k,2,0) = A71.at(2,0)*p3.at(k-1,0,0) + A71.at(2,1)*p3.at(k-1,1,0) + A71.at(2,2)*p3.at(k-1,2,0) + A72.at(2,0)*p3.at(k-2,0,0) + A72.at(2,1)*p3.at(k-2,1,0) + A72.at(2,2)*p3.at(k-2,2,0) + A73.at(2,0)*p3.at(k-3,0,0) + A73.at(2,1)*p3.at(k-3,1,0) + A73.at(2,2)*p3.at(k-3,2,0) + A74.at(2,0)*p3.at(k-4,0,0) + A74.at(2,1)*p3.at(k-4,1,0) + A74.at(2,2)*p3.at(k-4,2,0);
            p3.at(k,2,1) = A71.at(2,0)*p3.at(k-1,0,1) + A71.at(2,1)*p3.at(k-1,1,1) + A71.at(2,2)*p3.at(k-1,2,1) + A72.at(2,0)*p3.at(k-2,0,1) + A72.at(2,1)*p3.at(k-2,1,1) + A72.at(2,2)*p3.at(k-2,2,1) + A73.at(2,0)*p3.at(k-3,0,1) + A73.at(2,1)*p3.at(k-3,1,1) + A73.at(2,2)*p3.at(k-3,2,1) + A74.at(2,0)*p3.at(k-4,0,1) + A74.at(2,1)*p3.at(k-4,1,1) + A74.at(2,2)*p3.at(k-4,2,1);
            p3.at(k,2,2) = A71.at(2,0)*p3.at(k-1,0,2) + A71.at(2,1)*p3.at(k-1,1,2) + A71.at(2,2)*p3.at(k-1,2,2) + A72.at(2,0)*p3.at(k-2,0,2) + A72.at(2,1)*p3.at(k-2,1,2) + A72.at(2,2)*p3.at(k-2,2,2) + A73.at(2,0)*p3.at(k-3,0,2) + A73.at(2,1)*p3.at(k-3,1,2) + A73.at(2,2)*p3.at(k-3,2,2) + A74.at(2,0)*p3.at(k-4,0,2) + A74.at(2,1)*p3.at(k-4,1,2) + A74.at(2,2)*p3.at(k-4,2,2);

            //p2
            p2.at(k,0,0) = A71.at(0,0)*p2.at(k-1,0,0) + A71.at(0,1)*p2.at(k-1,1,0) + A71.at(0,2)*p2.at(k-1,2,0) + A72.at(0,0)*p2.at(k-2,0,0) + A72.at(0,1)*p2.at(k-2,1,0) + A72.at(0,2)*p2.at(k-2,2,0) + A73.at(0,0)*p2.at(k-3,0,0) + A73.at(0,1)*p2.at(k-3,1,0) + A73.at(0,2)*p2.at(k-3,2,0) + A74.at(0,0)*p2.at(k-4,0,0) + A74.at(0,1)*p2.at(k-4,1,0) + A74.at(0,2)*p2.at(k-4,2,0);
            p2.at(k,0,1) = A71.at(0,0)*p2.at(k-1,0,1) + A71.at(0,1)*p2.at(k-1,1,1) + A71.at(0,2)*p2.at(k-1,2,1) + A72.at(0,0)*p2.at(k-2,0,1) + A72.at(0,1)*p2.at(k-2,1,1) + A72.at(0,2)*p2.at(k-2,2,1) + A73.at(0,0)*p2.at(k-3,0,1) + A73.at(0,1)*p2.at(k-3,1,1) + A73.at(0,2)*p2.at(k-3,2,1) + A74.at(0,0)*p2.at(k-4,0,1) + A74.at(0,1)*p2.at(k-4,1,1) + A74.at(0,2)*p2.at(k-4,2,1);
            p2.at(k,0,2) = A71.at(0,0)*p2.at(k-1,0,2) + A71.at(0,1)*p2.at(k-1,1,2) + A71.at(0,2)*p2.at(k-1,2,2) + A72.at(0,0)*p2.at(k-2,0,2) + A72.at(0,1)*p2.at(k-2,1,2) + A72.at(0,2)*p2.at(k-2,2,2) + A73.at(0,0)*p2.at(k-3,0,2) + A73.at(0,1)*p2.at(k-3,1,2) + A73.at(0,2)*p2.at(k-3,2,2) + A74.at(0,0)*p2.at(k-4,0,2) + A74.at(0,1)*p2.at(k-4,1,2) + A74.at(0,2)*p2.at(k-4,2,2);

            p2.at(k,1,0) = A71.at(1,0)*p2.at(k-1,0,0) + A71.at(1,1)*p2.at(k-1,1,0) + A71.at(1,2)*p2.at(k-1,2,0) + A72.at(1,0)*p2.at(k-2,0,0) + A72.at(1,1)*p2.at(k-2,1,0) + A72.at(1,2)*p2.at(k-2,2,0) + A73.at(1,0)*p2.at(k-3,0,0) + A73.at(1,1)*p2.at(k-3,1,0) + A73.at(1,2)*p2.at(k-3,2,0) + A74.at(1,0)*p2.at(k-4,0,0) + A74.at(1,1)*p2.at(k-4,1,0) + A74.at(1,2)*p2.at(k-4,2,0);
            p2.at(k,1,1) = A71.at(1,0)*p2.at(k-1,0,1) + A71.at(1,1)*p2.at(k-1,1,1) + A71.at(1,2)*p2.at(k-1,2,1) + A72.at(1,0)*p2.at(k-2,0,1) + A72.at(1,1)*p2.at(k-2,1,1) + A72.at(1,2)*p2.at(k-2,2,1) + A73.at(1,0)*p2.at(k-3,0,1) + A73.at(1,1)*p2.at(k-3,1,1) + A73.at(1,2)*p2.at(k-3,2,1) + A74.at(1,0)*p2.at(k-4,0,1) + A74.at(1,1)*p2.at(k-4,1,1) + A74.at(1,2)*p2.at(k-4,2,1);
            p2.at(k,1,2) = A71.at(1,0)*p2.at(k-1,0,2) + A71.at(1,1)*p2.at(k-1,1,2) + A71.at(1,2)*p2.at(k-1,2,2) + A72.at(1,0)*p2.at(k-2,0,2) + A72.at(1,1)*p2.at(k-2,1,2) + A72.at(1,2)*p2.at(k-2,2,2) + A73.at(1,0)*p2.at(k-3,0,2) + A73.at(1,1)*p2.at(k-3,1,2) + A73.at(1,2)*p2.at(k-3,2,2) + A74.at(1,0)*p2.at(k-4,0,2) + A74.at(1,1)*p2.at(k-4,1,2) + A74.at(1,2)*p2.at(k-4,2,2);

            p2.at(k,2,0) = A71.at(2,0)*p2.at(k-1,0,0) + A71.at(2,1)*p2.at(k-1,1,0) + A71.at(2,2)*p2.at(k-1,2,0) + A72.at(2,0)*p2.at(k-2,0,0) + A72.at(2,1)*p2.at(k-2,1,0) + A72.at(2,2)*p2.at(k-2,2,0) + A73.at(2,0)*p2.at(k-3,0,0) + A73.at(2,1)*p2.at(k-3,1,0) + A73.at(2,2)*p2.at(k-3,2,0) + A74.at(2,0)*p2.at(k-4,0,0) + A74.at(2,1)*p2.at(k-4,1,0) + A74.at(2,2)*p2.at(k-4,2,0);
            p2.at(k,2,1) = A71.at(2,0)*p2.at(k-1,0,1) + A71.at(2,1)*p2.at(k-1,1,1) + A71.at(2,2)*p2.at(k-1,2,1) + A72.at(2,0)*p2.at(k-2,0,1) + A72.at(2,1)*p2.at(k-2,1,1) + A72.at(2,2)*p2.at(k-2,2,1) + A73.at(2,0)*p2.at(k-3,0,1) + A73.at(2,1)*p2.at(k-3,1,1) + A73.at(2,2)*p2.at(k-3,2,1) + A74.at(2,0)*p2.at(k-4,0,1) + A74.at(2,1)*p2.at(k-4,1,1) + A74.at(2,2)*p2.at(k-4,2,1);
            p2.at(k,2,2) = A71.at(2,0)*p2.at(k-1,0,2) + A71.at(2,1)*p2.at(k-1,1,2) + A71.at(2,2)*p2.at(k-1,2,2) + A72.at(2,0)*p2.at(k-2,0,2) + A72.at(2,1)*p2.at(k-2,1,2) + A72.at(2,2)*p2.at(k-2,2,2) + A73.at(2,0)*p2.at(k-3,0,2) + A73.at(2,1)*p2.at(k-3,1,2) + A73.at(2,2)*p2.at(k-3,2,2) + A74.at(2,0)*p2.at(k-4,0,2) + A74.at(2,1)*p2.at(k-4,1,2) + A74.at(2,2)*p2.at(k-4,2,2);

            //p1
            p1.at(k,0,0) = A71.at(0,0)*p1.at(k-1,0,0) + A71.at(0,1)*p1.at(k-1,1,0) + A71.at(0,2)*p1.at(k-1,2,0) + A72.at(0,0)*p1.at(k-2,0,0) + A72.at(0,1)*p1.at(k-2,1,0) + A72.at(0,2)*p1.at(k-2,2,0) + A73.at(0,0)*p1.at(k-3,0,0) + A73.at(0,1)*p1.at(k-3,1,0) + A73.at(0,2)*p1.at(k-3,2,0) + A74.at(0,0)*p1.at(k-4,0,0) + A74.at(0,1)*p1.at(k-4,1,0) + A74.at(0,2)*p1.at(k-4,2,0);
            p1.at(k,0,1) = A71.at(0,0)*p1.at(k-1,0,1) + A71.at(0,1)*p1.at(k-1,1,1) + A71.at(0,2)*p1.at(k-1,2,1) + A72.at(0,0)*p1.at(k-2,0,1) + A72.at(0,1)*p1.at(k-2,1,1) + A72.at(0,2)*p1.at(k-2,2,1) + A73.at(0,0)*p1.at(k-3,0,1) + A73.at(0,1)*p1.at(k-3,1,1) + A73.at(0,2)*p1.at(k-3,2,1) + A74.at(0,0)*p1.at(k-4,0,1) + A74.at(0,1)*p1.at(k-4,1,1) + A74.at(0,2)*p1.at(k-4,2,1);
            p1.at(k,0,2) = A71.at(0,0)*p1.at(k-1,0,2) + A71.at(0,1)*p1.at(k-1,1,2) + A71.at(0,2)*p1.at(k-1,2,2) + A72.at(0,0)*p1.at(k-2,0,2) + A72.at(0,1)*p1.at(k-2,1,2) + A72.at(0,2)*p1.at(k-2,2,2) + A73.at(0,0)*p1.at(k-3,0,2) + A73.at(0,1)*p1.at(k-3,1,2) + A73.at(0,2)*p1.at(k-3,2,2) + A74.at(0,0)*p1.at(k-4,0,2) + A74.at(0,1)*p1.at(k-4,1,2) + A74.at(0,2)*p1.at(k-4,2,2);

            p1.at(k,1,0) = A71.at(1,0)*p1.at(k-1,0,0) + A71.at(1,1)*p1.at(k-1,1,0) + A71.at(1,2)*p1.at(k-1,2,0) + A72.at(1,0)*p1.at(k-2,0,0) + A72.at(1,1)*p1.at(k-2,1,0) + A72.at(1,2)*p1.at(k-2,2,0) + A73.at(1,0)*p1.at(k-3,0,0) + A73.at(1,1)*p1.at(k-3,1,0) + A73.at(1,2)*p1.at(k-3,2,0) + A74.at(1,0)*p1.at(k-4,0,0) + A74.at(1,1)*p1.at(k-4,1,0) + A74.at(1,2)*p1.at(k-4,2,0);
            p1.at(k,1,1) = A71.at(1,0)*p1.at(k-1,0,1) + A71.at(1,1)*p1.at(k-1,1,1) + A71.at(1,2)*p1.at(k-1,2,1) + A72.at(1,0)*p1.at(k-2,0,1) + A72.at(1,1)*p1.at(k-2,1,1) + A72.at(1,2)*p1.at(k-2,2,1) + A73.at(1,0)*p1.at(k-3,0,1) + A73.at(1,1)*p1.at(k-3,1,1) + A73.at(1,2)*p1.at(k-3,2,1) + A74.at(1,0)*p1.at(k-4,0,1) + A74.at(1,1)*p1.at(k-4,1,1) + A74.at(1,2)*p1.at(k-4,2,1);
            p1.at(k,1,2) = A71.at(1,0)*p1.at(k-1,0,2) + A71.at(1,1)*p1.at(k-1,1,2) + A71.at(1,2)*p1.at(k-1,2,2) + A72.at(1,0)*p1.at(k-2,0,2) + A72.at(1,1)*p1.at(k-2,1,2) + A72.at(1,2)*p1.at(k-2,2,2) + A73.at(1,0)*p1.at(k-3,0,2) + A73.at(1,1)*p1.at(k-3,1,2) + A73.at(1,2)*p1.at(k-3,2,2) + A74.at(1,0)*p1.at(k-4,0,2) + A74.at(1,1)*p1.at(k-4,1,2) + A74.at(1,2)*p1.at(k-4,2,2);

            p1.at(k,2,0) = A71.at(2,0)*p1.at(k-1,0,0) + A71.at(2,1)*p1.at(k-1,1,0) + A71.at(2,2)*p1.at(k-1,2,0) + A72.at(2,0)*p1.at(k-2,0,0) + A72.at(2,1)*p1.at(k-2,1,0) + A72.at(2,2)*p1.at(k-2,2,0) + A73.at(2,0)*p1.at(k-3,0,0) + A73.at(2,1)*p1.at(k-3,1,0) + A73.at(2,2)*p1.at(k-3,2,0) + A74.at(2,0)*p1.at(k-4,0,0) + A74.at(2,1)*p1.at(k-4,1,0) + A74.at(2,2)*p1.at(k-4,2,0);
            p1.at(k,2,1) = A71.at(2,0)*p1.at(k-1,0,1) + A71.at(2,1)*p1.at(k-1,1,1) + A71.at(2,2)*p1.at(k-1,2,1) + A72.at(2,0)*p1.at(k-2,0,1) + A72.at(2,1)*p1.at(k-2,1,1) + A72.at(2,2)*p1.at(k-2,2,1) + A73.at(2,0)*p1.at(k-3,0,1) + A73.at(2,1)*p1.at(k-3,1,1) + A73.at(2,2)*p1.at(k-3,2,1) + A74.at(2,0)*p1.at(k-4,0,1) + A74.at(2,1)*p1.at(k-4,1,1) + A74.at(2,2)*p1.at(k-4,2,1);
            p1.at(k,2,2) = A71.at(2,0)*p1.at(k-1,0,2) + A71.at(2,1)*p1.at(k-1,1,2) + A71.at(2,2)*p1.at(k-1,2,2) + A72.at(2,0)*p1.at(k-2,0,2) + A72.at(2,1)*p1.at(k-2,1,2) + A72.at(2,2)*p1.at(k-2,2,2) + A73.at(2,0)*p1.at(k-3,0,2) + A73.at(2,1)*p1.at(k-3,1,2) + A73.at(2,2)*p1.at(k-3,2,2) + A74.at(2,0)*p1.at(k-4,0,2) + A74.at(2,1)*p1.at(k-4,1,2) + A74.at(2,2)*p1.at(k-4,2,2);

            //p0
            p0.at(k,0,0) = A71.at(0,0)*p0.at(k-1,0,0) + A71.at(0,1)*p0.at(k-1,1,0) + A71.at(0,2)*p0.at(k-1,2,0) + A72.at(0,0)*p0.at(k-2,0,0) + A72.at(0,1)*p0.at(k-2,1,0) + A72.at(0,2)*p0.at(k-2,2,0) + A73.at(0,0)*p0.at(k-3,0,0) + A73.at(0,1)*p0.at(k-3,1,0) + A73.at(0,2)*p0.at(k-3,2,0) + A74.at(0,0)*p0.at(k-4,0,0) + A74.at(0,1)*p0.at(k-4,1,0) + A74.at(0,2)*p0.at(k-4,2,0);
            p0.at(k,0,1) = A71.at(0,0)*p0.at(k-1,0,1) + A71.at(0,1)*p0.at(k-1,1,1) + A71.at(0,2)*p0.at(k-1,2,1) + A72.at(0,0)*p0.at(k-2,0,1) + A72.at(0,1)*p0.at(k-2,1,1) + A72.at(0,2)*p0.at(k-2,2,1) + A73.at(0,0)*p0.at(k-3,0,1) + A73.at(0,1)*p0.at(k-3,1,1) + A73.at(0,2)*p0.at(k-3,2,1) + A74.at(0,0)*p0.at(k-4,0,1) + A74.at(0,1)*p0.at(k-4,1,1) + A74.at(0,2)*p0.at(k-4,2,1);
            p0.at(k,0,2) = A71.at(0,0)*p0.at(k-1,0,2) + A71.at(0,1)*p0.at(k-1,1,2) + A71.at(0,2)*p0.at(k-1,2,2) + A72.at(0,0)*p0.at(k-2,0,2) + A72.at(0,1)*p0.at(k-2,1,2) + A72.at(0,2)*p0.at(k-2,2,2) + A73.at(0,0)*p0.at(k-3,0,2) + A73.at(0,1)*p0.at(k-3,1,2) + A73.at(0,2)*p0.at(k-3,2,2) + A74.at(0,0)*p0.at(k-4,0,2) + A74.at(0,1)*p0.at(k-4,1,2) + A74.at(0,2)*p0.at(k-4,2,2);

            p0.at(k,1,0) = A71.at(1,0)*p0.at(k-1,0,0) + A71.at(1,1)*p0.at(k-1,1,0) + A71.at(1,2)*p0.at(k-1,2,0) + A72.at(1,0)*p0.at(k-2,0,0) + A72.at(1,1)*p0.at(k-2,1,0) + A72.at(1,2)*p0.at(k-2,2,0) + A73.at(1,0)*p0.at(k-3,0,0) + A73.at(1,1)*p0.at(k-3,1,0) + A73.at(1,2)*p0.at(k-3,2,0) + A74.at(1,0)*p0.at(k-4,0,0) + A74.at(1,1)*p0.at(k-4,1,0) + A74.at(1,2)*p0.at(k-4,2,0);
            p0.at(k,1,1) = A71.at(1,0)*p0.at(k-1,0,1) + A71.at(1,1)*p0.at(k-1,1,1) + A71.at(1,2)*p0.at(k-1,2,1) + A72.at(1,0)*p0.at(k-2,0,1) + A72.at(1,1)*p0.at(k-2,1,1) + A72.at(1,2)*p0.at(k-2,2,1) + A73.at(1,0)*p0.at(k-3,0,1) + A73.at(1,1)*p0.at(k-3,1,1) + A73.at(1,2)*p0.at(k-3,2,1) + A74.at(1,0)*p0.at(k-4,0,1) + A74.at(1,1)*p0.at(k-4,1,1) + A74.at(1,2)*p0.at(k-4,2,1);
            p0.at(k,1,2) = A71.at(1,0)*p0.at(k-1,0,2) + A71.at(1,1)*p0.at(k-1,1,2) + A71.at(1,2)*p0.at(k-1,2,2) + A72.at(1,0)*p0.at(k-2,0,2) + A72.at(1,1)*p0.at(k-2,1,2) + A72.at(1,2)*p0.at(k-2,2,2) + A73.at(1,0)*p0.at(k-3,0,2) + A73.at(1,1)*p0.at(k-3,1,2) + A73.at(1,2)*p0.at(k-3,2,2) + A74.at(1,0)*p0.at(k-4,0,2) + A74.at(1,1)*p0.at(k-4,1,2) + A74.at(1,2)*p0.at(k-4,2,2);

            p0.at(k,2,0) = A71.at(2,0)*p0.at(k-1,0,0) + A71.at(2,1)*p0.at(k-1,1,0) + A71.at(2,2)*p0.at(k-1,2,0) + A72.at(2,0)*p0.at(k-2,0,0) + A72.at(2,1)*p0.at(k-2,1,0) + A72.at(2,2)*p0.at(k-2,2,0) + A73.at(2,0)*p0.at(k-3,0,0) + A73.at(2,1)*p0.at(k-3,1,0) + A73.at(2,2)*p0.at(k-3,2,0) + A74.at(2,0)*p0.at(k-4,0,0) + A74.at(2,1)*p0.at(k-4,1,0) + A74.at(2,2)*p0.at(k-4,2,0);
            p0.at(k,2,1) = A71.at(2,0)*p0.at(k-1,0,1) + A71.at(2,1)*p0.at(k-1,1,1) + A71.at(2,2)*p0.at(k-1,2,1) + A72.at(2,0)*p0.at(k-2,0,1) + A72.at(2,1)*p0.at(k-2,1,1) + A72.at(2,2)*p0.at(k-2,2,1) + A73.at(2,0)*p0.at(k-3,0,1) + A73.at(2,1)*p0.at(k-3,1,1) + A73.at(2,2)*p0.at(k-3,2,1) + A74.at(2,0)*p0.at(k-4,0,1) + A74.at(2,1)*p0.at(k-4,1,1) + A74.at(2,2)*p0.at(k-4,2,1);
            p0.at(k,2,2) = A71.at(2,0)*p0.at(k-1,0,2) + A71.at(2,1)*p0.at(k-1,1,2) + A71.at(2,2)*p0.at(k-1,2,2) + A72.at(2,0)*p0.at(k-2,0,2) + A72.at(2,1)*p0.at(k-2,1,2) + A72.at(2,2)*p0.at(k-2,2,2) + A73.at(2,0)*p0.at(k-3,0,2) + A73.at(2,1)*p0.at(k-3,1,2) + A73.at(2,2)*p0.at(k-3,2,2) + A74.at(2,0)*p0.at(k-4,0,2) + A74.at(2,1)*p0.at(k-4,1,2) + A74.at(2,2)*p0.at(k-4,2,2);

            //q0
            q0.at(k,0) = A71.at(0,0)*q0.at(k-1,0) + A71.at(0,1)*q0.at(k-1,1) + A71.at(0,2)*q0.at(k-1,2) + A72.at(0,0)*q0.at(k-2,0) + A72.at(0,1)*q0.at(k-2,1) + A72.at(0,2)*q0.at(k-2,2) + A73.at(0,0)*q0.at(k-3,0) + A73.at(0,1)*q0.at(k-3,1) + A73.at(0,2)*q0.at(k-3,2) + A74.at(0,0)*q0.at(k-4,0) + A74.at(0,1)*q0.at(k-4,1) + A74.at(0,2)*q0.at(k-4,2) + 0.48*h*b(1,k1);
            q0.at(k,1) = A71.at(1,0)*q0.at(k-1,0) + A71.at(1,1)*q0.at(k-1,1) + A71.at(1,2)*q0.at(k-1,2) + A72.at(1,0)*q0.at(k-2,0) + A72.at(1,1)*q0.at(k-2,1) + A72.at(1,2)*q0.at(k-2,2) + A73.at(1,0)*q0.at(k-3,0) + A73.at(1,1)*q0.at(k-3,1) + A73.at(1,2)*q0.at(k-3,2) + A74.at(1,0)*q0.at(k-4,0) + A74.at(1,1)*q0.at(k-4,1) + A74.at(1,2)*q0.at(k-4,2) + 0.48*h*b(2,k1);
            q0.at(k,2) = A71.at(2,0)*q0.at(k-1,0) + A71.at(2,1)*q0.at(k-1,1) + A71.at(2,2)*q0.at(k-1,2) + A72.at(2,0)*q0.at(k-2,0) + A72.at(2,1)*q0.at(k-2,1) + A72.at(2,2)*q0.at(k-2,2) + A73.at(2,0)*q0.at(k-3,0) + A73.at(2,1)*q0.at(k-3,1) + A73.at(2,2)*q0.at(k-3,2) + A74.at(2,0)*q0.at(k-4,0) + A74.at(2,1)*q0.at(k-4,1) + A74.at(2,2)*q0.at(k-4,2) + 0.48*h*b(3,k1);
        }

        x1.at(k) = (p3.at(k,0,0)*x1.at(3)+p3.at(k,0,1)*x2.at(3)+p3.at(k,0,2)*x3.at(3))
                + (p2.at(k,0,0)*x1.at(2)+p2.at(k,0,1)*x2.at(2)+p2.at(k,0,2)*x3.at(2))
                + (p1.at(k,0,0)*x1.at(1)+p1.at(k,0,1)*x2.at(1)+p1.at(k,0,2)*x3.at(1))
                + (p0.at(k,0,0)*x1.at(0)+p0.at(k,0,1)*x2.at(0)+p0.at(k,0,2)*x3.at(0)) + q0.at(k,0);
        x2.at(k) = (p3.at(k,1,0)*x1.at(3)+p3.at(k,1,1)*x2.at(3)+p3.at(k,1,2)*x3.at(3))
                + (p2.at(k,1,0)*x1.at(2)+p2.at(k,1,1)*x2.at(2)+p2.at(k,1,2)*x3.at(2))
                + (p1.at(k,1,0)*x1.at(1)+p1.at(k,1,1)*x2.at(1)+p1.at(k,1,2)*x3.at(1))
                + (p0.at(k,1,0)*x1.at(0)+p0.at(k,1,1)*x2.at(0)+p0.at(k,1,2)*x3.at(0)) + q0.at(k,1);
        x3.at(k) = (p3.at(k,2,0)*x1.at(3)+p3.at(k,2,1)*x2.at(3)+p3.at(k,2,2)*x3.at(3))
                + (p2.at(k,2,0)*x1.at(2)+p2.at(k,2,1)*x2.at(2)+p2.at(k,2,2)*x3.at(2))
                + (p1.at(k,2,0)*x1.at(1)+p1.at(k,2,1)*x2.at(1)+p1.at(k,2,2)*x3.at(1))
                + (p0.at(k,2,0)*x1.at(0)+p0.at(k,2,1)*x2.at(0)+p0.at(k,2,2)*x3.at(0)) + q0.at(k,2);

    }

    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
    puts("---");

    //DoubleMatrix eta(12,1);
    //DoubleVector gamma(12,12);
    //gamma[0,0] = 0.2; gamma[0,1] = 0.5; gamma[0,2] = 1.2; gamma[0,3] = 1.4; gamma[0,4] = 1.1; gamma[0,5] = 1.0; gamma[0,6] = 0.0; gamma[0,7] = 0.0; gamma[0,8] = 1.0; gamma[0,6] = 0.0; gamma[0,7] = 0.0; gamma[0,8] = 1.0;
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
