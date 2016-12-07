#include "example4.h"

double my_rand1()
{
    return ((rand() % 1000) + 1) / 1000.0;
}

void Example4::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Example4 e;
    //e.calculate1();
    e.calculate();
}

Example4::Example4()
{
    h = 0.001;
    N = 1000;
    n = 3;
    K = 4;
}

void Example4::calculate()
{
    DoubleMatrix x0(n,N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x0.at(0,i) = X1(i);
        x0.at(1,i) = X2(i);
        x0.at(2,i) = X3(i);
    }
    IPrinter::printVector(14,10,x0.row(0));
    IPrinter::printVector(14,10,x0.row(1));
    IPrinter::printVector(14,10,x0.row(2));
    puts("---------------------------------------------------------------------------------");

    DoubleMatrix x(n,N+1);

    x.at(0,0) = X1(0); x.at(0,1) = X1(1); x.at(0,2) = X1(2); x.at(0,3) = X1(3);
    x.at(1,0) = X2(0); x.at(1,1) = X2(1); x.at(1,2) = X2(2); x.at(1,3) = X2(3);
    x.at(2,0) = X3(0); x.at(2,1) = X3(1); x.at(2,2) = X3(2); x.at(2,3) = X3(3);

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

        if (k==K+1)
        {
            P3[k] = A[1]*P3[k-1] + A[2];
            P2[k] = A[1]*P2[k-1] + A[3];
            P1[k] = A[1]*P1[k-1] + A[4];
            P0[k] = A[1]*P0[k-1];
            Q[k]  = A[1]*Q[k-1] + A[0];
        }

        if (k==K+2)
        {
            P3[k] = A[1]*P3[k-1] + A[2]*P3[k-2] + A[3];
            P2[k] = A[1]*P2[k-1] + A[2]*P2[k-2] + A[4];
            P1[k] = A[1]*P1[k-1] + A[2]*P1[k-2];
            P0[k] = A[1]*P0[k-1] + A[2]*P0[k-2];
            Q[k]  = A[1]*DoubleMatrix(Q[k-1]) + A[2]*DoubleMatrix(Q[k-2]) + A[0];
        }

        if (k==K+3)
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

        DoubleVector xk = P3[k]*x.col(3) + P2[k]*x.col(2) + P1[k]*x.col(1) + P0[k]*x.col(0) + Q[k];

        x.at(0,k) = xk.at(0);
        x.at(1,k) = xk.at(1);
        x.at(2,k) = xk.at(2);
        xk.clear();

        //        x.at(0,k) = (P3[k].at(0,0)*x.at(0,3)+P3[k].at(0,1)*x.at(1,3)+P3[k].at(0,2)*x.at(2,3))
        //                  + (P2[k].at(0,0)*x.at(0,2)+P2[k].at(0,1)*x.at(1,2)+P2[k].at(0,2)*x.at(2,2))
        //                  + (P1[k].at(0,0)*x.at(0,1)+P1[k].at(0,1)*x.at(1,1)+P1[k].at(0,2)*x.at(2,1))
        //                  + (P0[k].at(0,0)*x.at(0,0)+P0[k].at(0,1)*x.at(1,0)+P0[k].at(0,2)*x.at(2,0)) + Q[k].at(0,0);
        //        x.at(1,k) = (P3[k].at(1,0)*x.at(0,3)+P3[k].at(1,1)*x.at(1,3)+P3[k].at(1,2)*x.at(2,3))
        //                  + (P2[k].at(1,0)*x.at(0,2)+P2[k].at(1,1)*x.at(1,2)+P2[k].at(1,2)*x.at(2,2))
        //                  + (P1[k].at(1,0)*x.at(0,1)+P1[k].at(1,1)*x.at(1,1)+P1[k].at(1,2)*x.at(2,1))
        //                  + (P0[k].at(1,0)*x.at(0,0)+P0[k].at(1,1)*x.at(1,0)+P0[k].at(1,2)*x.at(2,0)) + Q[k].at(1,0);
        //        x.at(2,k) = (P3[k].at(2,0)*x.at(0,3)+P3[k].at(2,1)*x.at(1,3)+P3[k].at(2,2)*x.at(2,3))
        //                  + (P2[k].at(2,0)*x.at(0,2)+P2[k].at(2,1)*x.at(1,2)+P2[k].at(2,2)*x.at(2,2))
        //                  + (P1[k].at(2,0)*x.at(0,1)+P1[k].at(2,1)*x.at(1,1)+P1[k].at(2,2)*x.at(2,1))
        //                  + (P0[k].at(2,0)*x.at(0,0)+P0[k].at(2,1)*x.at(1,0)+P0[k].at(2,2)*x.at(2,0)) + Q[k].at(2,0);

        //        x.at(0,k) = (P3[k].at(0,0)*x.at(0,3)+P3[k].at(0,1)*x.at(1,3)+P3[k].at(0,2)*x.at(2,3))
        //                  + (P2[k].at(0,0)*x.at(0,2)+P2[k].at(0,1)*x.at(1,2)+P2[k].at(0,2)*x.at(2,2))
        //                  + (P1[k].at(0,0)*x.at(0,1)+P1[k].at(0,1)*x.at(1,1)+P1[k].at(0,2)*x.at(2,1))
        //                  + (P0[k].at(0,0)*x.at(0,0)+P0[k].at(0,1)*x.at(1,0)+P0[k].at(0,2)*x.at(2,0)) + Q[k].at(0);
        //        x.at(1,k) = (P3[k].at(1,0)*x.at(0,3)+P3[k].at(1,1)*x.at(1,3)+P3[k].at(1,2)*x.at(2,3))
        //                  + (P2[k].at(1,0)*x.at(0,2)+P2[k].at(1,1)*x.at(1,2)+P2[k].at(1,2)*x.at(2,2))
        //                  + (P1[k].at(1,0)*x.at(0,1)+P1[k].at(1,1)*x.at(1,1)+P1[k].at(1,2)*x.at(2,1))
        //                  + (P0[k].at(1,0)*x.at(0,0)+P0[k].at(1,1)*x.at(1,0)+P0[k].at(1,2)*x.at(2,0)) + Q[k].at(1);
        //        x.at(2,k) = (P3[k].at(2,0)*x.at(0,3)+P3[k].at(2,1)*x.at(1,3)+P3[k].at(2,2)*x.at(2,3))
        //                  + (P2[k].at(2,0)*x.at(0,2)+P2[k].at(2,1)*x.at(1,2)+P2[k].at(2,2)*x.at(2,2))
        //                  + (P1[k].at(2,0)*x.at(0,1)+P1[k].at(2,1)*x.at(1,1)+P1[k].at(2,2)*x.at(2,1))
        //                  + (P0[k].at(2,0)*x.at(0,0)+P0[k].at(2,1)*x.at(1,0)+P0[k].at(2,2)*x.at(2,0)) + Q[k].at(2);
    }

    IPrinter::printVector(14,10,x.row(0));
    IPrinter::printVector(14,10,x.row(1));
    IPrinter::printVector(14,10,x.row(2));
    puts("--------------------------------------------------------------------------------------");


    A[4].clear();
    A[3].clear();
    A[2].clear();
    A[1].clear();
    A[0].clear();
    A.clear();

    DoubleMatrix M(K*n, K*n,0.0);
    DoubleVector b(K*n);

    DoubleVector x00(n); x00.at(0) = X1(0);  x00.at(1) = X2(0);  x00.at(2) = X3(0);
    DoubleVector x01(n); x01.at(0) = X1(1);  x01.at(1) = X2(1);  x01.at(2) = X3(1);
    DoubleVector x02(n); x02.at(0) = X1(2);  x02.at(1) = X2(2);  x02.at(2) = X3(2);
    DoubleVector x03(n); x03.at(0) = X1(3);  x03.at(1) = X2(3);  x03.at(2) = X3(3);

//    DoubleMatrix G000(n,n,0.0); G000.randomData();
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
        DoubleMatrix G100(n,n,0.0); G100.randomData();
        DoubleMatrix G900(n,n,0.0); G900.randomData();
        DoubleMatrix G00N(n,n,0.0); G00N.randomData();

        //DoubleVector x100(n,0.0); x100.at(0) = X1(100); x100.at(1) = X2(100); x100.at(2) = X3(100);
        //DoubleVector x900(n,0.0); x900.at(0) = X1(900); x900.at(1) = X2(900); x900.at(2) = X3(900);

        DoubleMatrix U3 = G100*P3[100] + G900*P3[900] + G00N*P3[N];
        DoubleMatrix U2 = G100*P2[100] + G900*P2[900] + G00N*P2[N];
        DoubleMatrix U1 = G100*P1[100] + G900*P1[900] + G00N*P1[N];
        DoubleMatrix U0 = G100*P0[100] + G900*P0[900] + G00N*P0[N] + G000;
        //DoubleMatrix V0 = G100*Q[100] + G900*Q[900];

        //DoubleVector B1 = G100*x100 + G900*x900;
        //printf("%14.10f %14.10f %14.10f\n", B1.at(0), B1.at(1), B1.at(2));
        DoubleVector B2 = U3*x03 + U2*x02 + U1*x01 + U0*x00;// + V0;
        //printf("%14.10f %14.10f %14.10f\n", B2.at(0), B2.at(1), B2.at(2));

        //puts("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        //IPrinter::print(G100,n,n);
        //puts("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        //IPrinter::print(U2,n,n);
        //puts("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

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
        //IPrinter::print(M,12,12);
        //puts("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    }

    {
        DoubleMatrix G000(n,n,0.0); G000.randomData();
        DoubleMatrix G200(n,n,0.0); G200.randomData();
        DoubleMatrix G300(n,n,0.0); G300.randomData();
        DoubleMatrix G00N(n,n,0.0); G00N.randomData();

        //DoubleVector x200(n); x200.at(0) = X1(200); x200.at(1) = X2(200); x200.at(2) = X3(200);
        //DoubleVector x300(n); x300.at(0) = X1(300); x300.at(1) = X2(300); x300.at(2) = X3(300);

        DoubleMatrix U3 = G200*P3[200] + G300*P3[300] + G00N*P3[N];
        DoubleMatrix U2 = G200*P2[200] + G300*P2[300] + G00N*P2[N];
        DoubleMatrix U1 = G200*P1[200] + G300*P1[300] + G00N*P1[N];
        DoubleMatrix U0 = G200*P0[200] + G300*P0[300] + G00N*P0[N] + G000;
        //DoubleMatrix V0 = G200*Q[200]  + G300*Q[300];

        //DoubleVector B1 = G200*x200 + G300*x300;
        //printf("%14.10f %14.10f %14.10f\n", B1.at(0), B1.at(1), B1.at(2));
        DoubleVector B2 = U3*x03 + U2*x02 + U1*x01 + U0*x00;// + V0;
        //printf("%14.10f %14.10f %14.10f\n", B2.at(0), B2.at(1), B2.at(2));

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

        //DoubleVector x800(n); x800.at(0) = X1(800); x800.at(1) = X2(800); x800.at(2) = X3(800);
        //DoubleVector x900(n); x900.at(0) = X1(900); x900.at(1) = X2(900); x900.at(2) = X3(900);

        DoubleMatrix U3 = G800*P3[800] + G900*P3[900] + G00N*P3[N];
        DoubleMatrix U2 = G800*P2[800] + G900*P2[900] + G00N*P2[N];
        DoubleMatrix U1 = G800*P1[800] + G900*P1[900] + G00N*P1[N];
        DoubleMatrix U0 = G800*P0[800] + G900*P0[900] + G00N*P0[N] + G000;
        //DoubleMatrix V0 = G800*Q[800]  + G900*Q[900];

        //DoubleVector B1 = G800*x800 + G900*x900;
        //printf("%14.10f %14.10f %14.10f\n", B1.at(0), B1.at(1), B1.at(2));
        DoubleVector B2 = U3*x03 + U2*x02 + U1*x01 + U0*x00;// + V0;
        //printf("%14.10f %14.10f %14.10f\n", B2.at(0), B2.at(1), B2.at(2));

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
        DoubleMatrix G250(n,n,0.0); G250.randomData();
        DoubleMatrix G450(n,n,0.0); G450.randomData();
        DoubleMatrix G00N(n,n,0.0); G00N.randomData();

        //DoubleVector xN1(n); xN1.at(0) = X1(N-1); xN1.at(1) = X2(N-1); xN1.at(2) = X3(N-1);
        //DoubleVector xN0(n); xN0.at(0) = X1(N-0); xN0.at(1) = X2(N-0); xN0.at(2) = X3(N-0);

        DoubleMatrix U3 = G250*P3[250] + G450*P3[450] + G00N*P3[N];
        DoubleMatrix U2 = G250*P2[250] + G450*P2[450] + G00N*P2[N];
        DoubleMatrix U1 = G250*P1[250] + G450*P1[450] + G00N*P1[N];
        DoubleMatrix U0 = G250*P0[250] + G450*P0[450] + G00N*P0[N] + G000;
        //DoubleMatrix V0 = G0N1*Q[N-1]  + G0N0*Q[N-0];

        //DoubleVector B1 = GN1*xN1 + GN0*xN0;
        //printf("%14.10f %14.10f %14.10f\n", B1.at(0), B1.at(1), B1.at(2));
        DoubleVector B2 = U3*x03 + U2*x02 + U1*x01 + U0*x00;// + V0;
        //printf("%14.10f %14.10f %14.10f\n", B2.at(0), B2.at(1), B2.at(2));

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
        //IPrinter::print(M,M.rows(),M.cols());
    }

    puts("+++");
    for (unsigned int i=0; i<M.rows(); i++)
    {
        for (unsigned int j=0; j<M.cols(); j++)
        {
            printf("%10.6f ",M[i][j]);
        }
        printf("%14.10f\n", b[i]);
    }
    puts("+++");

    DoubleVector xx(12);
    printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11]);
//    GaussianElimination(M, b, xx);

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

    DoubleMatrix M1 = M;
    DoubleVector b1 = b;
    //M1.randomData();
    puts("..............................................................");
    IPrinter::print(M1,12,12,18,14);
    for (unsigned int k=0; k<11; k++)
    {
        for (unsigned int j=k+1; j<12; j++)
        {
            printf("%d %d\n", k,j);
            double c = -M1.at(j,k)/M1.at(k,k);
            for (unsigned int i=k; i<12; i++) M1.at(j, i) = M1.at(j, i) + M1.at(k,i)*c;
            b1[j] = b1[j] + b1[k] *c;
            puts("..............................................................");
            system("pause");
            IPrinter::print(M1,12,12,18,14);
        }
    }
    puts("..............................................................");
    printf("%18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f\n", b1[0], b1[1], b1[2], b1[3], b1[4], b1[5], b1[6], b1[7], b1[8], b1[9], b1[10], b1[11]);

//    xx[11] = (b1[11]) / M1.at(11,11);
//    xx[10] = (b1[10] - M1.at(10,11)*xx[11]) / M1.at(10,10);
//    xx[9]  = (b1[9] - M1.at(9,11)*xx[11] - M1.at(9,10)*xx[10]) / M1.at(9,9);
//    xx[8]  = (b1[8] - /*M1.at(8,11)*xx[11] - M1.at(8,10)*xx[10] -*/ M1.at(8,9)*xx[9]) / M1.at(8,8);
//    xx[7]  = (b1[7] - M1.at(7,11)*xx[11] - M1.at(7,10)*xx[10] - M1.at(7,9)*xx[9] - M1.at(7,8)*xx[8]) / M1.at(7,7);


    for (unsigned int i=11; i!=UINT32_MAX; i--)
    {
        for (unsigned int j=11; j>i; j--) b1[i] -= (M1.at(i,j) * xx[j]);
        xx[i] = b1[i] / M1.at(i,i);
        printf("%d %14.10f\n", i, xx[i]);
    }

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

    puts("---------------------------------------------------");
    printf("%18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f\n",
           xx[0], xx[1], xx[2], xx[3], xx[4], xx[5], xx[6], xx[7], xx[8], xx[9], xx[10], xx[11]);
    printf("%18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f %18.14f\n",
           X1(3), X2(3), X3(3), X1(2), X2(2), X3(2), X1(1), X2(1), X3(1), X1(0), X2(0), X3(0));

    puts("---------------------------------------------------");

    printf("%14.10f %14.10f %14.10f\n", X1(3), X2(3), X3(3));
    printf("%14.10f %14.10f %14.10f\n", X1(2), X2(2), X3(2));
    printf("%14.10f %14.10f %14.10f\n", X1(1), X2(1), X3(1));
    printf("%14.10f %14.10f %14.10f\n", X1(0), X2(0), X3(0));

    puts("---------------------------------------------------");

    printf("%14.10f %14.10f %14.10f\n", xx[0], xx[1], xx[2]);
    printf("%14.10f %14.10f %14.10f\n", xx[3], xx[4], xx[5]);
    printf("%14.10f %14.10f %14.10f\n", xx[6], xx[7], xx[8]);
    printf("%14.10f %14.10f %14.10f\n", xx[9], xx[10], xx[11]);

    puts("---------------------------------------------------");

    for (unsigned int i=0; i<12; i++)
    {
        double a1 = M.at(i,0)*X1(3)+M.at(i,1)*X2(3)+M.at(i,2)*X3(3)+M.at(i,3)*X1(2)+M.at(i,4)*X2(2)+M.at(i,5)*X3(2)+M.at(i,6)*X1(1)+M.at(i,7)*X2(1)+M.at(i,8)*X3(1)+M.at(i,9)*X1(0)+M.at(i,10)*X2(0)+M.at(i,11)*X3(0);
        double a2 = M.at(i,0)*xx[0]+M.at(i,1)*xx[1]+M.at(i,2)*xx[2]+M.at(i,3)*xx[3]+M.at(i,4)*xx[4]+M.at(i,5)*xx[5]+M.at(i,6)*xx[6]+M.at(i,7)*xx[7]+M.at(i,8)*xx[8]+M.at(i,9)*xx[9]+M.at(i,10)*xx[10]+M.at(i,11)*xx[11];
        printf("%14.10f %14.10f\n", a1, a2);
    }

    puts("---------------------------------------------------");

    //    printf("%14.10f %14.10f %14.10f\n", x.at(0), x[1], x[2]);
    //    printf("%14.10f %14.10f %14.10f\n", x.at(3), x[4], x[5]);
    //    printf("%14.10f %14.10f %14.10f\n", x.at(6), x[7], x[8]);
    //    printf("%14.10f %14.10f %14.10f\n", x.at(9), x[10], x[11]);

    //x.clear();

    ///////////////////////////////////////////////////////////

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

void Example4::calculate1()
{
    DoubleVector x01(N+1);
    DoubleVector x02(N+1);
    DoubleVector x03(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x01.at(i) = X1(i);
        x02.at(i) = X2(i);
        x03.at(i) = X3(i);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    IPrinter::printVector(14,10,x03);

    DoubleVector x1(N+1);
    DoubleVector x2(N+2);
    DoubleVector x3(N+2);

    x1.at(0) = X1(0); x2.at(0) = X2(0); x3.at(0) = X3(0);
    x1.at(1) = X1(1); x2.at(1) = X2(1); x3.at(1) = X3(1);
    x1.at(2) = X1(2); x2.at(2) = X2(2); x3.at(2) = X3(2);
    x1.at(3) = X1(3); x2.at(3) = X2(3); x3.at(3) = X3(3);

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
        x01.at(i) = X1(i);
        x02.at(i) = X2(i);
        x03.at(i) = X3(i);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    IPrinter::printVector(14,10,x03);

    DoubleVector x1(N+1);
    DoubleVector x2(N+2);
    DoubleVector x3(N+2);

    x1.at(N-0) = X1(N-0); x2.at(N-0) = X2(N-0); x3.at(N-0) = X3(N-0);
    x1.at(N-1) = X1(N-1); x2.at(N-1) = X2(N-1); x3.at(N-1) = X3(N-1);
    x1.at(N-2) = X1(N-2); x2.at(N-2) = X2(N-2); x3.at(N-2) = X3(N-2);
    x1.at(N-3) = X1(N-3); x2.at(N-3) = X2(N-3); x3.at(N-3) = X3(N-3);

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
        x01.at(i) = X1(i);
        x02.at(i) = X2(i);
        x03.at(i) = X3(i);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    IPrinter::printVector(14,10,x03);
    puts("---");

    DoubleVector x1(N+1);
    DoubleVector x2(N+1);
    DoubleVector x3(N+1);

    x1.at(0) = X1(0); x1.at(1) = X1(1); x1.at(2) = X1(2); x1.at(3) = X1(3);
    x2.at(0) = X2(0); x2.at(1) = X2(1); x2.at(2) = X2(2); x2.at(3) = X2(3);
    x3.at(0) = X3(0); x3.at(1) = X3(1); x3.at(2) = X3(2); x3.at(3) = X3(3);

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

double Example4::X1(unsigned int k) const
{
    double t = k*h;
    return sin(2.0*t) + t*t;
}

double Example4::X2(unsigned int k) const
{
    double t = k*h;
    return 3.0*t;
}

double Example4::X3(unsigned int k) const
{
    double t = k*h;
    return cos(2.0*t) - sin(t);
}

double Example4::a(unsigned int i, unsigned int j, unsigned int k) const
{
    double t = k*h;

    if (i==1 && j==1) return +3.0;
    if (i==1 && j==2) return -t;
    if (i==1 && j==3) return +2.0;

    if (i==2 && j==1) return +1.0;
    if (i==2 && j==2) return +1.0;
    if (i==2 && j==3) return +3.0;

    if (i==3 && j==1) return -2.0;
    if (i==3 && j==2) return +t;
    if (i==3 && j==3) return +1.0;

    return 0.0;
}

double Example4::b(unsigned int i, unsigned int k) const
{
    double t = k*h;
    if (i==1) return 2.0*t + 2.0*sin(t) - 3.0*sin(2.0*t);
    if (i==2) return 3.0*sin(t) - sin(2.0*t) - 3.0*cos(2.0*t) - t*t - 3.0*t + 3.0;
    if (i==3) return sin(t) - cos(t) - cos(2.0*t) - t*t;
    return 0.0;
}
