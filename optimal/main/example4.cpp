#include "example4.h"

void Example4::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Example4 e;
    e.h = 0.001;
    e.N = 1000;
    e.F = e.N/10;
    e.n = 3;
    e.K = 4;
    e.w = 12;
    e.p = 8;
    //e.L = 5;

    std::vector<unsigned int> *s = new std::vector<unsigned int>[e.K];
    s[0].push_back(0); s[0].push_back(2*e.F); s[0].push_back(5*e.F); s[0].push_back(8*e.F); s[0].push_back(10*e.F);
    s[1].push_back(0); s[1].push_back(1);     s[1].push_back(2);     s[1].push_back(3);     s[1].push_back(4);
    s[2].push_back(0); s[2].push_back(1);     s[2].push_back(2);     s[2].push_back(3);     s[2].push_back(4);
    s[3].push_back(0); s[3].push_back(1);     s[3].push_back(2);     s[3].push_back(3);     s[3].push_back(4);

    //s[1].push_back(0); s[1].push_back(1*e.F); s[1].push_back(2*e.F); s[1].push_back(3*e.F); s[1].push_back(10*e.F);
    //s[2].push_back(0); s[2].push_back(3*e.F); s[2].push_back(4*e.F); s[2].push_back(8*e.F); s[2].push_back(10*e.F);
    //s[3].push_back(0); s[3].push_back(5*e.F); s[3].push_back(7*e.F); s[3].push_back(9*e.F); s[3].push_back(10*e.F);

    IPrinter::printSeperatorLine(NULL,'-',stdout);
    DoubleVector x(e.n*e.K);
    for (unsigned int i=0; i<e.K; i++)
    {
        for (unsigned int j=0; j<e.n; j++)
        {
            x.at(i*e.n+j) = e.fx(j+1,i);
        }
    }
    IPrinter::print(x,x.size(),e.w,e.p,stdout);
    IPrinter::printSeperatorLine(NULL,'-',stdout);

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

    puts("Method #1");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
    e.calculateM1(s, rx, nx);
    //--------------------------------------------------------------------------
    puts("Method #2");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
    e.calculateM2(s, rx, nx);
}

Example4::Example4()
{
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

    for (unsigned int i=0; i<K; i++)
    {
        nx.at(i) = rx.at(i);
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

void Example4::calculateM1(const std::vector<unsigned int> *s, const DoubleMatrix &rx UNUSED_PARAM, const DoubleMatrix &nx UNUSED_PARAM)
{
    DoubleMatrix M(K*n, K*n, 0.0);
    DoubleVector B(K*n,0.0);

    calculateM1BE(0,s[0],nx,M,B);
    calculateM1BE(1,s[1],nx,M,B);
    calculateM1BE(2,s[2],nx,M,B);
    calculateM1BE(3,s[3],nx,M,B);

    printf("det: %14.10f\n",M.determinant1());

    DoubleVector x(n*K);
    GaussianElimination(M,B,x);

    IPrinter::printSeperatorLine(NULL,'-',stdout);
    IPrinter::print(x,x.size(),w,p,stdout);
    IPrinter::printSeperatorLine(NULL,'-',stdout);

    DoubleMatrix cx(n,K);
    cx.at(0,0) = x[0]; cx.at(0,1) = x[3]; cx.at(0,2) = x[6]; cx.at(0,3) = x[9];
    cx.at(1,0) = x[1]; cx.at(1,1) = x[4]; cx.at(1,2) = x[7]; cx.at(1,3) = x[10];
    cx.at(2,0) = x[2]; cx.at(2,1) = x[5]; cx.at(2,2) = x[8]; cx.at(2,3) = x[11];

    DoubleMatrix nx1;
    calculateNX(cx,nx1);
    IPrinter::printVector(w,p,nx1.row(0),"x1: ");
    IPrinter::printVector(w,p,nx1.row(1),"x2: ");
    IPrinter::printVector(w,p,nx1.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
}

void Example4::calculateM2(const std::vector<unsigned int> *s, const DoubleMatrix &rx UNUSED_PARAM, const DoubleMatrix &nx UNUSED_PARAM)
{
    /* real solution vectors */
    std::vector<stdDoubleMatrixVector> P(K);
    for (unsigned int i=0; i<K; i++) P[i].resize(N+1);

//    stdDoubleMatrixVector P3(N+1);
//    stdDoubleMatrixVector P2(N+1);
//    stdDoubleMatrixVector P1(N+1);
//    stdDoubleMatrixVector P0(N+1);
    stdDoubleMatrixVector Q(N+1);
//    calculatePQ(P3,P2,P1,P0,Q);
    calculatePQ(P, Q);

    /* numerical solution vectors */
    DoubleMatrix nx1;
//    calculateNS(nx1,rx,P3,P2,P1,P0,Q);
    calculateNS(nx1,rx,P,Q);
    IPrinter::printVector(w,p,nx1.row(0),"x1: ");
    IPrinter::printVector(w,p,nx1.row(1),"x2: ");
    IPrinter::printVector(w,p,nx1.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
    nx1.clear();

    /* find x0, x1, x2, x3 */

    DoubleMatrix M(K*n, K*n, 0.0);
    DoubleVector B(K*n,0.0);

//    calculateM2BE(0,s[0],nx,M,B,P3,P2,P1,P0,Q);
//    calculateM2BE(1,s[1],nx,M,B,P3,P2,P1,P0,Q);
//    calculateM2BE(2,s[2],nx,M,B,P3,P2,P1,P0,Q);
//    calculateM2BE(3,s[3],nx,M,B,P3,P2,P1,P0,Q);
    calculateM2BE(0,s[0],nx,M,B,P,Q);
    calculateM2BE(1,s[1],nx,M,B,P,Q);
    calculateM2BE(2,s[2],nx,M,B,P,Q);
    calculateM2BE(3,s[3],nx,M,B,P,Q);

    DoubleVector x(n*K);
    GaussianElimination(M, B, x);

    IPrinter::printSeperatorLine(NULL,'-',stdout);
    IPrinter::print(x,x.size(),w,p,stdout);
    IPrinter::printSeperatorLine(NULL,'-',stdout);

    DoubleMatrix rx1(n,K);
    rx1.at(0,0) = x[0]; rx1.at(0,1) = x[3]; rx1.at(0,2) = x[6]; rx1.at(0,3) = x[9];
    rx1.at(1,0) = x[1]; rx1.at(1,1) = x[4]; rx1.at(1,2) = x[7]; rx1.at(1,3) = x[10];
    rx1.at(2,0) = x[2]; rx1.at(2,1) = x[5]; rx1.at(2,2) = x[8]; rx1.at(2,3) = x[11];

    DoubleMatrix cx;
//    calculateNS(cx,rx1,P3,P2,P1,P0,Q);
    calculateNS(cx,rx1,P,Q);
    IPrinter::printVector(w,p,cx.row(0),"x1: ");
    IPrinter::printVector(w,p,cx.row(1),"x2: ");
    IPrinter::printVector(w,p,cx.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
    cx.clear();

    for (unsigned int i=0; i<=N; i++)
    {
        P[3][i].clear();
        P[2][i].clear();
        P[1][i].clear();
        P[0][i].clear();
        Q[i].clear();
    }
    P[3].clear();
    P[2].clear();
    P[1].clear();
    P[0].clear();
    Q.clear();
}

void Example4::calculateM1BE(unsigned int c, const std::vector<unsigned int> s, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B)
{
    unsigned int L = s.size();
    std::vector<DoubleMatrix> GAMMA(L);
    DoubleVector ETA(n,0.0);

    std::vector<DoubleMatrix> betta(N+1);

    //    for (unsigned int i=0; i<L; i++)
    //    {
    //        GAMMA[i].resize(n,n,0.0);
    //        GAMMA[i].randomData();
    //        ETA = GAMMA[i]*nx.col(s[i]) + ETA;
    //    }

    fillGamma(GAMMA, ETA, c, K);
    if (c == 0) { for (unsigned int i=0; i<L; i++) ETA = GAMMA[i]*nx.col(s[i]) + ETA; }

    if (c==0)
    {
        betta[N-0] = GAMMA[L-1];
        betta[N-1].resize(n,n,0.0);
        betta[N-2].resize(n,n,0.0);
        betta[N-3].resize(n,n,0.0);

        std::vector<DoubleMatrix> A;
        initAMatrices(A);
        for (unsigned int k=N; k>=K; k--)
        {
            updateAMatrices(A,k);
            betta[k-1] = betta[k]*A[1] + betta[k-1];
            betta[k-2] = betta[k]*A[2] + betta[k-2];
            betta[k-3] = betta[k]*A[3] + betta[k-3];
            betta[k-4] = betta[k]*A[4];// + betta[k-4];
            ETA        = ETA - betta[k]*A[0];

            for (unsigned int i=0; i<L; i++)
            {
                if (k==(s[i]+K))
                {
                    betta[k-4] = betta[k-4] + GAMMA[i];
                }
            }
        }
        clearAMatrices(A);
    }
    else
    {
        betta[s[L-1]-0] = GAMMA[L-1];
        betta[s[L-1]-1] = GAMMA[L-2];
        betta[s[L-1]-2] = GAMMA[L-3];
        betta[s[L-1]-3] = GAMMA[L-4];
        //betta[s[L-1]-4] = GAMMA[L-5];

        std::vector<DoubleMatrix> A;
        initAMatrices(A);
        for (unsigned int k=s[L-1]; k>=K; k--)
        {
            updateAMatrices(A,k);
            betta[k-1] = betta[k]*A[1] + betta[k-1];
            betta[k-2] = betta[k]*A[2] + betta[k-2];
            betta[k-3] = betta[k]*A[3] + betta[k-3];
            betta[k-4] = betta[k]*A[4];// + betta[k-4];
            ETA        = ETA - betta[k]*A[0];

            for (unsigned int i=0; i<L; i++)
            {
                if (k==(s[i]+K))
                {
                    betta[k-4] = betta[k-4] + GAMMA[i];
                }
            }
        }
        clearAMatrices(A);
    }

    for (unsigned int i=0; i<n; i++)
    {
        for (unsigned int j=0; j<n; j++)
        {
            M[c*n+i][0*n+j] = betta[0][i][j];
            M[c*n+i][1*n+j] = betta[1][i][j];
            M[c*n+i][2*n+j] = betta[2][i][j];
            M[c*n+i][3*n+j] = betta[3][i][j];
        }
        B.at(c*n+i) = ETA.at(i);
    }

    ETA.clear();
    for (unsigned int i=0; i<L; i++) GAMMA[i].clear();
    GAMMA.clear();
    for (unsigned int i=0; i<=N; i++) betta[i].clear();
    betta.clear();
}

void Example4::calculateM2BE(unsigned int c, const std::vector<unsigned int> s, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B,
                             const std::vector<stdDoubleMatrixVector> &P, const stdDoubleMatrixVector &Q)
{
    unsigned int L = s.size();
    std::vector<DoubleMatrix> GAMMA(L);
    DoubleVector B1(n,0.0);

    //    for (unsigned int i=0; i<L; i++)
    //    {
    //        GAMMA[i].resize(n,n,0.0);
    //        GAMMA[i].randomData();
    //        B1 = GAMMA[i]*nx.col(s[i]) + B1;
    //    }

    fillGamma(GAMMA, B1, c, K);
    if (c == 0) { for (unsigned int i=0; i<L; i++) B1 = GAMMA[i]*nx.col(s[i]) + B1; }

    DoubleMatrix U3(n,n,0.0);
    DoubleMatrix U2(n,n,0.0);
    DoubleMatrix U1(n,n,0.0);
    DoubleMatrix U0(n,n,0.0);
    DoubleMatrix V0(n,1,0.0);

    for (unsigned int i=0; i<L; i++)
    {
        if (s[i] == 3) U3 = U3 + GAMMA[i];
        else if (s[i] == 2) U2 = U2 + GAMMA[i];
        else if (s[i] == 1) U1 = U1 + GAMMA[i];
        else if (s[i] == 0) U0 = U0 + GAMMA[i];
        else
        {
            U3 = U3 + GAMMA[i]*P[3][s[i]];
            U2 = U2 + GAMMA[i]*P[2][s[i]];
            U1 = U1 + GAMMA[i]*P[1][s[i]];
            U0 = U0 + GAMMA[i]*P[0][s[i]];
            V0 = V0 + GAMMA[i]*Q[s[i]];
        }
    }
    B1 = B1 - V0;

    for (unsigned int i=0; i<n; i++)
    {
        for (unsigned int j=0; j<n; j++)
        {
            M[c*n+i][0*n+j] = U0[i][j];
            M[c*n+i][1*n+j] = U1[i][j];
            M[c*n+i][2*n+j] = U2[i][j];
            M[c*n+i][3*n+j] = U3[i][j];
        }
        B.at(c*n+i) = B1.at(i);
    }
}

//void Example4::calculateM2BE(unsigned int c, const std::vector<unsigned int> s, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B,
//                             const stdDoubleMatrixVector &P3, const stdDoubleMatrixVector &P2,
//                             const stdDoubleMatrixVector &P1, const stdDoubleMatrixVector &P0,
//                             const stdDoubleMatrixVector &Q)
//{
//    unsigned int L = s.size();
//    std::vector<DoubleMatrix> GAMMA(L);
//    DoubleVector B1(n,0.0);

//    //    for (unsigned int i=0; i<L; i++)
//    //    {
//    //        GAMMA[i].resize(n,n,0.0);
//    //        GAMMA[i].randomData();
//    //        B1 = GAMMA[i]*nx.col(s[i]) + B1;
//    //    }

//    fillGamma(GAMMA, B1, c, K);
//    if (c == 0) { for (unsigned int i=0; i<L; i++) B1 = GAMMA[i]*nx.col(s[i]) + B1; }

//    DoubleMatrix U3(n,n,0.0);
//    DoubleMatrix U2(n,n,0.0);
//    DoubleMatrix U1(n,n,0.0);
//    DoubleMatrix U0(n,n,0.0);
//    DoubleMatrix V0(n,1,0.0);

//    for (unsigned int i=0; i<L; i++)
//    {
//        if (s[i] == 3) U3 = U3 + GAMMA[i];
//        else if (s[i] == 2) U2 = U2 + GAMMA[i];
//        else if (s[i] == 1) U1 = U1 + GAMMA[i];
//        else if (s[i] == 0) U0 = U0 + GAMMA[i];
//        else
//        {
//            U3 = U3 + GAMMA[i]*P3[s[i]];
//            U2 = U2 + GAMMA[i]*P2[s[i]];
//            U1 = U1 + GAMMA[i]*P1[s[i]];
//            U0 = U0 + GAMMA[i]*P0[s[i]];
//            V0 = V0 + GAMMA[i]*Q[s[i]];
//        }
//    }
//    B1 = B1 - V0;

//    for (unsigned int i=0; i<n; i++)
//    {
//        for (unsigned int j=0; j<n; j++)
//        {
//            M[c*n+i][0*n+j] = U0[i][j];
//            M[c*n+i][1*n+j] = U1[i][j];
//            M[c*n+i][2*n+j] = U2[i][j];
//            M[c*n+i][3*n+j] = U3[i][j];
//        }
//        B.at(c*n+i) = B1.at(i);
//    }
//}

void Example4::calculateRX(DoubleMatrix &rx)
{
    if (rx.empty())
    {
        rx.resize(n,N+1,0.0);
        for (unsigned int i=0; i<=N; i++)
        {
            rx.at(0,i) = fx(1,i);
            rx.at(1,i) = fx(2,i);
            rx.at(2,i) = fx(3,i);
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

void Example4::calculateNS(DoubleMatrix &nx, const DoubleMatrix &rx, const std::vector<stdDoubleMatrixVector> &P, const stdDoubleMatrixVector &Q)
{
    if (nx.empty())
    {
        nx.resize(n,N+1,0.0);

        nx.setColumn(0, rx.col(0));
        nx.setColumn(1, rx.col(1));
        nx.setColumn(2, rx.col(2));
        nx.setColumn(3, rx.col(3));

        for (unsigned int k=K; k<=N; k++)
        {
            DoubleVector ck = P[3][k]*nx.col(3) + P[2][k]*nx.col(2) + P[1][k]*nx.col(1) + P[0][k]*nx.col(0) + Q[k];
            nx.setColumn(k,ck);
        }
    }
}

//void Example4::calculateNS(DoubleMatrix &nx, const DoubleMatrix &rx, const stdDoubleMatrixVector &P3, const stdDoubleMatrixVector &P2,
//                           const stdDoubleMatrixVector &P1, const stdDoubleMatrixVector &P0, const stdDoubleMatrixVector &Q)
//{
//    if (nx.empty())
//    {
//        nx.resize(n,N+1,0.0);

//        nx.setColumn(0, rx.col(0));
//        nx.setColumn(1, rx.col(1));
//        nx.setColumn(2, rx.col(2));
//        nx.setColumn(3, rx.col(3));

//        for (unsigned int k=K; k<=N; k++)
//        {
//            DoubleVector ck = P3[k]*nx.col(3) + P2[k]*nx.col(2) + P1[k]*nx.col(1) + P0[k]*nx.col(0) + Q[k];
//            nx.setColumn(k,ck);
//        }
//    }
//}

//void Example4::calculatePQ(stdDoubleMatrixVector &P3, stdDoubleMatrixVector &P2, stdDoubleMatrixVector &P1, stdDoubleMatrixVector &P0, stdDoubleMatrixVector &Q)
//{
//    /* calculating P,Q matrices */
//    stdDoubleMatrixVector A;
//    initAMatrices(A);
//    for (unsigned int k=K; k<=N; k++)
//    {
//        updateAMatrices(A,k);

//        if (k==K)
//        {
//            P3[k] = A[1];
//            P2[k] = A[2];
//            P1[k] = A[3];
//            P0[k] = A[4];
//            Q[k]  = A[0];
//        }
//        else if (k==K+1)
//        {
//            P3[k] = A[1]*P3[k-1] + A[2];
//            P2[k] = A[1]*P2[k-1] + A[3];
//            P1[k] = A[1]*P1[k-1] + A[4];
//            P0[k] = A[1]*P0[k-1];
//            Q[k]  = A[1]*Q[k-1] + A[0];
//        }
//        else if (k==K+2)
//        {
//            P3[k] = A[1]*P3[k-1] + A[2]*P3[k-2] + A[3];
//            P2[k] = A[1]*P2[k-1] + A[2]*P2[k-2] + A[4];
//            P1[k] = A[1]*P1[k-1] + A[2]*P1[k-2];
//            P0[k] = A[1]*P0[k-1] + A[2]*P0[k-2];
//            Q[k]  = A[1]*Q[k-1] + A[2]*Q[k-2] + A[0];
//        }
//        else if (k==K+3)
//        {
//            P3[k] = A[1]*P3[k-1] + A[2]*P3[k-2] + A[3]*P3[k-3] + A[4];
//            P2[k] = A[1]*P2[k-1] + A[2]*P2[k-2] + A[3]*P2[k-3];
//            P1[k] = A[1]*P1[k-1] + A[2]*P1[k-2] + A[3]*P1[k-3];
//            P0[k] = A[1]*P0[k-1] + A[2]*P0[k-2] + A[3]*P0[k-3];
//            Q[k]  = A[1]*Q[k-1] + A[2]*Q[k-2] + A[3]*Q[k-3] + A[0];
//        }
//        if (k>=2*K)
//        {
//            P3[k] = A[1]*P3[k-1] + A[2]*P3[k-2] + A[3]*P3[k-3] + A[4]*P3[k-4];
//            P2[k] = A[1]*P2[k-1] + A[2]*P2[k-2] + A[3]*P2[k-3] + A[4]*P2[k-4];
//            P1[k] = A[1]*P1[k-1] + A[2]*P1[k-2] + A[3]*P1[k-3] + A[4]*P1[k-4];
//            P0[k] = A[1]*P0[k-1] + A[2]*P0[k-2] + A[3]*P0[k-3] + A[4]*P0[k-4];
//            Q[k]  = A[1]*Q[k-1] + A[2]*Q[k-2] + A[3]*Q[k-3] + A[4]*Q[k-4] + A[0];
//        }
//    }
//    A[4].clear();
//    A[3].clear();
//    A[2].clear();
//    A[1].clear();
//    A[0].clear();
//    A.clear();
//    /* calculating P,Q matrices */
//}

void Example4::calculatePQ(std::vector<stdDoubleMatrixVector> &P, stdDoubleMatrixVector &Q)
{
    /* calculating P,Q matrices */
    stdDoubleMatrixVector A;
    initAMatrices(A);
    for (unsigned int k=K; k<=N; k++)
    {
        updateAMatrices(A,k);

        if (k==K)
        {
            P[3][k] = A[1];
            P[2][k] = A[2];
            P[1][k] = A[3];
            P[0][k] = A[4];
            Q[k]  = A[0];
        }
        else if (k==K+1)
        {
            P[3][k] = A[1]*P[3][k-1] + A[2];
            P[2][k] = A[1]*P[2][k-1] + A[3];
            P[1][k] = A[1]*P[1][k-1] + A[4];
            P[0][k] = A[1]*P[0][k-1];
            Q[k]  = A[1]*Q[k-1] + A[0];
        }
        else if (k==K+2)
        {
            P[3][k] = A[1]*P[3][k-1] + A[2]*P[3][k-2] + A[3];
            P[2][k] = A[1]*P[2][k-1] + A[2]*P[2][k-2] + A[4];
            P[1][k] = A[1]*P[1][k-1] + A[2]*P[1][k-2];
            P[0][k] = A[1]*P[0][k-1] + A[2]*P[0][k-2];
            Q[k]  = A[1]*Q[k-1] + A[2]*Q[k-2] + A[0];
        }
        else if (k==K+3)
        {
            P[3][k] = A[1]*P[3][k-1] + A[2]*P[3][k-2] + A[3]*P[3][k-3] + A[4];
            P[2][k] = A[1]*P[2][k-1] + A[2]*P[2][k-2] + A[3]*P[2][k-3];
            P[1][k] = A[1]*P[1][k-1] + A[2]*P[1][k-2] + A[3]*P[1][k-3];
            P[0][k] = A[1]*P[0][k-1] + A[2]*P[0][k-2] + A[3]*P[0][k-3];
            Q[k]  = A[1]*Q[k-1] + A[2]*Q[k-2] + A[3]*Q[k-3] + A[0];
        }
        if (k>=2*K)
        {
            P[3][k] = A[1]*P[3][k-1] + A[2]*P[3][k-2] + A[3]*P[3][k-3] + A[4]*P[3][k-4];
            P[2][k] = A[1]*P[2][k-1] + A[2]*P[2][k-2] + A[3]*P[2][k-3] + A[4]*P[2][k-4];
            P[1][k] = A[1]*P[1][k-1] + A[2]*P[1][k-2] + A[3]*P[1][k-3] + A[4]*P[1][k-4];
            P[0][k] = A[1]*P[0][k-1] + A[2]*P[0][k-2] + A[3]*P[0][k-3] + A[4]*P[0][k-4];
            Q[k]  = A[1]*Q[k-1] + A[2]*Q[k-2] + A[3]*Q[k-3] + A[4]*Q[k-4] + A[0];
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

void Example4::fillGamma(std::vector<DoubleMatrix> &GAMMA, DoubleVector &ETA, unsigned int s, unsigned int k)
{
    unsigned int L = GAMMA.size();
    if (k == K)
    {
        if (s == 0)
        {
            for (unsigned int i=0; i<L; i++)
            {
                GAMMA[i].resize(n,n,1.0);
                //GAMMA[i].randomData();
                for (unsigned int i1=0; i1<n; i1++)
                {
                    for (unsigned int i2=0; i2<n; i2++)
                    {
                        GAMMA[i].at(i1,i2) = sin(i1+i)+cos(i2+i);
                    }
                }
                //eta = eta + DoubleVector(GAMMA[i]*nx.col(s[i]));
                //                ETA = ETA + DoubleVector(GAMMA[i]*nx.col(s[i]));
            }

            //            for (unsigned int i=0; i<L; i++)
            //            {
            //                ETA = ETA + DoubleVector(GAMMA[i]*nx.col(s[i]));
            //            }
        }
        if (s == 1)
        {
            GAMMA[0].resize(n,n,0.0);
            GAMMA[0].at(0,0) = GAMMA[0].at(1,1) = GAMMA[0].at(2,2) = -3.0;

            GAMMA[1].resize(n,n,0.0);
            GAMMA[1].at(0,0) = -12.0*h*a(1,1,1)-10.0; GAMMA[1].at(0,1) = -12.0*h*a(1,2,1);      GAMMA[1].at(0,2) = -12.0*h*a(1,3,1);
            GAMMA[1].at(1,0) = -12.0*h*a(2,1,1);      GAMMA[1].at(1,1) = -12.0*h*a(2,2,1)-10.0; GAMMA[1].at(1,2) = -12.0*h*a(2,3,1);
            GAMMA[1].at(2,0) = -12.0*h*a(3,1,1);      GAMMA[1].at(2,1) = -12.0*h*a(3,2,1);      GAMMA[1].at(2,2) = -12.0*h*a(3,3,1)-10.0;

            GAMMA[2].resize(n,n,0.0);
            GAMMA[2].at(0,0) = GAMMA[2].at(1,1) = GAMMA[2].at(2,2) = +18.0;

            GAMMA[3].resize(n,n,0.0);
            GAMMA[3].at(0,0) = GAMMA[3].at(1,1) = GAMMA[3].at(2,2) = -6.0;

            GAMMA[4].resize(n,n,0.0);
            GAMMA[4].at(0,0) = GAMMA[4].at(1,1) = GAMMA[4].at(2,2) = +1.0;

            ETA.at(0) = 12.0*h*b(1,1);
            ETA.at(1) = 12.0*h*b(2,1);
            ETA.at(2) = 12.0*h*b(3,1);
        }
        if (s == 2)
        {
            GAMMA[0].resize(n,n,0.0);
            GAMMA[0].at(0,0) = GAMMA[0].at(1,1) = GAMMA[0].at(2,2) = +1.0;

            GAMMA[1].resize(n,n,0.0);
            GAMMA[1].at(0,0) = GAMMA[1].at(1,1) = GAMMA[1].at(2,2) = -8.0;

            GAMMA[2].resize(n,n,0.0);
            GAMMA[2].at(0,0) = -12.0*h*a(1,1,2); GAMMA[2].at(0,1) = -12.0*h*a(1,2,2); GAMMA[2].at(0,2) = -12.0*h*a(1,3,2);
            GAMMA[2].at(1,0) = -12.0*h*a(2,1,2); GAMMA[2].at(1,1) = -12.0*h*a(2,2,2); GAMMA[2].at(1,2) = -12.0*h*a(2,3,2);
            GAMMA[2].at(2,0) = -12.0*h*a(3,1,2); GAMMA[2].at(2,1) = -12.0*h*a(3,2,2); GAMMA[2].at(2,2) = -12.0*h*a(3,3,2);

            GAMMA[3].resize(n,n,0.0);
            GAMMA[3].at(0,0) = GAMMA[3].at(1,1) = GAMMA[3].at(2,2) = +8.0;

            GAMMA[4].resize(n,n,0.0);
            GAMMA[4].at(0,0) = GAMMA[4].at(1,1) = GAMMA[4].at(2,2) = -1.0;

            ETA.at(0) = 12.0*h*b(1,2);
            ETA.at(1) = 12.0*h*b(2,2);
            ETA.at(2) = 12.0*h*b(3,2);
        }
        if (s == 3)
        {
            GAMMA[0].resize(n,n,0.0);
            GAMMA[0].at(0,0) = GAMMA[0].at(1,1) = GAMMA[0].at(2,2) = -1.0;

            GAMMA[1].resize(n,n,0.0);
            GAMMA[1].at(0,0) = GAMMA[1].at(1,1) = GAMMA[1].at(2,2) = +6.0;

            GAMMA[2].resize(n,n,0.0);
            GAMMA[2].at(0,0) = GAMMA[2].at(1,1) = GAMMA[2].at(2,2) = -18.0;

            GAMMA[3].resize(n,n,0.0);
            GAMMA[3].at(0,0) = -12.0*h*a(1,1,3)+10.0; GAMMA[3].at(0,1) = -12.0*h*a(1,2,3);      GAMMA[3].at(0,2) = -12.0*h*a(1,3,3);
            GAMMA[3].at(1,0) = -12.0*h*a(2,1,3);      GAMMA[3].at(1,1) = -12.0*h*a(2,2,3)+10.0; GAMMA[3].at(1,2) = -12.0*h*a(2,3,3);
            GAMMA[3].at(2,0) = -12.0*h*a(3,1,3);      GAMMA[3].at(2,1) = -12.0*h*a(3,2,3);      GAMMA[3].at(2,2) = -12.0*h*a(3,3,3)+10.0;

            GAMMA[4].resize(n,n,0.0);
            GAMMA[4].at(0,0) = GAMMA[4].at(1,1) = GAMMA[4].at(2,2) = +3.0;

            ETA.at(0) = 12.0*h*b(1,3);
            ETA.at(1) = 12.0*h*b(2,3);
            ETA.at(2) = 12.0*h*b(3,3);
        }
    }
}

double Example4::fx(unsigned int n, unsigned int i) const
{
    double t = i*h;

#ifdef SAMPLE_1
    if (n==1) return sin(2.0*t) + t*t;
    if (n==2) return 3.0*t;
    if (n==3) return cos(2.0*t) - sin(t);
#endif
#ifdef SAMPLE_2
    if (n==1) return t*t+t;
    if (n==2) return 2.0*t;
    if (n==3) return 3.0*t*t;
#endif
    return 0.0;
}

//double Example4::fx1(unsigned int k) const
//{
//    double t = k*h;
//#ifdef SAMPLE_1
//    return sin(2.0*t) + t*t;
//#endif
//#ifdef SAMPLE_2
//    return t*t+t;
//#endif
//}

//double Example4::fx2(unsigned int k) const
//{
//    double t = k*h;
//#ifdef SAMPLE_1
//    return 3.0*t;
//#endif
//#ifdef SAMPLE_2
//    return 2.0*t;
//#endif
//}

//double Example4::fx3(unsigned int k) const
//{
//    double t = k*h;
//#ifdef SAMPLE_1
//    return cos(2.0*t) - sin(t);
//#endif
//#ifdef SAMPLE_2
//    return 3.0*t*t;
//#endif
//}

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
