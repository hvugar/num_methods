#include "nonlocal.h"

INonLocal::~INonLocal() {}

void INonLocal::calculateAlphaSingle(double *alpha, unsigned int n, unsigned int k, double h) const
{
    if (k==2)
    {
        double t = n*h;
        double m = 1.0/(3.0+2.0*h*a(t));
        alpha[0] = -2.0*h*b(t)*m;
        alpha[1] = +4.0*m;
        alpha[2] = -1.0*m;
    }
    if (k==3)
    {
        double t = n*h;
        double m = 1.0/(11.0+6.0*h*a(t));
        alpha[0] = -6.0*h*b(t)*m;
        alpha[1] = +18.0*m;
        alpha[2] = -9.0*m;
        alpha[3] = +2.0*m;
    }
    if (k==4)
    {
        double t = n*h;
        double m = 1.0/(25.0+12.0*h*a(t));
        alpha[0] = -12.0*h*b(t)*m;
        alpha[1] = +48.0*m;
        alpha[2] = -36.0*m;
        alpha[3] = +16.0*m;
        alpha[4] = -3.0*m;
    }
    if (k==5)
    {
        double t = n*h;
        double m = 1.0/(137.0+60.0*h*a(t));
        alpha[0] = -60.0*h*b(t)*m;
        alpha[1] = +300.0*m;
        alpha[2] = -300.0*m;
        alpha[3] = +200.0*m;
        alpha[4] = -75.0*m;
        alpha[5] = +12.0*m;
    }
    if (k==6)
    {
        double t = n*h;
        double m = 1.0/(147.0+60.0*h*a(t));
        alpha[0] = -60.0*h*b(t)*m;
        alpha[1] = +360.0*m;
        alpha[2] = -450.0*m;
        alpha[3] = +400.0*m;
        alpha[4] = -225.0*m;
        alpha[5] = +72.0*m;
        alpha[6] = -10.0*m;
    }
}

void INonLocal::calculateOtherSingle(DoubleMatrix &c, DoubleVector &d, unsigned int N, unsigned int k, double h, double *betta) const
{
    c.clear(); c.resize(k+1, k+1);
    d.clear(); d.resize(k+1);

    if (k==2)
    {
        c[0][0] = 0.0;  c[0][1] = betta[1]; c[0][2] = betta[2]; d[0] = betta[0];
        //c[0][0] = -3.0; c[0][1] = +4.0;;    c[0][2] = -1.0;     d[0] = +2.0*h*b((N-2)*h); c[0][0] += -2.0*h*a((N-2)*h);
        c[1][0] = -1.0; c[1][1] = +0.0;;    c[1][2] = +1.0;     d[1] = +2.0*h*b((N-1)*h); c[1][1] += -2.0*h*a((N-1)*h);
        c[2][0] = +1.0; c[2][1] = -4.0;     c[2][2] = +3.0;     d[2] = +2.0*h*b((N-0)*h); c[2][2] += -2.0*h*a((N-0)*h);
    }
    if (k==3)
    {
        c[0][0] =  0.0;  c[0][1] = betta[1]; c[0][2] = betta[2]; c[0][3] = betta[3]; d[0] = betta[0];
        //c[0][0] = -11.0; c[0][1] = +18.0;    c[0][2] = -9.0;     c[0][3] = +2.0;     d[0] = +6.0*h*b((N-3)*h); c[0][1] += -6.0*h*a((N-3)*h);
        c[1][0] = -2.0;  c[1][1] = -3.0;     c[1][2] = +6.0;     c[1][3] = -1.0;     d[1] = +6.0*h*b((N-2)*h); c[1][1] += -6.0*h*a((N-2)*h);
        c[2][0] = +1.0;  c[2][1] = -6.0;     c[2][2] = +3.0;     c[2][3] = +2.0;     d[2] = +6.0*h*b((N-1)*h); c[2][2] += -6.0*h*a((N-1)*h);
        c[3][0] = -2.0;  c[3][1] = +9.0;     c[3][2] = -18.0;    c[3][3] = +11.0;    d[3] = +6.0*h*b((N-0)*h); c[3][3] += -6.0*h*a((N-0)*h);
    }
    if (k==4)
    {
        c[0][0] =  0.0;  c[0][1] = betta[1]; c[0][2] = betta[2]; c[0][3] = betta[3]; c[0][4] = betta[4]; d[0] = betta[0];
        //c[0][0] = -25.0; c[1][1] = +48.0;    c[1][2] = -36.0;    c[1][3] = +16.0;    c[1][4] = -3.0;     d[0] = +12.0*h*b((N-4)*h); c[0][0] += -12.0*h*a((N-4)*h);
        c[1][0] = -3.0;  c[1][1] = -10.0;    c[1][2] = +18.0;    c[1][3] = -6.0;     c[1][4] = +1.0;     d[1] = +12.0*h*b((N-3)*h); c[1][1] += -12.0*h*a((N-3)*h);
        c[2][0] = +1.0;  c[2][1] = -8.0;     c[2][2] = +0.0;     c[2][3] = +8.0;     c[2][4] = -1.0;     d[2] = +12.0*h*b((N-2)*h); c[2][2] += -12.0*h*a((N-2)*h);
        c[3][0] = -1.0;  c[3][1] = +6.0;     c[3][2] = -18.0;    c[3][3] = +10.0;    c[3][4] = +3.0;     d[3] = +12.0*h*b((N-1)*h); c[3][3] += -12.0*h*a((N-1)*h);
        c[4][0] = +3.0;  c[4][1] = -16.0;    c[4][2] = +36.0;    c[4][3] = -48.0;    c[4][4] = +25.0;    d[4] = +12.0*h*b((N-0)*h); c[4][4] += -12.0*h*a((N-0)*h);
    }
    if (k==5)
    {
        c[0][0] =  0.0;   c[0][1] = betta[1]; c[0][2] = betta[2]; c[0][3] = betta[3]; c[0][4] = betta[4]; c[0][5] = betta[5]; d[0] = betta[0];
        //c[0][0] = -137.0; c[0][1] = +300.0;   c[0][2] = -300.0;   c[0][3] = +200.0;   c[0][4] = -75.0;    c[0][5] = +12.0;    d[0] = +60.0*h*b((N-5)*h); c[0][0] += -60.0*h*a((N-5)*h);
        c[1][0] = -12.0;  c[1][1] = -65.0;    c[1][2] = +120.0;   c[1][3] = -60.0;    c[1][4] = +20.0;    c[1][5] = -3.0;     d[1] = +60.0*h*b((N-4)*h); c[1][1] += -60.0*h*a((N-4)*h);
        c[2][0] = +3.0;   c[2][1] = -30.0;    c[2][2] = -20.0;    c[2][3] = +60.0;    c[2][4] = -15.0;    c[2][5] = +2.0;     d[2] = +60.0*h*b((N-3)*h); c[2][2] += -60.0*h*a((N-3)*h);
        c[3][0] = -2.0;   c[3][1] = +15.0;    c[3][2] = -60.0;    c[3][3] = +20.0;    c[3][4] = +30.0;    c[3][5] = -3.0;     d[3] = +60.0*h*b((N-2)*h); c[3][3] += -60.0*h*a((N-2)*h);
        c[4][0] = +3.0;   c[4][1] = -20.0;    c[4][2] = +60.0;    c[4][3] = -120.0;   c[4][4] = +65.0;    c[4][5] = +12.0;    d[4] = +60.0*h*b((N-1)*h); c[4][4] += -60.0*h*a((N-1)*h);
        c[5][0] = -12.0;  c[5][1] = +75.0;    c[5][2] = -200.0;   c[5][3] = +300.0;   c[5][4] = -300.0;   c[5][5] = +137.0;   d[5] = +60.0*h*b((N-0)*h); c[5][5] += -60.0*h*a((N-0)*h);
    }
    if (k==6)
    {
        //        c[0][0] =  -147.0;  c[0][1] = 360.0; c[0][2] = -450.0; c[0][3] = +400.0; c[0][4] = -225.0; c[0][5] = +72.0; c[0][6] = -10.0; d[0] = +60.0*h*b((N-6)*h); c[0][0] += -60.0*h*a((N-6)*h);
        c[0][0] =  0.0;  c[0][1] = betta[1]; c[0][2] = betta[2]; c[0][3] = betta[3]; c[0][4] = betta[4]; c[0][5] = betta[5]; c[0][6] = betta[6]; d[0] = betta[0];
        c[1][0] = -10.0; c[1][1] = -77.0;    c[1][2] = +150.0;   c[1][3] = -100.0;   c[1][4] = +50.0;    c[1][5] = -15.0;    c[1][6] = +2.0;     d[1] = +60.0*h*b((N-5)*h); c[1][1] += -60.0*h*a((N-5)*h);
        c[2][0] = +2.0;  c[2][1] = -24.0;    c[2][2] = -35.0;    c[2][3] = +80.0;    c[2][4] = -30.0;    c[2][5] = +8.0;     c[2][6] = -1.0;     d[2] = +60.0*h*b((N-4)*h); c[2][2] += -60.0*h*a((N-4)*h);
        c[3][0] = -1.0;  c[3][1] = +9.0;     c[3][2] = -45.0;    c[3][3] = +0.0;     c[3][4] = +45.0;    c[3][5] = -9.0;     c[3][6] = +1.0;     d[3] = +60.0*h*b((N-3)*h); c[3][3] += -60.0*h*a((N-3)*h);
        c[4][0] = +1.0;  c[4][1] = -8.0;     c[4][2] = +30.0;    c[4][3] = -80.0;    c[4][4] = +35.0;    c[4][5] = +24.0;    c[4][6] = -2.0;     d[4] = +60.0*h*b((N-2)*h); c[4][4] += -60.0*h*a((N-2)*h);
        c[5][0] = -2.0;  c[5][1] = +15.0;    c[5][2] = -50.0;    c[5][3] = +100.0;   c[5][4] = -150.0;   c[5][5] = +77.0;    c[5][6] = +10.0;    d[5] = +60.0*h*b((N-1)*h); c[5][5] += -60.0*h*a((N-1)*h);
        c[6][0] = +10.0; c[6][1] = -72.0;    c[6][2] = +225.0;   c[6][3] = -400.0;   c[6][4] = +450.0;   c[6][5] = -360.0;   c[6][6] = +147.0;   d[6] = +60.0*h*b((N-0)*h); c[6][6] += -60.0*h*a((N-0)*h);
    }
}

void INonLocal::solveSingle(const DoubleVector &a, double b, double h) const
{
    const unsigned int k = 6;
    const unsigned int N = a.length()-1;
    const unsigned int end = N-k+1;

    double **betta = new double*[end+1];
    for (unsigned int n=0; n<=end; n++) betta[n] = new double[k+1];

    for (unsigned int n=0; n<=end; n++)
    {
        if (n==0)
        {
            betta[n][0] = b;
            for (unsigned int i=1; i<=k; i++) betta[n][i] = a[i-1];
        }
        else
        {
            double alpha[k+1];
            calculateAlphaSingle(alpha, n-1, k, h);

            betta[n][0] = -betta[n-1][1]*alpha[0] + betta[n-1][0];
            for (unsigned int i=1; i<=k-1; i++)
            {
                betta[n][i] = betta[n-1][1]*alpha[i] + betta[n-1][i+1];
            }
            betta[n][k] = betta[n-1][1]*alpha[k] + a[k+(n-1)];
        }
    }

    DoubleMatrix c;
    DoubleVector d;


    calculateOtherSingle(c, d, N, k, h, betta[end]);


    DoubleVector xf(k+1);
    LinearEquation::GaussianElimination(c, d, xf);

    double norm = 0.0;
    DoubleVector x0(N+1);
    for (unsigned int n=N; n>=N-k; n--)
    {
        x0[n] = xf[n-(N-k)];
        norm += (x0[n]-x(n*h))*(x0[n]-x(n*h));
    }
    xf.clear();

    c.clear();
    d.clear();

    unsigned int stop = static_cast<unsigned int>(0)-1;
    for (unsigned int n=N-k-1; n!=stop; n--)
    {
        x0[n] = betta[n][0];
        for (unsigned int j=2; j<=k; j++) x0[n] -= betta[n][j]*x0[n+j-1];
        for (unsigned int j=n+k; j<=N; j++) x0[n] -= a[j]*x0[j];
        x0[n] /= betta[n][1];
        //printf("%4d %.8f %.8f\n", n, x[n], NonLocal::x(n*h));
        printf("%.8f\n", x0[n]);
        norm += (x0[n]-x(n*h))*(x0[n]-x(n*h));
    }

    for (unsigned int n=0; n<=end; n++) delete [] (betta[n]);
    delete [] betta;

    printf("norm %.8f\n", sqrt(norm));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

INonLocalSystem::~INonLocalSystem() {}

void INonLocalSystem::getAlpha(DoubleMatrix *alpha, unsigned int n, unsigned int k, double h, unsigned int M) const
{
    if (k == 2)
    {
        double t = h*n;
        DoubleMatrix mx(M, M);
        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                mx[r][c] = 2.0*h*A(t, r+1, c+1);
            }
            mx[r][r] += 3.0;
        }
        mx.inverse();
        alpha[0].resize(M, 1, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[0])[m][0] = -2.0*h*B(t, m+1); } alpha[0] = mx*alpha[0];
        alpha[1].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[1])[m][m] = +4.0; }             alpha[1] = mx*alpha[1];
        alpha[2].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[2])[m][m] = -1.0; }             alpha[2] = mx*alpha[2];

//        double t = h*(n+2);
//        alpha[0].resize(M, 1, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[0])[m][0] = +2.0*h*B(t, m+1); }
//        alpha[1].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[1])[m][m] = +4.0; }
//        alpha[2].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[2])[m][m] = -3.0; }
//        for (unsigned int r=0; r<M; r++) for (unsigned int c=0; c<M; c++) (alpha[2])[r][c] += 2.0*h*A(t, r+1, c+1);

//        double t = h*(n+1);
//        alpha[0].resize(M, 1, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[0])[m][0] = -2.0*h*B(t, m+1); }
//        alpha[1].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[1])[m][m] = +0.0; }
//        alpha[2].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[2])[m][m] = +1.0; }
//        for (unsigned int r=0; r<M; r++) for (unsigned int c=0; c<M; c++) (alpha[1])[r][c] += -2.0*h*A(t, r+1, c+1);

    }
    if (k == 3)
    {
        //        for (unsigned int n=1; n<=N; n++)
        //        {
        //            double t = h*m;
        //            double mx = 1.0/(11.0+6.0*h*A(t, n, n));
        //            alpha[n][0] = -6.0*h*B(t, n)*mx;
        //            alpha[n][1] = +18.0*mx;
        //            alpha[n][2] = -9.0*mx;
        //            alpha[n][3] = +2.0*mx;
        //        }
    }
    if (k==4)
    {
        double t = h*n;
        DoubleMatrix mx(M, M);
        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                mx[r][c] = 12.0*h*A(t, r+1, c+1);
            }
            mx[r][r] += 25.0;
        }
        mx.inverse();

        alpha[0].resize(M, 1, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[0])[m][0] = -12.0*h*B(t, m+1); } alpha[0] = mx*alpha[0];
        alpha[1].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[1])[m][m] = +48.0; }             alpha[1] = mx*alpha[1];
        alpha[2].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[2])[m][m] = -36.0; }             alpha[2] = mx*alpha[2];
        alpha[3].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[3])[m][m] = +16.0; }             alpha[3] = mx*alpha[3];
        alpha[4].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[4])[m][m] = -3.0; }              alpha[4] = mx*alpha[4];
    }
    if (k==5)
    {
        //        for (unsigned int n=1; n<=N; n++)
        //        {
        //            double t = h*m;
        //            double mx = 1.0/(137.0+60.0*h*A(t, n, n));
        //            alpha[n][0] = -60.0*h*B(t, n)*mx;
        //            alpha[n][1] = +300.0*mx;
        //            alpha[n][2] = -300.0*mx;
        //            alpha[n][3] = +200.0*mx;
        //            alpha[n][4] = -75.0*mx;
        //            alpha[n][5] = +12.0*mx;
        //        }
    }
    if (k==6)
    {
        double t = h*n;
        DoubleMatrix mx(M, M);
        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                mx[r][c] = 60.0*h*A(t, r+1, c+1);
            }
            mx[r][r] += 147.0;
        }
        mx.inverse();

        alpha[0].resize(M, 1, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[0])[m][0] = -60.0*h*B(t, m+1); } alpha[0] = mx*alpha[0];
        alpha[1].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[1])[m][m] = +360.0; }            alpha[1] = mx*alpha[1];
        alpha[2].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[2])[m][m] = -450.0; }            alpha[2] = mx*alpha[2];
        alpha[3].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[3])[m][m] = +400.0; }            alpha[3] = mx*alpha[3];
        alpha[4].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[4])[m][m] = -225.0; }            alpha[4] = mx*alpha[4];
        alpha[5].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[5])[m][m] = +72.0; }             alpha[5] = mx*alpha[5];
        alpha[6].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[6])[m][m] = -10.0; }             alpha[6] = mx*alpha[6];


        //        for (unsigned int n=1; n<=N; n++)
        //        {
        //            double t = h*m;
        //            double mx = 1.0/(147.0+60.0*h*A(t, n, n));
        //            alpha[n][0] = -60.0*h*B(t, n)*mx;
        //            alpha[n][1] = +360.0*mx;
        //            alpha[n][2] = -450.0*mx;
        //            alpha[n][3] = +400.0*mx;
        //            alpha[n][4] = -225.0*mx;
        //            alpha[n][5] = +72.0*mx;
        //            alpha[n][6] = -10.0*mx;
        //        }
    }
}

void INonLocalSystem::getEquations(DoubleMatrix &C, DoubleVector &d, unsigned int N, unsigned int k, double h, DoubleMatrix *bt, unsigned int M) const
{
    C.clear(); C.resize((k+1)*M, (k+1)*M);
    d.clear(); d.resize((k+1)*M);

    if (k==2)
    {
        double t1 = (N-1)*h;
        double t0 = (N-0)*h;
        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                C[0*M+r][0*M+c] = 0.0;
                C[0*M+r][1*M+c] = bt[1][r][c];
                C[0*M+r][2*M+c] = bt[2][r][c];

                C[1*M+r][0*M+c] = 0.0;                  if (r==c) C[1*M+r][0*M+c] += -1.0;
                C[1*M+r][1*M+c] = -2.0*h*A(t1,r+1,c+1); if (r==c) C[1*M+r][1*M+c] += +0.0;
                C[1*M+r][2*M+c] = 0.0;                  if (r==c) C[1*M+r][2*M+c] += +1.0;

                C[2*M+r][0*M+c] = 0.0;                  if (r==c) C[2*M+r][0*M+c] += +1.0;
                C[2*M+r][1*M+c] = 0.0;                  if (r==c) C[2*M+r][1*M+c] += -4.0;
                C[2*M+r][2*M+c] = -2.0*h*A(t0,r+1,c+1); if (r==c) C[2*M+r][2*M+c] += +3.0;
            }
            d[0*M+r] = bt[0][r][0];
            d[1*M+r] = 2.0*h*B(t1,r+1);
            d[2*M+r] = 2.0*h*B(t0,r+1);
        }

//        double t2 = (N-2)*h;
//        double t1 = (N-1)*h;
//        for (unsigned int r=0; r<M; r++)
//        {
//            for (unsigned int c=0; c<M; c++)
//            {
//                C[0*M+r][0*M+c] = 0.0;
//                C[0*M+r][1*M+c] = bt[1][r][c];
//                C[0*M+r][2*M+c] = bt[2][r][c];

//                C[1*M+r][0*M+c] = -2.0*h*A(t2,r+1,c+1); if (r==c) C[1*M+r][0*M+c] += -3.0;
//                C[1*M+r][1*M+c] = 0.0;                  if (r==c) C[1*M+r][1*M+c] += +4.0;
//                C[1*M+r][2*M+c] = 0.0;                  if (r==c) C[1*M+r][2*M+c] += -1.0;

//                C[2*M+r][0*M+c] = 0.0;                  if (r==c) C[2*M+r][0*M+c] += -1.0;
//                C[2*M+r][1*M+c] = -2.0*h*A(t1,r+1,c+1); if (r==c) C[2*M+r][1*M+c] += +0.0;
//                C[2*M+r][2*M+c] = 0.0;                  if (r==c) C[2*M+r][2*M+c] += +1.0;
//            }
//            d[0*M+r] = bt[0][r][0];
//            d[1*M+r] = 2.0*h*B(t2,r+1);
//            d[2*M+r] = 2.0*h*B(t1,r+1);
//        }

//        double t2 = (N-2)*h;
//        double t0 = (N-0)*h;
//        for (unsigned int r=0; r<M; r++)
//        {
//            for (unsigned int c=0; c<M; c++)
//            {
//                C[0*M+r][0*M+c] = 0.0;
//                C[0*M+r][1*M+c] = bt[1][r][c];
//                C[0*M+r][2*M+c] = bt[2][r][c];

//                C[1*M+r][0*M+c] = -2.0*h*A(t2,r+1,c+1); if (r==c) C[1*M+r][0*M+c] += -3.0;
//                C[1*M+r][1*M+c] = 0.0;                  if (r==c) C[1*M+r][1*M+c] += +4.0;
//                C[1*M+r][2*M+c] = 0.0;                  if (r==c) C[1*M+r][2*M+c] += -1.0;

//                C[2*M+r][0*M+c] = 0.0;                  if (r==c) C[2*M+r][0*M+c] += +1.0;
//                C[2*M+r][1*M+c] = 0.0;                  if (r==c) C[2*M+r][1*M+c] += -4.0;
//                C[2*M+r][2*M+c] = -2.0*h*A(t0,r+1,c+1); if (r==c) C[2*M+r][2*M+c] += +3.0;
//            }
//            d[0*M+r] = bt[0][r][0];
//            d[1*M+r] = 2.0*h*B(t2,r+1);
//            d[2*M+r] = 2.0*h*B(t0,r+1);
//        }
    }
    if (k==3)
    {
        //        c[0][0] =  0.0;  c[0][1] = betta[1]; c[0][2] = betta[2]; c[0][3] = betta[3]; d[0] = betta[0];
        //        //c[0][0] = -11.0; c[0][1] = +18.0;    c[0][2] = -9.0;     c[0][3] = +2.0;     d[0] = +6.0*h*b((N-3)*h); c[0][1] += -6.0*h*a((N-3)*h);
        //        c[1][0] = -2.0;  c[1][1] = -3.0;     c[1][2] = +6.0;     c[1][3] = -1.0;     d[1] = +6.0*h*b((N-2)*h); c[1][1] += -6.0*h*a((N-2)*h);
        //        c[2][0] = +1.0;  c[2][1] = -6.0;     c[2][2] = +3.0;     c[2][3] = +2.0;     d[2] = +6.0*h*b((N-1)*h); c[2][2] += -6.0*h*a((N-1)*h);
        //        c[3][0] = -2.0;  c[3][1] = +9.0;     c[3][2] = -18.0;    c[3][3] = +11.0;    d[3] = +6.0*h*b((N-0)*h); c[3][3] += -6.0*h*a((N-0)*h);
    }
    if (k==4)
    {
        double t3 = (N-3)*h;
        double t2 = (N-2)*h;
        double t1 = (N-1)*h;
        double t0 = (N-0)*h;

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                C[0*M+r][0*M+c] = 0.0;
                C[0*M+r][1*M+c] = bt[1][r][c];
                C[0*M+r][2*M+c] = bt[2][r][c];
                C[0*M+r][3*M+c] = bt[3][r][c];
                C[0*M+r][4*M+c] = bt[4][r][c];

                C[1*M+r][0*M+c] = 0.0;                   if (r==c) C[1*M+r][0*M+c] += -3.0;
                C[1*M+r][1*M+c] = -12.0*h*A(t3,r+1,c+1); if (r==c) C[1*M+r][1*M+c] += -10.0;
                C[1*M+r][2*M+c] = 0.0;                   if (r==c) C[1*M+r][2*M+c] += +18.0;
                C[1*M+r][3*M+c] = 0.0;                   if (r==c) C[1*M+r][3*M+c] += -6.0;
                C[1*M+r][4*M+c] = 0.0;                   if (r==c) C[1*M+r][4*M+c] += +1.0;

                C[2*M+r][0*M+c] = 0.0;                   if (r==c) C[2*M+r][0*M+c] += +1.0;
                C[2*M+r][1*M+c] = 0.0;                   if (r==c) C[2*M+r][1*M+c] += -8.0;
                C[2*M+r][2*M+c] = -12.0*h*A(t2,r+1,c+1); if (r==c) C[2*M+r][2*M+c] += +0.0;
                C[2*M+r][3*M+c] = 0.0;                   if (r==c) C[2*M+r][3*M+c] += +8.0;
                C[2*M+r][4*M+c] = 0.0;                   if (r==c) C[2*M+r][4*M+c] += -1.0;

                C[3*M+r][0*M+c] = 0.0;                   if (r==c) C[3*M+r][0*M+c] += -1.0;
                C[3*M+r][1*M+c] = 0.0;                   if (r==c) C[3*M+r][1*M+c] += +6.0;
                C[3*M+r][2*M+c] = 0.0;                   if (r==c) C[3*M+r][2*M+c] += -18.0;
                C[3*M+r][3*M+c] = -12.0*h*A(t1,r+1,c+1); if (r==c) C[3*M+r][3*M+c] += +10.0;
                C[3*M+r][4*M+c] = 0.0;                   if (r==c) C[3*M+r][4*M+c] += +3.0;

                C[4*M+r][0*M+c] = 0.0;                   if (r==c) C[4*M+r][0*M+c] += +3.0;
                C[4*M+r][1*M+c] = 0.0;                   if (r==c) C[4*M+r][1*M+c] += -16.0;
                C[4*M+r][2*M+c] = 0.0;                   if (r==c) C[4*M+r][2*M+c] += +36.0;
                C[4*M+r][3*M+c] = 0.0;                   if (r==c) C[4*M+r][3*M+c] += -48.0;
                C[4*M+r][4*M+c] = -12.0*h*A(t0,r+1,c+1); if (r==c) C[4*M+r][4*M+c] += +25.0;
            }
            d[0*M+r] = bt[0][r][0];
            d[1*M+r] = +12.0*h*B(t3,r+1);
            d[2*M+r] = +12.0*h*B(t2,r+1);
            d[3*M+r] = +12.0*h*B(t1,r+1);
            d[4*M+r] = +12.0*h*B(t0,r+1);
        }

        //        c[0][0] =  0.0;  c[0][1] = betta[1]; c[0][2] = betta[2]; c[0][3] = betta[3]; c[0][4] = betta[4]; d[0] = betta[0];
        //        //c[0][0] = -25.0; c[1][1] = +48.0;    c[1][2] = -36.0;    c[1][3] = +16.0;    c[1][4] = -3.0;     d[0] = +12.0*h*b((N-4)*h); c[0][0] += -12.0*h*a((N-4)*h);
        //        c[1][0] = -3.0;  c[1][1] = -10.0;    c[1][2] = +18.0;    c[1][3] = -6.0;     c[1][4] = +1.0;     d[1] = +12.0*h*b((N-3)*h); c[1][1] += -12.0*h*a((N-3)*h);
        //        c[2][0] = +1.0;  c[2][1] = -8.0;     c[2][2] = +0.0;     c[2][3] = +8.0;     c[2][4] = -1.0;     d[2] = +12.0*h*b((N-2)*h); c[2][2] += -12.0*h*a((N-2)*h);
        //        c[3][0] = -1.0;  c[3][1] = +6.0;     c[3][2] = -18.0;    c[3][3] = +10.0;    c[3][4] = +3.0;     d[3] = +12.0*h*b((N-1)*h); c[3][3] += -12.0*h*a((N-1)*h);
        //        c[4][0] = +3.0;  c[4][1] = -16.0;    c[4][2] = +36.0;    c[4][3] = -48.0;    c[4][4] = +25.0;    d[4] = +12.0*h*b((N-0)*h); c[4][4] += -12.0*h*a((N-0)*h);
    }
    if (k==5)
    {
        //        c[0][0] =  0.0;   c[0][1] = betta[1]; c[0][2] = betta[2]; c[0][3] = betta[3]; c[0][4] = betta[4]; c[0][5] = betta[5]; d[0] = betta[0];
        //        //c[0][0] = -137.0; c[0][1] = +300.0;   c[0][2] = -300.0;   c[0][3] = +200.0;   c[0][4] = -75.0;    c[0][5] = +12.0;    d[0] = +60.0*h*b((N-5)*h); c[0][0] += -60.0*h*a((N-5)*h);
        //        c[1][0] = -12.0;  c[1][1] = -65.0;    c[1][2] = +120.0;   c[1][3] = -60.0;    c[1][4] = +20.0;    c[1][5] = -3.0;     d[1] = +60.0*h*b((N-4)*h); c[1][1] += -60.0*h*a((N-4)*h);
        //        c[2][0] = +3.0;   c[2][1] = -30.0;    c[2][2] = -20.0;    c[2][3] = +60.0;    c[2][4] = -15.0;    c[2][5] = +2.0;     d[2] = +60.0*h*b((N-3)*h); c[2][2] += -60.0*h*a((N-3)*h);
        //        c[3][0] = -2.0;   c[3][1] = +15.0;    c[3][2] = -60.0;    c[3][3] = +20.0;    c[3][4] = +30.0;    c[3][5] = -3.0;     d[3] = +60.0*h*b((N-2)*h); c[3][3] += -60.0*h*a((N-2)*h);
        //        c[4][0] = +3.0;   c[4][1] = -20.0;    c[4][2] = +60.0;    c[4][3] = -120.0;   c[4][4] = +65.0;    c[4][5] = +12.0;    d[4] = +60.0*h*b((N-1)*h); c[4][4] += -60.0*h*a((N-1)*h);
        //        c[5][0] = -12.0;  c[5][1] = +75.0;    c[5][2] = -200.0;   c[5][3] = +300.0;   c[5][4] = -300.0;   c[5][5] = +137.0;   d[5] = +60.0*h*b((N-0)*h); c[5][5] += -60.0*h*a((N-0)*h);
    }
    if (k==6)
    {
        double t5 = (N-5)*h;
        double t4 = (N-4)*h;
        double t3 = (N-3)*h;
        double t2 = (N-2)*h;
        double t1 = (N-1)*h;
        double t0 = (N-0)*h;

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                C[0*M+r][0*M+c] = 0.0;
                C[0*M+r][1*M+c] = bt[1][r][c];
                C[0*M+r][2*M+c] = bt[2][r][c];
                C[0*M+r][3*M+c] = bt[3][r][c];
                C[0*M+r][4*M+c] = bt[4][r][c];
                C[0*M+r][5*M+c] = bt[5][r][c];
                C[0*M+r][6*M+c] = bt[6][r][c];

                C[1*M+r][0*M+c] = 0.0;                   if (r==c) C[1*M+r][0*M+c] += -10.0;
                C[1*M+r][1*M+c] = -60.0*h*A(t5,r+1,c+1); if (r==c) C[1*M+r][1*M+c] += -77.0;
                C[1*M+r][2*M+c] = 0.0;                   if (r==c) C[1*M+r][2*M+c] += +150.0;
                C[1*M+r][3*M+c] = 0.0;                   if (r==c) C[1*M+r][3*M+c] += -100.0;
                C[1*M+r][4*M+c] = 0.0;                   if (r==c) C[1*M+r][4*M+c] += +50.0;
                C[1*M+r][5*M+c] = 0.0;                   if (r==c) C[1*M+r][5*M+c] += -15.0;
                C[1*M+r][6*M+c] = 0.0;                   if (r==c) C[1*M+r][6*M+c] += +2.0;

                C[2*M+r][0*M+c] = 0.0;                   if (r==c) C[2*M+r][0*M+c] += +2.0;
                C[2*M+r][1*M+c] = 0.0;                   if (r==c) C[2*M+r][1*M+c] += -24.0;
                C[2*M+r][2*M+c] = -60.0*h*A(t4,r+1,c+1); if (r==c) C[2*M+r][2*M+c] += -35.0;
                C[2*M+r][3*M+c] = 0.0;                   if (r==c) C[2*M+r][3*M+c] += +80.0;
                C[2*M+r][4*M+c] = 0.0;                   if (r==c) C[2*M+r][4*M+c] += -30.0;
                C[2*M+r][5*M+c] = 0.0;                   if (r==c) C[2*M+r][5*M+c] += +8.0;
                C[2*M+r][6*M+c] = 0.0;                   if (r==c) C[2*M+r][6*M+c] += -1.0;

                C[3*M+r][0*M+c] = 0.0;                   if (r==c) C[3*M+r][0*M+c] += -1.0;
                C[3*M+r][1*M+c] = 0.0;                   if (r==c) C[3*M+r][1*M+c] += +9.0;
                C[3*M+r][2*M+c] = 0.0;                   if (r==c) C[3*M+r][2*M+c] += -45.0;
                C[3*M+r][3*M+c] = -60.0*h*A(t3,r+1,c+1); if (r==c) C[3*M+r][3*M+c] += +0.0;
                C[3*M+r][4*M+c] = 0.0;                   if (r==c) C[3*M+r][4*M+c] += +45.0;
                C[3*M+r][5*M+c] = 0.0;                   if (r==c) C[3*M+r][5*M+c] += -9.0;
                C[3*M+r][6*M+c] = 0.0;                   if (r==c) C[3*M+r][6*M+c] += +1.0;

                C[4*M+r][0*M+c] = 0.0;                   if (r==c) C[4*M+r][0*M+c] += +1.0;
                C[4*M+r][1*M+c] = 0.0;                   if (r==c) C[4*M+r][1*M+c] += -8.0;
                C[4*M+r][2*M+c] = 0.0;                   if (r==c) C[4*M+r][2*M+c] += +30.0;
                C[4*M+r][3*M+c] = 0.0;                   if (r==c) C[4*M+r][3*M+c] += -80.0;
                C[4*M+r][4*M+c] = -60.0*h*A(t2,r+1,c+1); if (r==c) C[4*M+r][4*M+c] += +35.0;
                C[4*M+r][5*M+c] = 0.0;                   if (r==c) C[4*M+r][5*M+c] += +24.0;
                C[4*M+r][6*M+c] = 0.0;                   if (r==c) C[4*M+r][6*M+c] += -2.0;

                C[5*M+r][0*M+c] = 0.0;                   if (r==c) C[5*M+r][0*M+c] += -2.0;
                C[5*M+r][1*M+c] = 0.0;                   if (r==c) C[5*M+r][1*M+c] += +15.0;
                C[5*M+r][2*M+c] = 0.0;                   if (r==c) C[5*M+r][2*M+c] += -50.0;
                C[5*M+r][3*M+c] = 0.0;                   if (r==c) C[5*M+r][3*M+c] += +100.0;
                C[5*M+r][4*M+c] = 0.0;                   if (r==c) C[5*M+r][4*M+c] += -150.0;
                C[5*M+r][5*M+c] = -60.0*h*A(t1,r+1,c+1); if (r==c) C[5*M+r][5*M+c] += +77.0;
                C[5*M+r][6*M+c] = 0.0;                   if (r==c) C[5*M+r][6*M+c] += +10.0;

                C[6*M+r][0*M+c] = 0.0;                   if (r==c) C[6*M+r][0*M+c] += +10.0;
                C[6*M+r][1*M+c] = 0.0;                   if (r==c) C[6*M+r][1*M+c] += -72.0;
                C[6*M+r][2*M+c] = 0.0;                   if (r==c) C[6*M+r][2*M+c] += +225.0;
                C[6*M+r][3*M+c] = 0.0;                   if (r==c) C[6*M+r][3*M+c] += -400.0;
                C[6*M+r][4*M+c] = 0.0;                   if (r==c) C[6*M+r][4*M+c] += +450.0;
                C[6*M+r][5*M+c] = 0.0;                   if (r==c) C[6*M+r][5*M+c] += -360.0;
                C[6*M+r][6*M+c] = -60.0*h*A(t0,r+1,c+1); if (r==c) C[6*M+r][6*M+c] += +147.0;
            }
            d[0*M+r] = bt[0][r][0];
            d[1*M+r] = +60.0*h*B(t5,r+1);
            d[2*M+r] = +60.0*h*B(t4,r+1);
            d[3*M+r] = +60.0*h*B(t3,r+1);
            d[4*M+r] = +60.0*h*B(t2,r+1);
            d[5*M+r] = +60.0*h*B(t1,r+1);
            d[6*M+r] = +60.0*h*B(t0,r+1);
        }
        //        c[0][0] =  0.0;  c[0][1] = betta[1]; c[0][2] = betta[2]; c[0][3] = betta[3]; c[0][4] = betta[4]; c[0][5] = betta[5]; c[0][6] = betta[6]; d[0] = betta[0];
        //        //c[0][0] =  -147.0;  c[0][1] = 360.0; c[0][2] = -450.0; c[0][3] = +400.0; c[0][4] = -225.0; c[0][5] = +72.0; c[0][6] = -10.0; d[0] = +60.0*h*b((N-6)*h); c[0][0] += -60.0*h*a((N-6)*h);
        //        c[1][0] = -10.0; c[1][1] = -77.0;    c[1][2] = +150.0;   c[1][3] = -100.0;   c[1][4] = +50.0;    c[1][5] = -15.0;    c[1][6] = +2.0;     d[1] = +60.0*h*b((N-5)*h); c[1][1] += -60.0*h*a((N-5)*h);
        //        c[2][0] = +2.0;  c[2][1] = -24.0;    c[2][2] = -35.0;    c[2][3] = +80.0;    c[2][4] = -30.0;    c[2][5] = +8.0;     c[2][6] = -1.0;     d[2] = +60.0*h*b((N-4)*h); c[2][2] += -60.0*h*a((N-4)*h);
        //        c[3][0] = -1.0;  c[3][1] = +9.0;     c[3][2] = -45.0;    c[3][3] = +0.0;     c[3][4] = +45.0;    c[3][5] = -9.0;     c[3][6] = +1.0;     d[3] = +60.0*h*b((N-3)*h); c[3][3] += -60.0*h*a((N-3)*h);
        //        c[4][0] = +1.0;  c[4][1] = -8.0;     c[4][2] = +30.0;    c[4][3] = -80.0;    c[4][4] = +35.0;    c[4][5] = +24.0;    c[4][6] = -2.0;     d[4] = +60.0*h*b((N-2)*h); c[4][4] += -60.0*h*a((N-2)*h);
        //        c[5][0] = -2.0;  c[5][1] = +15.0;    c[5][2] = -50.0;    c[5][3] = +100.0;   c[5][4] = -150.0;   c[5][5] = +77.0;    c[5][6] = +10.0;    d[5] = +60.0*h*b((N-1)*h); c[5][5] += -60.0*h*a((N-1)*h);
        //        c[6][0] = +10.0; c[6][1] = -72.0;    c[6][2] = +225.0;   c[6][3] = -400.0;   c[6][4] = +450.0;   c[6][5] = -360.0;   c[6][6] = +147.0;   d[6] = +60.0*h*b((N-0)*h); c[6][6] += -60.0*h*a((N-0)*h);
    }
}

void INonLocalSystem::solve(const std::vector<DoubleMatrix> &C, const DoubleVector &d, std::vector<DoubleVector> &x, const Dimension &dim, unsigned int M) const
{
    const unsigned int k = 2;
    const unsigned int N = static_cast<unsigned int>(C.size()-1);
    const unsigned int end = N-k+1;
    const double h = dim.step();

    DoubleMatrix **betta = new DoubleMatrix*[end+1];
    for (unsigned int n=0; n<=end; n++)
    {
        betta[n] = new DoubleMatrix[k+1];
        betta[n][0].resize(N,1);
        for (unsigned int i=1; i<=k; i++) betta[n][k].resize(M,M);
    }

    for (unsigned int n=0; n<=end; n++)
    {
        if (n==0)
        {
            betta[n][0] = d;
            for (unsigned int i=1; i<=k; i++) betta[n][i] = C[i-1];
        }
        else
        {
            DoubleMatrix alpha[k+1];
            getAlpha(alpha, n-1, k, h, M);

            betta[n][0] = betta[n-1][0] - betta[n-1][1]*alpha[0];
            for (unsigned int i=1; i<=k-1; i++)
            {
                betta[n][i] = betta[n-1][i+1] + betta[n-1][1]*alpha[i];
            }
            betta[n][k] = C[k+(n-1)] + betta[n-1][1]*alpha[k];
        }
    }

    DoubleMatrix F;
    DoubleVector g;

    getEquations(F, g, N, k, h, betta[end], M);

//    IPrinter::printSeperatorLine("F");
//    IPrinter::print(F);
//    IPrinter::printSeperatorLine("g");
//    IPrinter::print(g, g.length());

//    IPrinter::printSeperatorLine("X");
//    DoubleVector X(9);
//    X[0] = x(0.98, 1); X[1] = x(0.98, 2); X[2] = x(0.98, 3);
//    X[3] = x(0.99, 1); X[4] = x(0.99, 2); X[5] = x(0.99, 3);
//    X[6] = x(1.00, 1); X[7] = x(1.00, 2); X[8] = x(1.00, 3);
//    IPrinter::print(X, X.length());
//    IPrinter::printSeperatorLine();

    //    double aaa = 0.0;
    //    for (unsigned int i=0; i<9; i++)
    //    {
    //        aaa += F[0][i]*X[i];
    //    }
    //    printf("%.14f\n", aaa);

    DoubleVector xf((k+1)*M);
    LinearEquation::GaussianElimination(F, g, xf);
//    IPrinter::print(xf, xf.length());
//    IPrinter::printSeperatorLine();

//    double aaa = 0.0;
//    for (unsigned int i=0; i<9; i++)
//    {
//        aaa += F[8][i]*xf[i];
//    }
//    printf("%.10f\n", aaa);

    //    //    double norm = 0.0;
    //    //    DoubleVector x(N+1);
    //    //    for (unsigned int n=N; n>=N-k; n--)
    //    //    {
    //    //        x[n] = xf[n-(N-k)];
    //    //        //printf("%4d %.8f %.8f\n", n, x[n], NonLocal::x(n*h));
    //    //        printf("%.8f\n", x[n]);
    //    //        norm += (x[n]-NonLocal::x(n*h))*(x[n]-NonLocal::x(n*h));
    //    //    }
    //    xf.clear();

    F.clear();
    g.clear();

    x.resize(N+1); for (unsigned int n=0; n<=N; n++) x[n].resize(M);

    unsigned int s = xf.length()-M;
    unsigned int e = xf.length()-1;
    for (unsigned int n=N; n>=N-k; n--)
    {
        x[n] = xf.mid(s, e);
        s -= M;
        e -= M;
    }
    xf.clear();

    unsigned int stop = static_cast<unsigned int>(0)-1;
    for (unsigned int n=(N-k)-1; n!=stop; n--)
    {
        x[n] = betta[n][0];
        for (unsigned int j=2; j<=k; j++) x[n] -= betta[n][j]*x[n+j-1];
        for (unsigned int j=n+k; j<=N; j++) x[n] -= C[j]*x[j];
        betta[n][1].inverse();
        x[n] = betta[n][1]*x[n];
    }

    for (unsigned int n=0; n<=end; n++)
    {
        betta[n][0].clear();
        for (unsigned int i=1; i<=k; i++) betta[n][k].clear();
        delete [] betta[n];
    }
    delete [] betta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

NonLocalSystem::~NonLocalSystem() {}

double NonLocalSystem::A(double t UNUSED_PARAM, unsigned int r UNUSED_PARAM, unsigned int c UNUSED_PARAM) const
{
    return 1.0;//cos(t*t)+0.5*sin(r*t);
}

double NonLocalSystem::B(double t UNUSED_PARAM, unsigned int n UNUSED_PARAM) const
{
    if (n == 1) return -(A(t,1,1)*x(t,1)+A(t,1,2)*x(t,2)+A(t,1,3)*x(t,3)) + 1.0;
    if (n == 2) return -(A(t,2,1)*x(t,1)+A(t,2,2)*x(t,2)+A(t,2,3)*x(t,3)) + 2.0*t;
    if (n == 3) return -(A(t,3,1)*x(t,1)+A(t,3,2)*x(t,2)+A(t,3,3)*x(t,3)) + 2.0*t+1.0;

//    if (n == 1) return -(A(t,1,1)*x(t,1)+A(t,1,2)*x(t,2)+A(t,1,3)*x(t,3)) + 6.0*M_PI*cos(6.0*M_PI*t);
//    if (n == 2) return -(A(t,2,1)*x(t,1)+A(t,2,2)*x(t,2)+A(t,2,3)*x(t,3)) - 8.0*M_PI*sin(8.0*M_PI*t);
//    if (n == 3) return -(A(t,3,1)*x(t,1)+A(t,3,2)*x(t,2)+A(t,3,3)*x(t,3)) - 6.0*exp(-6.0*t);
    return NAN;
}

double NonLocalSystem::x(double t UNUSED_PARAM, unsigned int n UNUSED_PARAM) const
{
    if (n == 1) return t;
    if (n == 2) return t*t;
    if (n == 3) return t*t+t;

//    if (n == 1) return sin(6.0*M_PI*t);
//    if (n == 2) return cos(8.0*M_PI*t);
//    if (n == 3) return exp(-6.0*t);


    return NAN;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

NonLocal::~NonLocal() {}

double NonLocal::a(double t UNUSED_PARAM) const
{
    return exp(t);
}

double NonLocal::b(double t UNUSED_PARAM) const
{
    //return 2.0*t-a(t)*x(t);
    //return 2.0*M_PI*cos(2.0*M_PI*t)-a(t)*x(t);
    //return 4.0*M_PI*t*cos(2.0*M_PI*t*t)-a(t)*x(t);
    return 8.0*M_PI*t*cos(4.0*M_PI*t*t)-a(t)*NonLocal::x(t);
    //return 10.0*M_PI*cos(10.0*M_PI*t)-a(t)*x(t);
}

double NonLocal::x(double t) const
{
    //return t*t;
    //return sin(2.0*M_PI*t);
    //return sin(2.0*M_PI*t*t);
    return sin(4.0*M_PI*t*t);
    //return sin(10.0*M_PI*t);
}
