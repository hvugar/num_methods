#include "slodenlcsv2.h"
#include <printer.h>
#include <math.h>

void SystemLinearODENonLocalContionsV2::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    SystemLinearODENonLocalContionsV2 slodenlcv;
    slodenlcv.setGrid(UniformODEGrid(0.1, 10, 0));

    unsigned int n0 = 3;
    unsigned int n1 = 0;
    unsigned int n2 = 0;
    unsigned int n = n0 + n1 + n2;

    /* Adding conditions */

    Condition nlc0;
    nlc0.time = 0.0;
    nlc0.nmbr = 0;
    nlc0.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc0.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nlc1;
    nlc1.time = 0.2;
    nlc1.nmbr = 2;
    nlc1.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc1.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nlc2;
    nlc2.time = 0.5;
    nlc2.nmbr = 5;
    nlc2.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc2.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nlc3;
    nlc3.time = 0.8;
    nlc3.nmbr = 8;
    nlc3.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc3.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nlc4;
    nlc4.time = 1.0;
    nlc4.nmbr = 10;
    nlc4.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc4.mtrx[row][col] = (rand() % 1000) / 1000.0;

    slodenlcv.addCondition(nlc0);
    slodenlcv.addCondition(nlc1);
    slodenlcv.addCondition(nlc2);
    slodenlcv.addCondition(nlc3);
    slodenlcv.addCondition(nlc4);

    /* Adding loading points */

    LoadPoint lpt0;
    lpt0.time = 0.3;
    lpt0.nmbr = 3;
    lpt0.mtrx.resize(n,n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) lpt0.mtrx[row][col] = (rand() % 1000) / 1000.0;

    LoadPoint lpt1;
    lpt1.time = 0.6;
    lpt1.nmbr = 6;
    lpt1.mtrx.resize(n,n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) lpt1.mtrx[row][col] = (rand() % 1000) / 1000.0;

    slodenlcv.addLoadPoint(lpt0);
    slodenlcv.addLoadPoint(lpt1);

    /* Calculating gamma */

    DoubleVector gamma;
    gamma.resize(n);

    for (unsigned int row=0; row<n; row++)
    {
        gamma[row] = 0.0;

        unsigned int k2 = slodenlcv.nlscs.size();
        for (unsigned int i=0; i<k2; i++)
        {
            Condition &cc = slodenlcv.nlscs[i];
            for (unsigned int col=0; col<n; col++)
            {
                DoubleMatrix &m = cc.mtrx;
                gamma[row] += m[row][col]*slodenlcv.X(cc.time,col);
            }
        }

        unsigned int k1 = slodenlcv.lpnts.size();
        for (unsigned int i=0; i<k1; i++)
        {
            LoadPoint &lp = slodenlcv.lpnts[i];
            for (unsigned int col=0; col<n; col++)
            {
                DoubleMatrix &m = lp.mtrx;
                gamma[row] += m[row][col]*slodenlcv.X(lp.time,col);
            }
        }
    }

    slodenlcv.setRightSize(gamma);

    slodenlcv.calculateForward2();
}

SystemLinearODENonLocalContionsV2::SystemLinearODENonLocalContionsV2() : ISystemLinearODENonLocalContionsV2()
{}

void SystemLinearODENonLocalContionsV2::fillMatrix(DoubleMatrix &M, DoubleVector &P UNUSED_PARAM, unsigned int st,
                                                   unsigned int n, unsigned int k1, unsigned int k2,
                                                   unsigned int row)
{
    for (unsigned int i=st; i<k2; i++)
    {
        Condition &cn = nlscs.at(i);
        for (unsigned col=0; col<n; col++)
        {
            M[row][i*n+col] = cn.mtrx.at(row, col);
        }
    }
    for (unsigned int s=0; s<k1; s++)
    {
        LoadPoint &lp = lpnts.at(s);
        for (unsigned col=0; col<n; col++)
        {
            M[row][k2*n+s*n+col] = lp.mtrx.at(row, col);
        }
    }
    //P[row] = gamma[row];
}

void SystemLinearODENonLocalContionsV2::calculateForward()
{
    unsigned int N = grid().sizeN();
    double h = grid().step();

    unsigned int k2 = nlscs.size();
    unsigned int k1 = lpnts.size();
    unsigned int n = 3;
    unsigned int n0 = nlscs.at(0).mtrx.rows();

    DoubleVector x0(n*(k1+1)+2);
    std::vector<DoubleVector> rx(n*(k1+1)+2);

    DoubleMatrix MX(n*(k1+k2),n*(k1+k2));
    DoubleVector RS(n*(k1+k2));

    fillMatrix(MX, RS, 0, n, k1, k2, 0); fflush(stdout);
    fillMatrix(MX, RS, 0, n, k1, k2, 1); fflush(stdout);
    fillMatrix(MX, RS, 0, n, k1, k2, 2); fflush(stdout);
    //IPrinter::print(M, M.rows(), M.cols(), 8, 6);

    /* Matrix */
    for (unsigned int row = 0; row < n; row++)
    {
        for (unsigned int i=0; i<k2; i++)
        {
            Condition &cn = nlscs.at(i);
            for (unsigned col=0; col<n; col++) MX[row][i*n+col] = cn.mtrx.at(row, col);
        }
        for (unsigned int s=0; s<k1; s++)
        {
            LoadPoint &lp = lpnts.at(s);
            for (unsigned col=0; col<n; col++) MX[row][k2*n+s*n+col] = lp.mtrx.at(row, col);
        }
        RS[row] = gamma[row];
    }

    unsigned int ROW = n;
    for (unsigned int i=0; i<k2-1; i++)
    {
        for (unsigned int row = 0; row<n0; row++)
        {
            /* Current segment and points */
            Condition &sc = nlscs.at(i);
            Condition &ec = nlscs.at(i+1);

            /* Forming arguments of right side of ODE */
            for (unsigned int j=0; j<n; j++) x0[j] = sc.mtrx[row][j];

            for (unsigned int s=1; s<=k1; s++)
            {
                LoadPoint &betta = lpnts.at(s-1);
                for (unsigned int j=0; j<n; j++)
                {
                    x0[s*n+j] = betta.mtrx[row][j];
                }
            }
            x0[n*(k1+1)] = gamma[row];
            x0[n*(k1+1)+1] = 1.0;

            /* Calculating differensial equation */
            calculateCauchyProblem(sc, ec, x0, rx, h);

            /* Copying results */
            unsigned int NUM = ec.nmbr - sc.nmbr;
            for (unsigned int j=0; j<n; j++) sc.mtrx[row][j] = rx[j][NUM];
            for (unsigned int s=1; s<=k1; s++)
            {
                LoadPoint &betta = lpnts.at(s-1);
                for (unsigned int j=0; j<n; j++)
                {
                    betta.mtrx[row][j] = rx[s*n+j][NUM];
                }
            }
            gamma[row] = rx[n*(k1+1)][NUM];
            double M = rx[n*(k1+1)+1][NUM];

            /* Calculating new coifficients */
            for (unsigned int k=i+1; k<k1; k++)
            {
                Condition &cc = nlscs.at(k);
                for (unsigned int j=0; j<n; j++) cc.mtrx[row][j] *= M;
            }

            for (unsigned int j=0; j<n; j++) ec.mtrx[row][j] += rx[j][NUM];
        }

        for (unsigned int p=0; p<n; p++)
        {
            for (unsigned int j=i+1; j<k2; j++)
            {
                Condition &cn = nlscs.at(j);
                for (unsigned col=0; col<n; col++) MX[ROW][j*n+col] = cn.mtrx.at(ROW%n, col);
            }
            for (unsigned int s=0; s<k1; s++)
            {
                LoadPoint &lp = lpnts.at(s);
                for (unsigned col=0; col<n; col++) MX[ROW][k2*n+s*n+col] = lp.mtrx.at(ROW%n, col);
            }

            RS[ROW] = gamma[ROW%n];
            ROW++;
        }
    }
    IPrinter::print(MX, MX.rows(), MX.cols(), 10, 6);
    x0.clear();
    rx.clear();
}

void SystemLinearODENonLocalContionsV2::calculateForward2()
{
    unsigned int N = grid().sizeN();
    double h = grid().step();
    unsigned int k2 = nlscs.size(); // number of condition points
    unsigned int k1 = lpnts.size(); // number of loaded points
    unsigned int n = nlscs.at(0).mtrx.rows();

    DoubleMatrix MX(n*(k1+k2),n*(k1+k2));
    DoubleVector GM(n*(k1+k2));

    std::vector<DoubleVector> mx((k1+k2)*n+2);
    for (unsigned int i=0; i<mx.size(); i++) mx[i].resize(N+1);

    unsigned int ROW = 0;
    unsigned int row = 0;
//    for (unsigned int row=0; row<n; row++)
    {
        /* Initializing alpha-s first elements */
        for (unsigned int i=0; i<k2; i++)
        {
            const Condition &cond = nlscs.at(i);
            for (unsigned int j=0; j<n; j++) mx[i*n+j][0] = cond.mtrx[row][j];
        }

        /* Initializing betta-s first elements */
        for (unsigned int s=0; s<k1; s++)
        {
            const LoadPoint &ldpnt = lpnts.at(s);
            for (unsigned int j=0; j<n; j++) mx[k2*n+s*n+j][0] = ldpnt.mtrx[row][j];
        }

        mx[(k1+k2)*n][0] = gamma[row];
        mx[(k1+k2)*n+1][0] = 1.0;

//        for (unsigned int i=0; i<(k1+k2)*n; i++)
//        {
//            MX[ROW][i] = mx[i][0];
//        }
//        ROW++;

        // Forming initial vector for Cauchy problem
        for (unsigned int cnd=0; cnd<k2-1; cnd++)
        {
            const Condition &sc = nlscs.at(cnd);
            const Condition &ec = nlscs.at(cnd+1);

            DoubleVector x0(n+n*k1+2);
            std::vector<DoubleVector> rx(n+n*k1+2);

            /* Forming arguments of right side of ODE */
            for (unsigned int i=0; i<n; i++) x0[i] = mx[cnd*n+i][sc.nmbr];

            for (unsigned int s=0; s<k1; s++)  //////// betta
            {
                for (unsigned int i=0; i<n; i++)
                {
                    x0[n+s*n+i] = mx[k2*n+s*n+i][sc.nmbr];
                }
            }
            x0[n+k1*n]   = mx[k2*n+k1*n][sc.nmbr]; // gamma
            x0[n+k1*n+1] = 1.0;                    // M

            /* Calculating differensial equation */
            calculateCauchyProblem(sc, ec, x0, rx, h);

            /* Current alpha */
            for (unsigned int i=0; i<n; i++)
            {
                for (unsigned int j=sc.nmbr; j<=ec.nmbr; j++) mx[cnd*n+i][j] = rx[i][j-sc.nmbr];
            }

            /* Bettas */
            for (unsigned int s=0; s<k1; s++)
            {
                for (unsigned int i=0; i<n; i++)
                {
                    for (unsigned int j=sc.nmbr; j<=ec.nmbr; j++) mx[k2*n+s*n+i][j] = rx[n+s*n+i][j-sc.nmbr];
                }
            }

            /* gamma */
            for (unsigned int j=sc.nmbr; j<=ec.nmbr; j++) mx[k2*n+k1*n][j] = rx[n+k1*n][j-sc.nmbr];
            /* M */
            for (unsigned int j=sc.nmbr; j<=ec.nmbr; j++) mx[k2*n+k1*n+1][j] = rx[n+k1*n+1][j-sc.nmbr];

            for (unsigned int s=cnd+1; s<k2; s++)
            {
                for (unsigned int i=0; i<n; i++)
                {
                    for (unsigned int j=sc.nmbr; j<=ec.nmbr; j++) mx[s*n+i][j] = mx[s*n+i][sc.nmbr]*rx[n+k1*n+1][j-sc.nmbr];
                }
            }
            for (unsigned int i=0; i<n; i++) mx[(cnd+1)*n+i][ec.nmbr] += mx[cnd*n+i][ec.nmbr];

            x0.clear();
            for (unsigned int i=0; i<rx.size(); i++) rx[i].clear();
            rx.clear();
        }

        IPrinter::printSeperatorLine("alpha_0");
        IPrinter::printVector(mx[0], "mx0 ");
        IPrinter::printVector(mx[1], "mx1 ");
        IPrinter::printVector(mx[2], "mx2 ");
        IPrinter::printSeperatorLine("alpha_1");
        IPrinter::printVector(mx[3], "mx3 ");
        IPrinter::printVector(mx[4], "mx4 ");
        IPrinter::printVector(mx[5], "mx5 ");
        IPrinter::printSeperatorLine("alpha_2");
        IPrinter::printVector(mx[6], "mx6 ");
        IPrinter::printVector(mx[7], "mx7 ");
        IPrinter::printVector(mx[8], "mx8 ");
        IPrinter::printSeperatorLine("alpha_3");
        IPrinter::printVector(mx[9], "mx9 ");
        IPrinter::printVector(mx[10], "mx10");
        IPrinter::printVector(mx[11], "mx11");
        IPrinter::printSeperatorLine("alpha_4");
        IPrinter::printVector(mx[12], "mx12");
        IPrinter::printVector(mx[13], "mx13");
        IPrinter::printVector(mx[14], "mx14");
        IPrinter::printSeperatorLine("alpha_5");
        IPrinter::printVector(mx[15], "mx15");
        IPrinter::printVector(mx[16], "mx16");
        IPrinter::printVector(mx[17], "mx17");
        IPrinter::printSeperatorLine("alpha_6");
        IPrinter::printVector(mx[18], "mx18");
        IPrinter::printVector(mx[19], "mx19");
        IPrinter::printVector(mx[20], "mx20");
        IPrinter::printSeperatorLine();


        for (unsigned int i=0; i<k2; i++)
        {
            const Condition &sc = nlscs.at(i);
            for (unsigned int j=i*n; j<(k1+k2)*n; j++) MX[ROW][j] = mx[j][sc.nmbr];
            GM[ROW] = mx[(k1+k2)*n+1][sc.nmbr];
            ROW++;
        }

//        for (unsigned int s=0; s<k1; s++)
//        {
//            const LoadPoint &lp = lpnts.at(s);
//            for (unsigned int j=0; j<(k1+k2)*n; j++) MX[ROW][j] = mx[j][lp.nmbr];
//            GM[ROW] = mx[(k1+k2)*n+1][lp.nmbr];
//            ROW++;
//        }



//        for (unsigned int i=0; i<mx.size(); i++)
//        {
//            IPrinter::printVector(mx[i],NULL);
//        }

        puts("---");
//        {
//            int kk = 0;
//            double aa = mx[0][kk] *X(0.0, 0) + mx[1][kk] *X(0.0, 1) + mx[2][kk]*X(0.0, 2) +
//                    mx[3][kk] *X(0.2, 0) + mx[4][kk] *X(0.2, 1) + mx[5][kk]*X(0.2, 2) +
//                    mx[6][kk] *X(0.5, 0) + mx[7][kk] *X(0.5, 1) + mx[8][kk]*X(0.5, 2) +
//                    mx[9][kk] *X(0.8, 0) + mx[10][kk]*X(0.8, 1) + mx[11][kk]*X(0.8, 2) +
//                    mx[12][kk]*X(1.0, 0) + mx[13][kk]*X(1.0, 1) + mx[14][kk]*X(1.0, 2) +
//                    mx[15][kk]*X(0.3, 0) + mx[16][kk]*X(0.3, 1) + mx[17][kk]*X(0.3, 2) +
//                    mx[18][kk]*X(0.6, 0) + mx[19][kk]*X(0.6, 1) + mx[20][kk]*X(0.6, 2);
//            printf("0   %14.10f %14.10f\n", aa, mx[21][kk]);
//        }
//        {
//            int kk = 2;
//            double aa =
//                    //mx[0][kk] *X(0.0, 0) + mx[1][kk] *X(0.0, 1) + mx[2][kk]*X(0.0, 2) +
//                    mx[3][kk] *X(0.2, 0) + mx[4][kk] *X(0.2, 1) + mx[5][kk]*X(0.2, 2) +
//                    mx[6][kk] *X(0.5, 0) + mx[7][kk] *X(0.5, 1) + mx[8][kk]*X(0.5, 2) +
//                    mx[9][kk] *X(0.8, 0) + mx[10][kk]*X(0.8, 1) + mx[11][kk]*X(0.8, 2) +
//                    mx[12][kk]*X(1.0, 0) + mx[13][kk]*X(1.0, 1) + mx[14][kk]*X(1.0, 2) +
//                    mx[15][kk]*X(0.3, 0) + mx[16][kk]*X(0.3, 1) + mx[17][kk]*X(0.3, 2) +
//                    mx[18][kk]*X(0.6, 0) + mx[19][kk]*X(0.6, 1) + mx[20][kk]*X(0.6, 2);
//            printf("20  %14.10f %14.10f\n", aa, mx[21][kk]);
//        }
//        {
//            int kk = 3;
//            double aa =
//                    mx[3][kk] *X(0.3, 0) + mx[4][kk] *X(0.3, 1) + mx[5][kk]*X(0.3, 2) +
//                    mx[6][kk] *X(0.5, 0) + mx[7][kk] *X(0.5, 1) + mx[8][kk]*X(0.5, 2) +
//                    mx[9][kk] *X(0.8, 0) + mx[10][kk]*X(0.8, 1) + mx[11][kk]*X(0.8, 2) +
//                    mx[12][kk]*X(1.0, 0) + mx[13][kk]*X(1.0, 1) + mx[14][kk]*X(1.0, 2) +
//                    mx[15][kk]*X(0.3, 0) + mx[16][kk]*X(0.3, 1) + mx[17][kk]*X(0.3, 2) +
//                    mx[18][kk]*X(0.6, 0) + mx[19][kk]*X(0.6, 1) + mx[20][kk]*X(0.6, 2);
//            printf("30  %14.10f %14.10f\n", aa, mx[21][kk]);
//        }
//        {
//            int kk = 5;
//            double aa =
//                    mx[6][kk] *X(0.5, 0) + mx[7][kk] *X(0.5, 1) + mx[8][kk]*X(0.5, 2) +
//                    mx[9][kk] *X(0.8, 0) + mx[10][kk]*X(0.8, 1) + mx[11][kk]*X(0.8, 2) +
//                    mx[12][kk]*X(1.0, 0) + mx[13][kk]*X(1.0, 1) + mx[14][kk]*X(1.0, 2) +
//                    mx[15][kk]*X(0.3, 0) + mx[16][kk]*X(0.3, 1) + mx[17][kk]*X(0.3, 2) +
//                    mx[18][kk]*X(0.6, 0) + mx[19][kk]*X(0.6, 1) + mx[20][kk]*X(0.6, 2);
//            printf("50  %14.10f %14.10f\n", aa, mx[21][kk]);
//        }
//        {
//            int kk = 6;
//            double aa =
//                    mx[6][kk] *X(0.6, 0) + mx[7][kk] *X(0.6, 1) + mx[8][kk]*X(0.6, 2) +
//                    mx[9][kk] *X(0.8, 0) + mx[10][kk]*X(0.8, 1) + mx[11][kk]*X(0.8, 2) +
//                    mx[12][kk]*X(1.0, 0) + mx[13][kk]*X(1.0, 1) + mx[14][kk]*X(1.0, 2) +
//                    mx[15][kk]*X(0.3, 0) + mx[16][kk]*X(0.3, 1) + mx[17][kk]*X(0.3, 2) +
//                    mx[18][kk]*X(0.6, 0) + mx[19][kk]*X(0.6, 1) + mx[20][kk]*X(0.6, 2);
//            printf("60  %14.10f %14.10f\n", aa, mx[21][kk]);
//        }
//        {
//            int kk = 8;
//            double aa =
//                    mx[9][kk] *X(0.8, 0) + mx[10][kk]*X(0.8, 1) + mx[11][kk]*X(0.8, 2) +
//                    mx[12][kk]*X(1.0, 0) + mx[13][kk]*X(1.0, 1) + mx[14][kk]*X(1.0, 2) +
//                    mx[15][kk]*X(0.3, 0) + mx[16][kk]*X(0.3, 1) + mx[17][kk]*X(0.3, 2) +
//                    mx[18][kk]*X(0.6, 0) + mx[19][kk]*X(0.6, 1) + mx[20][kk]*X(0.6, 2);
//            printf("80  %14.10f %14.10f\n", aa, mx[21][kk]);
//        }
//        {
//            int kk = 10;
//            double aa =
//                    mx[12][kk]*X(1.0, 0) + mx[13][kk]*X(1.0, 1) + mx[14][kk]*X(1.0, 2) +
//                    mx[15][kk]*X(0.3, 0) + mx[16][kk]*X(0.3, 1) + mx[17][kk]*X(0.3, 2) +
//                    mx[18][kk]*X(0.6, 0) + mx[19][kk]*X(0.6, 1) + mx[20][kk]*X(0.6, 2);
//            printf("100 %14.10f %14.10f\n", aa, mx[21][kk]);
//        }
    }


    IPrinter::print(MX, MX.rows(), MX.cols(), 10, 6);
}

double SystemLinearODENonLocalContionsV2::A(TimeNodePDE node UNUSED_PARAM, unsigned int row UNUSED_PARAM, unsigned int col UNUSED_PARAM) const
{
    double t = node.t;

    if (row==0) { if (col==0) { return t; }   if (col==1) { return 2.0; } if (col==2) { return 4.0; } }
    if (row==1) { if (col==0) { return 5.0; } if (col==1) { return 6.0; } if (col==2) { return 7.0; } }
    if (row==2) { if (col==0) { return 8.0; } if (col==1) { return t; }   if (col==2) { return 1.0; } }

    return NAN;
}

double SystemLinearODENonLocalContionsV2::B(TimeNodePDE node UNUSED_PARAM, unsigned int s UNUSED_PARAM, unsigned int row UNUSED_PARAM, unsigned int col UNUSED_PARAM) const
{
    if (s == 0)
    {
        if (row==0) { if (col==0) { return 0.001; } if (col==1) { return 0.008; } if (col==2) { return 0.005; } }
        if (row==1) { if (col==0) { return 0.004; } if (col==1) { return 0.005; } if (col==2) { return 0.007; } }
        if (row==2) { if (col==0) { return 0.009; } if (col==1) { return 0.004; } if (col==2) { return 0.009; } }
    }
    if (s == 1)
    {
        if (row==0) { if (col==0) { return 0.009; } if (col==1) { return 0.007; } if (col==2) { return 0.001; } }
        if (row==1) { if (col==0) { return 0.004; } if (col==1) { return 0.009; } if (col==2) { return 0.003; } }
        if (row==2) { if (col==0) { return 0.002; } if (col==1) { return 0.001; } if (col==2) { return 0.002; } }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsV2::C(TimeNodePDE node UNUSED_PARAM, unsigned int row) const
{
    double t = node.t;

    if (row==0) { return 6.0*cos(6.0*t) - t*sin(6.0*t) - 2.0*cos(4.0*t) - 4.0*t*t + 4.0*t - (-0.0065617303); }
    if (row==1) { return -4.0*sin(4.0*t) - 5.0*sin(6.0*t) - 6.0*cos(4.0*t) - 7.0*t*t + 7.0*t - (-0.0048894459);  }
    if (row==2) { return 2.0*t - 1.0 - 8.0*sin(6.0*t) - t*cos(4.0*t) - t*t + t - (+0.0062216251); }
    return NAN;
}

double SystemLinearODENonLocalContionsV2::X(double t, unsigned int row) const
{
    if (row==0) { return sin(6.0*t); }
    if (row==1) { return cos(4.0*t); }
    if (row==2) { return t*t-t; }
    return NAN;
}
