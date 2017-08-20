#include "slodenlcsv2.h"

void SystemLinearODENonLocalContionsV2::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Dimension dim(0.01, 100, 0);
    ODEGrid grid(dim);
    SystemLinearODENonLocalContionsV2 slodenlcv;

    unsigned int n0 = 3;
    unsigned int n1 = 0;
    unsigned int n2 = 0;
    unsigned int n = n0 + n1 + n2;

    Condition nlc0;
    nlc0.time = 0.0;
    nlc0.nmbr = 0;
    nlc0.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc0.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nlc1;
    nlc1.time = 0.2;
    nlc1.nmbr = 20;
    nlc1.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc1.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nlc2;
    nlc2.time = 0.5;
    nlc2.nmbr = 50;
    nlc2.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc2.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nlc3;
    nlc3.time = 0.8;
    nlc3.nmbr = 80;
    nlc3.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc3.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nlc4;
    nlc4.time = 1.0;
    nlc4.nmbr = 100;
    nlc4.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nlc4.mtrx[row][col] = (rand() % 1000) / 1000.0;

    slodenlcv.addCondition(nlc0);
    slodenlcv.addCondition(nlc1);
    slodenlcv.addCondition(nlc2);
    slodenlcv.addCondition(nlc3);
    slodenlcv.addCondition(nlc4);

    //    IPrinter::print(nlc0.mtrx, nlc0.mtrx.rows(), nlc0.mtrx.cols(), 8, 6);
    //    IPrinter::print(nlc1.mtrx, nlc1.mtrx.rows(), nlc1.mtrx.cols(), 8, 6);
    //    IPrinter::print(nlc2.mtrx, nlc2.mtrx.rows(), nlc2.mtrx.cols(), 8, 6);
    //    IPrinter::print(nlc3.mtrx, nlc3.mtrx.rows(), nlc3.mtrx.cols(), 8, 6);
    //    IPrinter::print(nlc4.mtrx, nlc4.mtrx.rows(), nlc4.mtrx.cols(), 8, 6);

    LoadPoint lpt0;
    lpt0.time = 0.3;
    lpt0.nmbr = 30;
    lpt0.mtrx.resize(n,n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) lpt0.mtrx[row][col] = (rand() % 1000) / 1000.0;

    LoadPoint lpt1;
    lpt1.time = 0.6;
    lpt1.nmbr = 60;
    lpt1.mtrx.resize(n,n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) lpt1.mtrx[row][col] = (rand() % 1000) / 1000.0;

    //    IPrinter::print(lpt0.mtrx, lpt0.mtrx.rows(), lpt0.mtrx.cols(), 8, 6);
    //    IPrinter::print(lpt1.mtrx, lpt1.mtrx.rows(), lpt1.mtrx.cols(), 8, 6)

    slodenlcv.addLoadPoint(lpt0);
    slodenlcv.addLoadPoint(lpt1);

    DoubleVector gamma;
    gamma.resize(n);

    slodenlcv.setRightSize(gamma);

    slodenlcv.calculateForward();
}

SystemLinearODENonLocalContionsV2::SystemLinearODENonLocalContionsV2() : ISystemLinearODENonLocalContionsV2()
{}

void SystemLinearODENonLocalContionsV2::fillMatrix(DoubleMatrix &M, DoubleVector &P, unsigned int st,
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
    unsigned int N = grid().dimension().sizeN();
    double h = grid().dimension().step();

    unsigned int k2 = nlscs.size();
    unsigned int k1 = lpnts.size();
    unsigned int n = 3;
    unsigned int n0 = nlscs.at(0).mtrx.rows();

    DoubleVector x0(n*(k1+1)+2);
    std::vector<DoubleVector> rx(n*(k1+1)+2);

    DoubleMatrix M(n*(k1+k2),n*(k1+k2));
    DoubleVector P(n*(k1+k2));

    fillMatrix(M, P, 0, n, k1, k2, 0);
    fillMatrix(M, P, 0, n, k1, k2, 1);
    fillMatrix(M, P, 0, n, k1, k2, 2);
    //IPrinter::print(M, M.rows(), M.cols(), 8, 6);

    for (unsigned int row = 0; row<n0; row++)
    {
        for (unsigned int i=0; i<k2-1; i++)
        {
            /* Current segment and points */
            Condition &sc = nlscs.at(i);
            Condition &ec = nlscs.at(i+1);

            /* Forming arguments of right side of ODE */
            for (unsigned int j=0; j<n; j++) x0[j] = sc.mtrx[row][j];

            for (unsigned int s=0; s<k1; s++)
            {
                LoadPoint &betta = lpnts.at(s);
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
            for (unsigned int j=0; j<n; j++) sc.mtrx[row][j] = rx[j][ec.nmbr-sc.nmbr];
            for (unsigned int s=1; s<=k1; s++)
            {
                LoadPoint &betta = lpnts.at(s);
                for (unsigned int j=0; j<n; j++)
                {
                    betta.mtrx[row][j] = rx[s*n+j][ec.nmbr-sc.nmbr];
                }
            }
            gamma[row] = rx[n*(k1+1)][ec.nmbr-sc.nmbr];
            double M = rx[n*(k1+1)+1][ec.nmbr-sc.nmbr];

            /* Calculating new coifficients */
            for (unsigned int k=i+1; k<k1; k++)
            {
                Condition &cc = nlscs.at(k);
                for (unsigned int j=0; j<n; j++) cc.mtrx[row][j] *= M;
            }

            for (unsigned int j=0; j<n; j++) ec.mtrx[row][j] += rx[j][ec.nmbr-sc.nmbr];

            //fillMatrix(M, P, i, n, k1, k2);
        }
    }
    IPrinter::print(M, M.rows(), M.cols(), 8, 6);
    x0.clear();
    rx.clear();
}

double SystemLinearODENonLocalContionsV2::A(TimeNode node UNUSED_PARAM, unsigned int row UNUSED_PARAM, unsigned int col UNUSED_PARAM) const
{
    double t = node.t;

    if (row==0) { if (col==0) { return t; }   if (col==1) { return 2.0; } if (col==2) { return 4.0; } }
    if (row==1) { if (col==0) { return 5.0; } if (col==1) { return 6.0; } if (col==2) { return 7.0; } }
    if (row==2) { if (col==0) { return 8.0; } if (col==1) { return t; }   if (col==2) { return 1.0; } }
    return NAN;
}

double SystemLinearODENonLocalContionsV2::B(TimeNode node UNUSED_PARAM, unsigned int s UNUSED_PARAM, unsigned int row UNUSED_PARAM, unsigned int col UNUSED_PARAM) const
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

double SystemLinearODENonLocalContionsV2::C(TimeNode node UNUSED_PARAM, unsigned int row) const
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
