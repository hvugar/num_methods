#include "iproblem2hbackward2d.h"
#include <imaging.h>

void IProblem2HBackward2D::calculateMVD(DoubleMatrix &p, vector<ExtendedSpaceNode2DH> &info, bool use) const
{
    Dimension dimX = spaceDimension(Dimension::DimensionX);
    Dimension dimY = spaceDimension(Dimension::DimensionY);
    Dimension time = timeDimension();

    unsigned int N = dimX.sizeN();
    unsigned int M = dimY.sizeN();
    unsigned int L = time.sizeN();

    double hx = dimX.step();
    double hy = dimY.step();
    double ht = time.step();

    double lambda = mEquParameter.lambda;
    //double lambda1 = mEquParameter.lambda1;
    double a = mEquParameter.a;
    unsigned int No = mEquParameter.No;
    unsigned int Nc = mEquParameter.No;
    unsigned int Ns = mEquParameter.Ns;

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p05(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    p.clear(); p.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<IProblem2H2D::ExtendedSpacePointNode> obsDeltaNodes;
    for (unsigned int j=0; j<No; j++) IProblem2H2D::distributeDelta(mOptParameter.xi[j], obsDeltaNodes, j, dimX, dimY);

    std::vector<IProblem2H2D::ExtendedSpacePointNode> cntPointNodes;
    for (unsigned int i=0; i<Nc; i++) IProblem2H2D::distributeDelta(mOptParameter.eta[i], cntPointNodes, i, dimX, dimY);

    SpaceNodePDE sn;

    //--------------------------------------------------------------------------------------------//
    vector<unsigned int> rows0;
    vector<unsigned int> rows1;
    vector<unsigned int> rows2;
    for (unsigned int ny=0; ny<=M; ny++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int j=0; j<obsDeltaNodes.size(); j++)
        {
            const IProblem2H2D::ExtendedSpacePointNode &odn = obsDeltaNodes.at(j);
            if (odn.j == ny)
            {
                found1 = true;
                for (unsigned int i=0; i<cntPointNodes.size(); i++)
                {
                    const IProblem2H2D::ExtendedSpacePointNode &cpn = cntPointNodes.at(i);
                    if (cpn.j == ny)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == true  && found2 == true)  rows2.push_back(ny);
        if (found1 == true)                     rows1.push_back(ny);
        if (found1 == false && found2 == false) rows0.push_back(ny);
    }

    //--------------------------------------------------------------------------------------------//
    vector<unsigned int> cols0;
    vector<unsigned int> cols1;
    vector<unsigned int> cols2;
    for (unsigned int nx=0; nx<=N; nx++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int j=0; j<obsDeltaNodes.size(); j++)
        {
            const IProblem2H2D::ExtendedSpacePointNode &odn = obsDeltaNodes.at(j);
            if (odn.i == nx)
            {
                found1 = true;
                for (unsigned int i=0; i<cntPointNodes.size(); i++)
                {
                    const IProblem2H2D::ExtendedSpacePointNode &cpn = cntPointNodes.at(i);
                    if (cpn.i == nx)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == true  && found2 == true)  cols2.push_back(nx);
        if (found1 == true)                     cols1.push_back(nx);
        if (found1 == false && found2 == false) cols0.push_back(nx);
    }
    //--------------------------------------------------------------------------------------------//

    //-------------------------------------------- info --------------------------------------------//
    if (use == true)
    {
        info.resize(mEquParameter.Nc);
        for (unsigned int i=0; i<mEquParameter.Nc; i++)
        {
            ExtendedSpaceNode2DH &e = info[i];
            e.setSpaceNode(mOptParameter.eta[i]);
            e.id = i;
            e.extendWeights(dimX, dimY);
            e.extendLayers(L+1);
        }
    }
    //-------------------------------------------- info --------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p00[m][n] = initial1(sn);
        }
    }
    if (use == true) add2Info(p00, info, L);
    layerInfo(p00, L);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;

            if (n==0)      sum += a*a*(p00[m][n]-2.0*p00[m][n+1]+p00[m][n+2])/(hx*hx);
            else if (n==N) sum += a*a*(p00[m][n-2]-2.0*p00[m][n-1]+p00[m][n])/(hx*hx);
            else           sum += a*a*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1])/(hx*hx);

            if (m==0)      sum += a*a*(p00[m][n]-2.0*p00[m+1][n]+p00[m+2][n])/(hy*hy);
            else if (m==M) sum += a*a*(p00[m-2][n]-2.0*p00[m-1][n]+p00[m][n])/(hy*hy);
            else           sum += a*a*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n])/(hy*hy);

            sum += lambda*initial2(sn);

            for (unsigned int odn=0; odn<obsDeltaNodes.size(); odn++)
            {
                IProblem2H2D::ExtendedSpacePointNode obsNode = obsDeltaNodes.at(odn);
                if (obsNode.i == n && obsNode.j == m)
                {
                    for (unsigned int cpn=0; cpn<cntPointNodes.size(); cpn++)
                    {
                        IProblem2H2D::ExtendedSpacePointNode cntNode = cntPointNodes.at(cpn);
                        if (cntNode.i == n && cntNode.j == m)
                        {
                            sum += mOptParameter.k[cntNode.id][obsNode.id] * p00[m][n]*(cntNode.w*(hx*hy)) * obsNode.w;
                        }
                    }
                }
            }

            p05[m][n] = p00[m][n] - initial2(sn)*ht*0.5 + 0.125*ht*ht*sum;
            p10[m][n] = p00[m][n] - initial2(sn)*ht     + 0.500*ht*ht*sum;
        }
    }

    if (use == true) add2Info(p10, info, L-1);
    layerInfo(p10, L-1);
    //------------------------------------- initial conditions -------------------------------------//

    double *a1X = (double *) malloc(sizeof(double)*(N+1));
    double *b1X = (double *) malloc(sizeof(double)*(N+1));
    double *c1X = (double *) malloc(sizeof(double)*(N+1));
    double *d1X = (double *) malloc(sizeof(double)*(N+1));
    double *x1X = (double *) malloc(sizeof(double)*(N+1));

    double *a1Y = (double *) malloc(sizeof(double)*(M+1));
    double *b1Y = (double *) malloc(sizeof(double)*(M+1));
    double *c1Y = (double *) malloc(sizeof(double)*(M+1));
    double *d1Y = (double *) malloc(sizeof(double)*(M+1));
    double *x1Y = (double *) malloc(sizeof(double)*(M+1));

    //TimeNodePDE tn;

    for (unsigned int l1=2; l1<=L; l1++)
    {
        unsigned int l = L-l1;

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        //tn.i = l; tn.t = tn.i*ht;

        //--------------------------------------------------------------------------//
        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;

                if (m == 0)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                if (m>0 && m<M) d1X[n] = ((a*a*ht*ht)/(hy*hy))*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                if (m == M)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                d1X[n] += 0.5*(p10[m][n]-p00[m][n]) + p10[m][n];
                d1X[n] += 0.5*lambda*ht*(4.0*p10[m][n]-p05[m][n]);

                if (n==0)
                {
                    a1X[0] = 0.0;
                    b1X[0] = +1.0 + (a*a*ht*ht)/(hx*hx) + 1.5*(lambda*ht);
                    c1X[0] = -(a*a*ht*ht)/(hx*hx);
                    //d1X[0] += 0.0;
                }
                else if(n == N)
                {
                    a1X[N] = -(a*a*ht*ht)/(hx*hx);
                    b1X[N] = +1.0 + (a*a*ht*ht)/(hx*hx) + 1.5*(lambda*ht);
                    c1X[N] = 0.0;
                    //d1X[N] += 0.0;
                }
                else
                {
                    a1X[n] = -0.5*(a*a*ht*ht)/(hx*hx);
                    b1X[n] = +1.0 + (a*a*ht*ht)/(hx*hx) + 1.5*(lambda*ht);
                    c1X[n] = -0.5*(a*a*ht*ht)/(hx*hx);
                    //d1X[n] += 0.0;
                }
            }
            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
            for (unsigned int n=0; n<=N; n++) p15[m][n] = x1X[n];
        }

        if (rows2.size() == 0)
        {
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (m == 0)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    if (m>0 && m<M) d1X[n] = ((a*a*ht*ht)/(hy*hy))*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    if (m == M)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    d1X[n] += 0.5*(p10[m][n]-p00[m][n]) + p10[m][n];
                    d1X[n] += 0.5*lambda*ht*(4.0*p10[m][n]-p05[m][n]);

                    if (n==0)
                    {
                        a1X[0] = 0.0;
                        b1X[0] = +1.0 + (a*a*ht*ht)/(hx*hx) + 1.5*(lambda*ht);
                        c1X[0] = -(a*a*ht*ht)/(hx*hx);
                        //d1X[0] += 0.0;
                    }
                    else if(n==N)
                    {
                        a1X[N] = -(a*a*ht*ht)/(hx*hx);
                        b1X[N] = +1.0 + (a*a*ht*ht)/(hx*hx) + 1.5*(lambda*ht);
                        c1X[N] = 0.0;
                        //d1X[N] += 0.0;
                    }
                    else
                    {
                        a1X[n] = -0.5*(a*a*ht*ht)/(hx*hx);
                        b1X[n] = +1.0 + (a*a*ht*ht)/(hx*hx) + 1.5*(lambda*ht);
                        c1X[n] = -0.5*(a*a*ht*ht)/(hx*hx);
                        //d1X[n] += 0.0;
                    }


                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &cpn = cntPointNodes[cni];
                                d1X[n] += 0.5*ht*ht * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                for (unsigned int n=0; n<=N; n++) p15[m][n] = x1X[n];
            }
        }
        else
        {
            /////////////////////////////////////
        }
        //--------------------------------------------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        //tn.i = l+1; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int col=0; col<cols0.size(); col++)
        {
            unsigned int n = cols0.at(col);
            sn.i = n; sn.x = n*hx;
            for (unsigned int m=0; m<=M; m++)
            {
                sn.j = m; sn.y = m*hy;

                if (n==0)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                if (n>0 && n<N) d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                if (n==N)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                d1Y[m] += 0.5*(p10[m][n]-p00[m][n]) + p15[m][n];
                d1Y[m] += 0.5*lambda*ht*(4.0*p15[m][n]-p10[m][n]);

                if (m == 0)
                {
                    a1Y[0] = 0.0;
                    b1Y[0] = +1.0 + (a*a*ht*ht)/(hy*hy) + 1.5*(lambda*ht);
                    c1Y[0] = -(a*a*ht*ht)/(hy*hy);
                    //d1Y[0] += 0.0;
                }
                else if (m == M)
                {
                    a1Y[M] = -(a*a*ht*ht)/(hy*hy);
                    b1Y[M] = +1.0 + (a*a*ht*ht)/(hy*hy) + 1.5*(lambda*ht);
                    c1Y[M] = 0.0;
                    //d1Y[M] += 0.0;
                }
                else
                {
                    a1Y[m] = -0.5*(a*a*ht*ht)/(hy*hy);
                    b1Y[m] = +1.0 + (a*a*ht*ht)/(hy*hy) + 1.5*(lambda*ht);
                    c1Y[m] = -0.5*(a*a*ht*ht)/(hy*hy);
                    //d1Y[m] += 0.0;
                }
            }
            tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
            for (unsigned int m=0; m<=M; m++) p[m][n] = x1Y[m];
        }

        if (cols2.size() == 0)
        {
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    if (n==0)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    if (n>0 && n<N) d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    if (n==N)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    d1Y[m] += 0.5*(p10[m][n]-p00[m][n]) + p15[m][n];
                    d1Y[m] += 0.5*lambda*ht*(4.0*p15[m][n]-p10[m][n]);

                    if (m == 0)
                    {
                        a1Y[0] = 0.0;
                        b1Y[0] = +1.0 + (a*a*ht*ht)/(hy*hy) + 1.5*(lambda*ht);
                        c1Y[0] = -(a*a*ht*ht)/(hy*hy);
                        //d1Y[0] += 0.0;
                    }
                    else if (m == M)
                    {
                        a1Y[M] = -(a*a*ht*ht)/(hy*hy);
                        b1Y[M] = +1.0 + (a*a*ht*ht)/(hy*hy) + 1.5*(lambda*ht);
                        c1Y[M] = 0.0;
                        //d1Y[M] += 0.0;
                    }
                    else
                    {
                        a1Y[m] = -0.5*(a*a*ht*ht)/(hy*hy);
                        b1Y[m] = +1.0 + (a*a*ht*ht)/(hy*hy) + 1.5*(lambda*ht);
                        c1Y[m] = -0.5*(a*a*ht*ht)/(hy*hy);
                        //d1Y[m] += 0.0;
                    }

                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &cpn = cntPointNodes.at(cni);
                                d1Y[m] += 0.5*ht*ht * mOptParameter.k[cpn.id][odn.id] * p[cpn.j][cpn.i]* (cpn.w * (hx*hy)) * odn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
                for (unsigned int m=0; m<=M; m++) p[m][n] = x1Y[m];
            }
        }
        else
        {
            ////////////////////////////////////
        }
        //--------------------------------------------------------------------------//
        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p[m][n];
                p05[m][n] = p15[m][n];
            }
        }

        if (use == true) add2Info(p, info, l);
        layerInfo(p, l);
    }

    free(x1X);
    free(d1X);
    free(c1X);
    free(b1X);
    free(a1X);

    free(x1Y);
    free(d1Y);
    free(c1Y);
    free(b1Y);
    free(a1Y);

    rows0.clear();
    rows1.clear();
    rows2.clear();

    cols0.clear();
    cols1.clear();
    cols2.clear();

    obsDeltaNodes.clear();
    cntPointNodes.clear();

    p00.clear();
    p05.clear();
    p10.clear();
    p15.clear();
}

void IProblem2HBackward2D::layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int layerNumber UNUSED_PARAM) const
{
    QPixmap px;
    visualizeMatrixHeat(u, u.min(), u.max(), px);
    char buffer[30] = {0};
    int c = sprintf(buffer, "images/b/%4d.png", layerNumber);
    buffer[c] = 0;
    px.save(buffer, "png", 0);
}

double IProblem2HBackward2D::initial1(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return -2.0*ifunc->alpha1*(UTt[sn.j][sn.i]-ifunc->V1[sn.j][sn.i]);
}

double IProblem2HBackward2D::initial2(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return 2.0*ifunc->alpha0*(UT[sn.j][sn.i]-ifunc->V0[sn.j][sn.i]) + mEquParameter.lambda * initial1(sn);
}

double IProblem2HBackward2D::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType) const
{
    return 0.0;
}

double IProblem2HBackward2D::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 0.0;
}

void IProblem2HBackward2D::add2Info(const DoubleMatrix &p, vector<ExtendedSpaceNode2DH> &info, unsigned int ln) const
{
    for (unsigned int i=0; i<mEquParameter.Nc; i++)
    {
        ExtendedSpaceNode2DH &pi = info[i];
        for (unsigned int r=0; r<4; r++)
        {
            for (unsigned int c=0; c<4; c++)
            {
                unsigned int x = pi.wi[r][c].i;
                unsigned int y = pi.wi[r][c].j;
                pi.wi[r][c].u[ln] = p[y][x];
            }
        }
    }
}
