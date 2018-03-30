#include "iproblem2hforward2d.h"
#include <imaging.h>

void IProblem2HForward2D::calculateMVD(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const
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
    //double lambda1 = mParameter.lambda1;
    double a = mEquParameter.a;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    u.clear(); u.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<IProblem2H2D::ExtendedSpacePointNode> obsPointNodes;
    for (unsigned int j=0; j<mEquParameter.No; j++) IProblem2H2D::distributeDelta(mOptParameter.xi[j], obsPointNodes, j, dimX, dimY);

    std::vector<IProblem2H2D::ExtendedSpacePointNode> cntDeltaNodes;
    for (unsigned int i=0; i<mEquParameter.Nc; i++) IProblem2H2D::distributeDelta(mOptParameter.eta[i], cntDeltaNodes, i, dimX, dimY);

    SpaceNodePDE sn;

    //--------------------------------------------------------------------------------------------//
    vector<unsigned int> rows0;
    vector<unsigned int> rows1;
    vector<unsigned int> rows2;
    for (unsigned int ny=0; ny<=M; ny++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const IProblem2H2D::ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);
            if (cdn.j == ny)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const IProblem2H2D::ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.j == ny)
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
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const IProblem2H2D::ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);
            if (cdn.i == nx)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const IProblem2H2D::ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.i == nx)
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
        info.resize(mEquParameter.No);
        for (unsigned int j=0; j<mEquParameter.No; j++)
        {
            ExtendedSpaceNode2DH &e = info[j];
            e.setSpaceNode(mOptParameter.xi[j]);
            e.id = j;
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
            u00[m][n] = initial1(sn);
        }
    }
    if (use == true) add2Info(u00, info, 0);
    layerInfo(u00, 0);

    std::vector<IProblem2H2D::ExtendedSpacePointNode> qPointNodes;
    for (unsigned int s=0; s<mEquParameter.Ns; s++) IProblem2H2D::distributeDelta(mEquParameter.theta[s], qPointNodes, s, dimX, dimY);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;

            if (n==0)      sum += a*a*(u00[m][n]-2.0*u00[m][n+1]+u00[m][n+2])/(hx*hx);
            else if (n==N) sum += a*a*(u00[m][n-2]-2.0*u00[m][n-1]+u00[m][n])/(hx*hx);
            else           sum += a*a*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1])/(hx*hx);

            if (m==0)      sum += a*a*(u00[m][n]-2.0*u00[m+1][n]+u00[m+2][n])/(hy*hy);
            else if (m==M) sum += a*a*(u00[m-2][n]-2.0*u00[m-1][n]+u00[m][n])/(hy*hy);
            else           sum += a*a*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n])/(hy*hy);

            sum -= mEquParameter.lambda*initial2(sn);

            for (unsigned int cdn=0; cdn<cntDeltaNodes.size(); cdn++)
            {
                IProblem2H2D::ExtendedSpacePointNode cntNode = cntDeltaNodes.at(cdn);
                if (cntNode.i == n && cntNode.j == m)
                {
                    for (unsigned int opn=0; opn<obsPointNodes.size(); opn++)
                    {
                        IProblem2H2D::ExtendedSpacePointNode obsNode = obsPointNodes.at(opn);
                        if (obsNode.i == n && obsNode.j == m)
                        {
                            sum += mOptParameter.k[obsNode.id][obsNode.id] * ((u00[m][n]*obsNode.w*(hx*hy)) - mOptParameter.z[obsNode.id][obsNode.id]) * cntNode.w;
                        }
                    }
                }
            }

            for (unsigned int si=0; si<qPointNodes.size(); si++)
            {
                IProblem2H2D::ExtendedSpacePointNode qNode = qPointNodes.at(si);
                if (qNode.i == n && qNode.j == m)
                {
                    sum += mEquParameter.q[qNode.id] * qNode.w * (1.0/ht);
                }
            }

            u05[m][n] = u00[m][n] + initial2(sn)*ht*0.5 + sum*ht*ht*0.125;
            u10[m][n] = u00[m][n] + initial2(sn)*ht     + sum*ht*ht*0.500;
        }
    }
    qPointNodes.clear();

    if (use == true) add2Info(u10, info, 1);
    layerInfo(u10, 1);
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

    for (unsigned int l=2; l<=L; l++)
    {
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

                if (m == 0)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                if (m>0 && m<M) d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                if (m == M)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                d1X[n] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                d1X[n] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

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
            for (unsigned int n=0; n<=N; n++) u15[m][n] = x1X[n];
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

                    if (m == 0)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    if (m>0 && m<M) d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    if (m == M)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1X[n] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                    d1X[n] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

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

                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &opn = obsPointNodes[odj];
                                d1X[n] += 0.5*ht*ht * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            }

                            for (unsigned int j=0; j<mEquParameter.No; j++)
                            {
                                d1X[n] -= 0.5*ht*ht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                for (unsigned int n=0; n<=N; n++) u15[m][n] = x1X[n];
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

                if (n==0)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                if (n>0 && n<N) d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                if (n==N)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                d1Y[m] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                d1Y[m] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

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
            for (unsigned int m=0; m<=M; m++) u[m][n] = x1Y[m];
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

                    if (n==0)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    if (n>0 && n<N) d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    if (n==N)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d1Y[m] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                    d1Y[m] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

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

                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &opn = obsPointNodes.at(onj);
                                d1Y[m] += 0.5*ht*ht * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            }

                            for (unsigned int j=0; j<mEquParameter.No; j++)
                            {
                                d1Y[m] -= 0.5*ht*ht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
                for (unsigned int m=0; m<=M; m++) u[m][n] = x1Y[m];
            }
        }
        else
        {
            ////////////////////////////////////
        }
        //--------------------------------------------------------------------------//
        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        if (l==L)
        {
            ut.clear();
            ut.resize(M+1, N+1);
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    ut[m][n] = (u[m][n]-u00[m][n])/(2.0*ht);
                }
            }
        }

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u[m][n];
                u05[m][n] = u15[m][n];
            }
        }

        if (use == true) add2Info(u, info, l);
        layerInfo(u, l);
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

    obsPointNodes.clear();
    cntDeltaNodes.clear();

    u00.clear();
    u05.clear();
    u10.clear();
    u15.clear();
}

void IProblem2HForward2D::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    QPixmap px;
    visualizeMatrixHeat(u, u.min(), u.max(), px);
    char buffer[30] = {0};
    int c = sprintf(buffer, "images/f/%4d.png", ln);
    buffer[c] = 0;
    px.save(buffer, "png", 0);
}

double IProblem2HForward2D::initial1(const SpaceNodePDE &) const
{
    return 0.0;
}

double IProblem2HForward2D::initial2(const SpaceNodePDE &) const
{
    return 0.0;
}

double IProblem2HForward2D::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const
{
    return 0.0;
}

double IProblem2HForward2D::f(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

void IProblem2HForward2D::add2Info(const DoubleMatrix &u, vector<ExtendedSpaceNode2DH> &info, unsigned int ln) const
{
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        ExtendedSpaceNode2DH &pi = info[j];
        for (unsigned int r=0; r<4; r++)
        {
            for (unsigned int c=0; c<4; c++)
            {
                unsigned int x = pi.wi[r][c].i;
                unsigned int y = pi.wi[r][c].j;
                pi.wi[r][c].u[ln] = u[y][x];
            }
        }
    }
}
