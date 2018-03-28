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

    double lambda = mParameter.lambda;
    //double lambda1 = mParameter.lambda1;
    double a = mParameter.a;

    DoubleMatrix u0(M+1, N+1);
    DoubleMatrix u1(M+1, N+1);
    DoubleMatrix uh(M+1, N+1);
    u.clear();
    u.resize(M+1, N+1);
    ut.clear();
    ut.resize(M+1, N+1);

    //double min = +1000.0;
    //double max = -1000.0;

    //--------------------------------------------------------------------------------------------//
    std::vector<IProblem2H2D::ExtendedSpacePointNode> obsPointNodes;
    for (unsigned int j=0; j<mParameter.No; j++) IProblem2H2D::distributeDelta(mParameter.xi[j], obsPointNodes, j, dimX, dimY);

    std::vector<IProblem2H2D::ExtendedSpacePointNode> cntDeltaNodes;
    for (unsigned int i=0; i<mParameter.Nc; i++) IProblem2H2D::distributeDelta(mParameter.eta[i], cntDeltaNodes, i, dimX, dimY);

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
        info.resize(mParameter.No);
        for (unsigned int j=0; j<mParameter.No; j++)
        {
            ExtendedSpaceNode2DH &e = info[j];
            e.setSpaceNode(mParameter.xi[j]);
            e.id = j;
            e.extendWeights(dimX, dimY);
            e.extendLayers(L/2+1);
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
            u0[m][n] = initial1(sn);
        }
    }
    if (use == true) add2Info(u0, info, 0);
    layerInfo(u0, 0);

    std::vector<IProblem2H2D::ExtendedSpacePointNode> qPointNodes;
    for (unsigned int s=0; s<mParameter.Ns; s++) IProblem2H2D::distributeDelta(mParameter.theta[s], qPointNodes, s, dimX, dimY);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u1[m][n] = u0[m][n] + 2.0*initial2(sn)*ht;

            double diff = 0.0;
            if (n==0 )     diff += a*a*(u0[m][n]-2.0*u0[m][n+1]+u0[m][n+2])/(hx*hx);
            else if (n==N) diff += a*a*(u0[m][n-2]-2.0*u0[m][n-1]+u0[m][n])/(hx*hx);
            else           diff += a*a*(u0[m][n-1]-2.0*u0[m][n]+u0[m][n+1])/(hx*hx);

            if (m==0 )     diff += a*a*(u0[m][n]-2.0*u0[m+1][n]+u0[m+2][n])/(hy*hy);
            else if (m==M) diff += a*a*(u0[m-2][n]-2.0*u0[m-1][n]+u0[m][n])/(hy*hy);
            else           diff += a*a*(u0[m-1][n]-2.0*u0[m][n]+u0[m+1][n])/(hy*hy);
            u1[m][n] += 4.0*ht*ht*0.5*diff;

            for (unsigned int si=0; si<qPointNodes.size(); si++)
            {
                IProblem2H2D::ExtendedSpacePointNode qNode = qPointNodes.at(si);
                if (qNode.i == n && qNode.j == m)
                {
                    u1[m][n] += 2.0*ht*0.5 * mParameter.q[qNode.id] * qNode.w;
                }
            }
        }
    }
    qPointNodes.clear();

    if (use == true) add2Info(u1, info, 1);
    layerInfo(u1, 1);
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

    for (unsigned int l=4; l<=L; l+=2)
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

                d1X[n] = (2.0+2.0*lambda*ht)*u1[m][n] - (1.0+(lambda*ht)/2.0)*u0[m][n];

                if (m == 0)     d1X[n] += ((a*a*ht*ht)/(hy*hy))*(u1[0][n]   - 2.0*u1[1][n]   + u1[2][n]);
                if (m>0 && m<M) d1X[n] += ((a*a*ht*ht)/(hy*hy))*(u1[m-1][n] - 2.0*u1[m][n]   + u1[m+1][n]);
                if (m == M)     d1X[n] += ((a*a*ht*ht)/(hy*hy))*(u1[M-2][n] - 2.0*u1[M-1][n] + u1[M][n]);

                if (n==0)
                {
                    a1X[0] = 0.0;
                    b1X[0] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                    c1X[0] = -2.0*(a*a*ht*ht)/(hx*hx);
                    //d1X[0] += 0.0;
                }
                else if(n == N)
                {
                    a1X[N] = -2.0*(a*a*ht*ht)/(hx*hx);
                    b1X[N] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                    c1X[N] = 0.0;
                    //d1X[N] += 0.0;
                }
                else
                {
                    a1X[n] = -(a*a*ht*ht)/(hx*hx);
                    b1X[n] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                    c1X[n] = -(a*a*ht*ht)/(hx*hx);
                }
            }
            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
            for (unsigned int n=0; n<=N; n++) uh[m][n] = x1X[n];
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

                    d1X[n] = (2.0+2.0*lambda*ht)*u1[m][n] - (1.0+(lambda*ht)/2.0)*u0[m][n];

                    if (m == 0)     d1X[n] += ((a*a*ht*ht)/(hy*hy))*(u1[0][n]  -2.0*u1[1][n]  +u1[2][n]);
                    if (m>0 && m<M) d1X[n] += ((a*a*ht*ht)/(hy*hy))*(u1[m-1][n]-2.0*u1[m][n]  +u1[m+1][n]);
                    if (m == M)     d1X[n] += ((a*a*ht*ht)/(hy*hy))*(u1[M-2][n]-2.0*u1[M-1][n]+u1[M][n]);

                    if (n==0)
                    {
                        a1X[0] = 0.0;
                        b1X[0] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                        c1X[0] = -2.0*(a*a*ht*ht)/(hx*hx);
                        d1X[0] += 0.0;
                    }
                    else if(n==N)
                    {
                        a1X[N] = -2.0*(a*a*ht*ht)/(hx*hx);
                        b1X[N] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                        c1X[N] = 0.0;
                        d1X[N] += 0.0;
                    }
                    else
                    {
                        a1X[n] = -(a*a*ht*ht)/(hx*hx);
                        b1X[n] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                        c1X[n] = -(a*a*ht*ht)/(hx*hx);
                    }

                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int s=0; s<obsPointNodes.size(); s++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &opn = obsPointNodes[s];
                                d1X[n] += ht*ht * mParameter.k[cdn.id][opn.id] * uh[opn.j][opn.i] * cdn.w * (opn.w * (hx*hy));
                            }

                            for (unsigned int j=0; j<mParameter.No; j++)
                            {
                                d1X[n] -= ht*ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                for (unsigned int n=0; n<=N; n++) uh[m][n] = x1X[n];
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

                d1Y[m] = (2.0+2.0*lambda*ht)*uh[m][n] - (1.0+(lambda*ht)/2.0)*u1[m][n];

                if (n==0)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(uh[m][0]   - 2.0*uh[m][1]   + uh[m][2]);
                if (n>0 && n<N) d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(uh[m][n-1] - 2.0*uh[m][n]   + uh[m][n+1]);
                if (n==N)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(uh[m][N-2] - 2.0*uh[m][N-1] + uh[m][N]);

                if (m == 0)
                {
                    a1Y[0] = 0.0;
                    b1Y[0] = +1.0 + 2.0*((a*a*ht*ht)/(hy*hy)) + 3.0*(lambda*ht)/2.0;
                    c1Y[0] = -2.0*(a*a*ht*ht)/(hy*hy);
                    d1Y[0] += 0.0;
                }
                else if (m == M)
                {
                    a1Y[M] = -2.0*(a*a*ht*ht)/(hy*hy);
                    b1Y[M] = +1.0 + 2.0*((a*a*ht*ht)/(hy*hy)) + 3.0*(lambda*ht)/2.0;
                    c1Y[M] = 0.0;
                    d1Y[M] += 0.0;
                }
                else
                {
                    a1Y[m] = -(a*a*ht*ht)/(hy*hy);
                    b1Y[m] = +1.0 + 2.0*(a*a*ht*ht)/(hy*hy) + 3.0*(lambda*ht)/2.0;
                    c1Y[m] = -(a*a*ht*ht)/(hy*hy);
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

                    d1Y[m] = (2.0+2.0*lambda*ht)*uh[m][n] - (1.0+(lambda*ht)/2.0)*u1[m][n];

                    if (n==0)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(uh[m][0]   - 2.0*uh[m][1]   + uh[m][2]);
                    if (n>0 && n<N) d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(uh[m][n-1] - 2.0*uh[m][n]   + uh[m][n+1]);
                    if (n==N)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(uh[m][N-2] - 2.0*uh[m][N-1] + uh[m][N]);

                    if (m == 0)
                    {
                        a1Y[0] = 0.0;
                        b1Y[0] = +1.0 + 2.0*((a*a*ht*ht)/(hy*hy)) + 3.0*(lambda*ht)/2.0;
                        c1Y[0] = -2.0*(a*a*ht*ht)/(hy*hy);
                        d1Y[0] += 0.0;
                    }
                    else if (m == M)
                    {
                        a1Y[M] = -2.0*(a*a*ht*ht)/(hy*hy);
                        b1Y[M] = +1.0 + 2.0*((a*a*ht*ht)/(hy*hy)) + 3.0*(lambda*ht)/2.0;
                        c1Y[M] = 0.0;
                        d1Y[M] += 0.0;
                    }
                    else
                    {
                        a1Y[m] = -(a*a*ht*ht)/(hy*hy);
                        b1Y[m] = +1.0 + 2.0*(a*a*ht*ht)/(hy*hy) + 3.0*(lambda*ht)/2.0;
                        c1Y[m] = -(a*a*ht*ht)/(hy*hy);
                    }

                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int s=0; s<obsPointNodes.size(); s++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &opn = obsPointNodes.at(s);
                                d1Y[m] += ht*ht * mParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * cdn.w * (opn.w * (hx*hy));
                            }

                            for (unsigned int j=0; j<mParameter.No; j++)
                            {
                                d1Y[m] -= ht*ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
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
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    ut[m][n] = (u[m][n]-u0[m][n])/(2.0*ht);
                }
            }
        }

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u0[m][n] = uh[m][n];
                u1[m][n] = u[m][n];
            }
        }

        if (use == true) add2Info(u, info, l/2);
        layerInfo(u, l/2);
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

    //u3.clear();
    uh.clear();
    u1.clear();
    u0.clear();
}

void IProblem2HForward2D::calculateMVD1(DoubleMatrix &u, DoubleMatrix &ut) const
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
    double hh = 0.5*ht;

    double lambda = mParameter.lambda;
    //double lambda1 = mParameter.lambda1;
    double a = mParameter.a;

    DoubleMatrix u0(M+1, N+1);
    DoubleMatrix u1(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);

    u.clear();
    u.resize(M+1, N+1);

//    double min = +1000.0;
//    double max = -1000.0;

    //--------------------------------------------------------------------------------------------//
    std::vector<IProblem2H2D::ExtendedSpacePointNode> obsPointNodes;
    for (unsigned int j=0; j<mParameter.No; j++) IProblem2H2D::distributeDelta(mParameter.xi[j], obsPointNodes, j, dimX, dimY);

    std::vector<IProblem2H2D::ExtendedSpacePointNode> cntdeltaNodes;
    for (unsigned int i=0; i<mParameter.Nc; i++) IProblem2H2D::distributeDelta(mParameter.eta[i], cntdeltaNodes, i, dimX, dimY);

    SpaceNodePDE sn;

    //--------------------------------------------------------------------------------------------//
    vector<unsigned int> rows0;
    vector<unsigned int> rows1;
    vector<unsigned int> rows2;
    for (unsigned int ny=0; ny<=M; ny++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntdeltaNodes.size(); i++)
        {
            const IProblem2H2D::ExtendedSpacePointNode &cdn = cntdeltaNodes.at(i);
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
        for (unsigned int i=0; i<cntdeltaNodes.size(); i++)
        {
            const IProblem2H2D::ExtendedSpacePointNode &cdn = cntdeltaNodes.at(i);
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
        info.resize(mParameter.No);
        for (unsigned int j=0; j<mParameter.No; j++)
        {
            ExtendedSpaceNode2DH &e = info[j];
            e.setSpaceNode(mParameter.xi[j]);
            e.id = j;
            e.extendWeights(dimX, dimY);
            e.extendLayers(L/2+1);
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
            u0[m][n] = initial1(sn);
        }
    }
    //if (use == true) add2Info(u0, info, 0);
    layerInfo(u0, 0);

    std::vector<IProblem2H2D::ExtendedSpacePointNode> qPointNodes;
    for (unsigned int s=0; s<mParameter.Ns; s++) IProblem2H2D::distributeDelta(mParameter.theta[s], qPointNodes, s, dimX, dimY);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u1[m][n]  = u0[m][n] + initial2(sn)*ht;
            u15[m][n] = u0[m][n] + initial2(sn)*hh;

            double diff = 0.0;
            if (n==0 )     diff += a*a*(u0[m][n]-2.0*u0[m][n+1]+u0[m][n+2])/(hx*hx);
            else if (n==N) diff += a*a*(u0[m][n-2]-2.0*u0[m][n-1]+u0[m][n])/(hx*hx);
            else           diff += a*a*(u0[m][n-1]-2.0*u0[m][n]+u0[m][n+1])/(hx*hx);

            if (m==0 )     diff += a*a*(u0[m][n]-2.0*u0[m+1][n]+u0[m+2][n])/(hy*hy);
            else if (m==M) diff += a*a*(u0[m-2][n]-2.0*u0[m-1][n]+u0[m][n])/(hy*hy);
            else           diff += a*a*(u0[m-1][n]-2.0*u0[m][n]+u0[m+1][n])/(hy*hy);

            u1[m][n]  += 0.5*ht*ht*diff;
            u15[m][n] += 0.5*hh*hh*diff;

            u1[m][n] -= mParameter.lambda*initial2(sn);

            for (unsigned int si=0; si<qPointNodes.size(); si++)
            {
                IProblem2H2D::ExtendedSpacePointNode qNode = qPointNodes.at(si);
                if (qNode.i == n && qNode.j == m)
                {
                    u1[m][n] += ht*0.5 * mParameter.q[qNode.id] * qNode.w;
                }
            }
        }
    }
    qPointNodes.clear();

    //if (use == true) add2Info(u1, info, 1);
    layerInfo(u1, 1);
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

                if (m == 0)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[0][n]   - 2.0*u1[1][n]   + u1[2][n]);
                if (m>0 && m<M) d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[m-1][n] - 2.0*u1[m][n]   + u1[m+1][n]);
                if (m == M)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[M-2][n] - 2.0*u1[M-1][n] + u1[M][n]);

                d1X[n] += 0.5*(u1[m][n]-u0[m][n]) + u1[m][n];
                d1X[n] += 0.5*lambda*ht*(4.0*u1[m][n]-u15[m][n]);

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
            for (unsigned int n=0; n<=N; n++) u05[m][n] = x1X[n];
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

                    if (m == 0)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[0][n]   - 2.0*u1[1][n]   + u1[2][n]);
                    if (m>0 && m<M) d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[m-1][n] - 2.0*u1[m][n]   + u1[m+1][n]);
                    if (m == M)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[M-2][n] - 2.0*u1[M-1][n] + u1[M][n]);

                    d1X[n] += 0.5*(u1[m][n]-u0[m][n]) + u1[m][n];
                    d1X[n] += 0.5*lambda*ht*(4.0*u1[m][n]-u15[m][n]);

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

                    for (unsigned int cni=0; cni<cntdeltaNodes.size(); cni++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &cdn = cntdeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int s=0; s<obsPointNodes.size(); s++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &opn = obsPointNodes[s];
                                d1X[n] += 0.5*ht*ht * mParameter.k[cdn.id][opn.id] * u05[opn.j][opn.i] * cdn.w * (opn.w * (hx*hy));
                            }

                            for (unsigned int j=0; j<mParameter.No; j++)
                            {
                                d1X[n] -= 0.5*ht*ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                for (unsigned int n=0; n<=N; n++) u05[m][n] = x1X[n];
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

                if (n==0)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u05[m][0]   - 2.0*u05[m][1]   + u05[m][2]);
                if (n>0 && n<N) d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u05[m][n-1] - 2.0*u05[m][n]   + u05[m][n+1]);
                if (n==N)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u05[m][N-2] - 2.0*u05[m][N-1] + u05[m][N]);

                d1Y[m] += 0.5*(u1[m][n]-u0[m][n]) + u05[m][n];
                d1Y[m] += 0.5*lambda*ht*(4.0*u05[m][n]-u1[m][n]);

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

                    if (n==0)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u05[m][0]   - 2.0*u05[m][1]   + u05[m][2]);
                    if (n>0 && n<N) d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u05[m][n-1] - 2.0*u05[m][n]   + u05[m][n+1]);
                    if (n==N)       d1Y[m] = ((a*a*ht*ht)/(hx*hx))*(u05[m][N-2] - 2.0*u05[m][N-1] + u05[m][N]);

                    d1Y[m] += 0.5*(u1[m][n]-u0[m][n]) + u05[m][n];
                    d1Y[m] += 0.5*lambda*ht*(4.0*u05[m][n]-u1[m][n]);

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

                    for (unsigned int cni=0; cni<cntdeltaNodes.size(); cni++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &cdn = cntdeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int s=0; s<obsPointNodes.size(); s++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &opn = obsPointNodes.at(s);
                                d1Y[m] += 0.5*ht*ht * mParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * cdn.w * (opn.w * (hx*hy));
                            }

                            for (unsigned int j=0; j<mParameter.No; j++)
                            {
                                d1Y[m] -= 0.5*ht*ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
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

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u0[m][n] = u1[m][n];
                u1[m][n] = u[m][n];
                u15[m][n] = u05[m][n];
            }
        }

        //if (use == true) add2Info(u1, info, l);
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
    cntdeltaNodes.clear();

    u1.clear();
    u0.clear();
}

void IProblem2HForward2D::layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const
{
    C_UNUSED(u);
    C_UNUSED(layerNumber);

    //    IPrinter::printSeperatorLine();
    //    printf("%d %f\n", l, l*ht);
    //    IPrinter::printMatrix(12, 6, u);

    {
        QPixmap px;
        visualizeMatrixHeat(u, u.min(), u.max(), px);
        char buffer[12] = {0};
        int c = sprintf(buffer, "d:/images1/%6d.png", layerNumber);
        buffer[c] = 0;
        px.save(buffer, "png", 0);
    }
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
    for (unsigned int j=0; j<mParameter.No; j++)
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
