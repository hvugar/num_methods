#include "iproblem2hforward2d.h"
#include <imaging.h>

using namespace IProblem2H;

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
    double a = mEquParameter.a;
    unsigned int No = mEquParameter.No;
    unsigned int Nc = mEquParameter.Nc;
    unsigned int Ns = mEquParameter.Ns;

    double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    double p_aa_htht__hxhx_h = +0.5*(a*a*ht*ht)/(hx*hx);
    double m_aa_htht__hxhx_h = -0.5*(a*a*ht*ht)/(hx*hx);
    double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    double p_aa_htht__hyhy_h = +0.5*(a*a*ht*ht)/(hy*hy);
    double m_aa_htht__hyhy_h = -0.5*(a*a*ht*ht)/(hy*hy);
    double p_aa_htht__hxhx___lambda_ht = +1.0 + (a*a*ht*ht)/(hx*hx) + 1.5*(lambda*ht);
    double p_aa_htht__hyhy___lambda_ht = +1.0 + (a*a*ht*ht)/(hy*hy) + 1.5*(lambda*ht);
    double aa__hxhx = ((a*a)/(hx*hx));
    double aa__hyhy = ((a*a)/(hy*hy));
    double htht_h = 0.5*ht*ht;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    u.clear(); u.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<ExtendedSpacePointNode> obsPointNodes;
    for (unsigned int j=0; j<No; j++) IProblem2H2D::distributeDelta(mOptParameter.xi[j], obsPointNodes, j, dimX, dimY);

    std::vector<ExtendedSpacePointNode> cntDeltaNodes;
    for (unsigned int i=0; i<Nc; i++) IProblem2H2D::distributeDelta(mOptParameter.eta[i], cntDeltaNodes, i, dimX, dimY);

    std::vector<ExtendedSpacePointNode> qPointNodes;
    for (unsigned int s=0; s<Ns; s++) IProblem2H2D::distributeDelta(mEquParameter.theta[s], qPointNodes, s, dimX, dimY);

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
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);
            if (cdn.j == ny)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
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
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);
            if (cdn.i == nx)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
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
        info.resize(No);
        for (unsigned int j=0; j<No; j++)
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

//            for (unsigned int si=0; si<qPointNodes.size(); si++)
//            {
//                const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
//                if (qNode.i == n && qNode.j == m)
//                {
//                    u00[m][n] += mEquParameter.q[qNode.id] * qNode.w * (1.0/ht);
//                }
//            }

            if (n==50 && m==50)
            {
                u00[m][n] += mEquParameter.q[0] * (1.0/(hx*hy)) * (1.0/ht);
            }
        }
    }
    if (use == true) add2Info(u00, info, 0);
    layerInfo(u00, 0);
    layerInfo(u00, u00, 0);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;

            if (n==0)      sum += aa__hxhx*(u00[m][n]-2.0*u00[m][n+1]+u00[m][n+2]);
            else if (n==N) sum += aa__hxhx*(u00[m][n-2]-2.0*u00[m][n-1]+u00[m][n]);
            else           sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);

            if (m==0)      sum += aa__hyhy*(u00[m][n]-2.0*u00[m+1][n]+u00[m+2][n]);
            else if (m==M) sum += aa__hyhy*(u00[m-2][n]-2.0*u00[m-1][n]+u00[m][n]);
            else           sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);

            sum -= lambda*initial2(sn);

            for (unsigned int cdn=0; cdn<cntDeltaNodes.size(); cdn++)
            {
                const ExtendedSpacePointNode &cntNode = cntDeltaNodes.at(cdn);
                if (cntNode.i == n && cntNode.j == m)
                {
                    for (unsigned int opn=0; opn<obsPointNodes.size(); opn++)
                    {
                        const ExtendedSpacePointNode &obsNode = obsPointNodes.at(opn);
                        if (obsNode.i == n && obsNode.j == m)
                        {
                            sum += mOptParameter.k[obsNode.id][obsNode.id] * ((u00[m][n]*obsNode.w*(hx*hy)) - mOptParameter.z[obsNode.id][obsNode.id]) * cntNode.w;
                        }
                    }
                }
            }

//            for (unsigned int si=0; si<qPointNodes.size(); si++)
//            {
//                const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
//                if (qNode.i == n && qNode.j == m)
//                {
//                    sum += mEquParameter.q[qNode.id] * qNode.w * (1.0/ht);
//                }
//            }

            u05[m][n] = u00[m][n] + initial2(sn)*ht*0.5 + sum*ht*ht*0.125;
            u10[m][n] = u00[m][n] + initial2(sn)*ht     + sum*ht*ht*0.500;
        }
    }
    qPointNodes.clear();

    if (use == true) add2Info(u10, info, 1);
    layerInfo(u10, u10, 1);
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

    TimeNodePDE tn;

    for (unsigned int l=2; l<=L; l++)
    {
        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht;

        //--------------------------------------------------------------------------//
        //for (unsigned int row=0; row<rows0.size(); row++)
        for (unsigned int m=0; m<=M; m++)
        {
            //unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;

                if (m == 0)     d1X[n] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                if (m>0 && m<M) d1X[n] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                if (m == M)     d1X[n] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                d1X[n] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                d1X[n] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);
                d1X[n] += htht_h * f(sn, tn);

                if (n==0)
                {
                    a1X[0] = 0.0;
                    b1X[0] = p_aa_htht__hxhx___lambda_ht;
                    c1X[0] = m_aa_htht__hxhx;
                }
                else if(n == N)
                {
                    a1X[N] = m_aa_htht__hxhx;
                    b1X[N] = p_aa_htht__hxhx___lambda_ht;
                    c1X[N] = 0.0;
                    d1X[n] -= m_aa_htht__hxhx * 2.0;
                }
                else
                {
                    a1X[n] = m_aa_htht__hxhx_h;
                    b1X[n] = p_aa_htht__hxhx___lambda_ht;
                    c1X[n] = m_aa_htht__hxhx_h;
                }
            }
            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
            for (unsigned int n=0; n<=N; n++) u15[m][n] = x1X[n];
        }

//        if (rows2.size() == 0)
//        {
//            for (unsigned int row=0; row<rows1.size(); row++)
//            {
//                unsigned int m = rows1.at(row);
//                sn.j = m; sn.y = m*hy;
//                for (unsigned int n=0; n<=N; n++)
//                {
//                    sn.i = n; sn.x = n*hx;

//                    if (m == 0)     d1X[n] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
//                    if (m>0 && m<M) d1X[n] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
//                    if (m == M)     d1X[n] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

//                    d1X[n] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
//                    d1X[n] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

//                    if (n==0)
//                    {
//                        a1X[0] = 0.0;
//                        b1X[0] = p_aa_htht__hxhx___lambda_ht;
//                        c1X[0] = m_aa_htht__hxhx;
//                    }
//                    else if(n==N)
//                    {
//                        a1X[N] = m_aa_htht__hxhx;
//                        b1X[N] = p_aa_htht__hxhx___lambda_ht;
//                        c1X[N] = 0.0;
//                    }
//                    else
//                    {
//                        a1X[n] = m_aa_htht__hxhx_h;
//                        b1X[n] = p_aa_htht__hxhx___lambda_ht;
//                        c1X[n] = m_aa_htht__hxhx_h;
//                    }

//                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
//                    {
//                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
//                        if (cdn.i == sn.i && cdn.j == sn.j)
//                        {
//                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
//                            {
//                                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
//                                d1X[n] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
//                            }

//                            for (unsigned int j=0; j<No; j++)
//                            {
//                                d1X[n] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
//                            }
//                        }
//                    }
//                }
//                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
//                for (unsigned int n=0; n<=N; n++) u15[m][n] = x1X[n];
//            }
//        }
//        else
//        {
//            double* a1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
//            double* b1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
//            double* c1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
//            double* d1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
//            double* x1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
//            DoubleMatrix w2(rows1.size()*(N+1), rows1.size()*(N+1), 0.0);

//            unsigned int offset = 0;
//            for (unsigned int row=0; row<rows1.size(); row++)
//            {
//                unsigned int m = rows1.at(row);
//                sn.j = m; sn.y = m*hy;

//                for (unsigned int n=0; n<=N; n++)
//                {
//                    sn.i = n; sn.x = n*hx;

//                    if (m==0)       d1[offset+n] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
//                    if (m>0 && m<M) d1[offset+n] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
//                    if (m==M)       d1[offset+n] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

//                    d1[offset+n] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
//                    d1[offset+n] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

//                    if (n == 0)
//                    {
//                        a1[offset+0] = 0.0;
//                        b1[offset+0] = p_aa_htht__hxhx___lambda_ht;
//                        c1[offset+0] = m_aa_htht__hxhx;
//                    }
//                    else if (n == N)
//                    {
//                        a1[offset+N] = m_aa_htht__hxhx;
//                        b1[offset+N] = p_aa_htht__hxhx___lambda_ht;
//                        c1[offset+N] = 0.0;
//                    }
//                    else
//                    {
//                        a1[offset+n] = m_aa_htht__hxhx_h;
//                        b1[offset+n] = p_aa_htht__hxhx___lambda_ht;
//                        c1[offset+n] = m_aa_htht__hxhx_h;
//                    }

//                    //------------------------------------- Adding delta part -------------------------------------//
//                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
//                    {
//                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
//                        if (cdn.i == sn.i && cdn.j == sn.j)
//                        {
//                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
//                            {
//                                const ExtendedSpacePointNode &opn = obsPointNodes.at(odj);

//                                bool found = false;
//                                for (unsigned int cs=0; cs<rows1.size(); cs++)
//                                {
//                                    if (opn.j == rows1[cs])
//                                    {
//                                        found = true;
//                                        w2[offset+n][cs*(N+1)+opn.i] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
//                                    }
//                                }

//                                if (!found)
//                                {
//                                    d1[offset+n] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
//                                }
//                            }

//                            for (unsigned int j=0; j<No; j++)
//                            {
//                                d1[offset+n] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
//                            }

//                        }
//                    }
//                }
//                offset += N+1;
//            }

//            LinearEquation::func1(a1, b1, c1, d1, w2.data(), x1, rows1.size()*(N+1));

//            offset = 0;
//            for (unsigned int row=0; row<rows1.size(); row++)
//            {
//                unsigned int m=rows1.at(row);
//                for (unsigned int n=0; n<=N; n++)
//                {
//                    u15[m][n] = x1[offset+n];
//                }
//                offset += N+1;
//            }

//            w2.clear();
//            free(x1);
//            free(d1);
//            free(c1);
//            free(b1);
//            free(a1);
//        }
        //--------------------------------------------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//
//        tn.i = l+1; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//
        //for (unsigned int col=0; col<cols0.size(); col++)
        for (unsigned int n=0; n<=N; n++)
        {
            //unsigned int n = cols0.at(col);
            sn.i = n; sn.x = n*hx;
            for (unsigned int m=0; m<=M; m++)
            {
                sn.j = m; sn.y = m*hy;

                if (n==0)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                if (n>0 && n<N) d1Y[m] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                if (n==N)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                d1Y[m] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                d1Y[m] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);
                d1Y[m] += htht_h * f(sn, tn);

                if (m == 0)
                {
                    a1Y[0] = 0.0;
                    b1Y[0] = p_aa_htht__hyhy___lambda_ht;
                    c1Y[0] = m_aa_htht__hyhy;
                }
                else if (m == M)
                {
                    a1Y[M] = m_aa_htht__hyhy;
                    b1Y[M] = p_aa_htht__hyhy___lambda_ht;
                    c1Y[M] = 0.0;
                    d1Y[M] -= m_aa_htht__hyhy * 2.0;
                }
                else
                {
                    a1Y[m] = m_aa_htht__hyhy_h;
                    b1Y[m] = p_aa_htht__hyhy___lambda_ht;
                    c1Y[m] = m_aa_htht__hyhy_h;
                }
            }
            tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
            for (unsigned int m=0; m<=M; m++) u[m][n] = x1Y[m];
        }

//        if (cols2.size() == 0)
//        {
//            for (unsigned int col=0; col<cols1.size(); col++)
//            {
//                unsigned int n = cols1.at(col);
//                sn.i = n; sn.x = n*hx;
//                for (unsigned int m=0; m<=M; m++)
//                {
//                    sn.j = m; sn.y = m*hy;

//                    if (n==0)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
//                    if (n>0 && n<N) d1Y[m] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
//                    if (n==N)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

//                    d1Y[m] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
//                    d1Y[m] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

//                    if (m == 0)
//                    {
//                        a1Y[0] = 0.0;
//                        b1Y[0] = p_aa_htht__hyhy___lambda_ht;
//                        c1Y[0] = m_aa_htht__hyhy;
//                    }
//                    else if (m == M)
//                    {
//                        a1Y[M] = m_aa_htht__hyhy;
//                        b1Y[M] = p_aa_htht__hyhy___lambda_ht;
//                        c1Y[M] = 0.0;
//                    }
//                    else
//                    {
//                        a1Y[m] = m_aa_htht__hyhy_h;
//                        b1Y[m] = p_aa_htht__hyhy___lambda_ht;
//                        c1Y[m] = m_aa_htht__hyhy_h;
//                    }

//                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
//                    {
//                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
//                        if (cdn.i == sn.i && cdn.j == sn.j)
//                        {
//                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
//                            {
//                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);
//                                d1Y[m] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
//                            }

//                            for (unsigned int j=0; j<No; j++)
//                            {
//                                d1Y[m] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
//                            }
//                        }
//                    }
//                }
//                tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
//                for (unsigned int m=0; m<=M; m++) u[m][n] = x1Y[m];
//            }
//        }
//        else
//        {
//            double* a2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
//            double* b2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
//            double* c2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
//            double* d2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
//            double* x2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
//            DoubleMatrix w2(cols1.size()*(M+1), cols1.size()*(M+1), 0.0);

//            unsigned int offset = 0;
//            for (unsigned int col=0; col<cols1.size(); col++)
//            {
//                unsigned int n = cols1.at(col);
//                sn.i = n; sn.x = n*hx;

//                for (unsigned int m=0; m<=M; m++)
//                {
//                    sn.j = m; sn.y = m*hy;

//                    if (n==0)       d2[offset+m] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
//                    if (n>0 && n<N) d2[offset+m] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
//                    if (n==N)       d2[offset+m] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

//                    d2[offset+m] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
//                    d2[offset+m] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

//                    if (m == 0)
//                    {
//                        a2[offset+0] = 0.0;
//                        b2[offset+0] = p_aa_htht__hyhy___lambda_ht;
//                        c2[offset+0] = m_aa_htht__hyhy;
//                        //d2[offset+0] += 0.0;
//                    }
//                    else if (m == M)
//                    {
//                        a2[offset+M] = m_aa_htht__hyhy;
//                        b2[offset+M] = p_aa_htht__hyhy___lambda_ht;
//                        c2[offset+M] = 0.0;
//                        //d2[offset+M] += 0.0;
//                    }
//                    else
//                    {
//                        a2[offset+m] = m_aa_htht__hyhy_h;
//                        b2[offset+m] = p_aa_htht__hyhy___lambda_ht;
//                        c2[offset+m] = m_aa_htht__hyhy_h;
//                    }

//                    //------------------------------------- Adding delta part -------------------------------------//
//                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
//                    {
//                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
//                        if (cdn.i == sn.i && cdn.j == sn.j)
//                        {
//                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
//                            {
//                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);

//                                bool found = false;
//                                for (unsigned int cs=0; cs<cols1.size(); cs++)
//                                {
//                                    if (opn.i == cols1[cs])
//                                    {
//                                        found = true;
//                                        w2[offset+m][cs*(M+1)+opn.j] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
//                                    }
//                                }

//                                if (!found)
//                                {
//                                    d2[offset+m] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
//                                }
//                            }

//                            for (unsigned int j=0; j<No; j++)
//                            {
//                                d2[offset+m] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
//                            }
//                        }
//                    }
//                    //------------------------------------- Adding delta part -------------------------------------//
//                }
//                offset += M+1;
//            }

//            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M+1));

//            offset = 0;
//            for (unsigned int col=0; col<cols1.size(); col++)
//            {
//                unsigned int n=cols1.at(col);
//                for (unsigned int m=0; m<=M; m++)
//                {
//                    u[m][n] = x2[offset+m];
//                }
//                offset += M+1;
//            }

//            w2.clear();
//            free(x2);
//            free(d2);
//            free(c2);
//            free(b2);
//            free(a2);
//        }
        //--------------------------------------------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        //if (l==L)
        //if (l>3)
        {
            ut.clear();
            ut.resize(M+1, N+1);
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    //ut[m][n] = (u[m][n]-u00[m][n])/(2.0*ht);
                    ut[m][n] = (3.0*u[m][n]-4.0*u10[m][n]+u00[m][n])/(2.0*ht);
                }
            }
            layerInfo(u, ut, l);
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

void IProblem2HForward2D::layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const
{
//    if (ln%1000==0)
//    {
//    QPixmap px;
//    visualizeMatrixHeat(u, -5200.0, 1500.0, px);
//    char buffer[30] = {0};
//    int c = sprintf(buffer, "images/s/%6d.png", ln);
//    buffer[c] = 0;
//    px.save(buffer, "png", 0);
//    }

//    if (ln==0)
//    {
//        FILE *file = fopen("images/0.txt", "w");
//        IPrinter::print(u, u.rows(), u.cols(), 14, 2, file);
//        fclose(file);
//    }
}

void IProblem2HForward2D::layerInfo(const DoubleMatrix &u, const DoubleMatrix &ut, unsigned int ln) const
{
    double hx = mspaceDimension[0].step();
    double hy = mspaceDimension[1].step();
    unsigned int N1 = mspaceDimension[0].sizeN();
    unsigned int N2 = mspaceDimension[1].sizeN();

    DoubleMatrix V0(101, 101, 0.0);
    DoubleMatrix V1(101, 101, 0.0);

    double sum0 = 0.0;
    double sum1 = 0.0;

    sum0 += 0.25*(u[0][0]   - V0[0][0])   * (u[0][0]   - V0[0][0]);
    sum0 += 0.25*(u[0][N1]  - V0[0][N1])  * (u[0][N1]  - V0[0][N1]);
    sum0 += 0.25*(u[N2][0]  - V0[N2][0])  * (u[N2][0]  - V0[N2][0]);
    sum0 += 0.25*(u[N2][N1] - V0[N2][N1]) * (u[N2][N1] - V0[N2][N1]);

    sum1 += 0.25*(ut[0][0]   - V1[0][0])   * (ut[0][0]   - V1[0][0]);
    sum1 += 0.25*(ut[0][N1]  - V1[0][N1])  * (ut[0][N1]  - V1[0][N1]);
    sum1 += 0.25*(ut[N2][0]  - V1[N2][0])  * (ut[N2][0]  - V1[N2][0]);
    sum1 += 0.25*(ut[N2][N1] - V1[N2][N1]) * (ut[N2][N1] - V1[N2][N1]);

    for (unsigned int n1=1; n1<=N1-1; n1++)
    {
        sum0 += 0.5*(u[0][n1]  - V0[0][n1]) *(u[0][n1]  - V0[0][n1]);
        sum0 += 0.5*(u[N2][n1] - V0[N2][n1])*(u[N2][n1] - V0[N2][n1]);
        sum1 += 0.5*(ut[0][n1]  - V1[0][n1]) *(ut[0][n1]  - V1[0][n1]);
        sum1 += 0.5*(ut[N2][n1] - V1[N2][n1])*(ut[N2][n1] - V1[N2][n1]);
    }

    for (unsigned int n2=1; n2<=N2-1; n2++)
    {
        sum0 += 0.5*(u[n2][0]  - V0[n2][0]) *(u[n2][0]  - V0[n2][0]);
        sum0 += 0.5*(u[n2][N1] - V0[n2][N1])*(u[n2][N1] - V0[n2][N1]);
        sum1 += 0.5*(ut[n2][0]  - V1[n2][0]) *(ut[n2][0]  - V1[n2][0]);
        sum1 += 0.5*(ut[n2][N1] - V1[n2][N1])*(ut[n2][N1] - V1[n2][N1]);
    }

    for (unsigned int n2 = 1; n2 <= N2-1; n2++)
    {
        for (unsigned int n1 = 1; n1 <= N1-1; n1++)
        {
            sum0 += (u[n2][n1] - V0[n2][n1])*(u[n2][n1] - V0[n2][n1]);
            sum1 += (ut[n2][n1] - V1[n2][n1])*(ut[n2][n1] - V1[n2][n1]);
        }
    }

    sum0 *= (hx*hy);
    sum1 *= (hx*hy);

    V0.clear();
    V1.clear();

    printf("%6d %20.8f %20.8f u:%20.8f %20.8f ut:%20.8f %20.8f\n", ln, sum0, sum1, u.min(), u.max(), ut.min(), ut.max());
}

double IProblem2HForward2D::initial1(const SpaceNodePDE & sn) const
{
    double hx = mspaceDimension[0].step();
    double hy = mspaceDimension[1].step();
    double x = hx*sn.i;
    double y = hy*sn.j;
    return x*x + y*y;
//    return 0.0;
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
    return -2.0;
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
