#include "iproblem2hforward2d.h"
#include <imaging.h>
#include <ctime>
#include <chrono>
#include <iostream>
#include <iomanip>

using namespace IProblem2H;

void IProblem2HForward2D::calculateMVD(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const
{
    //std::clock_t c_start = std::clock();
    calculateMVD_D3(u, ut, info, use);
    //std::clock_t c_end = std::clock();
    //std::cout << "CPU time used: " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";
}

void IProblem2HForward2D::calculateMVD_D(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const
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

    double p_aa_htht__hyhy_h = +0.5*(a*a*ht*ht)/(hy*hy);
    double m_aa_htht__hxhx_h = -0.5*(a*a*ht*ht)/(hx*hx);
    double p_aa_htht__hxhx___lambda_ht = +1.0 + (a*a*ht*ht)/(hx*hx) + 1.5*(lambda*ht);

    double p_aa_htht__hxhx_h = +0.5*(a*a*ht*ht)/(hx*hx);
    double m_aa_htht__hyhy_h = -0.5*(a*a*ht*ht)/(hy*hy);
    double p_aa_htht__hyhy___lambda_ht = +1.0 + (a*a*ht*ht)/(hy*hy) + 1.5*(lambda*ht);

    //double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    //double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    double aa__hxhx = (a*a)/(hx*hx);
    double aa__hyhy = (a*a)/(hy*hy);
    double htht_h = 0.5*ht*ht;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    u.clear(); u.resize(M+1, N+1);

    //----------------------------------------------------------------------------------------------//
    std::vector<ExtendedSpacePointNode> obsPointNodes;
    for (unsigned int j=0; j<No; j++) IProblem2H2D::distributeDelta(mOptParameter.xi[j], obsPointNodes, j, dimX, dimY, 4);

    std::vector<ExtendedSpacePointNode> cntDeltaNodes;
    for (unsigned int i=0; i<Nc; i++) IProblem2H2D::distributeDelta(mOptParameter.eta[i], cntDeltaNodes, i, dimX, dimY, 4);

    std::vector<ExtendedSpacePointNode> qPointNodes;
    for (unsigned int s=0; s<Ns; s++) IProblem2H2D::distributeDelta(mEquParameter.theta[s], qPointNodes, s, dimX, dimY, 4);

    SpaceNodePDE sn;

    //----------------------------------------------------------------------------------------------//
    std::vector<unsigned int> rows0, rows1, rows2;
    for (unsigned int m=1; m<=M-1; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);

            if (cdn.j == m)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.j == m)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
    }
    //printf("%d %d %d\n", rows0.size(), rows1.size(), rows2.size());

    //----------------------------------------------------------------------------------------------//
    vector<unsigned int> cols0, cols1, cols2;
    for (unsigned int n=1; n<=N-1; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);
            if (cdn.i == n)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.i == n)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
    }
    //printf("%d %d %d\n", cols0.size(), cols1.size(), cols2.size());

    //----------------------------------------------------------------------------------------------//
    //-------------------------------------------- info --------------------------------------------//
    if (use == true)
    {
        info.resize(No);
        for (unsigned int j=0; j<No; j++)
        {
            ExtendedSpaceNode2DH &e = info[j];
            e.id = j;
            e.setSpaceNode(mOptParameter.xi[j]);
            e.extendWeights(dimX, dimY, L+1, _INFO_ROWS_, _INFO_COLS_);
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
    if (use == true) add2Info(u00, info, 0, _INFO_ROWS_, _INFO_COLS_);
    layerInfo(u00, 0);
    layerInfo(u00, u00, 0);
    //------------------------------------- initial conditions -------------------------------------//
    double hh = 0.5*ht;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum1 = 0.0;
            double sum2 = 0.0;

            if (n==0)      sum1 += aa__hxhx*(u00[m][n]-2.0*u00[m][n+1]+u00[m][n+2]);
            else if (n==N) sum1 += aa__hxhx*(u00[m][n-2]-2.0*u00[m][n-1]+u00[m][n]);
            else           sum1 += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);

            if (m==0)      sum1 += aa__hyhy*(u00[m][n]-2.0*u00[m+1][n]+u00[m+2][n]);
            else if (m==M) sum1 += aa__hyhy*(u00[m-2][n]-2.0*u00[m-1][n]+u00[m][n]);
            else           sum1 += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);

            sum1 -= lambda*initial2(sn);

            for (unsigned int cdn=0; cdn<cntDeltaNodes.size(); cdn++)
            {
                const ExtendedSpacePointNode &cntNode = cntDeltaNodes.at(cdn);
                if (cntNode.i == n && cntNode.j == m)
                {
                    double kz = 0.0;
                    for (unsigned int j=0; j<No; j++)
                    {
                        kz += mOptParameter.k[cntNode.id][j]*mOptParameter.z[cntNode.id][j];
                    }
                    //sum1 -= kz * cntNode.w;
                }
            }

            sum2 = sum1;
            //            for (unsigned int si=0; si<qPointNodes.size(); si++)
            //            {
            //                const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
            //                if (qNode.i == n && qNode.j == m)
            //                {
            //                    sum1 += mEquParameter.q[qNode.id] * qNode.w * (1.0/ht);
            //                }
            //            }

            u05[m][n] = u00[m][n] + hh*initial2(sn) + 0.5*hh*hh*sum2;
            u10[m][n] = u00[m][n] + ht*initial2(sn) + 0.5*ht*ht*sum1;
        }
    }

    if (use == true) add2Info(u10, info, 1, _INFO_ROWS_, _INFO_COLS_);
    layerInfo(u10, u10, 1);
    layerInfo(u10, 1);
    //------------------------------------- initial conditions -------------------------------------//

    unsigned int sizeN = sizeof(double)*(N-1);
    unsigned int sizeM = sizeof(double)*(M-1);

    double *ax = (double *) malloc(sizeN); for (unsigned int n=1; n<=N-1; n++) ax[n-1] = m_aa_htht__hxhx_h;           ax[0] = 0.0;
    double *bx = (double *) malloc(sizeN); for (unsigned int n=1; n<=N-1; n++) bx[n-1] = p_aa_htht__hxhx___lambda_ht;
    double *cx = (double *) malloc(sizeN); for (unsigned int n=1; n<=N-1; n++) cx[n-1] = m_aa_htht__hxhx_h;           cx[N-2] = 0.0;
    double *dx = (double *) malloc(sizeN);
    double *rx = (double *) malloc(sizeN);

    double *ay = (double *) malloc(sizeM); for (unsigned int m=1; m<=M-1; m++) ay[m-1] = m_aa_htht__hyhy_h;           ay[0] = 0.0;
    double *by = (double *) malloc(sizeM); for (unsigned int m=1; m<=M-1; m++) by[m-1] = p_aa_htht__hyhy___lambda_ht;
    double *cy = (double *) malloc(sizeM); for (unsigned int m=1; m<=M-1; m++) cy[m-1] = m_aa_htht__hxhx_h;           cy[M-2] = 0.0;
    double *dy = (double *) malloc(sizeM);
    double *ry = (double *) malloc(sizeM);

    TimeNodePDE tn;

    for (unsigned int l=2; l<=L; l++)
    {
        //------------------------------------- border -------------------------------------//
        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy;
            sn1.j = m; sn1.y = m*hy;

            u15[m][0] = boundary(sn0, tn);
            u15[m][N] = boundary(sn1, tn);

            u[m][0] = boundary(sn0, tn);
            u[m][N] = boundary(sn1, tn);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx;
            sn1.i = n; sn1.x = n*hx;

            u15[0][n] = boundary(sn0, tn);
            u15[M][n] = boundary(sn1, tn);

            u[0][n] = boundary(sn0, tn);
            u[M][n] = boundary(sn1, tn);
        }
        //------------------------------------- border -------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht;

        //--------------------------------------------------------------------------//

        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = n; sn.x = n*hx;

                if (m == 0)     dx[n-1] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                if (m>0 && m<M) dx[n-1] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                if (m == M)     dx[n-1] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                dx[n-1] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                dx[n-1] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

                //ax[n-1] = m_aa_htht__hxhx_h;
                //bx[n-1] = p_aa_htht__hxhx___lambda_ht;
                //cx[n-1] = m_aa_htht__hxhx_h;
            }

            //ax[0]   = 0.0;
            //cx[N-2] = 0.0;

            dx[0]   -= m_aa_htht__hxhx_h * u15[m][0];
            dx[N-2] -= m_aa_htht__hxhx_h * u15[m][N];

            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }

        if (rows2.size() == 0)
        {
            double U15[No]; for (unsigned int j=0; j<No; j++) U15[j] = 0.0;
            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                U15[opn.id] += u15[opn.j][opn.i] * (opn.w * (hx*hy));
            }

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (m == 0)     dx[n-1] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    if (m>0 && m<M) dx[n-1] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    if (m == M)     dx[n-1] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    dx[n-1] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                    dx[n-1] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

                    //ax[n-1] = m_aa_htht__hxhx_h;
                    //bx[n-1] = p_aa_htht__hxhx___lambda_ht;
                    //cx[n-1] = m_aa_htht__hxhx_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int j=0; j<No; j++)
                            {
                                dx[n-1] += htht_h * mOptParameter.k[cdn.id][j] * (U15[j]-mOptParameter.z[cdn.id][j]) * cdn.w;
                            }
                            //                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            //                            {
                            //                                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                            //                                d1X[n-1] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            //                            }
                            //                            for (unsigned int j=0; j<No; j++)
                            //                            {
                            //                                d1X[n-1] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            //                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                //ax[0]   = 0.0;
                //cx[N-2] = 0.0;

                dx[0]   -= m_aa_htht__hxhx_h * u15[m][0];
                dx[N-2] -= m_aa_htht__hxhx_h * u15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }
        }

        if (rows2.size() != 0)
        {
            //throw std::exception();

            unsigned int row1_size = rows1.size()*(N-1);
            double* a1 = (double*) malloc(sizeof(double)*row1_size);
            double* b1 = (double*) malloc(sizeof(double)*row1_size);
            double* c1 = (double*) malloc(sizeof(double)*row1_size);
            double* d1 = (double*) malloc(sizeof(double)*row1_size);
            double* x1 = (double*) malloc(sizeof(double)*row1_size);
            DoubleMatrix w2(row1_size, row1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (m==0)       d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    if (m>0 && m<M) d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    if (m==M)       d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1[offset+(n-1)] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                    d1[offset+(n-1)] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

                    a1[offset+(n-1)] = m_aa_htht__hxhx_h;
                    b1[offset+(n-1)] = p_aa_htht__hxhx___lambda_ht;
                    c1[offset+(n-1)] = m_aa_htht__hxhx_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(odj);

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (opn.j == rows1[rs])
                                    {
                                        found = true;
                                        w2[offset+(n-1)][rs*(N-1)+(opn.i-1)] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                        break;
                                    }
                                }

                                if (found == false)
                                {
                                    d1[offset+(n-1)] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }
                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[offset+(n-1)] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx_h * u15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx_h * u15[m][N];

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w2.data(), x1, row1_size);

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=1; n<=N-1; n++)
                {
                    u15[m][n] = x1[offset+(n-1)];
                }
                offset += N-1;
            }

            w2.clear();
            free(x1);
            free(d1);
            free(c1);
            free(b1);
            free(a1);
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
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = m; sn.y = m*hy;

                if (n==0)       dy[m-1] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                if (n==N)       dy[m-1] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                dy[m-1] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                dy[m-1] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

                //ay[m-1] = m_aa_htht__hyhy_h;
                //by[m-1] = p_aa_htht__hyhy___lambda_ht;
                //cy[m-1] = m_aa_htht__hyhy_h;
            }

            //ay[0]   = 0.0;
            //cy[M-2] = 0.0;

            dy[0]   -= m_aa_htht__hyhy_h * u[0][n];
            dy[M-2] -= m_aa_htht__hyhy_h * u[M][n];

            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u[m][n] = ry[m-1];
        }

        if (cols2.size() == 0)
        {
            double U[No]; for (unsigned int j=0; j<No; j++) U[j] = 0.0;
            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                U[opn.id] += u[opn.j][opn.i] * (opn.w * (hx*hy));
            }

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    if (n==0)       dy[m-1] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    if (n==N)       dy[m-1] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    dy[m-1] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                    dy[m-1] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

                    //ay[m-1] = m_aa_htht__hyhy_h;
                    //by[m-1] = p_aa_htht__hyhy___lambda_ht;
                    //cy[m-1] = m_aa_htht__hyhy_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int j=0; j<No; j++)
                            {
                                dy[m-1] += htht_h * mOptParameter.k[cdn.id][j] * (U[j]-mOptParameter.z[cdn.id][j]) * cdn.w;
                            }
                            //                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            //                            {
                            //                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);
                            //                                d1Y[m-1] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            //                            }
                            //                            for (unsigned int j=0; j<No; j++)
                            //                            {
                            //                                d1Y[m-1] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            //                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                //ay[0]   = 0.0;
                //cy[M-2] = 0.0;

                dy[0]   -= m_aa_htht__hyhy_h * u[0][n];
                dy[M-2] -= m_aa_htht__hyhy_h * u[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u[m][n] = ry[m-1];
            }
        }

        if (cols2.size() != 0)
        {
            //throw std::exception();

            double* a2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* b2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* c2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* d2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* x2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            DoubleMatrix w2(cols1.size()*(M-1), cols1.size()*(M-1), 0.0);

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    if (n==0)       d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    if (n>0 && n<N) d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    if (n==N)       d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d2[offset+(m-1)] += 0.5*(u10[m][n] - u00[m][n]) + u15[m][n];
                    d2[offset+(m-1)] += 0.5*lambda*ht*(4.0*u15[m][n] - u10[m][n]);

                    a2[offset+(m-1)] = m_aa_htht__hyhy_h;
                    b2[offset+(m-1)] = p_aa_htht__hyhy___lambda_ht;
                    c2[offset+(m-1)] = m_aa_htht__hyhy_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (opn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+(m-1)][cs*(M-1)+(opn.j-0)] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+(m-1)] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }
                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[offset+(m-1)] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy_h * u[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy_h * u[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M-1));

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    u[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }

            w2.clear();
            free(x2);
            free(d2);
            free(c2);
            free(b2);
            free(a2);
        }
        //--------------------------------------------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        if (l==2)
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    for (unsigned int si=0; si<qPointNodes.size(); si++)
                    {
                        const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                        if (qNode.i == n && qNode.j == m)
                        {
                            u[m][n] += (mEquParameter.q[qNode.id] * qNode.w * (1.0/ht))*ht*ht;
                        }
                    }
                }
            }
        }

        if (l==L)
        {
            ut.clear();
            ut.resize(M+1, N+1);
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
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

        if (use == true) add2Info(u, info, l, _INFO_ROWS_, _INFO_COLS_);
        layerInfo(u, l);
    }

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    rows0.clear();
    rows1.clear();
    rows2.clear();

    cols0.clear();
    cols1.clear();
    cols2.clear();

    qPointNodes.clear();
    obsPointNodes.clear();
    cntDeltaNodes.clear();

    u00.clear();
    u05.clear();
    u10.clear();
    u15.clear();
}

void IProblem2HForward2D::calculateMVD_D1(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const
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
    //unsigned int Nc = mEquParameter.Nc;
    //unsigned int Ns = mEquParameter.Ns;

    //double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    double p_aa_htht__hxhx_h = +0.5*(a*a*ht*ht)/(hx*hx);
    double m_aa_htht__hxhx_h = -0.5*(a*a*ht*ht)/(hx*hx);
    //double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    double p_aa_htht__hyhy_h = +0.5*(a*a*ht*ht)/(hy*hy);
    double m_aa_htht__hyhy_h = -0.5*(a*a*ht*ht)/(hy*hy);
    double p_aa_htht__hxhx___lambda_ht = +1.0 + (a*a*ht*ht)/(hx*hx) + 1.5*(lambda*ht);
    double p_aa_htht__hyhy___lambda_ht = +1.0 + (a*a*ht*ht)/(hy*hy) + 1.5*(lambda*ht);
    double aa__hxhx = (a*a)/(hx*hx);
    double aa__hyhy = (a*a)/(hy*hy);
    double htht_h = 0.5*ht*ht;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    u.clear(); u.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<ExtendedSpacePointNode> obsPointNodes, cntDeltaNodes, qPointNodes;
    extendPoints(obsPointNodes, cntDeltaNodes, qPointNodes, dimX, dimY);
    //--------------------------------------------------------------------------------------------//
    std::vector<unsigned int> rows0, rows1, rows2, cols0, cols1, cols2;
    findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, obsPointNodes, cntDeltaNodes, N, M);
    //--------------------------------------------------------------------------------------------//

    //-------------------------------------------- info --------------------------------------------//
    if (use == true) prepareInfo(No, mOptParameter.xi, info, L, dimX, dimY);
    //-------------------------------------------- info --------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    initialLayers(u00, u05, u10, info, use, obsPointNodes, cntDeltaNodes, qPointNodes, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, lambda);
    //------------------------------------- initial conditions -------------------------------------//

    double *a1X, *b1X, *c1X, *d1X, *x1X, *a1Y, *b1Y, *c1Y, *d1Y, *x1Y;
    mallocVectors(a1X, b1X, c1X, d1X, x1X, N-1, a1Y, b1Y, c1Y, d1Y, x1Y, M-1);

    SpaceNodePDE sn;
    TimeNodePDE tn;
    for (unsigned int l=2; l<=L; l++)
    {
        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht;
        initBorders(N, M, hx, hy, ht, l, u15, u);

        //--------------------------------------------------------------------------//
        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = n; sn.x = n*hx;

                //                if (m == 0)     d1X[n-1] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                //                if (m>0 && m<M) d1X[n-1] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                //                if (m == M)     d1X[n-1] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                d1X[n-1] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);

                d1X[n-1] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                d1X[n-1] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

                a1X[n-1] = m_aa_htht__hxhx_h;
                b1X[n-1] = p_aa_htht__hxhx___lambda_ht;
                c1X[n-1] = m_aa_htht__hxhx_h;
            }

            a1X[0]   = 0.0;
            c1X[N-2] = 0.0;

            d1X[0]   -= m_aa_htht__hxhx_h * u15[m][0];
            d1X[N-2] -= m_aa_htht__hxhx_h * u15[m][N];

            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = x1X[n-1];
        }

        if (rows2.size() == 0)
        {
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    //                    if (m == 0)     d1X[n-1] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    //                    if (m>0 && m<M) d1X[n-1] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    //                    if (m == M)     d1X[n-1] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1X[n-1] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);

                    d1X[n-1] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                    d1X[n-1] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

                    a1X[n-1] = m_aa_htht__hxhx_h;
                    b1X[n-1] = p_aa_htht__hxhx___lambda_ht;
                    c1X[n-1] = m_aa_htht__hxhx_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                                d1X[n-1] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            }
                            for (unsigned int j=0; j<No; j++)
                            {
                                d1X[n-1] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1X[0]   = 0.0;
                c1X[N-2] = 0.0;

                d1X[0]   -= m_aa_htht__hxhx_h * u15[m][0];
                d1X[N-2] -= m_aa_htht__hxhx_h * u15[m][N];

                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = x1X[n-1];
            }
        }

        if (rows2.size() != 0)
        {
            double* a1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            double* b1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            double* c1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            double* d1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            double* x1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            DoubleMatrix w2(rows1.size()*(N-1), rows1.size()*(N-1), 0.0);

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    //                    if (m==0)       d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    //                    if (m>0 && m<M) d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    //                    if (m==M)       d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);

                    d1[offset+(n-1)] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                    d1[offset+(n-1)] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

                    a1[offset+(n-1)] = m_aa_htht__hxhx_h;
                    b1[offset+(n-1)] = p_aa_htht__hxhx___lambda_ht;
                    c1[offset+(n-1)] = m_aa_htht__hxhx_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(odj);

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (opn.j == rows1[rs])
                                    {
                                        found = true;
                                        w2[offset+(n-1)][rs*(N-1)+(opn.i-1)] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                    }
                                }

                                if (found == false)
                                {
                                    d1[offset+(n-1)] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[offset+(n-1)] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx_h * u15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx_h * u15[m][N];

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w2.data(), x1, rows1.size()*(N-1));

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=1; n<=N-1; n++)
                {
                    u15[m][n] = x1[offset+(n-1)];
                }
                offset += N-1;
            }

            w2.clear();
            free(x1);
            free(d1);
            free(c1);
            free(b1);
            free(a1);
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
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = m; sn.y = m*hy;

                //                if (n==0)       d1Y[m-1] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                //                if (n>0 && n<N) d1Y[m-1] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                //                if (n==N)       d1Y[m-1] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                d1Y[m-1] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);

                d1Y[m-1] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                d1Y[m-1] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

                a1Y[m-1] = m_aa_htht__hyhy_h;
                b1Y[m-1] = p_aa_htht__hyhy___lambda_ht;
                c1Y[m-1] = m_aa_htht__hyhy_h;
            }

            a1Y[0]   = 0.0;
            c1Y[M-2] = 0.0;

            d1Y[0]   -= m_aa_htht__hyhy_h * u[0][n];
            d1Y[M-2] -= m_aa_htht__hyhy_h * u[M][n];

            tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M-1);
            for (unsigned int m=1; m<=M-1; m++) u[m][n] = x1Y[m-1];
        }

        if (cols2.size() == 0)
        {
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    //                    if (n==0)       d1Y[m-1] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    //                    if (n>0 && n<N) d1Y[m-1] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    //                    if (n==N)       d1Y[m-1] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d1Y[m-1] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);

                    d1Y[m-1] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                    d1Y[m-1] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

                    a1Y[m-1] = m_aa_htht__hyhy_h;
                    b1Y[m-1] = p_aa_htht__hyhy___lambda_ht;
                    c1Y[m-1] = m_aa_htht__hyhy_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);
                                d1Y[m-1] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1Y[m-1] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1Y[0]   = 0.0;
                c1Y[M-2] = 0.0;

                d1Y[0]   -= m_aa_htht__hyhy_h * u[0][n];
                d1Y[M-2] -= m_aa_htht__hyhy_h * u[M][n];

                tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M-1);
                for (unsigned int m=1; m<=M-1; m++) u[m][n] = x1Y[m-1];
            }
        }

        if (cols2.size() != 0)
        {
            double* a2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* b2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* c2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* d2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* x2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            DoubleMatrix w2(cols1.size()*(M-1), cols1.size()*(M-1), 0.0);

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    //                    if (n==0)       d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    //                    if (n>0 && n<N) d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    //                    if (n==N)       d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);

                    d2[offset+(m-1)] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                    d2[offset+(m-1)] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

                    a2[offset+(m-1)] = m_aa_htht__hyhy_h;
                    b2[offset+(m-1)] = p_aa_htht__hyhy___lambda_ht;
                    c2[offset+(m-1)] = m_aa_htht__hyhy_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (opn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+(m-1)][cs*(M-1)+(opn.j-0)] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+(m-1)] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[offset+(m-1)] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy_h * u[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy_h * u[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M-1));

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    u[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }

            w2.clear();
            free(x2);
            free(d2);
            free(c2);
            free(b2);
            free(a2);
        }
        //--------------------------------------------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        // if (l==L)
        {
            ut.clear();
            ut.resize(M+1, N+1);
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
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

        if (use == true) add2Info(u, info, l, _INFO_ROWS_, _INFO_COLS_);
        layerInfo(u, l);
    }

    freeVectors(a1X, b1X, c1X, d1X, x1X, a1Y, b1Y, c1Y, d1Y, x1Y);

    rows0.clear();
    rows1.clear();
    rows2.clear();

    cols0.clear();
    cols1.clear();
    cols2.clear();

    qPointNodes.clear();
    obsPointNodes.clear();
    cntDeltaNodes.clear();

    u00.clear();
    u05.clear();
    u10.clear();
    u15.clear();
}

void IProblem2HForward2D::calculateMVD_D2(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const
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
    //unsigned int Nc = mEquParameter.Nc;
    //unsigned int Ns = mEquParameter.Ns;

    double aa__hxhx = (a*a)/(hx*hx);
    double aa__hyhy = (a*a)/(hy*hy);
    double htht = ht*ht;
    double lambda_ht = lambda*ht;

    double p_aa_htht__hyhy = +(a*a*ht*ht)/(hy*hy);
    double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    double p_aa_htht__hxhx___lambda_ht = +8.0 + (2.0*a*a*ht*ht)/(hx*hx) + 3.0*lambda*ht;

    double p_aa_htht__hxhx = +(a*a*ht*ht)/(hx*hx);
    double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    double p_aa_htht__hyhy___lambda_ht = +8.0 + (2.0*a*a*ht*ht)/(hy*hy) + 3.0*lambda*ht;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    u.clear(); u.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<ExtendedSpacePointNode> obsPointNodes, cntDeltaNodes, qPointNodes;
    extendPoints(obsPointNodes, cntDeltaNodes, qPointNodes, dimX, dimY);
    //--------------------------------------------------------------------------------------------//
    std::vector<unsigned int> rows0, rows1, rows2, cols0, cols1, cols2;
    findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, obsPointNodes, cntDeltaNodes, N, M);
    //--------------------------------------------------------------------------------------------//

    //-------------------------------------------- info --------------------------------------------//
    if (use == true)
    {
        info.resize(No);
        for (unsigned int j=0; j<No; j++)
        {
            ExtendedSpaceNode2DH &e = info[j];
            e.id = j;
            e.setSpaceNode(mOptParameter.xi[j]);
            e.extendWeights(dimX, dimY, L+1, _INFO_ROWS_, _INFO_COLS_);
        }
    }
    //-------------------------------------------- info --------------------------------------------//

    SpaceNodePDE sn;
    TimeNodePDE tn;
    //------------------------------------- initial conditions -------------------------------------//
    initialLayers(u00, u05, u10, info, use, obsPointNodes, cntDeltaNodes, qPointNodes, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, lambda);
    //------------------------------------- initial conditions -------------------------------------//
    double *a1X, *b1X, *c1X, *d1X, *x1X, *a1Y, *b1Y, *c1Y, *d1Y, *x1Y;
    mallocVectors(a1X, b1X, c1X, d1X, x1X, N-1, a1Y, b1Y, c1Y, d1Y, x1Y, M-1);

    for (unsigned int l=2; l<=L; l++)
    {
        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht;
        initBorders(N, M, hx, hy, ht, l, u15, u);

        //--------------------------------------------------------------------------//
        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = n; sn.x = n*hx;

                d1X[n-1] = 0.0;
                if (m == 0)     d1X[n-1] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                if (m>0 && m<M) d1X[n-1] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                if (m == M)     d1X[n-1] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                d1X[n-1] += 4.0*(5.0*u10[m][n] - 4.0*u05[m][n] + u00[m][n]);
                //d1X[n-1] += 4.0*(2.0*u10[m][n] - u05[m][n]);
                d1X[n-1] += lambda_ht*(4.0*u10[m][n] - u05[m][n]);

                a1X[n-1] = m_aa_htht__hxhx;
                b1X[n-1] = p_aa_htht__hxhx___lambda_ht;
                c1X[n-1] = m_aa_htht__hxhx;
            }

            a1X[0]   = 0.0;
            c1X[N-2] = 0.0;

            d1X[0]   -= m_aa_htht__hxhx * u15[m][0];
            d1X[N-2] -= m_aa_htht__hxhx * u15[m][N];

            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = x1X[n-1];
        }

        if (rows2.size() == 0)
        {
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    d1X[n-1] = 0.0;
                    if (m == 0)     d1X[n-1] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    if (m>0 && m<M) d1X[n-1] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    if (m == M)     d1X[n-1] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1X[n-1] += 4.0*(5.0*u10[m][n] - 4.0*u05[m][n] + u00[m][n]);
                    //d1X[n-1] += 4.0*(2.0*u10[m][n] - u05[m][n]);
                    d1X[n-1] += lambda_ht*(4.0*u10[m][n] - u05[m][n]);

                    a1X[n-1] = m_aa_htht__hxhx;
                    b1X[n-1] = p_aa_htht__hxhx___lambda_ht;
                    c1X[n-1] = m_aa_htht__hxhx;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                                d1X[n-1] += htht * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1X[n-1] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1X[0]   = 0.0;
                c1X[N-2] = 0.0;

                d1X[0]   -= m_aa_htht__hxhx * u15[m][0];
                d1X[N-2] -= m_aa_htht__hxhx * u15[m][N];

                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = x1X[n-1];
            }
        }

        //        if (rows2.size() != 0)
        //        {
        //            double* a1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
        //            double* b1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
        //            double* c1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
        //            double* d1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
        //            double* x1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
        //            DoubleMatrix w2(rows1.size()*(N-1), rows1.size()*(N-1), 0.0);

        //            unsigned int offset = 0;
        //            for (unsigned int row=0; row<rows1.size(); row++)
        //            {
        //                unsigned int m = rows1.at(row);
        //                sn.j = m; sn.y = m*hy;

        //                for (unsigned int n=1; n<=N-1; n++)
        //                {
        //                    sn.i = n; sn.x = n*hx;

        //                    if (m==0)       d1[offset+(n-1)] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
        //                    if (m>0 && m<M) d1[offset+(n-1)] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
        //                    if (m==M)       d1[offset+(n-1)] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

        //                    d1[offset+(n-1)] += 4.0*(5.0*u10[m][n] - 4.0*u05[m][n] + u00[m][n]);
        //                    //d1[offset+(n-1)] += 4.0*(2.0*u10[m][n] - u05[m][n]);
        //                    d1[offset+(n-1)] += lambda_ht*(4.0*u10[m][n] - u05[m][n]);

        //                    a1[offset+(n-1)] = m_aa_htht__hxhx;
        //                    b1[offset+(n-1)] = p_aa_htht__hxhx___lambda_ht;
        //                    c1[offset+(n-1)] = m_aa_htht__hxhx;

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
        //                                for (unsigned int rs=0; rs<rows1.size(); rs++)
        //                                {
        //                                    if (opn.j == rows1[rs])
        //                                    {
        //                                        found = true;
        //                                        w2[offset+(n-1)][rs*(N-1)+(opn.i-1)] -= htht * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
        //                                    }
        //                                }

        //                                if (found == false)
        //                                {
        //                                    d1[offset+(n-1)] += htht * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
        //                                }
        //                            }

        //                            for (unsigned int j=0; j<No; j++)
        //                            {
        //                                d1[offset+(n-1)] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
        //                            }
        //                        }
        //                    }
        //                    //------------------------------------- Adding delta part -------------------------------------//
        //                }

        //                a1[offset+0]   = 0.0;
        //                c1[offset+N-2] = 0.0;

        //                d1[offset+0]   -= m_aa_htht__hxhx * u15[m][0];
        //                d1[offset+N-2] -= m_aa_htht__hxhx * u15[m][N];

        //                offset += N-1;
        //            }

        //            LinearEquation::func1(a1, b1, c1, d1, w2.data(), x1, rows1.size()*(N-1));

        //            offset = 0;
        //            for (unsigned int row=0; row<rows1.size(); row++)
        //            {
        //                unsigned int m=rows1.at(row);
        //                for (unsigned int n=1; n<=N-1; n++)
        //                {
        //                    u15[m][n] = x1[offset+(n-1)];
        //                }
        //                offset += N-1;
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
        //tn.i = l+1; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//

        for (unsigned int col=0; col<cols0.size(); col++)
        {
            unsigned int n = cols0.at(col);
            sn.i = n; sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = m; sn.y = m*hy;

                d1Y[m-1] = 0.0;

                if (n==0)       d1Y[m-1] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                if (n>0 && n<N) d1Y[m-1] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                if (n==N)       d1Y[m-1] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                d1Y[m-1] += 4.0*(5.0*u15[m][n] - 4.0*u10[m][n] + u05[m][n]);
                //d1Y[m-1] += 4.0*(2.0*u15[m][n] - u10[m][n]);
                d1Y[m-1] += lambda_ht*(4.0*u15[m][n] - u10[m][n]);

                a1Y[m-1] = m_aa_htht__hyhy;
                b1Y[m-1] = p_aa_htht__hyhy___lambda_ht;
                c1Y[m-1] = m_aa_htht__hyhy;
            }

            a1Y[0]   = 0.0;
            c1Y[M-2] = 0.0;

            d1Y[0]   -= m_aa_htht__hyhy * u[0][n];
            d1Y[M-2] -= m_aa_htht__hyhy * u[M][n];

            tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M-1);
            for (unsigned int m=1; m<=M-1; m++) u[m][n] = x1Y[m-1];
        }

        if (cols2.size() == 0)
        {
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    d1Y[m-1] = 0.0;

                    if (n==0)       d1Y[m-1] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    if (n>0 && n<N) d1Y[m-1] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    if (n==N)       d1Y[m-1] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d1Y[m-1] += 4.0*(5.0*u15[m][n] - 4.0*u10[m][n] + u05[m][n]);
                    //d1Y[m-1] += 4.0*(2.0*u15[m][n] - u10[m][n]);
                    d1Y[m-1] += lambda_ht*(4.0*u15[m][n] - u10[m][n]);

                    a1Y[m-1] = m_aa_htht__hyhy;
                    b1Y[m-1] = p_aa_htht__hyhy___lambda_ht;
                    c1Y[m-1] = m_aa_htht__hyhy;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);
                                d1Y[m-1] += htht * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1Y[m-1] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1Y[0]   = 0.0;
                c1Y[M-2] = 0.0;

                d1Y[0]   -= m_aa_htht__hyhy * u[0][n];
                d1Y[M-2] -= m_aa_htht__hyhy * u[M][n];

                tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M-1);
                for (unsigned int m=1; m<=M-1; m++) u[m][n] = x1Y[m-1];
            }
        }

        //        if (cols2.size() != 0)
        //        {
        //            double* a2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
        //            double* b2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
        //            double* c2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
        //            double* d2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
        //            double* x2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
        //            DoubleMatrix w2(cols1.size()*(M-1), cols1.size()*(M-1), 0.0);

        //            unsigned int offset = 0;
        //            for (unsigned int col=0; col<cols1.size(); col++)
        //            {
        //                unsigned int n = cols1.at(col);
        //                sn.i = n; sn.x = n*hx;

        //                for (unsigned int m=1; m<=M-1; m++)
        //                {
        //                    sn.j = m; sn.y = m*hy;

        //                    if (n==0)       d2[offset+(m-1)] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
        //                    if (n>0 && n<N) d2[offset+(m-1)] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
        //                    if (n==N)       d2[offset+(m-1)] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

        //                    d2[offset+(m-1)] += 4.0*(5.0*u15[m][n] - 4.0*u10[m][n] + u05[m][n]);
        //                    //d2[offset+(m-1)] += 4.0*(2.0*u15[m][n] - u10[m][n]);
        //                    d2[offset+(m-1)] += lambda_ht*(4.0*u15[m][n] - u10[m][n]);

        //                    a2[offset+(m-1)] = m_aa_htht__hyhy;
        //                    b2[offset+(m-1)] = p_aa_htht__hyhy___lambda_ht;
        //                    c2[offset+(m-1)] = m_aa_htht__hyhy;

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
        //                                        w2[offset+(m-1)][cs*(M-1)+(opn.j-0)] -= htht * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
        //                                    }
        //                                }

        //                                if (!found)
        //                                {
        //                                    d2[offset+(m-1)] += htht * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
        //                                }
        //                            }

        //                            for (unsigned int j=0; j<No; j++)
        //                            {
        //                                d2[offset+(m-1)] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
        //                            }
        //                        }
        //                    }
        //                    //------------------------------------- Adding delta part -------------------------------------//
        //                }

        //                a2[offset+0]   = 0.0;
        //                c2[offset+M-2] = 0.0;

        //                d2[offset+0]   -= m_aa_htht__hyhy * u[0][n];
        //                d2[offset+M-2] -= m_aa_htht__hyhy * u[M][n];

        //                offset += M-1;
        //            }

        //            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M-1));

        //            offset = 0;
        //            for (unsigned int col=0; col<cols1.size(); col++)
        //            {
        //                unsigned int n=cols1.at(col);
        //                for (unsigned int m=1; m<=M-1; m++)
        //                {
        //                    u[m][n] = x2[offset+(m-1)];
        //                }
        //                offset += M-1;
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

        // if (l==L)
        {
            ut.clear();
            ut.resize(M+1, N+1);
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
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

        if (use == true) add2Info(u, info, l, _INFO_ROWS_, _INFO_COLS_);
        layerInfo(u, l);
    }

    freeVectors(a1X, b1X, c1X, d1X, x1X, a1Y, b1Y, c1Y, d1Y, x1Y);

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

void IProblem2HForward2D::calculateMVD_D3(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const
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

    double lambda = mEquParameter.lambda;
    double a = mEquParameter.a;
    unsigned int No = mEquParameter.No;
    unsigned int Nc = mEquParameter.Nc;
    unsigned int Ns = mEquParameter.Ns;

    //----------------------------------------------------------------------------------------------//

    double m_aa_htht__hxhx_h = -(a*a*hh*hh)/(hx*hx);
    double p_aa_htht__hxhx___lambda_ht = +2.0 + 2.0*(a*a*hh*hh)/(hx*hx) + (11.0/6.0)*(lambda*hh);
    double p_aa_htht__hyhy_h = +(a*a*hh*hh)/(hy*hy);

    double m_aa_htht__hyhy_h = -(a*a*hh*hh)/(hy*hy);
    double p_aa_htht__hyhy___lambda_ht = +2.0 + 2.0*(a*a*hh*hh)/(hy*hy) + (11.0/6.0)*(lambda*hh);
    double p_aa_htht__hxhx_h = +(a*a*hh*hh)/(hx*hx);

    double lambda_hh__2 = 0.5*lambda*hh;
    double htht_h = hh*hh;

    double aa__hxhx = (a*a)/(hx*hx);
    double aa__hyhy = (a*a)/(hy*hy);

    //----------------------------------------------------------------------------------------------//

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    u.clear(); u.resize(M+1, N+1);

    //----------------------------------------------------------------------------------------------//
    std::vector<ExtendedSpacePointNode> obsPointNodes;
    for (unsigned int j=0; j<No; j++) IProblem2H2D::distributeDelta(mOptParameter.xi[j], obsPointNodes, j, dimX, dimY, 4);

    std::vector<ExtendedSpacePointNode> cntDeltaNodes;
    for (unsigned int i=0; i<Nc; i++) IProblem2H2D::distributeDelta(mOptParameter.eta[i], cntDeltaNodes, i, dimX, dimY, 4);

    std::vector<ExtendedSpacePointNode> qPointNodes;
    for (unsigned int s=0; s<Ns; s++) IProblem2H2D::distributeDelta(mEquParameter.theta[s], qPointNodes, s, dimX, dimY, 4);

    SpaceNodePDE sn;

    //----------------------------------------------------------------------------------------------//
    std::vector<unsigned int> rows0, rows1, rows2;
    for (unsigned int m=1; m<=M-1; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);

            if (cdn.j == m)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.j == m)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
    }

    //----------------------------------------------------------------------------------------------//
    vector<unsigned int> cols0, cols1, cols2;
    for (unsigned int n=1; n<=N-1; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);
            if (cdn.i == n)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.i == n)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
    }

    //----------------------------------------------------------------------------------------------//

    //-------------------------------------------- info --------------------------------------------//
    if (use == true)
    {
        info.resize(No);
        for (unsigned int j=0; j<No; j++)
        {
            ExtendedSpaceNode2DH &e = info[j];
            e.id = j;
            e.setSpaceNode(mOptParameter.xi[j]);
            e.extendWeights(dimX, dimY, L+1, _INFO_ROWS_, _INFO_COLS_);
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
    if (use == true) add2Info(u00, info, 0, _INFO_ROWS_, _INFO_COLS_);
    layerInfo(u00, 0);
    layerInfo(u00, u00, 0);
    //------------------------------------- initial conditions -------------------------------------//

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum1 = 0.0;

            if (n==0)      sum1 += aa__hxhx*(u00[m][n]-2.0*u00[m][n+1]+u00[m][n+2]);
            else if (n==N) sum1 += aa__hxhx*(u00[m][n-2]-2.0*u00[m][n-1]+u00[m][n]);
            else           sum1 += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);

            if (m==0)      sum1 += aa__hyhy*(u00[m][n]-2.0*u00[m+1][n]+u00[m+2][n]);
            else if (m==M) sum1 += aa__hyhy*(u00[m-2][n]-2.0*u00[m-1][n]+u00[m][n]);
            else           sum1 += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);

            sum1 -= lambda*initial2(sn);

            for (unsigned int cdn=0; cdn<cntDeltaNodes.size(); cdn++)
            {
                const ExtendedSpacePointNode &cntNode = cntDeltaNodes.at(cdn);
                if (cntNode.i == n && cntNode.j == m)
                {
                    double kz = 0.0;
                    for (unsigned int j=0; j<No; j++)
                    {
                        kz += mOptParameter.k[cntNode.id][j]*mOptParameter.z[cntNode.id][j];
                    }
                    //sum1 -= kz * cntNode.w;
                }
            }

            u05[m][n] = u00[m][n] + hh*initial2(sn) + 0.5*hh*hh*sum1;
            u10[m][n] = u00[m][n] + ht*initial2(sn) + 0.5*ht*ht*sum1;
        }
    }

    if (use == true) add2Info(u10, info, 1, _INFO_ROWS_, _INFO_COLS_);
    layerInfo(u10, u10, 1);
    layerInfo(u10, 1);
    //------------------------------------- initial conditions -------------------------------------//

    unsigned int sizeN = sizeof(double)*(N-1);
    unsigned int sizeM = sizeof(double)*(M-1);

    double *ax = (double *) malloc(sizeN); for (unsigned int n=1; n<=N-1; n++) ax[n-1] = m_aa_htht__hxhx_h;           ax[0] = 0.0;
    double *bx = (double *) malloc(sizeN); for (unsigned int n=1; n<=N-1; n++) bx[n-1] = p_aa_htht__hxhx___lambda_ht;
    double *cx = (double *) malloc(sizeN); for (unsigned int n=1; n<=N-1; n++) cx[n-1] = m_aa_htht__hxhx_h;           cx[N-2] = 0.0;
    double *dx = (double *) malloc(sizeN);
    double *rx = (double *) malloc(sizeN);

    double *ay = (double *) malloc(sizeM); for (unsigned int m=1; m<=M-1; m++) ay[m-1] = m_aa_htht__hyhy_h;           ay[0] = 0.0;
    double *by = (double *) malloc(sizeM); for (unsigned int m=1; m<=M-1; m++) by[m-1] = p_aa_htht__hyhy___lambda_ht;
    double *cy = (double *) malloc(sizeM); for (unsigned int m=1; m<=M-1; m++) cy[m-1] = m_aa_htht__hxhx_h;           cy[M-2] = 0.0;
    double *dy = (double *) malloc(sizeM);
    double *ry = (double *) malloc(sizeM);

    TimeNodePDE tn;

    for (unsigned int l=2; l<=L; l++)
    {
        //------------------------------------- border -------------------------------------//
        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy;
            sn1.j = m; sn1.y = m*hy;

            u15[m][0] = boundary(sn0, tn);
            u15[m][N] = boundary(sn1, tn);

            u[m][0] = boundary(sn0, tn);
            u[m][N] = boundary(sn1, tn);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx;
            sn1.i = n; sn1.x = n*hx;

            u15[0][n] = boundary(sn0, tn);
            u15[M][n] = boundary(sn1, tn);

            u[0][n] = boundary(sn0, tn);
            u[M][n] = boundary(sn1, tn);
        }
        //------------------------------------- border -------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht;

        //--------------------------------------------------------------------------//

        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = n; sn.x = n*hx;

                if (m == 0)     dx[n-1] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                if (m>0 && m<M) dx[n-1] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                if (m == M)     dx[n-1] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                dx[n-1] += (5.0*u10[m][n] - 4.0*u05[m][n] + u00[m][n]);
                dx[n-1] += lambda_hh__2*(6.0*u10[m][n] - 3.0*u05[m][n] + 8.0*u00[m][n]);
            }

            dx[0]   -= m_aa_htht__hxhx_h * u15[m][0];
            dx[N-2] -= m_aa_htht__hxhx_h * u15[m][N];

            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }

        if (rows2.size() == 0)
        {
            double U15[No]; for (unsigned int j=0; j<No; j++) U15[j] = 0.0;
            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                U15[opn.id] += u15[opn.j][opn.i] * (opn.w * (hx*hy));
            }

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (m == 0)     dx[n-1] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    if (m>0 && m<M) dx[n-1] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    if (m == M)     dx[n-1] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    dx[n-1] += (5.0*u10[m][n] - 4.0*u05[m][n] + u00[m][n]);
                    dx[n-1] += lambda_hh__2*(6.0*u10[m][n] - 3.0*u05[m][n] + 8.0*u00[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int j=0; j<No; j++)
                            {
                                dx[n-1] += htht_h * mOptParameter.k[cdn.id][j] * (U15[j]-mOptParameter.z[cdn.id][j]) * cdn.w;
                            }
                            //                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            //                            {
                            //                                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                            //                                d1X[n-1] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            //                            }
                            //                            for (unsigned int j=0; j<No; j++)
                            //                            {
                            //                                d1X[n-1] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            //                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                dx[0]   -= m_aa_htht__hxhx_h * u15[m][0];
                dx[N-2] -= m_aa_htht__hxhx_h * u15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }
        }

        if (rows2.size() != 0)
        {
            //throw std::exception();

            unsigned int row1_size = rows1.size()*(N-1);
            double* a1 = (double*) malloc(sizeof(double)*row1_size);
            double* b1 = (double*) malloc(sizeof(double)*row1_size);
            double* c1 = (double*) malloc(sizeof(double)*row1_size);
            double* d1 = (double*) malloc(sizeof(double)*row1_size);
            double* x1 = (double*) malloc(sizeof(double)*row1_size);
            DoubleMatrix w2(row1_size, row1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (m==0)       d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    if (m>0 && m<M) d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    if (m==M)       d1[offset+(n-1)] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1[offset+(n-1)] += (5.0*u10[m][n] - 4.0*u05[m][n] + u00[m][n]);
                    d1[offset+(n-1)] += lambda_hh__2*(6.0*u10[m][n] - 3.0*u05[m][n] + 8.0*u00[m][n]);

                    a1[offset+(n-1)] = m_aa_htht__hxhx_h;
                    b1[offset+(n-1)] = p_aa_htht__hxhx___lambda_ht;
                    c1[offset+(n-1)] = m_aa_htht__hxhx_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(odj);

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (opn.j == rows1[rs])
                                    {
                                        found = true;
                                        w2[offset+(n-1)][rs*(N-1)+(opn.i-1)] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                        break;
                                    }
                                }

                                if (found == false)
                                {
                                    d1[offset+(n-1)] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }
                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[offset+(n-1)] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx_h * u15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx_h * u15[m][N];

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w2.data(), x1, row1_size);

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=1; n<=N-1; n++)
                {
                    u15[m][n] = x1[offset+(n-1)];
                }
                offset += N-1;
            }

            w2.clear();
            free(x1);
            free(d1);
            free(c1);
            free(b1);
            free(a1);
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
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = m; sn.y = m*hy;

                if (n==0)       dy[m-1] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                if (n==N)       dy[m-1] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                dy[m-1] += (5.0*u15[m][n] - 4.0*u10[m][n] + u05[m][n]);
                dy[m-1] += lambda_hh__2*(6.0*u15[m][n] - 3.0*u10[m][n] + 8.0*u05[m][n]);
            }

            dy[0]   -= m_aa_htht__hyhy_h * u[0][n];
            dy[M-2] -= m_aa_htht__hyhy_h * u[M][n];

            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u[m][n] = ry[m-1];
        }

        if (cols2.size() == 0)
        {
            double U[No]; for (unsigned int j=0; j<No; j++) U[j] = 0.0;
            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                U[opn.id] += u[opn.j][opn.i] * (opn.w * (hx*hy));
            }

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    if (n==0)       dy[m-1] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    if (n==N)       dy[m-1] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    dy[m-1] += (5.0*u15[m][n] - 4.0*u10[m][n] + u05[m][n]);
                    dy[m-1] += lambda_hh__2*(6.0*u15[m][n] - 3.0*u10[m][n] + 8.0*u05[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int j=0; j<No; j++)
                            {
                                dy[m-1] += htht_h * mOptParameter.k[cdn.id][j] * (U[j]-mOptParameter.z[cdn.id][j]) * cdn.w;
                            }
                            //                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            //                            {
                            //                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);
                            //                                d1Y[m-1] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            //                            }
                            //                            for (unsigned int j=0; j<No; j++)
                            //                            {
                            //                                d1Y[m-1] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            //                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                dy[0]   -= m_aa_htht__hyhy_h * u[0][n];
                dy[M-2] -= m_aa_htht__hyhy_h * u[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u[m][n] = ry[m-1];
            }
        }

        if (cols2.size() != 0)
        {
            //throw std::exception();

            double* a2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* b2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* c2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* d2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            double* x2 = (double*) malloc(sizeof(double)*cols1.size()*(M-1));
            DoubleMatrix w2(cols1.size()*(M-1), cols1.size()*(M-1), 0.0);

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    if (n==0)       d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    if (n>0 && n<N) d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    if (n==N)       d2[offset+(m-1)] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d2[offset+(m-1)] += (5.0*u15[m][n] - 4.0*u10[m][n] + u05[m][n]);
                    d2[offset+(m-1)] += lambda_hh__2*(6.0*u15[m][n] - 3.0*u10[m][n] + 8.0*u05[m][n]);

                    a2[offset+(m-1)] = m_aa_htht__hyhy_h;
                    b2[offset+(m-1)] = p_aa_htht__hyhy___lambda_ht;
                    c2[offset+(m-1)] = m_aa_htht__hyhy_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (opn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+(m-1)][cs*(M-1)+(opn.j-0)] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+(m-1)] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }
                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[offset+(m-1)] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy_h * u[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy_h * u[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M-1));

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    u[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }

            w2.clear();
            free(x2);
            free(d2);
            free(c2);
            free(b2);
            free(a2);
        }
        //--------------------------------------------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        if (l==2)
        {
            puts("------------------------------------------");
            for (unsigned int n=0; n<=N; n++)
            {
                for (unsigned int m=0; m<=M; m++)
                {
                    for (unsigned int si=0; si<qPointNodes.size(); si++)
                    {
                        const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                        if (qNode.i == n && qNode.j == m)
                        {
                            u[m][n] += (mEquParameter.q[qNode.id] * qNode.w * (1.0/ht))*ht*ht;
                        }
                    }
                }
            }
        }

        if (l==L)
        {
            ut.clear();
            ut.resize(M+1, N+1);
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
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

        if (use == true) add2Info(u, info, l, _INFO_ROWS_, _INFO_COLS_);
        layerInfo(u, l);
    }

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    rows0.clear();
    rows1.clear();
    rows2.clear();

    cols0.clear();
    cols1.clear();
    cols2.clear();

    qPointNodes.clear();
    obsPointNodes.clear();
    cntDeltaNodes.clear();

    u00.clear();
    u05.clear();
    u10.clear();
    u15.clear();
}

void IProblem2HForward2D::layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const
{
    //    IPrinter::printSeperatorLine();
    //    IPrinter::printMatrix(u);
    //    IPrinter::printSeperatorLine();

    //    if (ln%1000==0)
    //    {
    printf("%4d %18.6f %18.6f\n", ln, u.min(), u.max());
    QPixmap px;
    visualizeMatrixHeat(u, u.min(), u.max(), px);
    char buffer[30] = {0};
    int c = sprintf(buffer, "images/s/%6d.png", ln);
    buffer[c] = 0;
    px.save(buffer, "png", 0);
    //    }

    //    if (ln==0)
    //    {
    //        FILE *file = fopen("images/0.txt", "w");
    //        IPrinter::print(u, u.rows(), u.cols(), 14, 2, file);
    //        fclose(file);
    //    }
}

void IProblem2HForward2D::layerInfo(const DoubleMatrix &u UNUSED_PARAM, const DoubleMatrix &ut UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const
{
    //    double hx = mspaceDimension[0].step();
    //    double hy = mspaceDimension[1].step();
    //    unsigned int N1 = mspaceDimension[0].sizeN();
    //    unsigned int N2 = mspaceDimension[1].sizeN();

    //    DoubleMatrix V0(101, 101, 0.0);
    //    DoubleMatrix V1(101, 101, 0.0);

    //    double sum0 = 0.0;
    //    double sum1 = 0.0;

    //    sum0 += 0.25*(u[0][0]   - V0[0][0])   * (u[0][0]   - V0[0][0]);
    //    sum0 += 0.25*(u[0][N1]  - V0[0][N1])  * (u[0][N1]  - V0[0][N1]);
    //    sum0 += 0.25*(u[N2][0]  - V0[N2][0])  * (u[N2][0]  - V0[N2][0]);
    //    sum0 += 0.25*(u[N2][N1] - V0[N2][N1]) * (u[N2][N1] - V0[N2][N1]);

    //    sum1 += 0.25*(ut[0][0]   - V1[0][0])   * (ut[0][0]   - V1[0][0]);
    //    sum1 += 0.25*(ut[0][N1]  - V1[0][N1])  * (ut[0][N1]  - V1[0][N1]);
    //    sum1 += 0.25*(ut[N2][0]  - V1[N2][0])  * (ut[N2][0]  - V1[N2][0]);
    //    sum1 += 0.25*(ut[N2][N1] - V1[N2][N1]) * (ut[N2][N1] - V1[N2][N1]);

    //    for (unsigned int n1=1; n1<=N1-1; n1++)
    //    {
    //        sum0 += 0.5*(u[0][n1]  - V0[0][n1]) *(u[0][n1]  - V0[0][n1]);
    //        sum0 += 0.5*(u[N2][n1] - V0[N2][n1])*(u[N2][n1] - V0[N2][n1]);
    //        sum1 += 0.5*(ut[0][n1]  - V1[0][n1]) *(ut[0][n1]  - V1[0][n1]);
    //        sum1 += 0.5*(ut[N2][n1] - V1[N2][n1])*(ut[N2][n1] - V1[N2][n1]);
    //    }

    //    for (unsigned int n2=1; n2<=N2-1; n2++)
    //    {
    //        sum0 += 0.5*(u[n2][0]  - V0[n2][0]) *(u[n2][0]  - V0[n2][0]);
    //        sum0 += 0.5*(u[n2][N1] - V0[n2][N1])*(u[n2][N1] - V0[n2][N1]);
    //        sum1 += 0.5*(ut[n2][0]  - V1[n2][0]) *(ut[n2][0]  - V1[n2][0]);
    //        sum1 += 0.5*(ut[n2][N1] - V1[n2][N1])*(ut[n2][N1] - V1[n2][N1]);
    //    }

    //    for (unsigned int n2 = 1; n2 <= N2-1; n2++)
    //    {
    //        for (unsigned int n1 = 1; n1 <= N1-1; n1++)
    //        {
    //            sum0 += (u[n2][n1] - V0[n2][n1])*(u[n2][n1] - V0[n2][n1]);
    //            sum1 += (ut[n2][n1] - V1[n2][n1])*(ut[n2][n1] - V1[n2][n1]);
    //        }
    //    }

    //    sum0 *= (hx*hy);
    //    sum1 *= (hx*hy);

    //    V0.clear();
    //    V1.clear();

    //    printf("%6d %20.8f %20.8f u:%20.8f %20.8f ut:%20.8f %20.8f\n", ln, sum0, sum1, u.min(), u.max(), ut.min(), ut.max());
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

void IProblem2HForward2D::add2Info(const DoubleMatrix &u, vector<ExtendedSpaceNode2DH> &info, unsigned int ln, unsigned int rows, unsigned int cols) const
{
    double hx = 0.01;
    double hy = 0.01;
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        ExtendedSpaceNode2DH &ui = info[j];
        ui.u[ln] = u[ui.j][ui.i];
        //ui.ux[ln] = (u[ui.j][ui.i+1] - u[ui.j][ui.i-1])/(2.0*hx);
        //ui.uy[ln] = (u[ui.j+1][ui.i] - u[ui.j-1][ui.i])/(2.0*hy);
        ui.ux[ln] = (u[ui.j][ui.i-2]-8.0*u[ui.j][ui.i-1]+8.0*u[ui.j][ui.i+1]-u[ui.j][ui.i+2])/(12.0*hx);
        ui.uy[ln] = (u[ui.j-2][ui.i]-8.0*u[ui.j-1][ui.i]+8.0*u[ui.j+1][ui.i]-u[ui.j+2][ui.i])/(12.0*hy);
    }
    //    for (unsigned int j=0; j<mEquParameter.No; j++)
    //    {
    //        ExtendedSpaceNode2DH &ui = info[j];
    //        for (unsigned int r=0; r<rows; r++)
    //        {
    //            for (unsigned int c=0; c<cols; c++)
    //            {
    //                unsigned int x = ui.wi[r][c].i;
    //                unsigned int y = ui.wi[r][c].j;
    //                ui.wi[r][c].u[ln] = u[y][x];
    //            }
    //        }
    //    }
}

void IProblem2HForward2D::calculateMVD_N(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const
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
    for (unsigned int j=0; j<No; j++) IProblem2H2D::distributeDelta(mOptParameter.xi[j], obsPointNodes, j, dimX, dimY, 4);

    std::vector<ExtendedSpacePointNode> cntDeltaNodes;
    for (unsigned int i=0; i<Nc; i++) IProblem2H2D::distributeDelta(mOptParameter.eta[i], cntDeltaNodes, i, dimX, dimY, 4);

    std::vector<ExtendedSpacePointNode> qPointNodes;
    for (unsigned int s=0; s<Ns; s++) IProblem2H2D::distributeDelta(mEquParameter.theta[s], qPointNodes, s, dimX, dimY, 4);

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
            e.extendWeights(dimX, dimY, L+1, _INFO_ROWS_, _INFO_COLS_);
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
        }
    }
    if (use == true) add2Info(u00, info, 0, _INFO_ROWS_, _INFO_COLS_);
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

            u05[m][n] = u00[m][n] + initial2(sn)*ht*0.5 + sum*ht*ht*0.125;
            u10[m][n] = u00[m][n] + initial2(sn)*ht     + sum*ht*ht*0.500;
        }
    }
    qPointNodes.clear();

    if (use == true) add2Info(u10, info, 1, _INFO_ROWS_, _INFO_COLS_);
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

    //TimeNodePDE tn;

    for (unsigned int l=2; l<=L; l++)
    {
        if (l == 2)
        {
            for (unsigned int m=0; m<=M; m++)
            {
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;
                    for (unsigned int si=0; si<qPointNodes.size(); si++)
                    {
                        const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                        if (qNode.i == n && qNode.j == m)
                        {
                            u10[m][n] += mEquParameter.q[qNode.id] * qNode.w * (1.0/ht);
                        }
                    }
                }
            }
        }

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

                if (m == 0)     d1X[n] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                if (m>0 && m<M) d1X[n] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                if (m == M)     d1X[n] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                d1X[n] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                d1X[n] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

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

        if (rows2.size() == 0)
        {
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (m == 0)     d1X[n] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    if (m>0 && m<M) d1X[n] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    if (m == M)     d1X[n] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1X[n] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                    d1X[n] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

                    if (n==0)
                    {
                        a1X[0] = 0.0;
                        b1X[0] = p_aa_htht__hxhx___lambda_ht;
                        c1X[0] = m_aa_htht__hxhx;
                    }
                    else if(n==N)
                    {
                        a1X[N] = m_aa_htht__hxhx;
                        b1X[N] = p_aa_htht__hxhx___lambda_ht;
                        c1X[N] = 0.0;
                    }
                    else
                    {
                        a1X[n] = m_aa_htht__hxhx_h;
                        b1X[n] = p_aa_htht__hxhx___lambda_ht;
                        c1X[n] = m_aa_htht__hxhx_h;
                    }

                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                                d1X[n] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1X[n] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
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
            double* a1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
            double* b1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
            double* c1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
            double* d1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
            double* x1 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
            DoubleMatrix w2(rows1.size()*(N+1), rows1.size()*(N+1), 0.0);

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (m==0)       d1[offset+n] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    if (m>0 && m<M) d1[offset+n] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    if (m==M)       d1[offset+n] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1[offset+n] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                    d1[offset+n] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);

                    if (n == 0)
                    {
                        a1[offset+0] = 0.0;
                        b1[offset+0] = p_aa_htht__hxhx___lambda_ht;
                        c1[offset+0] = m_aa_htht__hxhx;
                    }
                    else if (n == N)
                    {
                        a1[offset+N] = m_aa_htht__hxhx;
                        b1[offset+N] = p_aa_htht__hxhx___lambda_ht;
                        c1[offset+N] = 0.0;
                    }
                    else
                    {
                        a1[offset+n] = m_aa_htht__hxhx_h;
                        b1[offset+n] = p_aa_htht__hxhx___lambda_ht;
                        c1[offset+n] = m_aa_htht__hxhx_h;
                    }

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(odj);

                                bool found = false;
                                for (unsigned int cs=0; cs<rows1.size(); cs++)
                                {
                                    if (opn.j == rows1[cs])
                                    {
                                        found = true;
                                        w2[offset+n][cs*(N+1)+opn.i] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                    }
                                }

                                if (!found)
                                {
                                    d1[offset+n] += htht_h * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[offset+n] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }

                        }
                    }
                }
                offset += N+1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w2.data(), x1, rows1.size()*(N+1));

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=0; n<=N; n++)
                {
                    u15[m][n] = x1[offset+n];
                }
                offset += N+1;
            }

            w2.clear();
            free(x1);
            free(d1);
            free(c1);
            free(b1);
            free(a1);
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

                if (n==0)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                if (n>0 && n<N) d1Y[m] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                if (n==N)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                d1Y[m] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                d1Y[m] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

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

        if (cols2.size() == 0)
        {
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    if (n==0)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    if (n>0 && n<N) d1Y[m] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    if (n==N)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d1Y[m] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                    d1Y[m] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

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
                    }
                    else
                    {
                        a1Y[m] = m_aa_htht__hyhy_h;
                        b1Y[m] = p_aa_htht__hyhy___lambda_ht;
                        c1Y[m] = m_aa_htht__hyhy_h;
                    }

                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);
                                d1Y[m] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1Y[m] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
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
            double* a2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
            double* b2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
            double* c2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
            double* d2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
            double* x2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
            DoubleMatrix w2(cols1.size()*(M+1), cols1.size()*(M+1), 0.0);

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    if (n==0)       d2[offset+m] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    if (n>0 && n<N) d2[offset+m] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    if (n==N)       d2[offset+m] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d2[offset+m] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                    d2[offset+m] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);

                    if (m == 0)
                    {
                        a2[offset+0] = 0.0;
                        b2[offset+0] = p_aa_htht__hyhy___lambda_ht;
                        c2[offset+0] = m_aa_htht__hyhy;
                    }
                    else if (m == M)
                    {
                        a2[offset+M] = m_aa_htht__hyhy;
                        b2[offset+M] = p_aa_htht__hyhy___lambda_ht;
                        c2[offset+M] = 0.0;
                    }
                    else
                    {
                        a2[offset+m] = m_aa_htht__hyhy_h;
                        b2[offset+m] = p_aa_htht__hyhy___lambda_ht;
                        c2[offset+m] = m_aa_htht__hyhy_h;
                    }

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (opn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+m][cs*(M+1)+opn.j] -= htht_h * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+m] += htht_h * mOptParameter.k[cdn.id][opn.id] * u[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[offset+m] -= htht_h * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }
                offset += M+1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M+1));

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=0; m<=M; m++)
                {
                    u[m][n] = x2[offset+m];
                }
                offset += M+1;
            }

            w2.clear();
            free(x2);
            free(d2);
            free(c2);
            free(b2);
            free(a2);
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

        if (use == true) add2Info(u, info, l, _INFO_ROWS_, _INFO_COLS_);
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

void IProblem2HForward2D::calculateMVD_N(DoubleMatrix &u, DoubleMatrix &ut) const
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
    //double htht_h = 0.5*ht*ht;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    u.clear(); u.resize(M+1, N+1);

    SpaceNodePDE sn;

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
    layerInfo(u00, 0);
    layerInfo(u00, u00, 0);

    TimeNodePDE tn; tn.i = 0; tn.t = 0.0;

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
            sum += f(sn, tn);

            u05[m][n] = u00[m][n] + initial2(sn)*ht*0.5 + sum*ht*ht*0.125;
            u10[m][n] = u00[m][n] + initial2(sn)*ht     + sum*ht*ht*0.500;
        }
    }
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

    for (unsigned int l=2; l<=L; l++)
    {
        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int m=0; m<=M; m++)
        {
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;

                if (m == 0)     d1X[n] = p_aa_htht__hyhy_h*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                if (m>0 && m<M) d1X[n] = p_aa_htht__hyhy_h*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                if (m == M)     d1X[n] = p_aa_htht__hyhy_h*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                d1X[n] += 0.5*(u10[m][n]-u00[m][n]) + u10[m][n];
                d1X[n] += 0.5*lambda*ht*(4.0*u10[m][n]-u05[m][n]);
                //d1X[n] += htht_h * (f(sn, tn)+2.0*lambda*tn.t);

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
                    //d1X[n] += (a*a*ht*ht)/hx * 2.0;
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

        //--------------------------------------------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht-0.5*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            for (unsigned int m=0; m<=M; m++)
            {
                sn.j = m; sn.y = m*hy;

                if (n==0)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                if (n>0 && n<N) d1Y[m] = p_aa_htht__hxhx_h*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                if (n==N)       d1Y[m] = p_aa_htht__hxhx_h*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                d1Y[m] += 0.5*(u10[m][n]-u00[m][n]) + u15[m][n];
                d1Y[m] += 0.5*lambda*ht*(4.0*u15[m][n]-u10[m][n]);
                //d1Y[m] += htht_h * (f(sn, tn)+2.0*lambda*tn.t);

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
                    //d1Y[M] += (a*a*ht*ht)/hy * 2.0;
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

    u00.clear();
    u05.clear();
    u10.clear();
    u15.clear();
}

void IProblem2HForward2D::extendPoints(std::vector<ExtendedSpacePointNode> &obsPointNodes, std::vector<ExtendedSpacePointNode> &cntDeltaNodes,
                                       std::vector<ExtendedSpacePointNode> &qPointNodes, const Dimension &dimX, const Dimension &dimY) const
{
    unsigned int No = mEquParameter.No;
    unsigned int Nc = mEquParameter.Nc;
    unsigned int Ns = mEquParameter.Ns;

    for (unsigned int j=0; j<No; j++) IProblem2H2D::distributeDelta(mOptParameter.xi[j], obsPointNodes, j, dimX, dimY, 4);
    for (unsigned int i=0; i<Nc; i++) IProblem2H2D::distributeDelta(mOptParameter.eta[i], cntDeltaNodes, i, dimX, dimY, 4);
    for (unsigned int s=0; s<Ns; s++) IProblem2H2D::distributeDelta(mEquParameter.theta[s], qPointNodes, s, dimX, dimY, 4);
}

void IProblem2HForward2D::findRowsCols(std::vector<unsigned int> &rows0, std::vector<unsigned int> &rows1, std::vector<unsigned int> &rows2,
                                       std::vector<unsigned int> &cols0, std::vector<unsigned int> &cols1, std::vector<unsigned int> &cols2,
                                       std::vector<ExtendedSpacePointNode> &obsPointNodes, std::vector<ExtendedSpacePointNode> &cntDeltaNodes,
                                       unsigned int N, unsigned int M) const
{
    for (unsigned int m=1; m<=M-1; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);

            if (cdn.j == m)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.j == m)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
    }

    for (unsigned int n=1; n<=N-1; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);
            if (cdn.i == n)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.i == n)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
    }
}

void IProblem2HForward2D::initialLayers(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10, vector<ExtendedSpaceNode2DH> &info, bool use,
                                        std::vector<ExtendedSpacePointNode> &obsPointNodes, std::vector<ExtendedSpacePointNode> &cntDeltaNodes,
                                        std::vector<ExtendedSpacePointNode> &qPointNodes, unsigned int N, unsigned int M,
                                        double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const
{
    unsigned int No = mEquParameter.No;
    //unsigned int Nc = mEquParameter.Nc;

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u00[m][n] = initial1(sn);
        }
    }
    if (use == true) add2Info(u00, info, 0, _INFO_ROWS_, _INFO_COLS_);
    layerInfo(u00, 0);
    layerInfo(u00, u00, 0);

    double hh = 0.5*ht;
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
                double kz = 0.0;
                const ExtendedSpacePointNode &cntNode = cntDeltaNodes.at(cdn);
                if (cntNode.i == n && cntNode.j == m)
                {
                    for (unsigned int j=0; j<No; j++)
                    {
                        kz += mOptParameter.k[cntNode.id][j]*mOptParameter.z[cntNode.id][j];
                    }
                    sum -= kz*cntNode.w;
                }
            }

            for (unsigned int cdn=0; cdn<cntDeltaNodes.size(); cdn++)
            {
                double sm2 = 0.0;
                const ExtendedSpacePointNode &cntNode = cntDeltaNodes.at(cdn);
                if (cntNode.i == n && cntNode.j == m)
                {
                    for (unsigned int opn=0; opn<obsPointNodes.size(); opn++)
                    {
                        const ExtendedSpacePointNode &obsNode = obsPointNodes.at(opn);
                        if (obsNode.i == n && obsNode.j == m)
                        {
                            sm2 += mOptParameter.k[cntNode.id][obsNode.id] * (u00[m][n]*(obsNode.w*hx*hy));
                        }
                    }
                    sum += sm2*cntNode.w;
                }
            }

            for (unsigned int si=0; si<qPointNodes.size(); si++)
            {
                const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                if (qNode.i == n && qNode.j == m)
                {
                    sum += mEquParameter.q[qNode.id] * qNode.w * (1.0/ht);
                }
            }

            u05[m][n] = u00[m][n] + hh*initial2(sn);// + hh*hh*sum;
            u10[m][n] = u00[m][n] + ht*initial2(sn) + 0.5*ht*ht*sum;
        }
    }
    qPointNodes.clear();

    if (use == true) add2Info(u10, info, 1, _INFO_ROWS_, _INFO_COLS_);
    layerInfo(u10, u10, 1);
    layerInfo(u10, 1);
}

void IProblem2HForward2D::mallocVectors(double *&a1X, double *&b1X, double *&c1X, double *&d1X, double *&x1X, unsigned int sizeX,
                                        double *&a1Y, double *&b1Y, double *&c1Y, double *&d1Y, double *&x1Y, unsigned int sizeY) const
{
    a1X = (double *) malloc(sizeof(double)*sizeX);
    b1X = (double *) malloc(sizeof(double)*sizeX);
    c1X = (double *) malloc(sizeof(double)*sizeX);
    d1X = (double *) malloc(sizeof(double)*sizeX);
    x1X = (double *) malloc(sizeof(double)*sizeX);

    a1Y = (double *) malloc(sizeof(double)*sizeY);
    b1Y = (double *) malloc(sizeof(double)*sizeY);
    c1Y = (double *) malloc(sizeof(double)*sizeY);
    d1Y = (double *) malloc(sizeof(double)*sizeY);
    x1Y = (double *) malloc(sizeof(double)*sizeY);
}

void IProblem2HForward2D::initBorders(unsigned int N, unsigned int M, double hx, double hy, double ht, unsigned int l,
                                      DoubleMatrix &u15, DoubleMatrix &u) const
{
    SpaceNodePDE sn0;
    SpaceNodePDE sn1;
    TimeNodePDE tn;

    tn.i = l; tn.t = tn.i*ht;

    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = hx*N;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy;
        sn1.j = m; sn1.y = m*hy;

        u15[m][0] = boundary(sn0, tn);
        u15[m][N] = boundary(sn1, tn);

        u[m][0] = boundary(sn0, tn);
        u[m][N] = boundary(sn1, tn);
    }

    sn0.j = 0; sn0.y = 0.0;
    sn1.j = M; sn1.y = hy*M;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx;
        sn1.i = n; sn1.x = n*hx;

        u15[0][n] = boundary(sn0, tn);
        u15[M][n] = boundary(sn1, tn);

        u[0][n] = boundary(sn0, tn);
        u[M][n] = boundary(sn1, tn);
    }
}

void IProblem2HForward2D::freeVectors(double *a1X, double *b1X, double *c1X, double *d1X, double *x1X,
                                      double *a1Y, double *b1Y, double *c1Y, double *d1Y, double *x1Y) const
{
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
}

void IProblem2HForward2D::prepareInfo(unsigned int N, std::vector<SpacePoint> points, std::vector<ExtendedSpaceNode2DH> &info,
                                      unsigned int L, const Dimension &dimX, const Dimension &dimY) const
{
    info.resize(N);
    for (unsigned int j=0; j<N; j++)
    {
        ExtendedSpaceNode2DH &e = info[j];
        e.id = j;
        e.setSpaceNode(points[j]);
        e.extendWeights(dimX, dimY, L+1, _INFO_ROWS_, _INFO_COLS_);
    }
}
