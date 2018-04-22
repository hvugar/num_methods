#include "iproblem2hbackward2d.h"
#include <imaging.h>

using namespace IProblem2H;

void IProblem2HBackward2D::calculateMVD(DoubleMatrix &p, vector<ExtendedSpaceNode2DH> &info, bool use, const vector<ExtendedSpaceNode2DH> &u_info) const
{
    calculateMVD_D(p, info, use, u_info);
}

void IProblem2HBackward2D::calculateMVD_D(DoubleMatrix &p, vector<ExtendedSpaceNode2DH> &info, bool use, const vector<ExtendedSpaceNode2DH> &u_info) const
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

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p05(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    p.clear(); p.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<ExtendedSpacePointNode> obsDeltaNodes;
    for (unsigned int j=0; j<No; j++) IProblem2H2D::distributeDelta(mOptParameter.xi[j], obsDeltaNodes, j, dimX, dimY, 4);

    std::vector<ExtendedSpacePointNode> cntPointNodes;
    for (unsigned int i=0; i<Nc; i++) IProblem2H2D::distributeDelta(mOptParameter.eta[i], cntPointNodes, i, dimX, dimY, 4);

    SpaceNodePDE sn;

    //--------------------------------------------------------------------------------------------//
    vector<unsigned int> rows0, rows1, rows2;
    for (unsigned int m=1; m<=M-1; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int j=0; j<obsDeltaNodes.size(); j++)
        {
            const ExtendedSpacePointNode &odn = obsDeltaNodes.at(j);
            if (odn.j == m)
            {
                found1 = true;
                for (unsigned int i=0; i<cntPointNodes.size(); i++)
                {
                    const ExtendedSpacePointNode &cpn = cntPointNodes.at(i);
                    if (cpn.j == m)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
    }

    //--------------------------------------------------------------------------------------------//
    vector<unsigned int> cols0, cols1, cols2;
    for (unsigned int n=1; n<=N-1; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int j=0; j<obsDeltaNodes.size(); j++)
        {
            const ExtendedSpacePointNode &odn = obsDeltaNodes.at(j);
            if (odn.i == n)
            {
                found1 = true;
                for (unsigned int i=0; i<cntPointNodes.size(); i++)
                {
                    const ExtendedSpacePointNode &cpn = cntPointNodes.at(i);
                    if (cpn.i == n)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
    }
    //----------------------------------------------------------------------------------------------//

    //-------------------------------------------- info --------------------------------------------//
    if (use == true)
    {
        info.resize(Nc);
        for (unsigned int i=0; i<Nc; i++)
        {
            ExtendedSpaceNode2DH &e = info[i];
            e.id = i;
            e.setSpaceNode(mOptParameter.eta[i]);
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
            p00[m][n] = initial1(sn);
        }
    }
    if (use == true) add2Info(p00, info, L, _INFO_ROWS_, _INFO_COLS_);
    layerInfo(p00, L);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;

            if (n==0)      sum += aa__hxhx*(p00[m][n]-2.0*p00[m][n+1]+p00[m][n+2]);
            else if (n==N) sum += aa__hxhx*(p00[m][n-2]-2.0*p00[m][n-1]+p00[m][n]);
            else           sum += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);

            if (m==0)      sum += aa__hyhy*(p00[m][n]-2.0*p00[m+1][n]+p00[m+2][n]);
            else if (m==M) sum += aa__hyhy*(p00[m-2][n]-2.0*p00[m-1][n]+p00[m][n]);
            else           sum += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);

            sum += lambda*initial2(sn);

            for (unsigned int odn=0; odn<obsDeltaNodes.size(); odn++)
            {
                const ExtendedSpacePointNode &obsNode = obsDeltaNodes.at(odn);
                if (obsNode.i == n && obsNode.j == m)
                {
                    for (unsigned int cpn=0; cpn<cntPointNodes.size(); cpn++)
                    {
                        const ExtendedSpacePointNode &cntNode = cntPointNodes.at(cpn);
                        if (cntNode.i == n && cntNode.j == m)
                        {
                            sum += mOptParameter.k[cntNode.id][obsNode.id] * p00[m][n]*(cntNode.w*(hx*hy)) * obsNode.w;
                        }
                    }
                }
            }

            p05[m][n] = p00[m][n] - initial2(sn)*ht*0.5;// + 0.125*ht*ht*sum;
            p10[m][n] = p00[m][n] - initial2(sn)*ht    ;// + 0.500*ht*ht*sum;
        }
    }

    if (use == true) add2Info(p10, info, L-1, _INFO_ROWS_, _INFO_COLS_);
    layerInfo(p10, L-1);
    //------------------------------------- initial conditions -------------------------------------//

    double *a1X = (double *) malloc(sizeof(double)*(N-1));
    double *b1X = (double *) malloc(sizeof(double)*(N-1));
    double *c1X = (double *) malloc(sizeof(double)*(N-1));
    double *d1X = (double *) malloc(sizeof(double)*(N-1));
    double *x1X = (double *) malloc(sizeof(double)*(N-1));

    double *a1Y = (double *) malloc(sizeof(double)*(M-1));
    double *b1Y = (double *) malloc(sizeof(double)*(M-1));
    double *c1Y = (double *) malloc(sizeof(double)*(M-1));
    double *d1Y = (double *) malloc(sizeof(double)*(M-1));
    double *x1Y = (double *) malloc(sizeof(double)*(M-1));

    TimeNodePDE tn;

    for (unsigned int l1=2; l1<=L; l1++)
    {
        unsigned int l = L-l1;

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht;

        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy;
            sn1.j = m; sn1.y = m*hy;

            p15[m][0] = boundary(sn0, tn);
            p15[m][N] = boundary(sn1, tn);

            p[m][0] = boundary(sn0, tn);
            p[m][N] = boundary(sn1, tn);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx;
            sn1.i = n; sn1.x = n*hx;

            p15[0][n] = boundary(sn0, tn);
            p15[M][n] = boundary(sn1, tn);

            p[0][n] = boundary(sn0, tn);
            p[M][n] = boundary(sn1, tn);
        }

        //--------------------------------------------------------------------------//
        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = n; sn.x = n*hx;

                if (m == 0)     d1X[n-1] = p_aa_htht__hyhy_h*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                if (m>0 && m<M) d1X[n-1] = p_aa_htht__hyhy_h*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                if (m == M)     d1X[n-1] = p_aa_htht__hyhy_h*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                d1X[n-1] += 0.5*(p10[m][n]-p00[m][n]) + p10[m][n];
                d1X[n-1] += 0.5*lambda*ht*(4.0*p10[m][n]-p05[m][n]);

                a1X[n-1] = m_aa_htht__hxhx_h;
                b1X[n-1] = p_aa_htht__hxhx___lambda_ht;
                c1X[n-1] = m_aa_htht__hxhx_h;
            }

            a1X[0]   = 0.0;
            c1X[N-2] = 0.0;

            d1X[0]   -= m_aa_htht__hxhx_h * p15[m][0];
            d1X[N-2] -= m_aa_htht__hxhx_h * p15[m][N];

            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N-1);
            for (unsigned int n=1; n<=N-1; n++) p15[m][n] = x1X[n-1];
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

                    if (m == 0)     d1X[n-1] = p_aa_htht__hyhy_h*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    if (m>0 && m<M) d1X[n-1] = p_aa_htht__hyhy_h*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    if (m == M)     d1X[n-1] = p_aa_htht__hyhy_h*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    d1X[n-1] += 0.5*(p10[m][n]-p00[m][n]) + p10[m][n];
                    d1X[n-1] += 0.5*lambda*ht*(4.0*p10[m][n]-p05[m][n]);

                    a1X[n-1] = m_aa_htht__hxhx_h;
                    b1X[n-1] = p_aa_htht__hxhx___lambda_ht;
                    c1X[n-1] = m_aa_htht__hxhx_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];
                                d1X[n-1] += htht_h * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1X[n-1] += 2.0*ifunc->r * htht_h * mOptParameter.k[i][odn.id] * ifunc->gpi(i, l, u_info, mOptParameter)*sgn(ifunc->g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1X[0]   = 0.0;
                c1X[N-2] = 0.0;

                d1X[0]   -= m_aa_htht__hxhx_h * p15[m][0];
                d1X[N-2] -= m_aa_htht__hxhx_h * p15[m][N];

                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N-1);
                for (unsigned int n=1; n<=N-1; n++) p15[m][n] = x1X[n-1];
            }
        }

        if (rows2.size() != 0)
        {
            double* a1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            double* b1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            double* c1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            double* d1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            double* x1 = (double*) malloc(sizeof(double)*rows1.size()*(N-1));
            DoubleMatrix w1(rows1.size()*(N-1), rows1.size()*(N-1), 0.0);

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (m == 0)     d1[offset+(n-1)] = p_aa_htht__hyhy_h*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    if (m>0 && m<M) d1[offset+(n-1)] = p_aa_htht__hyhy_h*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    if (m == M)     d1[offset+(n-1)] = p_aa_htht__hyhy_h*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    d1[offset+(n-1)] += 0.5*(p10[m][n]-p00[m][n]) + p10[m][n];
                    d1[offset+(n-1)] += 0.5*lambda*ht*(4.0*p10[m][n]-p05[m][n]);

                    a1[offset+(n-1)] = m_aa_htht__hxhx_h;
                    b1[offset+(n-1)] = p_aa_htht__hxhx___lambda_ht;
                    c1[offset+(n-1)] = m_aa_htht__hxhx_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (cpn.j == rows1[rs])
                                    {
                                        found = true;
                                        w1[offset+(n-1)][rs*(N-1)+(cpn.i-1)] -= htht_h * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                    }
                                }

                                if (!found)
                                {
                                    d1[offset+(n-1)] += htht_h * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1[offset+(n-1)] += 2.0*ifunc->r * ht *  mOptParameter.k[i][odn.id] * ifunc->gpi(i, l, u_info, mOptParameter)*sgn(ifunc->g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx_h * p15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx_h * p15[m][N];

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w1.data(), x1, rows1.size()*(N-1));

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=1; n<=N-1; n++)
                {
                    p15[m][n] = x1[offset+(n-1)];
                }
                offset += N-1;
            }

            w1.clear();
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

                if (n==0)       d1Y[m-1] = p_aa_htht__hxhx_h*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                if (n>0 && n<N) d1Y[m-1] = p_aa_htht__hxhx_h*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                if (n==N)       d1Y[m-1] = p_aa_htht__hxhx_h*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                d1Y[m-1] += 0.5*(p10[m][n]-p00[m][n]) + p15[m][n];
                d1Y[m-1] += 0.5*lambda*ht*(4.0*p15[m][n]-p10[m][n]);

                a1Y[m-1] = m_aa_htht__hyhy_h;
                b1Y[m-1] = p_aa_htht__hyhy___lambda_ht;
                c1Y[m-1] = m_aa_htht__hyhy_h;
            }

            a1Y[0]   = 0.0;
            c1Y[M-2] = 0.0;

            d1Y[0]   -= m_aa_htht__hyhy_h * p[0][n];
            d1Y[M-2] -= m_aa_htht__hyhy_h * p[M][n];

            tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M-1);
            for (unsigned int m=1; m<=M-1; m++) p[m][n] = x1Y[m-1];
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

                    if (n==0)       d1Y[m-1] = p_aa_htht__hxhx_h*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    if (n>0 && n<N) d1Y[m-1] = p_aa_htht__hxhx_h*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    if (n==N)       d1Y[m-1] = p_aa_htht__hxhx_h*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    d1Y[m-1] += 0.5*(p10[m][n]-p00[m][n]) + p15[m][n];
                    d1Y[m-1] += 0.5*lambda*ht*(4.0*p15[m][n]-p10[m][n]);

                    a1Y[m-1] = m_aa_htht__hyhy_h;
                    b1Y[m-1] = p_aa_htht__hyhy___lambda_ht;
                    c1Y[m-1] = m_aa_htht__hyhy_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes.at(cni);
                                d1Y[m-1] += htht_h * mOptParameter.k[cpn.id][odn.id] * p[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                            }

                            //                            for (unsigned int i=0; i<Nc; i++)
                            //                            {
                            //                                d1Y[m-1] += 2.0*ifunc->r * htht_h * mOptParameter.k[i][odn.id] * ifunc->gpi(i, l, u_info, mOptParameter)*sgn(ifunc->g0i(i, l, u_info, mOptParameter)) * odn.w;
                            //                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1Y[0]   = 0.0;
                c1Y[M-2] = 0.0;

                d1Y[0]   -= m_aa_htht__hyhy_h * p[0][n];
                d1Y[M-2] -= m_aa_htht__hyhy_h * p[M][n];

                tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M-1);
                for (unsigned int m=1; m<=M-1; m++) p[m][n] = x1Y[m-1];
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

                    if (n==0)       d2[offset+(m-1)] = p_aa_htht__hxhx_h*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    if (n>0 && n<N) d2[offset+(m-1)] = p_aa_htht__hxhx_h*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    if (n==N)       d2[offset+(m-1)] = p_aa_htht__hxhx_h*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    d2[offset+(m-1)] += 0.5*(p10[m][n]-p00[m][n]) + p15[m][n];
                    d2[offset+(m-1)] += 0.5*lambda*ht*(4.0*p15[m][n]-p10[m][n]);

                    a2[offset+(m-1)] = m_aa_htht__hyhy_h;
                    b2[offset+(m-1)] = p_aa_htht__hyhy___lambda_ht;
                    c2[offset+(m-1)] = m_aa_htht__hyhy_h;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (cpn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+(m-1)][cs*(M-1)+(cpn.j-1)] -= htht_h * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+(m-1)] += htht_h * mOptParameter.k[cpn.id][odn.id] * p[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d2[offset+(m-1)] += 2.0*ifunc->r * htht_h *  mOptParameter.k[i][odn.id] * ifunc->gpi(i, l, u_info, mOptParameter)*sgn(ifunc->g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy_h * p[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy_h * p[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M-1));

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    p[m][n] = x2[offset+(m-1)];
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

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p[m][n];
                p05[m][n] = p15[m][n];
            }
        }

        if (use == true) add2Info(p, info, l, _INFO_ROWS_, _INFO_COLS_);
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
    //QPixmap px;
    //visualizeMatrixHeat(u, u.min(), u.max(), px);
    //char buffer[30] = {0};
    //int c = sprintf(buffer, "images/b/%4d.png", layerNumber);
    //buffer[c] = 0;
    //px.save(buffer, "png", 0);
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

void IProblem2HBackward2D::add2Info(const DoubleMatrix &p, vector<ExtendedSpaceNode2DH> &info, unsigned int ln, unsigned int rows, unsigned int cols) const
{
    for (unsigned int i=0; i<mEquParameter.Nc; i++)
    {
        ExtendedSpaceNode2DH &pi = info[i];
        for (unsigned int r=0; r<rows; r++)
        {
            for (unsigned int c=0; c<cols; c++)
            {
                unsigned int x = pi.wi[r][c].i;
                unsigned int y = pi.wi[r][c].j;
                pi.wi[r][c].u[ln] = p[y][x];
            }
        }
    }
}

void IProblem2HBackward2D::calculateMVD_N(DoubleMatrix &p, vector<ExtendedSpaceNode2DH> &info, bool use, const vector<ExtendedSpaceNode2DH> &u_info) const
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
    //unsigned int Ns = mEquParameter.Ns;

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

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p05(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    p.clear(); p.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<ExtendedSpacePointNode> obsDeltaNodes;
    for (unsigned int j=0; j<No; j++) IProblem2H2D::distributeDelta(mOptParameter.xi[j], obsDeltaNodes, j, dimX, dimY, 4);

    std::vector<ExtendedSpacePointNode> cntPointNodes;
    for (unsigned int i=0; i<Nc; i++) IProblem2H2D::distributeDelta(mOptParameter.eta[i], cntPointNodes, i, dimX, dimY, 4);

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
            const ExtendedSpacePointNode &odn = obsDeltaNodes.at(j);
            if (odn.j == ny)
            {
                found1 = true;
                for (unsigned int i=0; i<cntPointNodes.size(); i++)
                {
                    const ExtendedSpacePointNode &cpn = cntPointNodes.at(i);
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
            const ExtendedSpacePointNode &odn = obsDeltaNodes.at(j);
            if (odn.i == nx)
            {
                found1 = true;
                for (unsigned int i=0; i<cntPointNodes.size(); i++)
                {
                    const ExtendedSpacePointNode &cpn = cntPointNodes.at(i);
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
        info.resize(Nc);
        for (unsigned int i=0; i<Nc; i++)
        {
            ExtendedSpaceNode2DH &e = info[i];
            e.id = i;
            e.setSpaceNode(mOptParameter.eta[i]);
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
            p00[m][n] = initial1(sn);
        }
    }
    if (use == true) add2Info(p00, info, L, _INFO_ROWS_, _INFO_COLS_);
    layerInfo(p00, L);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;

            if (n==0)      sum += aa__hxhx*(p00[m][n]-2.0*p00[m][n+1]+p00[m][n+2]);
            else if (n==N) sum += aa__hxhx*(p00[m][n-2]-2.0*p00[m][n-1]+p00[m][n]);
            else           sum += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);

            if (m==0)      sum += aa__hyhy*(p00[m][n]-2.0*p00[m+1][n]+p00[m+2][n]);
            else if (m==M) sum += aa__hyhy*(p00[m-2][n]-2.0*p00[m-1][n]+p00[m][n]);
            else           sum += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);

            sum += lambda*initial2(sn);

            for (unsigned int odn=0; odn<obsDeltaNodes.size(); odn++)
            {
                const ExtendedSpacePointNode &obsNode = obsDeltaNodes.at(odn);
                if (obsNode.i == n && obsNode.j == m)
                {
                    for (unsigned int cpn=0; cpn<cntPointNodes.size(); cpn++)
                    {
                        const ExtendedSpacePointNode &cntNode = cntPointNodes.at(cpn);
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

    if (use == true) add2Info(p10, info, L-1, _INFO_ROWS_, _INFO_COLS_);
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

    TimeNodePDE tn;

    for (unsigned int l1=2; l1<=L; l1++)
    {
        unsigned int l = L-l1;

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht;

        //--------------------------------------------------------------------------//
        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;

                if (m == 0)     d1X[n] = p_aa_htht__hyhy_h*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                if (m>0 && m<M) d1X[n] = p_aa_htht__hyhy_h*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                if (m == M)     d1X[n] = p_aa_htht__hyhy_h*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                d1X[n] += 0.5*(p10[m][n]-p00[m][n]) + p10[m][n];
                d1X[n] += 0.5*lambda*ht*(4.0*p10[m][n]-p05[m][n]);

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

                    if (m == 0)     d1X[n] = p_aa_htht__hyhy_h*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    if (m>0 && m<M) d1X[n] = p_aa_htht__hyhy_h*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    if (m == M)     d1X[n] = p_aa_htht__hyhy_h*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    d1X[n] += 0.5*(p10[m][n]-p00[m][n]) + p10[m][n];
                    d1X[n] += 0.5*lambda*ht*(4.0*p10[m][n]-p05[m][n]);

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

                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];
                                d1X[n] += htht_h * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1X[n] += 2.0*ifunc->r * htht_h * mOptParameter.k[i][odn.id] * ifunc->gpi(i, l, u_info, mOptParameter)*sgn(ifunc->g0i(i, l, u_info, mOptParameter)) * odn.w;
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

                    if (m == 0)     d1[offset+n] = p_aa_htht__hyhy_h*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    if (m>0 && m<M) d1[offset+n] = p_aa_htht__hyhy_h*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    if (m == M)     d1[offset+n] = p_aa_htht__hyhy_h*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    d1[offset+n] += 0.5*(p10[m][n]-p00[m][n]) + p10[m][n];
                    d1[offset+n] += 0.5*lambda*ht*(4.0*p10[m][n]-p05[m][n]);

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
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];

                                bool found = false;
                                for (unsigned int cs=0; cs<rows1.size(); cs++)
                                {
                                    if (cpn.j == rows1[cs])
                                    {
                                        found = true;
                                        w2[offset+n][cs*(N+1)+cpn.i] -= htht_h * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                    }
                                }

                                if (!found)
                                {
                                    d1[offset+n] += htht_h * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1[offset+n] += 2.0*ifunc->r * ht *  mOptParameter.k[i][odn.id] * ifunc->gpi(i, l, u_info, mOptParameter)*sgn(ifunc->g0i(i, l, u_info, mOptParameter)) * odn.w;
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
                    p15[m][n] = x1[offset+n];
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

                if (n==0)       d1Y[m] = p_aa_htht__hxhx_h*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                if (n>0 && n<N) d1Y[m] = p_aa_htht__hxhx_h*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                if (n==N)       d1Y[m] = p_aa_htht__hxhx_h*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                d1Y[m] += 0.5*(p10[m][n]-p00[m][n]) + p15[m][n];
                d1Y[m] += 0.5*lambda*ht*(4.0*p15[m][n]-p10[m][n]);

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

                    if (n==0)       d1Y[m] = p_aa_htht__hxhx_h*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    if (n>0 && n<N) d1Y[m] = p_aa_htht__hxhx_h*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    if (n==N)       d1Y[m] = p_aa_htht__hxhx_h*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    d1Y[m] += 0.5*(p10[m][n]-p00[m][n]) + p15[m][n];
                    d1Y[m] += 0.5*lambda*ht*(4.0*p15[m][n]-p10[m][n]);

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

                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes.at(cni);
                                d1Y[m] += htht_h * mOptParameter.k[cpn.id][odn.id] * p[cpn.j][cpn.i]* (cpn.w * (hx*hy)) * odn.w;
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1Y[m] += 2.0*ifunc->r * htht_h * mOptParameter.k[i][odn.id] * ifunc->gpi(i, l, u_info, mOptParameter)*sgn(ifunc->g0i(i, l, u_info, mOptParameter)) * odn.w;
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

                    if (n==0)       d2[offset+m] = p_aa_htht__hxhx_h*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    if (n>0 && n<N) d2[offset+m] = p_aa_htht__hxhx_h*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    if (n==N)       d2[offset+m] = p_aa_htht__hxhx_h*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    d2[offset+m] += 0.5*(p10[m][n]-p00[m][n]) + p15[m][n];
                    d2[offset+m] += 0.5*lambda*ht*(4.0*p15[m][n]-p10[m][n]);

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
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (cpn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+m][cs*(M+1)+cpn.j] -= htht_h * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+m] += htht_h * mOptParameter.k[cpn.id][odn.id] * p[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d2[offset+m] += 2.0*ifunc->r * htht_h *  mOptParameter.k[i][odn.id] * ifunc->gpi(i, l, u_info, mOptParameter)*sgn(ifunc->g0i(i, l, u_info, mOptParameter)) * odn.w;
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
                    p[m][n] = x2[offset+m];
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

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p[m][n];
                p05[m][n] = p15[m][n];
            }
        }

        if (use == true) add2Info(p, info, l, _INFO_ROWS_, _INFO_COLS_);
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
