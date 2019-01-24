#include "problem2h_solver1.h"

Problem2HDirichlet1::Problem2HDirichlet1() : Problem2HDirichletBase () {}

Problem2HDirichlet1::~Problem2HDirichlet1()
{}

double Problem2HDirichlet1::integral(const std::vector<DoubleMatrix> &vu) const
{
    const double ht = timeDimension().step();
    double sum = 0.0;
    sum += 0.5*integralU(vu[0]);
    for (unsigned int ln=2; ln<=2*(LD-1); ln+=2)
    {
        sum += integralU(vu[ln]);
    }
    sum += 0.5*integralU(vu[2*LD]);
    return sum*ht;
}

auto Problem2HDirichlet1::penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const -> double
{
    const double ht = mtimeDimension.step();
    const unsigned int L = static_cast<const unsigned int> ( mtimeDimension.size() );

    double pnlt = 0.0;
    for (unsigned int i=0; i<mEquParameter.Nc; i++)
    {
        double pnlt_i = 0.0;
        double _gpi_0 = gpi(i, 0, info, o_prm);
        pnlt_i += 0.5*_gpi_0*_gpi_0;
        for (unsigned int l=1; l<=L+LD-1; l++)
        {
            double _gpi_l = gpi(i, 2*l, info, o_prm);
            pnlt_i += _gpi_l*_gpi_l;
        }
        double _gpi_L = gpi(i, 2*(L+LD), info, o_prm);
        pnlt_i += 0.5*_gpi_L*_gpi_L;

        pnlt += pnlt_i*ht;
    }

    return pnlt;
}

auto Problem2HDirichlet1::gradient(const DoubleVector &pv, DoubleVector &g) const -> void
{
    const double ht        = mtimeDimension.step();
    const unsigned int Nc  = mEquParameter.Nc;
    const unsigned int No  = mEquParameter.No;

#if defined(DISCRETE_DELTA_TIME)
    const unsigned int Nt  = mEquParameter.Nt;
#else
    const unsigned int L   = static_cast<unsigned int>(mtimeDimension.size());
    const unsigned int LLD = L + LD;
#endif
    Problem2HDirichlet1* prob = const_cast<Problem2HDirichlet1*>(this);

    g.clear();
    g.resize(pv.length(), 0.0);

    OptimizeParameterH o_prm;
    VectorToPrm(pv, o_prm);
    prob->mOptParameter = o_prm;

    const DoubleVector &Q1 = mEquParameter.Q1;
    const DoubleVector &Q2 = mEquParameter.Q2;

    for (unsigned int q1=0; q1<Q1.length(); q1++)
    {
        prob->mEquParameter.pulses[0].q = Q1[q1];
        for (unsigned int q2=0; q2<Q2.length(); q2++)
        {
            prob->mEquParameter.pulses[1].q = Q2[q2];

            std::vector<DoubleMatrix> u;
            spif_vectorH u_info;
            solveForwardIBVP(u, u_info, true, pv);
            spif_vectorH p_info;
            solveBackwardIBVP(u, p_info, true, u_info, pv);

            unsigned int gi = 0;

            // k
            if (optimizeK)
            {
#if defined(DISCRETE_DELTA_TIME)
                for (unsigned int s=0; s<Nt; s++)
                {
                    unsigned int ln = 2*static_cast<unsigned int>(mEquParameter.timeMoments[s]/ht);
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const SpacePointInfoH &pi = p_info[i];
                        for (unsigned int j=0; j<No; j++)
                        {
                            const SpacePointInfoH &uj = u_info[j];
                            double zij = o_prm.z[s][i][j];
                            double grad_Kij = 0.0;
                            grad_Kij += -(uj.vl[ln] - zij) * pi.vl[ln];
                            //grad_Kij += -(uj.vl[ln] - zij) * 2.0*r*gpi(i,ln,u_info,o_prm)*sgn(g0i(i,ln,u_info,o_prm));
                            //grad_Kij += +2.0*regEpsilon*(o_prm.k[s][i][j] - mRegParameter.k[s][i][j]);
                            g[gi++] += grad_Kij * (1.0/(double(Q1.length())*double(Q2.length())));
                        }
                    }
                }
#else
                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfoH &pi = p_info[i];

                    for (unsigned int j=0; j<No; j++)
                    {
                        const SpacePointInfoH &uj = u_info[j];

                        double grad_Kij = 0.0;
                        double zij = o_prm.z[i][j];

                        grad_Kij += 0.5 * (pi.vl[0] + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * (uj.vl[0] - zij);
                        for (unsigned int ln=1; ln<=LLD-1; ln++)
                        {
                            grad_Kij += (pi.vl[2*ln] + 2.0*r*gpi(i,2*ln,u_info,o_prm)*sgn(g0i(i,2*ln,u_info,o_prm))) * (uj.vl[2*ln] - zij);
                        }
                        grad_Kij += 0.5 * (pi.vl[2*LLD] + 2.0*r*gpi(i,2*LLD,u_info,o_prm)*sgn(g0i(i,2*LLD,u_info,o_prm))) * (uj.vl[2*LLD] - zij);

                        grad_Kij *= -ht;

                        g[gi++] = grad_Kij + 2.0*regEpsilon*(o_prm.k[i][j] - mRegParameter.k[i][j]);
                    }
                }
#endif
            }
            else
            {
                for (unsigned int i=0; i<Nc; i++)
                {
                    for (unsigned int j=0; j<No; j++)
                    {
                        g[gi++] = 0.0;
                    }
                }
            }

            // z
            if (optimizeZ)
            {
#if defined(DISCRETE_DELTA_TIME)
                for (unsigned int s=0; s<Nt; s++)
                {
                    unsigned int ln = 2*static_cast<unsigned int>(mEquParameter.timeMoments[s]/ht);
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const SpacePointInfoH &pi = p_info[i];
                        for (unsigned int j=0; j<No; j++)
                        {
                            double kij = o_prm.k[s][i][j];
                            double grad_Zij = 0.0;
                            grad_Zij += kij * pi.vl[ln];
                            //grad_Zij += kij * 2.0*r*gpi(i,ln,u_info,o_prm)*sgn(g0i(i,ln,u_info,o_prm));
                            //grad_Zij += +2.0*regEpsilon*(o_prm.z[s][i][j] - mRegParameter.z[s][i][j]);
                            g[gi++] += grad_Zij * (1.0/(double(Q1.length())*double(Q2.length())));
                        }
                    }
                }
#else
                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfoH &pi = p_info[i];

                    for (unsigned int j=0; j<No; j++)
                    {
                        double grad_Zij = 0.0;
                        double kij = o_prm.k[i][j];

                        grad_Zij += 0.5 * (pi.vl[0] + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * kij;
                        for (unsigned int ln=1; ln<=LLD-1; ln++)
                        {
                            grad_Zij += (pi.vl[2*ln] + 2.0*r*gpi(i,2*ln,u_info,o_prm)*sgn(g0i(i,2*ln,u_info,o_prm))) * kij;
                        }
                        grad_Zij += 0.5 * (pi.vl[2*LLD] + 2.0*r*gpi(i,2*LLD,u_info,o_prm)*sgn(g0i(i,2*LLD,u_info,o_prm))) * kij;
                        grad_Zij *= ht;

                        g[gi++] = grad_Zij + 2.0*regEpsilon*(o_prm.z[i][j] - mRegParameter.z[i][j]);
                    }
                }
#endif
            }
            else
            {
                for (unsigned int i=0; i<Nc; i++)
                {
                    for (unsigned int j=0; j<No; j++)
                    {
                        g[gi++] = 0.0;
                    }
                }
            }

            // xi
            if (optimizeO)
            {
#if defined(DISCRETE_DELTA_TIME)
                for (unsigned int j=0; j<No; j++)
                {
                    const SpacePointInfoH &uj = u_info[j];

                    double gradXijX = 0.0;
                    double gradXijY = 0.0;

                    for (unsigned int s=0; s<Nt; s++)
                    {
                        unsigned int ln = 2*static_cast<unsigned int>(mEquParameter.timeMoments[s]/ht);

                        for (unsigned int i=0; i<Nc; i++)
                        {
                            gradXijX += -o_prm.k[s][i][j] * uj.dx[ln] * p_info[i].vl[ln];
                            gradXijY += -o_prm.k[s][i][j] * uj.dy[ln] * p_info[i].vl[ln];
                            //gradXijX += -o_prm.k[s][i][j] * uj.dx[ln] * 2.0*r*gpi(i,ln,u_info,o_prm)*sgn(g0i(i,ln,u_info,o_prm));
                            //gradXijY += -o_prm.k[s][i][j] * uj.dy[ln] * 2.0*r*gpi(i,ln,u_info,o_prm)*sgn(g0i(i,ln,u_info,o_prm));
                        }
                    }

                    //gradXijX += 2.0*regEpsilon*(o_prm.xi[j].x - mRegParameter.xi[j].x);
                    //gradXijY += 2.0*regEpsilon*(o_prm.xi[j].y - mRegParameter.xi[j].y);

                    g[gi++] += gradXijX * (1.0/(double(Q1.length())*double(Q2.length())));
                    g[gi++] += gradXijY * (1.0/(double(Q1.length())*double(Q2.length())));
                }
#else
                for (unsigned int j=0; j<No; j++)
                {
                    const SpacePointInfoH &uj = u_info[j];

                    double gradXijX = 0.0;
                    double gradXijY = 0.0;
                    double vi = 0.0;

                    vi = 0.0;
                    for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j] * (p_info[i].vl[0] + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm)));
                    gradXijX += 0.5 * uj.dx[0] * vi;
                    gradXijY += 0.5 * uj.dy[0] * vi;

                    for (unsigned int ln=1; ln<=LLD-1; ln++)
                    {
                        vi = 0.0;
                        for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].vl[2*ln] + 2.0*r*gpi(i,2*ln,u_info,o_prm)*sgn(g0i(i,2*ln,u_info,o_prm)));
                        gradXijX += uj.dx[2*ln] * vi;
                        gradXijY += uj.dy[2*ln] * vi;
                    }

                    vi = 0.0;
                    for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].vl[2*LLD] + 2.0*r*gpi(i,2*LLD,u_info,o_prm)*sgn(g0i(i,2*LLD,u_info,o_prm)));
                    gradXijX += 0.5 * uj.dx[2*LLD] * vi;
                    gradXijY += 0.5 * uj.dy[2*LLD] * vi;

                    gradXijX *= -ht;
                    gradXijY *= -ht;

                    g[gi++] = gradXijX + 2.0*regEpsilon*(o_prm.xi[j].x - mRegParameter.xi[j].x);
                    g[gi++] = gradXijY + 2.0*regEpsilon*(o_prm.xi[j].y - mRegParameter.xi[j].y);
                }
#endif
            }
            else
            {
                for (unsigned int j=0; j<No; j++)
                {
                    g[gi++] = 0.0;
                    g[gi++] = 0.0;
                }
            }

            // eta
            if (optimizeC)
            {
#if defined(DISCRETE_DELTA_TIME)
                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfoH &pi = p_info[i];

                    double gradEtaiX = 0.0;
                    double gradEtaiY = 0.0;

                    for (unsigned int s=0; s<Nt; s++)
                    {
                        unsigned int ln = 2*static_cast<unsigned int>(mEquParameter.timeMoments[s]/ht);

                        for (unsigned int j=0; j<No; j++)
                        {
                            gradEtaiX += -pi.dx[ln] * o_prm.k[s][i][j] * (u_info[j].vl[ln] - o_prm.z[s][i][j]);
                            gradEtaiY += -pi.dy[ln] * o_prm.k[s][i][j] * (u_info[j].vl[ln] - o_prm.z[s][i][j]);
                        }
                    }

                    //gradEtaiX += 2.0*regEpsilon*(o_prm.eta[i].x - mRegParameter.eta[i].x);
                    //gradEtaiY += 2.0*regEpsilon*(o_prm.eta[i].y - mRegParameter.eta[i].y);

                    g[gi++] += gradEtaiX * (1.0/(double(Q1.length())*double(Q2.length())));
                    g[gi++] += gradEtaiY * (1.0/(double(Q1.length())*double(Q2.length())));
                }
#else
                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfoH &pi = p_info[i];

                    double gradEtaiX = 0.0;
                    double gradEtaiY = 0.0;
                    double vi = 0.0;

                    vi = 0.0;
                    for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].vl[0] - o_prm.z[i][j]);
                    gradEtaiX += 0.5 * pi.dx[0] * vi;
                    gradEtaiY += 0.5 * pi.dy[0] * vi;

                    for (unsigned int ln=1; ln<=LLD-1; ln++)
                    {
                        vi = 0.0;
                        for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].vl[2*ln] - o_prm.z[i][j]);
                        gradEtaiX += pi.dx[2*ln] * vi;
                        gradEtaiY += pi.dy[2*ln] * vi;
                    }

                    vi = 0.0;
                    for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].vl[2*LLD] - o_prm.z[i][j]);
                    gradEtaiX += 0.5 * pi.dx[2*LLD] * vi;
                    gradEtaiY += 0.5 * pi.dy[2*LLD] * vi;

                    gradEtaiX *= -ht;
                    gradEtaiY *= -ht;

                    g[gi++] = gradEtaiX + 2.0*regEpsilon*(o_prm.eta[i].x - mRegParameter.eta[i].x);
                    g[gi++] = gradEtaiY + 2.0*regEpsilon*(o_prm.eta[i].y - mRegParameter.eta[i].y);
                }
#endif
            }
            else
            {
                for (unsigned int i=0; i<Nc; i++)
                {
                    g[gi++] = 0.0;
                    g[gi++] = 0.0;
                }
            }

            for (unsigned int i=0; i<u_info.size(); i++) u_info[i].clear(); u_info.clear();
            for (unsigned int i=0; i<p_info.size(); i++) p_info[i].clear(); p_info.clear();
            for (unsigned int i=0; i<u.size(); i++) u[i].clear(); u.clear();
        }
    }
}

auto Problem2HDirichlet1::solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use, const DoubleVector &pv, double lambda) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );
    const unsigned int L = static_cast<const unsigned int> ( time.size() );
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a        = mEquParameter.a;
    const double alpha    = mEquParameter.alpha;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;

    const double ht_ht_025 = ht*ht*0.25;
    const double alpha_ht_025 = alpha*ht*0.25;

    const double m_aa_htht__hxhx_025_lambda = -(0.25*a*a)*((ht*ht)/(hx*hx))*lambda;
    const double b_aa_htht__hxhx = +(1.0 + 0.5*(a*a)*((ht*ht)/(hx*hx))*lambda + alpha_ht_025);
    const double p_aa_htht__hyhy_025 = +(0.25*a*a)*((ht*ht)/(hy*hy));
    const double p_aa_htht__hxhx_025_lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*lambda;
    const double p_aa_htht__hxhx_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*(1.0-2.0*lambda);

    const double m_aa_htht__hyhy_025_lambda = -(0.25*a*a)*((ht*ht)/(hy*hy))*lambda;
    const double b_aa_htht__hyhy = +(1.0 + 0.5*(a*a)*((ht*ht)/(hy*hy))*lambda + alpha_ht_025);
    const double p_aa_htht__hxhx_025 = +(0.25*a*a)*((ht*ht)/(hx*hx));
    const double p_aa_htht__hyhy_025_lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*lambda;
    const double p_aa_htht__hyhy_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*(1.0-2.0*lambda);

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    for (unsigned int ln=0; ln<u.size(); ln++) u[ln].clear(); u.clear();
    unsigned int u_size = 2*LD + 1;
    u.resize(u_size); for (unsigned int ln=0; ln<u_size; ln++) u[ln].resize(M+1, N+1);

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    for (unsigned int n=1; n<=N-1; n++)
    {
        ax[n-1] = m_aa_htht__hxhx_025_lambda;
        bx[n-1] = b_aa_htht__hxhx;
        cx[n-1] = m_aa_htht__hxhx_025_lambda;
    }
    ax[0] = cx[N-2] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    for (unsigned int m=1; m<=M-1; m++)
    {
        ay[m-1] = m_aa_htht__hyhy_025_lambda;
        by[m-1] = b_aa_htht__hyhy;
        cy[m-1] = m_aa_htht__hyhy_025_lambda;
    }
    ay[0] = cy[M-2] = 0.0;

    //----------------------------------------------------------------------------------------------//
    std::vector<DeltaGrid2D> measuremntGirdList(No);
    std::vector<DeltaGrid2D> cntrlDeltaGridList(Nc);
    for (unsigned int j=0; j<No; j++)
    {
        measuremntGirdList[j].initGrid(N, hx, M, hy);
        measuremntGirdList[j].distributeGauss(mOptParameter.xi[j], 1, 1);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        cntrlDeltaGridList[i].initGrid(N, hx, M, hy);
        cntrlDeltaGridList[i].distributeGauss(mOptParameter.eta[i], 1, 1);
    }
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) prepareInfo(No,  mOptParameter.xi, u_info, 2*LLD+1);
    //----------------------------------------------------------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    initPulseWeightMatrix(mEquParameter.pulses);
    /************************************************************************/
    SpaceNodePDE sn00;
    for (unsigned int m=0; m<=M; m++)
    {
        sn00.j = static_cast<int>(m); sn00.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn00.i = static_cast<int>(n); sn00.x = n*hx;
            u00[m][n] = f_initial1(sn00);
        }
    }
    if (use == true) add2Info(u00, u_info, 0, hx, hy, measuremntGirdList); f_layerInfo(u00, 0);
    /************************************************************************/
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*ht;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = 1.0*ht;
    SpaceNodePDE sn05, sn10;
    sn05.i = 0; sn05.x = 0.0; sn10.i = static_cast<int>(N); sn10.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn05.j = static_cast<int>(m); sn05.y = m*hy; u10[m][0] = f_boundary(sn05, tn10); u05[m][0] = f_boundary(sn00, tn05);
        sn10.j = static_cast<int>(m); sn10.y = m*hy; u10[m][N] = f_boundary(sn10, tn10); u05[m][N] = f_boundary(sn10, tn05);
    }
    sn05.j = 0; sn05.y = 0.0; sn10.j = static_cast<int>(M); sn10.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn05.i = static_cast<int>(n); sn05.x = n*hx; u10[0][n] = f_boundary(sn05, tn10); u05[0][n] = f_boundary(sn05, tn05);
        sn10.i = static_cast<int>(n); sn10.x = n*hx; u10[M][n] = f_boundary(sn10, tn10); u05[M][n] = f_boundary(sn10, tn05);
    }
    for (unsigned int m=1; m<=M-1; m++)
    {
        sn10.j = static_cast<int>(m); sn10.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn10.i = static_cast<int>(n); sn10.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= alpha*(f_initial2(sn10));

            u05[m][n] = u00[m][n] + 0.5*ht*f_initial2(sn10) + (0.125*ht*ht) * sum;
            u10[m][n] = u00[m][n] + 1.0*ht*f_initial2(sn10) + (0.500*ht*ht) * sum;
        }
    }
    if (use == true) add2Info(u05, u_info, 1, hx, hy, measuremntGirdList); f_layerInfo(u05, 1);
    if (use == true) add2Info(u10, u_info, 2, hx, hy, measuremntGirdList); f_layerInfo(u10, 2);
    /************************************************************************/
    //------------------------------------- initial conditions -------------------------------------//

    SpaceNodePDE sn;
    for (unsigned int ln=2; ln<=LLD; ln++)
    {
        TimeNodePDE tn05; tn05.i = ln-1; tn05.t = tn05.i*ht-0.5*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;
        TimeNodePDE tn15; tn15.i = ln;   tn15.t = tn15.i*ht-0.5*ht;
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;

        /**************************************************** border conditions ***************************************************/
        SpaceNodePDE sn0, sn1;
        sn0.i = static_cast<int>(0); sn0.x = 0*hx; sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = static_cast<int>(m); sn0.y = m*hy; u15[m][0] = f_boundary(sn0, tn15); u20[m][0] = f_boundary(sn0, tn20);
            sn1.j = static_cast<int>(m); sn1.y = m*hy; u15[m][N] = f_boundary(sn1, tn15); u20[m][N] = f_boundary(sn1, tn20);
        }
        sn0.j = static_cast<int>(0); sn0.y = 0*hy; sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = static_cast<int>(n); sn0.x = n*hx; u15[0][n] = f_boundary(sn0, tn15); u20[0][n] = f_boundary(sn0, tn20);
            sn1.i = static_cast<int>(n); sn1.x = n*hx; u15[M][n] = f_boundary(sn1, tn15); u20[M][n] = f_boundary(sn1, tn20);
        }
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        currentLayerFGrid(u10, cntrlDeltaGridList, measuremntGirdList, 2*ln-2);
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1] = 0.0;
                dx[n-1] = (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025;
                dx[n-1] += 2.0*u10[m][n] - u05[m][n] + alpha_ht_025*u05[m][n];
                dx[n-1] += (u10[m][n-1] - 2.0*u10[m][n] + u10[m][n+1])*p_aa_htht__hxhx_025_1m2lambda;
                dx[n-1] += (u05[m][n-1] - 2.0*u05[m][n] + u05[m][n+1])*p_aa_htht__hxhx_025_lambda;
                dx[n-1] += ht_ht_025 * mCrFfxWeightMatrix[m][n];
            }
            dx[0]   -= u15[m][0]*m_aa_htht__hxhx_025_lambda;
            dx[N-2] -= u15[m][N]*m_aa_htht__hxhx_025_lambda;
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }
        if (use == true) add2Info(u15, u_info, 2*ln-1, hx, hy, measuremntGirdList); f_layerInfo(u15, 2*ln-1);
        /**************************************************** x direction apprx ***************************************************/
        /**************************************************** y direction apprx ***************************************************/
        currentLayerFGrid(u15, cntrlDeltaGridList, measuremntGirdList, 2*ln-1);
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = static_cast<int>(m); sn.y = m*hy;
                dy[m-1] = 0.0;
                dy[m-1] += (u15[m][n-1] - 2.0*u15[m][n] + u15[m][n+1])*p_aa_htht__hxhx_025;
                dy[m-1] += 2.0*u15[m][n] - u10[m][n] + alpha_ht_025*u10[m][n];
                dy[m-1] += (u15[m-1][n] - 2.0*u15[m][n] + u15[m+1][n])*p_aa_htht__hyhy_025_1m2lambda;
                dy[m-1] += (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025_lambda;
                dy[m-1] += ht_ht_025 * mCrFfxWeightMatrix[m][n];
            }
            dy[0]   -= u20[0][n]*m_aa_htht__hyhy_025_lambda;
            dy[M-2] -= u20[M][n]*m_aa_htht__hyhy_025_lambda;
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
        }
        if (use == true) add2Info(u20, u_info, 2*ln+0, hx, hy, measuremntGirdList); f_layerInfo(u20, 2*ln+0);
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u05[m][n] = u15[m][n];
                u10[m][n] = u20[m][n];
            }
        }
        /**************************************************** saving last LD layers ***********************************************/
        if ( L == ln ) u[0] = u20;
        if ( L+1 <= ln && ln <= LLD ) { u[2*(ln-L)-1] = u15; u[2*(ln-L)+0] = u20; }
        /**************************************************** saving last LD layers ***********************************************/
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

    for (unsigned int j=0; j<No; j++) measuremntGirdList[j].cleanGrid(); measuremntGirdList.clear();
    for (unsigned int i=0; i<Nc; i++) cntrlDeltaGridList[i].cleanGrid(); cntrlDeltaGridList.clear();

    u00.clear();
    u05.clear();
    u10.clear();
    u15.clear();
    u20.clear();
}



auto Problem2HDirichlet1::solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info, const DoubleVector &pv, double lambda) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<const unsigned int>( dimX.size() );
    const unsigned int M = static_cast<const unsigned int>( dimY.size() );
    const unsigned int L = static_cast<const unsigned int>( time.size() );
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a        = mEquParameter.a;
    const double alpha    = mEquParameter.alpha;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;

    const double ht_ht_025 = ht*ht*0.25;
    const double alpha_ht_025 = alpha*ht*0.25;

    const double m_aa_htht__hxhx_025_lambda = -(0.25*a*a)*((ht*ht)/(hx*hx))*lambda;
    const double b_aa_htht__hxhx = +(1.0 + 0.5*(a*a)*((ht*ht)/(hx*hx))*lambda + alpha_ht_025);
    const double p_aa_htht__hyhy_025 = +(0.25*a*a)*((ht*ht)/(hy*hy));
    const double p_aa_htht__hxhx_025_lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*lambda;
    const double p_aa_htht__hxhx_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*(1.0-2.0*lambda);

    const double m_aa_htht__hyhy_025_lambda = -(0.25*a*a)*((ht*ht)/(hy*hy))*lambda;
    const double b_aa_htht__hyhy = +(1.0 + 0.5*(a*a)*((ht*ht)/(hy*hy))*lambda + alpha_ht_025);
    const double p_aa_htht__hxhx_025 = +(0.25*a*a)*((ht*ht)/(hx*hx));
    const double p_aa_htht__hyhy_025_lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*lambda;
    const double p_aa_htht__hyhy_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*(1.0-2.0*lambda);

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p05(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    double *ax = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *bx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *cx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *dx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *rx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    for (unsigned int n=1; n<=N-1; n++)
    {
        ax[n-1] = m_aa_htht__hxhx_025_lambda;
        bx[n-1] = b_aa_htht__hxhx;
        cx[n-1] = m_aa_htht__hxhx_025_lambda;
    }
    ax[0] = cx[N-2] = 0.0;

    double *ay = static_cast<double*>( malloc(sizeof(double)*(M-1)) );
    double *by = static_cast<double*>( malloc(sizeof(double)*(M-1)) );
    double *cy = static_cast<double*>( malloc(sizeof(double)*(M-1)) );
    double *dy = static_cast<double*>( malloc(sizeof(double)*(M-1)) );
    double *ry = static_cast<double*>( malloc(sizeof(double)*(M-1)) );
    for (unsigned int m=1; m<=M-1; m++)
    {
        ay[m-1] = m_aa_htht__hyhy_025_lambda;
        by[m-1] = b_aa_htht__hyhy;
        cy[m-1] = m_aa_htht__hyhy_025_lambda;
    }
    ay[0] = cy[M-2] = 0.0;

    OptimizeParameterH mOptParameter;
    VectorToPrm(pv, mOptParameter);

    //--------------------------------------------------------------------------------------------//
    std::vector<DeltaGrid2D> measuremntGirdList(No);
    std::vector<DeltaGrid2D> cntrlDeltaGridList(Nc);
    for (unsigned int j=0; j<No; j++)
    {
        measuremntGirdList[j].initGrid(N, hx, M, hy);
        measuremntGirdList[j].distributeGauss(mOptParameter.xi[j], 1, 1);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        cntrlDeltaGridList[i].initGrid(N, hx, M, hy);
        cntrlDeltaGridList[i].distributeGauss(mOptParameter.eta[i], 1, 1);
    }
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) prepareInfo(Nc,  mOptParameter.eta, p_info, 2*LLD+1);
    //-------------------------------------------- info --------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    /************************************************************************/
    SpaceNodePDE sn00;
    for (unsigned int m=0; m<=M; m++)
    {
        sn00.j = static_cast<int>(m); sn00.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn00.i = static_cast<int>(n); sn00.x = n*hx;
            p00[m][n] = b_initial1(sn00);
        }
    }
    if (use == true) add2Info(p00, p_info, 2*LLD, hx, hy, cntrlDeltaGridList); b_layerInfo(p00, 2*LLD);
    /************************************************************************/
    TimeNodePDE tn05; tn05.i = LLD-1; tn05.t = LLD*ht - 0.5*ht;
    TimeNodePDE tn10; tn10.i = LLD-1; tn10.t = LLD*ht - 1.0*ht;
    SpaceNodePDE sn05, sn10;
    sn05.i = 0; sn05.x = 0.0; sn10.i = static_cast<int>(N); sn10.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn05.j = static_cast<int>(m); sn05.y = m*hy; p10[m][0] = b_boundary(sn05, tn10); p05[m][0] = b_boundary(sn05, tn05);
        sn10.j = static_cast<int>(m); sn10.y = m*hy; p10[m][N] = b_boundary(sn10, tn10); p05[m][N] = b_boundary(sn10, tn05);
    }
    sn05.j = 0;  sn05.y = 0.0; sn10.j = static_cast<int>(M); sn10.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn05.i = static_cast<int>(n); sn05.x = n*hx; p10[0][n] = b_boundary(sn05, tn10); p05[0][n] = b_boundary(sn05, tn05);
        sn10.i = static_cast<int>(n); sn10.x = n*hx; p10[M][n] = b_boundary(sn10, tn10); p05[M][n] = b_boundary(sn10, tn05);
    }
    /************************************************************************/
    double *_w = new double[No];
    double* _p00 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p00[i] = 0.0;
    for (unsigned int i=0; i<Nc; i++)
    {
        const DeltaGrid2D &dg = cntrlDeltaGridList[i];
        for (unsigned int m=dg.minY(); m<=dg.maxY(); m++)
        {
            for (unsigned int n=dg.minX(); n<=dg.maxX(); n++)
            {
                _p00[i] += p00[m][n] * (dg.weight(n,m) * (hx*hy));
            }
        }
    }
    for (unsigned int j=0; j<No; j++)
    {
        _w[j] = 0.0;
        for (unsigned int i=0; i<Nc; i++)
        {
            //_w[j] += mOptParameter.k[s][i][j] * (_p00[i] + 2.0*r*gpi(i, 2*LLD, u_info, mOptParameter)*sgn(g0i(i,2*LLD,u_info,mOptParameter)));
        }
    }
    delete [] _p00;
    /************************************************************************/
    for (unsigned int m=1; m<=M-1; m++)
    {
        sn10.j = static_cast<int>(m); sn10.y = static_cast<int>(m)*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn10.i = static_cast<int>(n); sn10.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);
            sum += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);
            sum += alpha*b_initial2(sn10);
            for (unsigned int j=0; j<No; j++)
            {
                sum += _w[j] * measuremntGirdList[j].weight(n,m);
            }
            sum -= 2.0*mu(n,m)*(u.at(2*LD)[m][n]);

            p05[m][n] = p00[m][n] - (ht*0.5) * b_initial2(sn10) + 0.125*ht*ht*sum;
            p10[m][n] = p00[m][n] - (ht*1.0) * b_initial2(sn10) + 0.500*ht*ht*sum;
        }
    }
    if (use == true) add2Info(p05, p_info, 2*LLD-1, hx, hy, cntrlDeltaGridList); b_layerInfo(p05, 2*LLD-1);
    if (use == true) add2Info(p10, p_info, 2*LLD-2, hx, hy, cntrlDeltaGridList); b_layerInfo(p10, 2*LLD-2);
    /************************************************************************/
    delete [] _w;
    SpaceNodePDE sn;
    const unsigned int stop = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=LLD-2; ln != stop; ln--)
    {
        TimeNodePDE tn05; tn05.i = ln+1; tn05.t = tn05.i*ht+0.5*ht;
        TimeNodePDE tn10; tn10.i = ln+1; tn10.t = tn10.i*ht;
        TimeNodePDE tn15; tn15.i = ln;   tn15.t = tn15.i*ht+0.5*ht;
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;

        /**************************************************** border conditions ***************************************************/
        SpaceNodePDE sn0, sn1;
        sn0.i = static_cast<int>(0); sn0.x = 0.0; sn1.i = static_cast<int>(N); sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = static_cast<int>(m); sn0.y = m*hy; p15[m][0] = b_boundary(sn0, tn15); p20[m][0] = b_boundary(sn0, tn20);
            sn1.j = static_cast<int>(m); sn1.y = m*hy; p15[m][N] = b_boundary(sn1, tn15); p20[m][N] = b_boundary(sn1, tn20);
        }
        sn0.j = static_cast<int>(0); sn0.y = 0.0; sn1.j = static_cast<int>(M); sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = static_cast<int>(n); sn0.x = n*hx; p15[0][n] = b_boundary(sn0, tn15); p20[0][n] = b_boundary(sn0, tn20);
            sn1.i = static_cast<int>(n); sn1.x = n*hx; p15[M][n] = b_boundary(sn1, tn15); p20[M][n] = b_boundary(sn1, tn20);
        }
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        currentLayerBGrid(p10, cntrlDeltaGridList, measuremntGirdList, 2*ln+2, u_info);
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1] = 0.0;
                dx[n-1] = (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n]) * p_aa_htht__hyhy_025;
                dx[n-1] += 2.0*p10[m][n] - p05[m][n] + alpha_ht_025*p05[m][n];
                dx[n-1] += (p10[m][n-1] - 2.0*p10[m][n] + p10[m][n+1]) * p_aa_htht__hxhx_025_1m2lambda;
                dx[n-1] += (p05[m][n-1] - 2.0*p05[m][n] + p05[m][n+1]) * p_aa_htht__hxhx_025_lambda;
                dx[n-1] += ht_ht_025 * mCrBfxWeightMatrix[m][n];
                //------------------------------------- Adding functional part --------------------------------//
                if (L <= ln && ln <= LLD) dx[n-1] += -2.0*mu(n,m)*(u.at(2*(ln-L)+2)[m][n]) * ht_ht_025;
                //------------------------------------- Adding functional part --------------------------------//
            }
            dx[0]   -= p15[m][0] * m_aa_htht__hxhx_025_lambda;
            dx[N-2] -= p15[m][N] * m_aa_htht__hxhx_025_lambda;
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
        }
        if (use == true) add2Info(p15, p_info, 2*ln+1, hx, hy, cntrlDeltaGridList); b_layerInfo(p15, 2*ln+1);
        /**************************************************** x direction apprx ***************************************************/
        /**************************************************** y direction apprx ***************************************************/
        currentLayerBGrid(p15, cntrlDeltaGridList, measuremntGirdList, 2*ln+1, u_info);
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = sn.i = static_cast<int>(m); sn.y = m*hy;

                dy[m-1] = 0.0;
                dy[m-1] += (p15[m][n-1] - 2.0*p15[m][n] + p15[m][n+1]) * p_aa_htht__hxhx_025;
                dy[m-1] += 2.0*p15[m][n] - p10[m][n] + alpha_ht_025*p10[m][n];
                dy[m-1] += (p15[m-1][n] - 2.0*p15[m][n] + p15[m+1][n]) * p_aa_htht__hyhy_025_1m2lambda;
                dy[m-1] += (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n]) * p_aa_htht__hyhy_025_lambda;
                dy[m-1] += ht_ht_025 * mCrBfxWeightMatrix[m][n];
                //------------------------------------- Adding functional part --------------------------------//
                if (L <= ln && ln <= LLD) dy[m-1] += -2.0*mu(n,m)*(u.at(2*(ln-L)+1)[m][n]) * ht_ht_025;
                //------------------------------------- Adding functional part --------------------------------//
            }
            dy[0]   -= p20[0][n] * m_aa_htht__hyhy_025_lambda;
            dy[M-2] -= p20[M][n] * m_aa_htht__hyhy_025_lambda;
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
        }
        if (use == true) add2Info(p20, p_info, 2*ln+0, hx, hy, cntrlDeltaGridList); b_layerInfo(p20, 2*ln+0);
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p05[m][n] = p15[m][n];
                p10[m][n] = p20[m][n];
            }
        }
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

    p00.clear();
    p05.clear();
    p10.clear();
    p15.clear();
    p20.clear();
}


