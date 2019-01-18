#include "problem1h_solver1.h"

Problem1HDirichlet1::Problem1HDirichlet1() : Problem1HDirichletBase () {}

Problem1HDirichlet1::~Problem1HDirichlet1() {}

auto Problem1HDirichlet1::integral(const std::vector<DoubleVector> &vu) const -> double
{
    const double ht = timeDimension().step();
    double sum = 0.0;
    sum += 0.5*integralU(vu[0]);
    for (unsigned int ln=1; ln<=(LD-1); ln++)
    {
        sum += integralU(vu[ln]);
    }
    sum += 0.5*integralU(vu[LD]);
    return sum*ht;
}

auto Problem1HDirichlet1::penalty(const spif_vector1H &info, const OptimizeParameter1H &o_prm) const -> double
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
            double _gpi_l = gpi(i, l, info, o_prm);
            pnlt_i += _gpi_l*_gpi_l;
        }
        double _gpi_L = gpi(i, (L+LD), info, o_prm);
        pnlt_i += 0.5*_gpi_L*_gpi_L;

        pnlt += pnlt_i*ht;
    }

    return pnlt;
}

auto Problem1HDirichlet1::gradient(const DoubleVector &pv, DoubleVector &g) const -> void
{
    const unsigned int L   = static_cast<unsigned int>(mtimeDimension.size());
    const double ht        = mtimeDimension.step();
    const unsigned int Nc  = mEquParameter.Nc;
    const unsigned int No  = mEquParameter.No;
    const unsigned int LLD = L + LD;

    OptimizeParameter1H o_prm;
    VectorToPrm(pv, o_prm);

    Problem1HDirichlet1* prob = const_cast<Problem1HDirichlet1*>(this);
    prob->mOptParameter = o_prm;

    std::vector<DoubleVector> u;
    spif_vector1H u_info;
    solveForwardIBVP(u, u_info, true);
    spif_vector1H p_info;
    solveBackwardIBVP(u, p_info, true, u_info);

    g.clear();
    g.resize(pv.length(), 0.0);
    unsigned int gi = 0;

    // k
    if (optimizeK)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo1H &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointInfo1H &uj = u_info[j];

                double grad_Kij = 0.0;
                double zij = o_prm.z[i][j];

                grad_Kij += 0.5 * (pi.vl[0] + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * (uj.vl[0] - zij);
                for (unsigned int ln=1; ln<=LLD-1; ln++)
                {
                    grad_Kij += (pi.vl[ln] + 2.0*r*gpi(i,ln,u_info,o_prm)*sgn(g0i(i,ln,u_info,o_prm))) * (uj.vl[ln] - zij);
                }
                grad_Kij += 0.5 * (pi.vl[LLD] + 2.0*r*gpi(i, LLD,u_info,o_prm)*sgn(g0i(i,LLD,u_info,o_prm))) * (uj.vl[LLD] - zij);

                grad_Kij *= -ht;

                g[gi++] = grad_Kij + 2.0*regEpsilon*(o_prm.k[i][j] - mRegParameter.k[i][j]);
            }
        }
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
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo1H &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                double grad_Zij = 0.0;
                double kij = o_prm.k[i][j];

                grad_Zij += 0.5 * (pi.vl[0] + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * kij;
                for (unsigned int ln=1; ln<=LLD-1; ln++)
                {
                    grad_Zij += (pi.vl[ln] + 2.0*r*gpi(i,ln,u_info,o_prm)*sgn(g0i(i,ln,u_info,o_prm))) * kij;
                }
                grad_Zij += 0.5 * (pi.vl[LLD] + 2.0*r*gpi(i,LLD,u_info,o_prm)*sgn(g0i(i,LLD,u_info,o_prm))) * kij;
                grad_Zij *= ht;

                g[gi++] = grad_Zij + 2.0*regEpsilon*(o_prm.z[i][j] - mRegParameter.z[i][j]);
            }
        }
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
        for (unsigned int j=0; j<No; j++)
        {
            const SpacePointInfo1H &uj = u_info[j];

            double gradXijX = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j] * (p_info[i].vl[0] + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm)));
            gradXijX += 0.5 * uj.dx[0] * vi;

            for (unsigned int ln=1; ln<=LLD-1; ln++)
            {
                vi = 0.0;
                for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].vl[ln] + 2.0*r*gpi(i,ln,u_info,o_prm)*sgn(g0i(i,ln,u_info,o_prm)));
                gradXijX += uj.dx[ln] * vi;
            }

            vi = 0.0;
            for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].vl[LLD] + 2.0*r*gpi(i,LLD,u_info,o_prm)*sgn(g0i(i,LLD,u_info,o_prm)));
            gradXijX += 0.5 * uj.dx[LLD] * vi;

            gradXijX *= -ht;

            g[gi++] = gradXijX + 2.0*regEpsilon*(o_prm.ksi[j].x - mRegParameter.ksi[j].x);
        }
    }
    else
    {
        for (unsigned int j=0; j<No; j++)
        {
            g[gi++] = 0.0;
        }
    }

    // eta
    if (optimizeC)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo1H &pi = p_info[i];

            double gradEtaiX = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].vl[0] - o_prm.z[i][j]);
            gradEtaiX += 0.5 * pi.dx[0] * vi;

            for (unsigned int ln=1; ln<=LLD-1; ln++)
            {
                vi = 0.0;
                for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].vl[ln] - o_prm.z[i][j]);
                gradEtaiX += pi.dx[ln] * vi;
            }

            vi = 0.0;
            for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].vl[LLD] - o_prm.z[i][j]);
            gradEtaiX += 0.5 * pi.dx[LLD] * vi;

            gradEtaiX *= -ht;

            g[gi++] = gradEtaiX + 2.0*regEpsilon*(o_prm.eta[i].x - mRegParameter.eta[i].x);
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            g[gi++] = 0.0;
        }
    }

    for (unsigned int i=0; i<u_info.size(); i++)
    {
        u_info[i].clear();
    }

    for (unsigned int i=0; i<p_info.size(); i++)
    {
        p_info[i].clear();
    }

    u_info.clear();
    p_info.clear();
}

auto Problem1HDirichlet1::solveForwardIBVP(std::vector<DoubleVector> &u, spif_vector1H &u_info, bool use, double lambda) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int L = static_cast<const unsigned int> ( time.size() );
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double ht = time.step();

    const double a        = mEquParameter.a;
    const double alpha    = mEquParameter.alpha;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;

    const double ht_ht = ht*ht;
    const double alpha_ht_05 = alpha*ht*0.5;

    const double m_aa_htht__hxhx_lambda = -(a*a)*((ht*ht)/(hx*hx))*lambda;
    const double b_aa_htht__hxhx = +(1.0 + 2.0*(a*a)*((ht*ht)/(hx*hx))*lambda + alpha_ht_05);
    const double p_aa_htht__hxhx_1m2lambda = +(a*a)*((ht*ht)/(hx*hx))*(1.0-2.0*lambda);
    const double p_aa_htht__hxhx_lambda = +(a*a)*((ht*ht)/(hx*hx))*lambda;

    const double aa__hxhx = (a*a)/(hx*hx);

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);
    DoubleVector u20(N+1);

    for (unsigned int ln=0; ln<u.size(); ln++) u[ln].clear(); u.clear();
    unsigned int u_size = LD + 1;
    u.resize(u_size); for (unsigned int ln=0; ln<u_size; ln++) u[ln].resize(N+1);

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    for (unsigned int n=1; n<=N-1; n++)
    {
        ax[n-1] = m_aa_htht__hxhx_lambda;
        bx[n-1] = b_aa_htht__hxhx;
        cx[n-1] = m_aa_htht__hxhx_lambda;
    }
    ax[0] = cx[N-2] = 0.0;

    //----------------------------------------------------------------------------------------------//
    std::vector<DeltaGrid1D> measuremntGirdList(No);
    std::vector<DeltaGrid1D> cntrlDeltaGridList(Nc);
    for (unsigned int j=0; j<No; j++)
    {
        measuremntGirdList[j].initGrid(N, hx);
        measuremntGirdList[j].distributeGauss(mOptParameter.ksi[j], 1);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        cntrlDeltaGridList[i].initGrid(N, hx);
        cntrlDeltaGridList[i].distributeGauss(mOptParameter.eta[i], 1);
    }
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) prepareInfo(No, mOptParameter.ksi, u_info, LLD+1);
    //----------------------------------------------------------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    initPulseWeightVector(mEquParameter.theta);
    /************************************************************************/
    SpaceNodePDE sn00;
    for (unsigned int n=0; n<=N; n++)
    {
        sn00.i = static_cast<int>(n); sn00.x = n*hx;
        u00[n] = f_initial1(sn00);
    }
    if (use == true) add2Info(u00, u_info, 0, hx, measuremntGirdList); f_layerInfo(u00, 0);
    /************************************************************************/
    TimeNodePDE tn10; tn10.i = 1; tn10.t = 1.0*ht;
    SpaceNodePDE sn10;
    sn10.i = static_cast<int>(0); sn10.x = 0*hx; u10[0] = f_boundary(sn10, tn10);
    sn10.i = static_cast<int>(N); sn10.x = N*hx; u10[N] = f_boundary(sn10, tn10);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sn10.i = static_cast<int>(n); sn10.x = n*hx;

        double sum = 0.0;
        sum += aa__hxhx*(u00[n-1]-2.0*u00[n]+u00[n+1]);
        sum -= alpha*(f_initial2(sn10));
        u10[n] = u00[n] + ht*f_initial2(sn10) + (0.5*ht*ht) * sum;
    }
    if (use == true) add2Info(u10, u_info, 1, hx, measuremntGirdList); f_layerInfo(u10, 2);
    /************************************************************************/
    //------------------------------------- initial conditions -------------------------------------//

    SpaceNodePDE sn;
    for (unsigned int ln=2; ln<=LLD; ln++)
    {
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;

        /**************************************************** border conditions ***************************************************/
        SpaceNodePDE sn0;
        sn0.i = static_cast<int>(0); sn0.x = 0*hx; u20[0] = f_boundary(sn0, tn20);
        sn0.i = static_cast<int>(N); sn0.x = N*hx; u20[N] = f_boundary(sn0, tn20);
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        currentLayerFGrid(u10, cntrlDeltaGridList, measuremntGirdList);
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            dx[n-1] = 0.0;
            dx[n-1] += 2.0*u10[n] - u00[n] + alpha_ht_05*u00[n];
            dx[n-1] += (u10[n-1] - 2.0*u10[n] + u10[n+1])*p_aa_htht__hxhx_1m2lambda;
            dx[n-1] += (u00[n-1] - 2.0*u00[n] + u00[n+1])*p_aa_htht__hxhx_lambda;
            dx[n-1] += ht_ht * mCrFfxWeightMatrix[n];
        }
        dx[0]   -= u20[0]*m_aa_htht__hxhx_lambda;
        dx[N-2] -= u20[N]*m_aa_htht__hxhx_lambda;
        tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
        for (unsigned int n=1; n<=N-1; n++) u20[n] = rx[n-1];
        if (use == true) add2Info(u20, u_info, ln, hx, measuremntGirdList); f_layerInfo(u20, ln);
        /**************************************************** x direction apprx ***************************************************/

        for (unsigned int n=0; n<=N; n++)
        {
            u00[n] = u10[n];
            u10[n] = u20[n];
        }

        /**************************************************** saving last LD layers ***********************************************/
        if ( L == ln ) u[0] = u20;
        if ( L+1 <= ln && ln <= LLD ) { u[(ln-L)] = u20; }
        /**************************************************** saving last LD layers ***********************************************/
    }

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);

    for (unsigned int j=0; j<No; j++) measuremntGirdList[j].cleanGrid(); measuremntGirdList.clear();
    for (unsigned int i=0; i<Nc; i++) cntrlDeltaGridList[i].cleanGrid(); cntrlDeltaGridList.clear();

    u00.clear();
    u10.clear();
    u20.clear();
}



auto Problem1HDirichlet1::solveBackwardIBVP(const std::vector<DoubleVector> &u, spif_vector1H &p_info, bool use, const spif_vector1H &u_info, double lambda) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<const unsigned int>( dimX.size() );
    const unsigned int L = static_cast<const unsigned int>( time.size() );
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double ht = time.step();

    const double a        = mEquParameter.a;
    const double alpha    = mEquParameter.alpha;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;

    const double ht_ht = ht*ht;
    const double alpha_ht_05 = alpha*ht*0.5;

    const double m_aa_htht__hxhx_lambda = -(a*a)*((ht*ht)/(hx*hx))*lambda;
    const double b_aa_htht__hxhx = +(1.0 + 2.0*(a*a)*((ht*ht)/(hx*hx))*lambda + alpha_ht_05);
    const double p_aa_htht__hxhx_1m2lambda = +(a*a)*((ht*ht)/(hx*hx))*(1.0-2.0*lambda);
    const double p_aa_htht__hxhx_lambda = +(a*a)*((ht*ht)/(hx*hx))*lambda;

    const double aa__hxhx = (a*a)/(hx*hx);

    DoubleVector p00(N+1);
    DoubleVector p10(N+1);
    DoubleVector p20(N+1);

    double *ax = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *bx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *cx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *dx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *rx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    for (unsigned int n=1; n<=N-1; n++)
    {
        ax[n-1] = m_aa_htht__hxhx_lambda;
        bx[n-1] = b_aa_htht__hxhx;
        cx[n-1] = m_aa_htht__hxhx_lambda;
    }
    ax[0] = cx[N-2] = 0.0;

    //--------------------------------------------------------------------------------------------//
    std::vector<DeltaGrid1D> measuremntGirdList(No);
    std::vector<DeltaGrid1D> cntrlDeltaGridList(Nc);
    for (unsigned int j=0; j<No; j++)
    {
        measuremntGirdList[j].initGrid(N, hx);
        measuremntGirdList[j].distributeGauss(mOptParameter.ksi[j], 1);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        cntrlDeltaGridList[i].initGrid(N, hx);
        cntrlDeltaGridList[i].distributeGauss(mOptParameter.eta[i], 1);
    }
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) prepareInfo(Nc,  mOptParameter.eta, p_info, LLD+1);
    //-------------------------------------------- info --------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    /************************************************************************/
    SpaceNodePDE sn00;
    for (unsigned int n=0; n<=N; n++)
    {
        sn00.i = static_cast<int>(n); sn00.x = n*hx;
        p00[n] = b_initial1(sn00);
    }
    if (use == true) add2Info(p00, p_info, LLD, hx, cntrlDeltaGridList); b_layerInfo(p00, LLD);
    /************************************************************************/
    TimeNodePDE tn10; tn10.i = LLD-1; tn10.t = LLD*ht - 1.0*ht;
    SpaceNodePDE sn10;
    sn10.i = static_cast<int>(0); sn10.x = 0*hx; p10[0] = b_boundary(sn10, tn10);
    sn10.i = static_cast<int>(N); sn10.x = N*hx; p10[N] = b_boundary(sn10, tn10);
    /************************************************************************/
    double *_w = new double[No];
    double* _p00 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p00[i] = 0.0;
    for (unsigned int i=0; i<Nc; i++)
    {
        const DeltaGrid1D &dg = cntrlDeltaGridList[i];
        for (unsigned int n=dg.minX(); n<=dg.maxX(); n++)
        {
            _p00[i] += p00[n] * dg.weight(n) * hx;
        }
    }
    for (unsigned int j=0; j<No; j++)
    {
        _w[j] = 0.0;
        for (unsigned int i=0; i<Nc; i++)
        {
            _w[j] += mOptParameter.k[i][j] * (_p00[i] + 2.0*r*gpi(i, LLD, u_info, mOptParameter)*sgn(g0i(i,LLD,u_info,mOptParameter)));
        }
    }
    delete [] _p00;
    /************************************************************************/
    for (unsigned int n=1; n<=N-1; n++)
    {
        sn10.i = static_cast<int>(n); sn10.x = n*hx;

        double sum = 0.0;
        sum += aa__hxhx*(p00[n-1]-2.0*p00[n]+p00[n+1]);
        sum += alpha*b_initial2(sn10);
        for (unsigned int j=0; j<No; j++)
        {
            sum += _w[j] * measuremntGirdList[j].weight(n);
        }
        sum -= 2.0*mu(n)*(u.at(LD)[n]);
        p10[n] = p00[n] - ht*b_initial2(sn10) + 0.5*ht*ht*sum;
    }
    if (use == true) add2Info(p10, p_info, LLD-1, hx, cntrlDeltaGridList); b_layerInfo(p10, LLD-1);
    /************************************************************************/
    delete [] _w;
    SpaceNodePDE sn;
    const unsigned int stop = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=LLD-2; ln != stop; ln--)
    {
        TimeNodePDE tn10; tn10.i = ln+1; tn10.t = tn10.i*ht;
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;

        /**************************************************** border conditions ***************************************************/
        SpaceNodePDE sn0;
        sn0.i = static_cast<int>(0); sn0.x = 0*hx; p20[0] = b_boundary(sn0, tn20);
        sn0.i = static_cast<int>(N); sn0.x = N*hx; p20[N] = b_boundary(sn0, tn20);
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        currentLayerBGrid(p10, cntrlDeltaGridList, measuremntGirdList, ln+1, u_info);
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            dx[n-1] = 0.0;
            dx[n-1] += 2.0*p10[n] - p00[n] + alpha_ht_05*p00[n];
            dx[n-1] += (p10[n-1] - 2.0*p10[n] + p10[n+1]) * p_aa_htht__hxhx_1m2lambda;
            dx[n-1] += (p00[n-1] - 2.0*p00[n] + p00[n+1]) * p_aa_htht__hxhx_lambda;
            dx[n-1] += ht_ht * mCrBfxWeightMatrix[n];
            //------------------------------------- Adding functional part --------------------------------//
            if (L <= ln && ln <= LLD) dx[n-1] += -2.0*mu(n)*(u.at(ln-L+1)[n]) * ht_ht;
            //------------------------------------- Adding functional part --------------------------------//
        }
        dx[0]   -= p20[0] * m_aa_htht__hxhx_lambda;
        dx[N-2] -= p20[N] * m_aa_htht__hxhx_lambda;
        tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
        for (unsigned int n=1; n<=N-1; n++) p20[n] = rx[n-1];
        if (use == true) add2Info(p20, p_info, ln, hx, cntrlDeltaGridList); b_layerInfo(p20, ln);
        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int n=0; n<=N; n++)
        {
            p00[n] = p10[n];
            p10[n] = p20[n];
        }
    }
    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);

    for (unsigned int j=0; j<No; j++) measuremntGirdList[j].cleanGrid(); measuremntGirdList.clear();
    for (unsigned int i=0; i<Nc; i++) cntrlDeltaGridList[i].cleanGrid(); cntrlDeltaGridList.clear();

    p00.clear();
    p10.clear();
    p20.clear();
}


