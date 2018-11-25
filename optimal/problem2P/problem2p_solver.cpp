#include "problem2p_solver.h"

Problem2PNeumann::Problem2PNeumann() {}

Problem2PNeumann::~Problem2PNeumann() {}

auto Problem2PNeumann::gradient(const DoubleVector &pv, DoubleVector &gv) const -> void
{
    const unsigned int L   = static_cast<const unsigned int>( mtimeDimension.sizeN() );
    const double ht        = mtimeDimension.step();
    const unsigned int Nc  = mEquParameter.Nc;
    const unsigned int No  = mEquParameter.No;
#ifdef TIME_DISCRETE
    const unsigned int Nt = mEquParameter.Nt;
#endif

    OptimizeParameterP mOptParameter;
    VectorToPrm(pv, mOptParameter);

    DoubleMatrix u;

    spif_vector u_info;
    solveForwardIBVP(u, u_info, true, mOptParameter);
    spif_vector p_info;
    solveBackwardIBVP(u, p_info, true, u_info, mOptParameter);

    gv.clear();
    gv.resize(pv.length(), 0.0);
    unsigned int gi = 0;

#ifdef TIME_DISCRETE
    std::vector<unsigned int> discrete_times;
    discrete_times.push_back(0);
    for (unsigned int s=0; s<Nt; s++)
    {
        const double tau = mOptParameter.tau[s];
        discrete_times.push_back(static_cast<unsigned int>(round(tau*L)));
    }
    discrete_times.push_back(L);
#endif

    // k
    if (optimizeK)
    {
#ifdef TIME_DISCRETE
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfoP &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointInfoP &uj = u_info[j];

                const double zij = mOptParameter.z[i][j];

                double grad_Kij = 0.0;
                for (unsigned int s=0; s<=Nt; s++)
                {
                    const unsigned int s0 = discrete_times[s+0];
                    const unsigned int s1 = discrete_times[s+1];

                    double pi_in = 0.0;
                    pi_in += 0.5 * ht * (pi.vl[2*s0] /*+ 2.0*r*gpi(i,(2*s0),u_info,mOptParameter)*sgn(g0i(i,(2*s0),u_info,mOptParameter))*/);
                    for (unsigned int ln=s0+1; ln<=s1-1; ln++)
                    {
                        pi_in += ht * (pi.vl[2*ln] /*+ 2.0*r*gpi(i,(2*s0),u_info,mOptParameter)*sgn(g0i(i,(2*s0),u_info,mOptParameter))*/);
                    }
                    pi_in += 0.5 * ht * (pi.vl[2*s1] /*+ 2.0*r*gpi(i,(2*s0),u_info,mOptParameter)*sgn(g0i(i,(2*s0),u_info,mOptParameter))*/);

                    const double us = (uj.vl[2*s0] - zij);

                    grad_Kij += pi_in*us;
                }

                gv[gi++] = -grad_Kij /*+ 2.0*regEpsilon*(mOptParameter.k[i][j] - mRegParameter.k[i][j])*/;
            }
        }
#else
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfoP &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointInfoP &uj = u_info[j];

                double grad_Kij = 0.0;
                double zij = mOptParameter.z[i][j];

                grad_Kij += 0.5 * (pi.vl[0] + 2.0*r1*gpi(i,0,u_info,mOptParameter)*sgn(g0i(i,0,u_info,mOptParameter))) * (uj.vl[0] - zij);
                for (unsigned int ln=1; ln<=L-1; ln++)
                {
                    grad_Kij += (pi.vl[2*ln] + 2.0*r1*gpi(i,2*ln,u_info,mOptParameter)*sgn(g0i(i,2*ln,u_info,mOptParameter))) * (uj.vl[2*ln] - zij);
                }
                grad_Kij += 0.5 * (pi.vl[2*L] + 2.0*r1*gpi(i,2*L,u_info,mOptParameter)*sgn(g0i(i,2*L,u_info,mOptParameter))) * (uj.vl[2*L] - zij);

                grad_Kij *= -ht;

                gv[gi++] = grad_Kij + 2.0*regEpsilon*(mOptParameter.k[i][j] - mRegParameter.k[i][j]);
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
                gv[gi++] = 0.0;
            }
        }
    }

    // z
    if (optimizeZ)
    {
#ifdef TIME_DISCRETE
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfoP &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                const double kij = mOptParameter.k[i][j];

                double grad_Zij = 0.0;
                for (unsigned int s=0; s<=Nt; s++)
                {
                    const unsigned int s0 = discrete_times[s+0];
                    const unsigned int s1 = discrete_times[s+1];

                    double pi_in = 0.0;
                    pi_in += 0.5 * ht * (pi.vl[2*s0] /*+ 2.0*r*gpi(i,(2*s0),u_info,mOptParameter)*sgn(g0i(i,(2*s0),u_info,mOptParameter))*/);
                    for (unsigned int ln=s0+1; ln<=s1-1; ln++)
                    {
                        pi_in += ht * (pi.vl[2*ln] /*+ 2.0*r*gpi(i,(2*s0),u_info,mOptParameter)*sgn(g0i(i,(2*s0),u_info,mOptParameter)*/);
                    }
                    pi_in += 0.5 * ht * (pi.vl[2*s1] /*+ 2.0*r*gpi(i,(2*s1),u_info,mOptParameter)*sgn(g0i(i,(2*s1),u_info,mOptParameter))*/);

                    grad_Zij += pi_in*kij;
                }

                gv[gi++] = grad_Zij/* + 2.0*regEpsilon*(mOptParameter.z[i][j] - mRegParameter.z[i][j])*/;
            }
        }
#else
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfoP &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                double grad_Zij = 0.0;
                double kij = mOptParameter.k[i][j];

                grad_Zij += 0.5 * (pi.vl[0] + 2.0*r1*gpi(i,0,u_info,mOptParameter)*sgn(g0i(i,0,u_info,mOptParameter))) * kij;
                for (unsigned int ln=1; ln<=L-1; ln++)
                {
                    grad_Zij += (pi.vl[2*ln] + 2.0*r1*gpi(i,2*ln,u_info,mOptParameter)*sgn(g0i(i,2*ln,u_info,mOptParameter))) * kij;
                }
                grad_Zij += 0.5 * (pi.vl[2*L] + 2.0*r1*gpi(i,2*L,u_info,mOptParameter)*sgn(g0i(i,2*L,u_info,mOptParameter))) * kij;
                grad_Zij *= ht;

                gv [gi++] = grad_Zij + 2.0*regEpsilon*(mOptParameter.z[i][j] - mRegParameter.z[i][j]);
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
                gv[gi++] = 0.0;
            }
        }
    }

    // xi
    if (optimizeO)
    {
#ifdef TIME_DISCRETE
        for (unsigned int j=0; j<No; j++)
        {
            const SpacePointInfoP &uj = u_info[j];

            double gradXijX = 0.0;
            double gradXijY = 0.0;
            double vi = 0.0;

            for (unsigned int s=0; s<=Nt; s++)
            {
                const unsigned int s0 = discrete_times[s+0];
                const unsigned int s1 = discrete_times[s+1];

                vi = 0.0;
                for (unsigned int i=0; i<Nc; i++) vi += mOptParameter.k[i][j] * (0.5 * ht * p_info[i].vl[2*s0] /*+ 2.0*r*gpi(i,2*s0,u_info,mOptParameter)*sgn(g0i(i,2*s0,u_info,mOptParameter))*/);
                gradXijX += uj.dx[2*s0] * vi;
                gradXijY += uj.dy[2*s0] * vi;

                for (unsigned int ln=s0+1; ln<=s1-1; ln++)
                {
                    vi = 0.0;
                    for (unsigned int i=0; i<Nc; i++) vi += mOptParameter.k[i][j]*(ht * p_info[i].vl[2*ln] /*+ 2.0*r*gpi(i,2*s0,u_info,mOptParameter)*sgn(g0i(i,2*s0,u_info,mOptParameter))*/);
                    gradXijX += uj.dx[2*s0] * vi;
                    gradXijY += uj.dy[2*s0] * vi;
                }

                vi = 0.0;
                for (unsigned int i=0; i<Nc; i++) vi += mOptParameter.k[i][j]*(0.5 * ht * p_info[i].vl[2*s1] /*+ 2.0*r*gpi(i,2*s0,u_info,mOptParameter)*sgn(g0i(i,2*s0,u_info,mOptParameter))*/);
                gradXijX += uj.dx[2*s0] * vi;
                gradXijY += uj.dy[2*s0] * vi;
            }

            //gradXijX *= -ht;
            //gradXijY *= -ht;

            gv[gi++] = -gradXijX /*+ 2.0*regEpsilon*(mOptParameter.xi[j].x - mRegParameter.xi[j].x)*/;
            gv[gi++] = -gradXijY /*+ 2.0*regEpsilon*(mOptParameter.xi[j].y - mRegParameter.xi[j].y)*/;
        }
#else
        for (unsigned int j=0; j<No; j++)
        {
            const SpacePointInfoP &uj = u_info[j];

            double gradXijX = 0.0;
            double gradXijY = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int i=0; i<Nc; i++) vi += mOptParameter.k[i][j] * (p_info[i].vl[0] + 2.0*r1*gpi(i,0,u_info,mOptParameter)*sgn(g0i(i,0,u_info,mOptParameter)));
            gradXijX += 0.5 * uj.dx[0] * vi;
            gradXijY += 0.5 * uj.dy[0] * vi;

            for (unsigned int ln=1; ln<=L-1; ln++)
            {
                vi = 0.0;
                for (unsigned int i=0; i<Nc; i++) vi += mOptParameter.k[i][j]*(p_info[i].vl[2*ln] + 2.0*r1*gpi(i,2*ln,u_info,mOptParameter)*sgn(g0i(i,2*ln,u_info,mOptParameter)));
                gradXijX += uj.dx[2*ln] * vi;
                gradXijY += uj.dy[2*ln] * vi;
            }

            vi = 0.0;
            for (unsigned int i=0; i<Nc; i++) vi += mOptParameter.k[i][j]*(p_info[i].vl[2*L] + 2.0*r1*gpi(i,2*L,u_info,mOptParameter)*sgn(g0i(i,2*L,u_info,mOptParameter)));
            gradXijX += 0.5 * uj.dx[2*L] * vi;
            gradXijY += 0.5 * uj.dy[2*L] * vi;

            gradXijX *= -ht;
            gradXijY *= -ht;

            double p2x = 0.0;
            double p2y = 0.0;
            const SpacePoint &eta = mOptParameter.eta[j];
            for (unsigned int i=0; i<Nc; i++)
            {
                const SpacePoint &xi = mOptParameter.xi[i];
                p2x += (eta.x - xi.x)*hpij(eta, xi);
                p2y += (eta.y - xi.y)*hpij(eta, xi);
            }

            gv[gi++] = gradXijX + 2.0*regEpsilon*(mOptParameter.xi[j].x - mRegParameter.xi[j].x) - 4.0*r2*p2x;
            gv[gi++] = gradXijY + 2.0*regEpsilon*(mOptParameter.xi[j].y - mRegParameter.xi[j].y) - 4.0*r2*p2y;

        }
#endif
    }
    else
    {
        for (unsigned int j=0; j<No; j++)
        {
            gv[gi++] = 0.0;
            gv[gi++] = 0.0;
        }
    }

    // eta
    if (optimizeC)
    {
#ifdef TIME_DISCRETE
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfoP &pi = p_info[i];

            double gradEtaiX = 0.0;
            double gradEtaiY = 0.0;

            for (unsigned int s=0; s<=Nt; s++)
            {
                const unsigned int s0 = discrete_times[s+0];
                const unsigned int s1 = discrete_times[s+1];

                double vi = 0.0;
                for (unsigned int j=0; j<No; j++) vi += mOptParameter.k[i][j] * (u_info[j].vl[2*s0] - mOptParameter.z[i][j]);

                gradEtaiX += 0.5 * ht * pi.dx[2*s0] * vi;
                gradEtaiY += 0.5 * ht * pi.dy[2*s0] * vi;

                for (unsigned int ln=s0+1; ln<=s1-1; ln++)
                {
                    gradEtaiX += ht * pi.dx[2*ln] * vi;
                    gradEtaiY += ht * pi.dy[2*ln] * vi;
                }

                gradEtaiX += 0.5 * ht * pi.dx[2*s1] * vi;
                gradEtaiY += 0.5 * ht * pi.dy[2*s1] * vi;
            }
            gv[gi++] = -gradEtaiX /*+ 2.0*regEpsilon*(mOptParameter.eta[i].x - mRegParameter.eta[i].x)*/;
            gv[gi++] = -gradEtaiY /*+ 2.0*regEpsilon*(mOptParameter.eta[i].y - mRegParameter.eta[i].y)*/;
        }
#else
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfoP &pi = p_info[i];

            double gradEtaiX = 0.0;
            double gradEtaiY = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int j=0; j<No; j++) vi += mOptParameter.k[i][j] * (u_info[j].vl[0] - mOptParameter.z[i][j]);
            gradEtaiX += 0.5 * pi.dx[0] * vi;
            gradEtaiY += 0.5 * pi.dy[0] * vi;

            for (unsigned int ln=1; ln<=L-1; ln++)
            {
                vi = 0.0;
                for (unsigned int j=0; j<No; j++) vi += mOptParameter.k[i][j] * (u_info[j].vl[2*ln] - mOptParameter.z[i][j]);
                gradEtaiX += pi.dx[2*ln] * vi;
                gradEtaiY += pi.dy[2*ln] * vi;
            }

            vi = 0.0;
            for (unsigned int j=0; j<No; j++) vi += mOptParameter.k[i][j] * (u_info[j].vl[2*L] - mOptParameter.z[i][j]);
            gradEtaiX += 0.5 * pi.dx[2*L] * vi;
            gradEtaiY += 0.5 * pi.dy[2*L] * vi;

            gradEtaiX *= -ht;
            gradEtaiY *= -ht;

            double p2x = 0.0;
            double p2y = 0.0;
            const SpacePoint &xi = mOptParameter.xi[i];
            for (unsigned int j=0; j<No; j++)
            {
                const SpacePoint &eta = mOptParameter.eta[j];
                p2x += (eta.x - xi.x)*hpij(eta, xi);
                p2y += (eta.y - xi.y)*hpij(eta, xi);
            }

            gv[gi++] = gradEtaiX + 2.0*regEpsilon*(mOptParameter.eta[i].x - mRegParameter.eta[i].x) + 4.0*r2*p2x;
            gv[gi++] = gradEtaiY + 2.0*regEpsilon*(mOptParameter.eta[i].y - mRegParameter.eta[i].y) + 4.0*r2*p2x;
        }
#endif
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            gv[gi++] = 0.0;
            gv[gi++] = 0.0;
        }
    }

#ifdef TIME_DISCRETE
    // tau
    if (optimezeT)
    {
        double a = mEquParameter.a;
        double alpha = mEquParameter.alpha;
        double theta = mEquParameter.theta;

        for (unsigned int s=0; s<Nt; s++)
        {
            const unsigned int taum = discrete_times[s+0];
            const unsigned int taus = discrete_times[s+1];
            const unsigned int taup = discrete_times[s+2];

            printf("%4u %4u %4u %4u\n", s, taum, taus, taup);
            //continue;

            double grad = 0.0;
            for (unsigned int i=0; i<Nc; i++)
            {
                const SpacePointInfoP &pi = p_info[i];
                const double pi_vl = pi.vl[taus*2];

                double pi_in = 0.0;
                pi_in += 0.5 * ht * pi.vl[2*taus];
                for (unsigned ln=taus+1; ln<=taup-1; ln++)
                {
                    pi_in += ht * pi.vl[2*ln];
                }
                pi_in += 0.5 * ht * pi.vl[2*taup];

                for (unsigned int j=0; j<No; j++)
                {
                    const SpacePointInfoP &uj = u_info[j];

                    double kij = mOptParameter.k[i][j];

                    const double uj_vlm = uj.vl[2*taum];
                    const double uj_vls = uj.vl[2*taus];

                    grad += pi_vl * kij * (uj_vlm - uj_vls);
                    //grad += pi_in * kij * (uj.vl[2*taus+2] - uj.vl[2*taus-2])/(2*ht);
                    grad += pi_in * kij * a*a*(uj.dxx[2*taus] + uj.dyy[2*taus] - alpha*(uj.vl[2*taus]-theta));
                }
            }
            //printf("---- %u %u %u %u\n", s, taum, taus, taup);
            gv[gi++] = -grad;
        }
    }
    else
    {
        for (unsigned int nt=0; nt<Nt; nt++) gv[gi++] = 0.0;
    }
#endif

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

auto Problem2PNeumann::fx(const DoubleVector &pv) const -> double
{
    OptimizeParameterP mOptParameter;

    VectorToPrm(pv, mOptParameter);

    DoubleMatrix u;
    spif_vector u_info;
    solveForwardIBVP(u, u_info, true, mOptParameter);

    double intgrl = fx_integral(u);

    double nrm = fx_norm(mEquParameter, mOptParameter, mRegParameter);

    double pnt1 = fx_penalty1(u_info, mOptParameter);
    double pnt2 = fx_penalty2(mOptParameter);

    double sum = intgrl + regEpsilon*nrm + r1*pnt1 + r2*pnt2;

    for (unsigned int j=0; j<u_info.size(); j++)
    {
        u_info[j].clear();
    }
    u_info.clear();

    return sum;

}

auto Problem2PNeumann::mu(double, double) const -> double
{
    return 1.0;
}

auto Problem2PNeumann::fx_integral(const DoubleMatrix &u) const -> double
{
    const double hx      = spaceDimension(Dimension::DimensionX).step();
    const double hy      = spaceDimension(Dimension::DimensionY).step();
    const unsigned int N = static_cast<const unsigned int> ( spaceDimension(Dimension::DimensionX).sizeN() );
    const unsigned int M = static_cast<const unsigned int> ( spaceDimension(Dimension::DimensionY).sizeN() );

    double udiff = 0.0;
    double usum = 0.0;

    udiff = (u[0][0]-U[0][0]); usum += 0.25 * udiff * udiff;
    udiff = (u[0][N]-U[0][N]); usum += 0.25 * udiff * udiff;
    udiff = (u[M][0]-U[M][0]); usum += 0.25 * udiff * udiff;
    udiff = (u[M][N]-U[M][N]); usum += 0.25 * udiff * udiff;

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = (u[0][n]-U[0][n]); usum += 0.5 * udiff * udiff;
        udiff = (u[M][n]-U[M][n]); usum += 0.5 * udiff * udiff;
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = (u[m][0]-U[m][0]); usum += 0.5 * udiff * udiff;
        udiff = (u[m][N]-U[m][N]); usum += 0.5 * udiff * udiff;
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = (u[m][n]-U[m][n]); usum += udiff * udiff;
        }
    }

    return usum*(hx*hy);
}

auto Problem2PNeumann::fx_norm(const EquationParameterP& e_prm, const OptimizeParameterP &o_prm, const OptimizeParameterP &r_prm) const -> double
{
    double _norm = 0.0;
    const unsigned int Nc = e_prm.Nc;
    const unsigned int No = e_prm.No;

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            _norm += (o_prm.k[i][j] - r_prm.k[i][j])*(o_prm.k[i][j] - r_prm.k[i][j]);
            _norm += (o_prm.z[i][j] - r_prm.z[i][j])*(o_prm.z[i][j] - r_prm.z[i][j]);
        }
    }

    for (unsigned int j=0; j<No; j++)
    {
        _norm += (o_prm.xi[j].x - r_prm.xi[j].x)*(o_prm.xi[j].x - r_prm.xi[j].x);
        _norm += (o_prm.xi[j].y - r_prm.xi[j].y)*(o_prm.xi[j].y - r_prm.xi[j].y);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        _norm += (o_prm.eta[i].x - r_prm.eta[i].x)*(o_prm.eta[i].x - r_prm.eta[i].x);
        _norm += (o_prm.eta[i].y - r_prm.eta[i].y)*(o_prm.eta[i].y - r_prm.eta[i].y);
    }

    return _norm;
}

auto Problem2PNeumann::fx_penalty1(const spif_vector &info, const OptimizeParameterP &o_prm) const -> double
{
    const double ht = mtimeDimension.step();

    const unsigned int L = static_cast<const unsigned int> ( mtimeDimension.sizeN() );

    double _penalty = 0.0;
    for (unsigned int i=0; i<mEquParameter.Nc; i++)
    {
        double pnlt_i = 0.0;
        double _gpi_0 = gpi(i, 0, info, o_prm);
        pnlt_i += 0.5*_gpi_0*_gpi_0;

        for (unsigned int l=1; l<=L-1; l++)
        {
            double _gpi_l = gpi(i, 2*l, info, o_prm);
            pnlt_i += _gpi_l*_gpi_l;
        }

        double _gpi_L = gpi(i, 2*L, info, o_prm);
        pnlt_i += 0.5*_gpi_L*_gpi_L;

        _penalty += pnlt_i*ht;
    }

    return _penalty;
}

auto Problem2PNeumann::gpi(unsigned int i, unsigned int layer, const spif_vector &u_info, const OptimizeParameterP &o_prm) const -> double
{
    double gpi_ln = fabs( g0i(i, layer, u_info, o_prm) ) - ( vmax.at(i) - vmin.at(i) )/2.0;
    return gpi_ln > 0.0 ? gpi_ln : 0.0;
}

auto Problem2PNeumann::g0i(unsigned int i, unsigned int layer, const spif_vector &u_info, const OptimizeParameterP &o_prm) const -> double
{
    double vi = 0.0;
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        const SpacePointInfoP &u_xij = u_info[j];
        vi += o_prm.k[i][j] * ( u_xij.vl[layer] - o_prm.z[i][j] );
    }
    return ( vmax.at(i) + vmin.at(i) )/2.0 - vi;
}

auto Problem2PNeumann::sign(double x) const -> double
{
    if (x < 0.0)       return -1.0;
    else if (x > 0.0)  return +1.0;
    else               return  0.0;
}

auto Problem2PNeumann::fx_penalty2(const OptimizeParameterP &o_prm) const -> double
{
    double pntl = 0.0;
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    for (unsigned int i=0; i<Nc; i++)
    {
        const SpacePoint &eta = o_prm.eta[i];
        for (unsigned int j=0; j<No; j++)
        {
            const SpacePoint &xi = o_prm.xi[j];
            pntl += hpij(eta, xi);
        }
    }

    return pntl*pntl;
}

auto Problem2PNeumann::hpij(const SpacePoint &eta, const SpacePoint &xi) const -> double
{
    double pij = d*d - (eta.x-xi.x)*(eta.x-xi.x) - (eta.y-xi.y)*(eta.y-xi.y);
    //printf("%f %f %f %f %f %f\n", d, pij, eta.x, eta.y, xi.x, xi.y);
    return (pij < 0.0) ? 0.0 : pij;
}

auto Problem2PNeumann::print(unsigned int i, const DoubleVector &x, const DoubleVector &g,
                             double f, double alpha, GradientMethod::MethodResult result) const -> void
{
    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(f); C_UNUSED(alpha); C_UNUSED(result);
    const char* msg = ""; C_UNUSED(msg);
    if (result == GradientMethod::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION   ";
    if (result == GradientMethod::FIRST_ITERATION)          msg = "FIRST_ITERATION         ";
    if (result == GradientMethod::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    if (result == GradientMethod::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS     ";
    if (result == GradientMethod::NEXT_ITERATION)           msg = "NEXT_ITERATION          ";

    Problem2PNeumann* prob = const_cast<Problem2PNeumann*>(this);
    OptimizeParameterP mOptParameter;
    VectorToPrm(x, mOptParameter);
    DoubleMatrix u;
    spif_vector u_info;
    solveForwardIBVP(u, u_info, true, mOptParameter);
    double ing = fx_integral(u);
    double pnt1 = fx_penalty1(u_info, mOptParameter);
    double pnt2 = fx_penalty2(mOptParameter);
    double nrm = fx_norm(mEquParameter, mOptParameter, mRegParameter);

    const unsigned int v_length = static_cast<const unsigned int>( timeDimension().sizeN() );
    DoubleVector v1(v_length+1);
    DoubleVector v2(v_length+1);

    for (unsigned int ln=0; ln<=v_length; ln++)
    {
        v1[ln] = mOptParameter.k[0][0] * (u_info.at(0).vl[2*ln] - mOptParameter.z[0][0]) + mOptParameter.k[0][1] * (u_info.at(1).vl[2*ln] - mOptParameter.z[0][1]);
        v2[ln] = mOptParameter.k[1][0] * (u_info.at(0).vl[2*ln] - mOptParameter.z[1][0]) + mOptParameter.k[1][1] * (u_info.at(1).vl[2*ln] - mOptParameter.z[1][1]);
    }

    IPrinter::printVector(v1, "v1: ", 10);
    IPrinter::printVector(v2, "v2: ", 10);

    printf("I[%3d]: F:%10.6f I:%10.6f P1:%12.6f P2:%12.6f N:%10.6f R:%7.3f e:%5.3f a:%10.6f min:%10.6f max:%10.6f \n", i, f, ing, pnt1, pnt2, nrm, r1, regEpsilon, alpha, u.min(), u.max());
    printf("k:%10.4f %10.4f %10.4f %10.4f z:%10.4f %10.4f %10.4f %10.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    printf("k:%10.4f %10.4f %10.4f %10.4f z:%10.4f %10.4f %10.4f %10.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);

    u.clear();
    u_info.clear();

    C_UNUSED(prob);
    IPrinter::printSeperatorLine();

//        prob->optimizeK = (i%2==1);
//        prob->optimizeZ = (i%2==0);
//        prob->optimizeO = (i%4==2);
//        prob->optimizeC = (i%4==3);
}

auto Problem2PNeumann::norm(const DoubleVector &v) const -> double
{
    return EuclideanNorm(v);
}

auto Problem2PNeumann::normalize(DoubleVector &v) const -> void
{
    DoubleVector kv = v.mid(0, 3);   IVectorNormalizer::EuclideanNormalize(kv);
    DoubleVector zv = v.mid(4, 7);   IVectorNormalizer::EuclideanNormalize(zv);
    DoubleVector ov = v.mid(8, 11);  IVectorNormalizer::EuclideanNormalize(ov);
    DoubleVector cv = v.mid(12, 15); IVectorNormalizer::EuclideanNormalize(cv);

    v[0]  = kv[0]; v[1]  = kv[1];  v[2]  = kv[2]; v[3]  = kv[3];
    v[4]  = zv[0]; v[5]  = zv[1];  v[6]  = zv[2]; v[7]  = zv[3];
    v[8]  = ov[0]; v[9]  = ov[1];  v[10] = ov[2]; v[11] = ov[3];
    v[12] = cv[0]; v[13] = cv[1];  v[14] = cv[2]; v[15] = cv[3];

    kv.clear(); zv.clear(); ov.clear(); cv.clear();
}

auto Problem2PNeumann::project(DoubleVector &, unsigned int) -> void
{}

auto Problem2PNeumann::project(DoubleVector &pv) const -> void
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    unsigned int start = 2*Nc*No;
    unsigned int end  =  2*Nc*No + 2*No + 2*Nc - 1;

    for (unsigned int index = start; index <= end; index++)
    {
        if (pv[index] < 0.05) pv[index] = 0.05;
        if (pv[index] > 0.95) pv[index] = 0.95;
    }

    //IPrinter::print(pv.mid(start, end));
    for (unsigned int index = start; index <=end; index++)
    {
        //projectControlPoints(pv, index);
        projectMeasurePoints(pv, index);
    }
    //IPrinter::print(pv.mid(start, end));
}

auto Problem2PNeumann::projectControlPoints(DoubleVector &pv, unsigned int index) const -> void
{
    double dist = 0.10;

    if (index == 12)
    {
        if (fabs(pv[12] - pv[ 8])<dist) pv[12] = pv[ 8] + sign(pv[12] - pv[ 8])*dist;
        if (fabs(pv[12] - pv[10])<dist) pv[12] = pv[10] + sign(pv[12] - pv[10])*dist;
        if (pv[12] < 0.05) pv[12] = 0.05;
        if (pv[12] > 0.95) pv[12] = 0.95;
    }

    if (index == 13)
    {
        if (fabs(pv[13] - pv[ 9])<dist) pv[13] = pv[ 9] + sign(pv[13] - pv[ 9])*dist;
        if (fabs(pv[13] - pv[11])<dist) pv[13] = pv[11] + sign(pv[13] - pv[11])*dist;
        if (pv[13] < 0.05) pv[13] = 0.05;
        if (pv[13] > 0.95) pv[13] = 0.95;
    }

    if (index == 14)
    {
        if (fabs(pv[14] - pv[ 8])<dist) pv[14] = pv[ 8] + sign(pv[14] - pv[ 8])*dist;
        if (fabs(pv[14] - pv[10])<dist) pv[14] = pv[10] + sign(pv[14] - pv[10])*dist;
        if (pv[14] < 0.05) pv[14] = 0.05;
        if (pv[14] > 0.95) pv[14] = 0.95;
    }

    if (index == 15)
    {
        if (fabs(pv[15] - pv[ 9])<dist) pv[15] = pv[ 9] + sign(pv[15] - pv[ 9])*dist;
        if (fabs(pv[15] - pv[11])<dist) pv[15] = pv[11] + sign(pv[15] - pv[11])*dist;
        if (pv[15] < 0.05) pv[15] = 0.05;
        if (pv[15] > 0.95) pv[15] = 0.95;
    }
}

auto Problem2PNeumann::projectMeasurePoints(DoubleVector &pv, unsigned int index) const -> void
{
    double dist = 0.10;

    if (index == 8)
    {
        if (fabs(pv[ 8] - pv[12])<dist) pv[8] = pv[12] + sign(pv[ 8] - pv[12])*dist;
        if (fabs(pv[ 8] - pv[14])<dist) pv[8] = pv[14] + sign(pv[ 8] - pv[14])*dist;
        if (pv[8] < 0.05) pv[8] = 0.05;
        if (pv[8] > 0.95) pv[8] = 0.95;
    }

    if (index == 9)
    {
        if (fabs(pv[ 9] - pv[13])<dist) pv[9] = pv[13] + sign(pv[ 9] - pv[13])*dist;
        if (fabs(pv[ 9] - pv[15])<dist) pv[9] = pv[15] + sign(pv[ 9] - pv[15])*dist;
        if (pv[9] < 0.05) pv[9] = 0.05;
        if (pv[9] > 0.95) pv[9] = 0.95;
    }

    if (index == 10)
    {
        if (fabs(pv[10] - pv[12])<dist) pv[10] = pv[12] + sign(pv[10] - pv[12])*dist;
        if (fabs(pv[10] - pv[14])<dist) pv[10] = pv[14] + sign(pv[10] - pv[14])*dist;
        if (pv[10] < 0.05) pv[10] = 0.05;
        if (pv[10] > 0.95) pv[10] = 0.95;
    }

    if (index == 11)
    {
        if (fabs(pv[11] - pv[13])<dist) pv[11] = pv[13] + sign(pv[11] - pv[13])*dist;
        if (fabs(pv[11] - pv[15])<dist) pv[11] = pv[15] + sign(pv[11] - pv[15])*dist;
        if (pv[11] < 0.05) pv[11] = 0.05;
        if (pv[11] > 0.95) pv[11] = 0.95;
    }
}

auto Problem2PNeumann::solveForwardIBVP(DoubleMatrix &u, spif_vector &u_info, bool use, const OptimizeParameterP &mOptParameter) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<const unsigned int> ( dimX.sizeN() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.sizeN() );
    const unsigned int L = static_cast<const unsigned int> ( time.sizeN() );

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a        = mEquParameter.a;
    const double alpha    = mEquParameter.alpha;
    const double lambda   = mEquParameter.lambda;
    const double theta    = mEquParameter.theta;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;
#ifdef TIME_DISCRETE
    const unsigned int Nt = mEquParameter.Nt;
#endif

    const double m_aa_ht05__hxhx = -(0.5*a*a*ht)/(hx*hx);
    const double p_aa_ht__hxhx___alpha_ht05 = +1.0 + (a*a*ht)/(hx*hx) + 0.5*alpha*ht;
    const double p_aa_ht05__hyhy = +(0.5*a*a*ht)/(hy*hy);

    const double m_aa_ht05__hyhy = -(0.5*a*a*ht)/(hy*hy);
    const double p_aa_ht__hyhy___alpha_ht05 = +1.0 + (a*a*ht)/(hy*hy) + 0.5*alpha*ht;
    const double p_aa_ht05__hxhx = +(0.5*a*a*ht)/(hx*hx);

    const double ht05             = 0.5*ht;
    const double alpha_ht05       = 0.5*alpha*ht;
    const double aa_lambda_ht__hx = (a*a*lambda*ht)/hx;
    const double aa_lambda_ht__hy = (a*a*lambda*ht)/hy;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);

    u.clear();
    u.resize(M+1, N+1);

    //----------------------------------------------------------------------------------------------//
    std::vector<ExtendedSpacePointP> msnExtSpacePoints, cntExtSpacePoints;
    newDistributeDeltaGaussCntrl(mOptParameter.eta, cntExtSpacePoints, dimX, dimY);
    newDistributeDeltaGaussMsmnt(mOptParameter.xi,  msnExtSpacePoints, dimX, dimY);
    //----------------------------------------------------------------------------------------------//
    uint_vector rows0, rows1, rows2, cols0, cols1, cols2;
    f_findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, N, M, cntExtSpacePoints, msnExtSpacePoints);
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) f_prepareInfo(No, mOptParameter.xi, u_info, L);
    //----------------------------------------------------------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    f_initialLayers(u00, u_info, use, N, M, hx, hy, msnExtSpacePoints);
    //------------------------------------- initial conditions -------------------------------------//

    double *ax = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    double *bx = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    double *cx = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    double *dx = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    double *rx = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = m_aa_ht05__hxhx;
        bx[n] = p_aa_ht__hxhx___alpha_ht05;
        cx[n] = m_aa_ht05__hxhx;
    }
    ax[0] = 0.0; bx[0] += aa_lambda_ht__hx; cx[0] += m_aa_ht05__hxhx;
    cx[N] = 0.0; bx[N] += aa_lambda_ht__hx; ax[N] += m_aa_ht05__hxhx;

    double *ay = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    double *by = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    double *cy = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    double *dy = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    double *ry = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    for (unsigned int m=0; m<=M; m++)
    {
        ay[m] = m_aa_ht05__hyhy;
        by[m] = p_aa_ht__hyhy___alpha_ht05;
        cy[m] = m_aa_ht05__hyhy;
    }
    ay[0] = 0.0; by[0] += aa_lambda_ht__hy; cy[0] += m_aa_ht05__hyhy;
    cy[M] = 0.0; by[M] += aa_lambda_ht__hy; ay[M] += m_aa_ht05__hyhy;

    const unsigned int rows1_size = static_cast<const unsigned int>( rows1.size()*(N+1) );
    double *a1=nullptr, *b1=nullptr, *c1=nullptr, *d1=nullptr, *x1=nullptr, **w1=nullptr;
    if (rows1.size() != 0 && rows2.size() != 0)
    {
        a1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        b1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        c1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        d1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        x1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        w1 = static_cast<double**> ( malloc(sizeof(double*)*rows1_size) );
        for (unsigned int row=0; row < rows1_size; row++) w1[row] = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
    }

    const unsigned int cols1_size = static_cast<const unsigned int>( cols1.size()*(M+1) );
    double *a2=nullptr, *b2=nullptr, *c2=nullptr, *d2=nullptr, *x2=nullptr, **w2=nullptr;
    if (cols1.size() != 0 && cols2.size() != 0)
    {
        a2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        b2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        c2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        d2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        x2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        w2 = static_cast<double**> ( malloc(sizeof(double*)*cols1_size) );
        for (unsigned int col=0; col < cols1_size; col++) w2[col] = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
    }

    SpaceNodePDE sn;

#ifdef TIME_DISCRETE
    std::vector<unsigned int> discrete_times;
    discrete_times.push_back(0);
    for (unsigned int s=0; s<Nt; s++)
    {
        const double tau = mOptParameter.tau[s];
        discrete_times.push_back(static_cast<unsigned int>(round(tau*L)));
    }
    discrete_times.push_back(L);
    //printf("%4u %4u %4u %4u %4u %4u\n", discrete_times.size(), discrete_times[0], discrete_times[1], discrete_times[2], discrete_times[3], discrete_times[4]);
#endif

    for (unsigned int l=1; l<=L; l++)
    {

#ifdef TIME_DISCRETE
        unsigned int current_time = 0;
        for (unsigned int i=0; i<discrete_times.size(); i++)
        {
            if (discrete_times[i] < l) current_time = discrete_times[i];
        }
        //printf("F %4u %4u\n", l, current_time);
#endif

        /**************************************************** x direction apprx ***************************************************/

        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = rows0.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n] = 0.0;
#ifdef OH1
                    if (m>0 && m<M)  dx[n] += p_aa_ht05__hyhy*(u00[m-1][n] - 2.0*u00[m][n]   + u00[m+1][n]);
                    else if (m == 0) dx[n] += p_aa_ht05__hyhy*(u00[0][n]   - 2.0*u00[1][n]   + u00[2][n]);
                    else if (m == M) dx[n] += p_aa_ht05__hyhy*(u00[M-2][n] - 2.0*u00[M-1][n] + u00[M][n]);
#else
                    if (m>0 && m<M)  dx[n] += p_aa_ht05__hyhy*(u00[m-1][n] - 2.0*u00[m][n]   + u00[m+1][n]);
                    else if (m == 0) dx[n] += p_aa_ht05__hyhy*(2.0*u00[0][n] - 5.0*u00[1][n]   + 4.0*u00[2][n]   - u00[3][n]);
                    else if (m == M) dx[n] += p_aa_ht05__hyhy*(2.0*u00[M][n] - 5.0*u00[M-1][n] + 4.0*u00[M-2][n] - u00[M-3][n]);
#endif

                    dx[n] += u00[m][n] + alpha_ht05*theta;
                }

                dx[0] += aa_lambda_ht__hx*theta;
                dx[N] += aa_lambda_ht__hx*theta;

                tomasAlgorithm(ax, bx, cx, dx, rx, N+1);
                for (unsigned int n=0; n<=N; n++) u05[m][n] = rx[n];
            }
        }

        if (rows1.size() != 0 && rows2.size() == 0)
        {
            //throw grid_exception("forward x1");

            double* _u05 = new double[No]; for (unsigned int j=0; j<No; j++) _u05[j] = 0.0;

#ifdef TIME_DISCRETE
            for (unsigned int j=0; j<No; j++)
            {
                _u05[j] = u_info[j].vl[2*current_time];
                _u05[j] *= (1.0 + noise);
            }
#else
            for (unsigned int j=0; j<No; j++)
            {
                const ExtendedSpacePointP &extendedSpacePoint = msnExtSpacePoints.at(j);
                const std::vector<ExtendedSpacePointNodeP> &nodes = extendedSpacePoint.nodes;
                const unsigned int nodes_size = static_cast<const unsigned int>( nodes.size() );
                for (unsigned int nj=0; nj<nodes_size; nj++)
                {
                    const ExtendedSpacePointNodeP &node = nodes.at(nj);
                    unsigned int node_nx = static_cast<const unsigned int>(node.nx);
                    unsigned int node_ny = static_cast<const unsigned int>(node.ny);
                    _u05[j] += u05[node_ny][node_nx] * (node.w * (hx*hy));
                }
                _u05[j] *= (1.0 + noise);
            }
#endif

            double *_v05 = new double[Nc];

            for (unsigned int i=0; i<Nc; i++)
            {
                _v05[i] = 0.0;
                for (unsigned int j=0; j<No; j++)
                {
                    _v05[i] += mOptParameter.k[i][j] * ( _u05[j] - mOptParameter.z[i][j] );
                }
            }

            delete [] _u05;

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n] = 0.0;
#ifdef OH1
                    if (m>0 && m<M)  dx[n] += p_aa_ht05__hyhy*(u00[m-1][n] - 2.0*u00[m][n]   + u00[m+1][n]);
                    else if (m == 0) dx[n] += p_aa_ht05__hyhy*(u00[0][n]   - 2.0*u00[1][n]   + u00[2][n]);
                    else if (m == M) dx[n] += p_aa_ht05__hyhy*(u00[M-2][n] - 2.0*u00[M-1][n] + u00[M][n]);
#else
                    if (m>0 && m<M)  dx[n] += p_aa_ht05__hyhy*(u00[m-1][n] - 2.0*u00[m][n]   + u00[m+1][n]);
                    else if (m == 0) dx[n] += p_aa_ht05__hyhy*(2.0*u00[0][n] - 5.0*u00[1][n]   + 4.0*u00[2][n]   - u00[3][n]);
                    else if (m == M) dx[n] += p_aa_ht05__hyhy*(2.0*u00[M][n] - 5.0*u00[M-1][n] + 4.0*u00[M-2][n] - u00[M-3][n]);
#endif

                    dx[n] += u00[m][n] + alpha_ht05*theta;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePointP &extendedSpacePoint = cntExtSpacePoints.at(i);
                        if (extendedSpacePoint.contains(sn))
                        {
                            const std::vector<ExtendedSpacePointNodeP> &nodes = extendedSpacePoint.nodes;
                            const unsigned int nodes_size = static_cast<const unsigned int>( nodes.size() );
                            for (unsigned int ni=0; ni<nodes_size; ni++)
                            {
                                const ExtendedSpacePointNodeP &node = nodes.at(ni);
                                if (node.equals(sn)) dx[n] += ht05 * _v05[i] * node.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                dx[0] += aa_lambda_ht__hx*theta;
                dx[N] += aa_lambda_ht__hx*theta;

                tomasAlgorithm(ax, bx, cx, dx, rx, N+1);
                for (unsigned int n=0; n<=N; n++) u05[m][n] = rx[n];
            }

            delete [] _v05;
        }

        if (rows1.size() != 0 && rows2.size() != 0)
        {
            //throw grid_exception("forward x2");

            for (unsigned int m=0; m<rows1_size; m++) for (unsigned int n=0; n<rows1_size; n++) w1[m][n] = 0.0;

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    const unsigned int index = offset+n;
                    d1[index] = 0.0;
#ifdef OH1
                    if (m>0 && m<M) d1[index] = p_aa_ht05__hyhy*(u00[m-1][n] - 2.0*u00[m][n]   + u00[m+1][n]);
                    else if (m==0)  d1[index] = p_aa_ht05__hyhy*(u00[0][n]   - 2.0*u00[1][n]   + u00[2][n]);
                    else if (m==M)  d1[index] = p_aa_ht05__hyhy*(u00[M-2][n] - 2.0*u00[M-1][n] + u00[M][n]);
#else
                    if (m>0 && m<M)  d1[index] += p_aa_ht05__hyhy*(u00[m-1][n] - 2.0*u00[m][n]   + u00[m+1][n]);
                    else if (m == 0) d1[index] += p_aa_ht05__hyhy*(2.0*u00[0][n] - 5.0*u00[1][n]   + 4.0*u00[2][n]   - u00[3][n]);
                    else if (m == M) d1[index] += p_aa_ht05__hyhy*(2.0*u00[M][n] - 5.0*u00[M-1][n] + 4.0*u00[M-2][n] - u00[M-3][n]);
#endif

                    d1[index] += u00[m][n] + alpha_ht05*theta;

                    a1[index] = m_aa_ht05__hxhx;
                    b1[index] = p_aa_ht__hxhx___alpha_ht05;
                    c1[index] = m_aa_ht05__hxhx;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePointP &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                        if (cExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNodeP> &nodes1 = cExtendedSpacePoint.nodes;
                            for (unsigned int ni=0; ni<nodes1.size(); ni++)
                            {
                                const ExtendedSpacePointNodeP &node1 = nodes1.at(ni);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                const ExtendedSpacePointP &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                                const std::vector<ExtendedSpacePointNodeP> &nodes2 = mExtendedSpacePoint.nodes;
                                for (unsigned int nj=0; nj<nodes2.size(); nj++)
                                {
                                    const ExtendedSpacePointNodeP &node2 = nodes2.at(nj);
                                    unsigned int node2_nx = static_cast<unsigned int>(node2.nx);
                                    unsigned int node2_ny = static_cast<unsigned int>(node2.ny);

                                    bool found = false;
                                    for (unsigned int rs=0; rs<rows1.size(); rs++)
                                    {
                                        if (node2_ny == rows1[rs])
                                        {
                                            found = true;
                                            w1[index][rs*(N+1)+node2_nx] -= ht05 * mOptParameter.k[i][j] * (node2.w * (hx*hy)) * w * (1.0 + noise);
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d1[index] += ht05 * mOptParameter.k[i][j] * u00[node2_ny][node2_nx] * (node2.w * (hx*hy)) * w * (1.0 + noise);
                                    }
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[index] -= ht05 * mOptParameter.k[i][j] * mOptParameter.z[i][j] * w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1[offset+0] = 0.0; b1[offset+0] += aa_lambda_ht__hx; c1[offset+0] += m_aa_ht05__hxhx;
                c1[offset+N] = 0.0; b1[offset+N] += aa_lambda_ht__hx; a1[offset+N] += m_aa_ht05__hxhx;
                d1[offset+0] += aa_lambda_ht__hx*theta;
                d1[offset+N] += aa_lambda_ht__hx*theta;

                offset += N+1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w1, x1, rows1_size);

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=0; n<=N; n++)
                {
                    u05[m][n] = x1[offset+n];
                }
                offset += N+1;
            }
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        if (cols0.size() != 0)
        {
            for (unsigned int col=0; col<cols0.size(); col++)
            {
                unsigned int n = cols0.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m] = 0.0;
#ifdef OH1
                    if (n>0 && n<N)  dy[m] += p_aa_ht05__hxhx*(u05[m][n-1] - 2.0*u05[m][n]   + u05[m][n+1]);
                    else if (n == 0) dy[m] += p_aa_ht05__hxhx*(u05[m][0]   - 2.0*u05[m][1]   + u05[m][2]);
                    else if (n == N) dy[m] += p_aa_ht05__hxhx*(u05[m][N-2] - 2.0*u05[m][N-1] + u05[m][N]);
#else
                    if (n>0 && n<N)  dy[m] += p_aa_ht05__hxhx*(u05[m][n-1] - 2.0*u05[m][n]   + u05[m][n+1]);
                    else if (n == 0) dy[m] += p_aa_ht05__hxhx*(2.0*u05[m][0] - 5.0*u05[m][1]   + 4.0*u05[m][2]   - u05[m][3]);
                    else if (n == N) dy[m] += p_aa_ht05__hxhx*(2.0*u05[m][N] - 5.0*u05[m][N-1] + 4.0*u05[m][N-2] - u05[m][N-3]);
#endif

                    dy[m] += u05[m][n] + alpha_ht05*theta;
                }

                dy[0] += aa_lambda_ht__hy*theta;
                dy[M] += aa_lambda_ht__hy*theta;

                tomasAlgorithm(ay, by, cy, dy, ry, M+1);
                for (unsigned int m=0; m<=M; m++) u10[m][n] = ry[m];
            }
        }

        if (cols1.size() != 0 && cols2.size() == 0)
        {
            //throw grid_exception("forward y1");

            double* _u10 = new double[No]; for (unsigned int j=0; j<No; j++) _u10[j] = 0.0;

#ifdef TIME_DISCRETE
            for (unsigned int j=0; j<No; j++)
            {
                _u10[j] = u_info[j].vl[2*current_time];
                _u10[j] *= (1.0+noise);
            }
#else
            for (unsigned int j=0; j<No; j++)
            {
                const ExtendedSpacePointP &extendedSpacePoint = msnExtSpacePoints.at(j);
                const std::vector<ExtendedSpacePointNodeP> &nodes = extendedSpacePoint.nodes;
                const unsigned int nodes_size = static_cast<const unsigned int>( nodes.size() );
                for (unsigned int nj=0; nj<nodes_size; nj++)
                {
                    const ExtendedSpacePointNodeP &node = nodes.at(nj);
                    unsigned int node_nx = static_cast<unsigned int>(node.nx);
                    unsigned int node_ny = static_cast<unsigned int>(node.ny);
                    _u10[j] += u10[node_ny][node_nx] * (node.w * (hx*hy));
                }
                _u10[j] *= (1.0+noise);
            }
#endif

            double *_v10 = new double[Nc];

            for (unsigned int i=0; i<Nc; i++)
            {
                _v10[i] = 0.0;
                for (unsigned int j=0; j<No; j++)
                {
                    _v10[i] += mOptParameter.k[i][j] * ( _u10[j] - mOptParameter.z[i][j] );
                }
            }

            delete [] _u10;

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m] = 0.0;
#ifdef OH1
                    if (n>0 && n<N) dy[m] += p_aa_ht05__hxhx*(u05[m][n-1] - 2.0*u05[m][n]   + u05[m][n+1]);
                    else if (n==0)  dy[m] += p_aa_ht05__hxhx*(u05[m][0]   - 2.0*u05[m][1]   + u05[m][2]);
                    else if (n==N)  dy[m] += p_aa_ht05__hxhx*(u05[m][N-2] - 2.0*u05[m][N-1] + u05[m][N]);
#else
                    if (n>0 && n<N)  dy[m] += p_aa_ht05__hxhx*(u05[m][n-1] - 2.0*u05[m][n]   + u05[m][n+1]);
                    else if (n == 0) dy[m] += p_aa_ht05__hxhx*(2.0*u05[m][0] - 5.0*u05[m][1]   + 4.0*u05[m][2]   - u05[m][3]);
                    else if (n == N) dy[m] += p_aa_ht05__hxhx*(2.0*u05[m][N] - 5.0*u05[m][N-1] + 4.0*u05[m][N-2] - u05[m][N-3]);
#endif

                    dy[m] += u05[m][n] + alpha_ht05*theta;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePointP &extendedSpacePoint = cntExtSpacePoints.at(i);
                        if (extendedSpacePoint.contains(sn))
                        {
                            const std::vector<ExtendedSpacePointNodeP> &nodes = extendedSpacePoint.nodes;
                            const unsigned int nodes_size = static_cast<const unsigned int>( nodes.size() );
                            for (unsigned int ni=0; ni<nodes_size; ni++)
                            {
                                const ExtendedSpacePointNodeP &node = nodes.at(ni);
                                if (node.equals(sn)) dy[m] += ht05 * _v10[i] * node.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                dy[0] += aa_lambda_ht__hy*theta;
                dy[M] += aa_lambda_ht__hy*theta;

                tomasAlgorithm(ay, by, cy, dy, ry, M+1);
                for (unsigned int m=0; m<=M; m++) u10[m][n] = ry[m];
            }

            delete [] _v10;
        }

        if (cols1.size() != 0 && cols2.size() != 0)
        {
            //throw grid_exception("forward y2");

            for (unsigned int m=0; m<cols1_size; m++) for (unsigned int n=0; n<cols1_size; n++) w2[m][n] = 0.0;

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    const unsigned int index = offset+m;
                    d2[index] = 0.0;
#ifdef OH1
                    if (n>0 && n<N) d2[index] = p_aa_ht05__hxhx*(u05[m][n-1] - 2.0*u05[m][n]   + u05[m][n+1]);
                    else if (n==0)  d2[index] = p_aa_ht05__hxhx*(u05[m][0]   - 2.0*u05[m][1]   + u05[m][2]);
                    else if (n==N)  d2[index] = p_aa_ht05__hxhx*(u05[m][N-2] - 2.0*u05[m][N-1] + u05[m][N]);
#else
                    if (n>0 && n<N)  d2[index] += p_aa_ht05__hxhx*(u05[m][n-1] - 2.0*u05[m][n]   + u05[m][n+1]);
                    else if (n == 0) d2[index] += p_aa_ht05__hxhx*(2.0*u05[m][0] - 5.0*u05[m][1]   + 4.0*u05[m][2]   - u05[m][3]);
                    else if (n == N) d2[index] += p_aa_ht05__hxhx*(2.0*u05[m][N] - 5.0*u05[m][N-1] + 4.0*u05[m][N-2] - u05[m][N-3]);
#endif

                    d2[index] += u05[m][n] + alpha_ht05*theta;

                    a2[index] = m_aa_ht05__hyhy;
                    b2[index] = p_aa_ht__hyhy___alpha_ht05;
                    c2[index] = m_aa_ht05__hyhy;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePointP &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                        if (cExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNodeP> &nodes1 = cExtendedSpacePoint.nodes;
                            for (unsigned int ni=0; ni<nodes1.size(); ni++)
                            {
                                const ExtendedSpacePointNodeP &node1 = nodes1.at(ni);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                const ExtendedSpacePointP &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                                const std::vector<ExtendedSpacePointNodeP> &nodes2 = mExtendedSpacePoint.nodes;
                                for (unsigned int nj=0; nj<nodes2.size(); nj++)
                                {
                                    const ExtendedSpacePointNodeP &node2 = nodes2.at(nj);
                                    unsigned int node2_nx = static_cast<unsigned int>(node2.nx);
                                    unsigned int node2_ny = static_cast<unsigned int>(node2.ny);

                                    bool found = false;
                                    for (unsigned int cs=0; cs<cols1.size(); cs++)
                                    {
                                        if (node2_nx == cols1[cs])
                                        {
                                            found = true;
                                            w2[index][cs*(M+1)+node2_ny] -= ht05 * mOptParameter.k[i][j] * (node2.w * (hx*hy)) * w * (1.0 + noise);
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d2[index] += ht05 * mOptParameter.k[i][j] * u10[node2_ny][node2_nx] * (node2.w * (hx*hy)) * w * (1.0 + noise);
                                    }
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[index] -= ht05 * mOptParameter.k[i][j] * mOptParameter.z[i][j] * w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0] = 0.0; b2[offset+0] += aa_lambda_ht__hy; c2[offset+0] += m_aa_ht05__hyhy;
                c2[offset+M] = 0.0; b2[offset+M] += aa_lambda_ht__hy; a2[offset+M] += m_aa_ht05__hyhy;
                d2[offset+0] += aa_lambda_ht__hy*theta;
                d2[offset+M] += aa_lambda_ht__hy*theta;

                offset += M+1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2, x2, cols1_size);

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=0; m<=M; m++)
                {
                    u10[m][n] = x2[offset+m];
                }
                offset += M+1;
            }
        }

        /**************************************************** y direction apprx ***************************************************/

        if (use == true) f_add2Info(u05, u_info, 2*l-1, hx, hy, msnExtSpacePoints); f_layerInfo(u05, 2*l-1);
        if (use == true) f_add2Info(u10, u_info, 2*l+0, hx, hy, msnExtSpacePoints); f_layerInfo(u10, 2*l+0);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
            }
        }
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            u[m][n] = u00[m][n];
        }
    }

    if (rows1.size() != 0 && rows2.size() != 0)
    {
        for (unsigned int row=0; row < rows1_size; row++) free(w1[row]); free(w1);
        free(x1);
        free(d1);
        free(c1);
        free(b1);
        free(a1);
    }

    if (cols1.size() != 0 && cols2.size() != 0)
    {
        for (unsigned int col=0; col < cols1_size; col++) free(w2[col]); free(w2);
        free(x2);
        free(d2);
        free(c2);
        free(b2);
        free(a2);
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

    msnExtSpacePoints.clear();
    cntExtSpacePoints.clear();

    u00.clear();
    u10.clear();
}

auto Problem2PNeumann::solveBackwardIBVP(const DoubleMatrix &u, spif_vector &p_info, bool use, const spif_vector &u_info, const OptimizeParameterP &mOptParameter) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<const unsigned int>( dimX.sizeN() );
    const unsigned int M = static_cast<const unsigned int>( dimY.sizeN() );
    const unsigned int L = static_cast<const unsigned int>( time.sizeN() );

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a        = mEquParameter.a;
    const double alpha    = mEquParameter.alpha;
    const double lambda   = mEquParameter.lambda;
    //const double theta    = mEquParameter.theta;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;
#ifdef TIME_DISCRETE
    const unsigned int Nt = mEquParameter.Nt;
#endif

    const double m_aa_ht05__hxhx = -(0.5*a*a*ht)/(hx*hx);
    const double p_aa_ht__hxhx___alpha_ht05 = +1.0 + (a*a*ht)/(hx*hx) + 0.5*alpha*ht;
    const double p_aa_ht05__hyhy = +(0.5*a*a*ht)/(hy*hy);

    const double m_aa_ht05__hyhy = -(0.5*a*a*ht)/(hy*hy);
    const double p_aa_ht__hyhy___alpha_ht05 = +1.0 + (a*a*ht)/(hy*hy) + 0.5*alpha*ht;
    const double p_aa_ht05__hxhx = +(0.5*a*a*ht)/(hx*hx);

    const double ht05             = 0.5*ht;
    //const double alpha_ht05       = 0.5*alpha*ht;
    const double aa_lambda_ht__hx = (a*a*lambda*ht)/hx;
    const double aa_lambda_ht__hy = (a*a*lambda*ht)/hy;

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p05(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);




    //--------------------------------------------------------------------------------------------//
    std::vector<ExtendedSpacePointP> cntExtSpacePoints, msnExtSpacePoints;
    newDistributeDeltaGaussCntrl(mOptParameter.eta, cntExtSpacePoints, dimX, dimY);
    newDistributeDeltaGaussMsmnt(mOptParameter.xi,  msnExtSpacePoints, dimX, dimY);
    //----------------------------------------------------------------------------------------------//
    uint_vector rows0, rows1, rows2, cols0, cols1, cols2;
    b_findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, N, M, msnExtSpacePoints, cntExtSpacePoints);
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) b_prepareInfo(Nc, mOptParameter.eta, p_info, L);
    //-------------------------------------------- info --------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    b_initialLayers(p00, p_info, use, N, M, hx, hy, cntExtSpacePoints, u);
    //------------------------------------- initial conditions -------------------------------------//

    double *ax = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    double *bx = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    double *cx = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    double *dx = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    double *rx = static_cast<double*>( malloc(sizeof(double)*(N+1)) );
    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = m_aa_ht05__hxhx;
        bx[n] = p_aa_ht__hxhx___alpha_ht05;
        cx[n] = m_aa_ht05__hxhx;
    }
    ax[0] = 0.0; bx[0] += aa_lambda_ht__hx; cx[0] += m_aa_ht05__hxhx;
    cx[N] = 0.0; bx[N] += aa_lambda_ht__hx; ax[N] += m_aa_ht05__hxhx;

    double *ay = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    double *by = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    double *cy = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    double *dy = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    double *ry = static_cast<double*>( malloc(sizeof(double)*(M+1)) );
    for (unsigned int m=0; m<=M; m++)
    {
        ay[m] = m_aa_ht05__hyhy;
        by[m] = p_aa_ht__hyhy___alpha_ht05;
        cy[m] = m_aa_ht05__hyhy;
    }
    ay[0] = 0.0; by[0] += aa_lambda_ht__hy; cy[0] += m_aa_ht05__hyhy;
    cy[M] = 0.0; by[M] += aa_lambda_ht__hy; ay[M] += m_aa_ht05__hyhy;

    const unsigned int rows1_size = static_cast<const unsigned int>( rows1.size()*(N+1) );
    double *a1=nullptr, *b1=nullptr, *c1=nullptr, *d1=nullptr, *x1=nullptr, **w1=nullptr;
    if (rows1.size() != 0 && rows2.size() != 0)
    {
        a1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        b1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        c1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        d1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        x1 = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
        w1 = static_cast<double**> ( malloc(sizeof(double*)*rows1_size) );
        for (unsigned int row=0; row < rows1_size; row++) w1[row] = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
    }

    const unsigned int cols1_size = static_cast<const unsigned int>( cols1.size()*(M+1) );
    double *a2=nullptr, *b2=nullptr, *c2=nullptr, *d2=nullptr, *x2=nullptr, **w2=nullptr;
    if (cols1.size() != 0 && cols2.size() != 0)
    {
        a2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        b2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        c2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        d2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        x2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        w2 = static_cast<double**> ( malloc(sizeof(double*)*cols1_size) );
        for (unsigned int col=0; col < cols1_size; col++) w2[col] = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
    }

    SpaceNodePDE sn;

#ifdef TIME_DISCRETE
    std::vector<unsigned int> discrete_times;
    discrete_times.push_back(0);
    for (unsigned int s=0; s<Nt; s++)
    {
        const double tau = mOptParameter.tau[s];
        discrete_times.push_back(static_cast<unsigned int>(round(tau*L)));
    }
    discrete_times.push_back(L);
#endif

    for (unsigned int l=L-1; l != static_cast<unsigned int>(0)-1; l--)
    {

#ifdef TIME_DISCRETE
        unsigned int current_time = 0;  unsigned int current_indx = 0;
        for (unsigned int i=0; i<discrete_times.size(); i++)
        {
            if (discrete_times[i] <= l) current_indx = i;
            current_time = discrete_times[current_indx];
        }
        //printf("B %4u %4u\n", l, current_time);
#endif
        /**************************************************** x direction apprx ***************************************************/

        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = rows0.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n] = 0.0;
#ifdef OH1
                    if (m>0 && m<M)  dx[n] += p_aa_ht05__hyhy*(p00[m-1][n] - 2.0*p00[m][n]   + p00[m+1][n]);
                    else if (m == 0) dx[n] += p_aa_ht05__hyhy*(p00[0][n]   - 2.0*p00[1][n]   + p00[2][n]);
                    else if (m == M) dx[n] += p_aa_ht05__hyhy*(p00[M-2][n] - 2.0*p00[M-1][n] + p00[M][n]);
#else
                    if (m>0 && m<M)  dx[n] += p_aa_ht05__hyhy*(p00[m-1][n] - 2.0*p00[m][n]   + p00[m+1][n]);
                    else if (m == 0) dx[n] += p_aa_ht05__hyhy*(2.0*p00[0][n] - 5.0*p00[1][n]   + 4.0*p00[2][n]   - p00[3][n]);
                    else if (m == M) dx[n] += p_aa_ht05__hyhy*(2.0*p00[M][n] - 5.0*p00[M-1][n] + 4.0*p00[M-2][n] - p00[M-3][n]);
#endif


                    dx[n] += p00[m][n];
                }




                tomasAlgorithm(ax, bx, cx, dx, rx, N+1);
                for (unsigned int n=0; n<=N; n++) p05[m][n] = rx[n];
            }
        }

        if (rows1.size() != 0 && rows2.size() == 0)
        {
            //throw std::exception();
            //throw grid_exception("backward x1");

            double *_w05 = new double[No];

            double* _p05 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p05[i] = 0.0;

#ifdef TIME_DISCRETE
            if (l == current_time)
            {
                unsigned int taus = discrete_times[current_indx+0];
                unsigned int taup = discrete_times[current_indx+1];
                //                if (current_indx == discrete_times.size()-1)
                //                {
                //                    taup = L;
                //                }
                //                else
                //                {
                //                    taup = discrete_times[current_indx+1];
                //                }

                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfoP &pi = p_info[i];
                    double pi_in = 0.0;
                    pi_in += 0.5 * ht * pi.vl[2*taus];
                    for (unsigned int ln=taus+1; ln<=taup-1; ln++) pi_in += ht * pi.vl[2*ln];
                    pi_in += 0.5 * ht * pi.vl[2*taup];

                    _p05[i] = pi_in * (1.0/ht);
                }
            }
#else

            for (unsigned int i=0; i<Nc; i++)
            {
                const ExtendedSpacePointP &extendedSpacePoint = cntExtSpacePoints.at(i);
                const std::vector<ExtendedSpacePointNodeP> &nodes = extendedSpacePoint.nodes;
                const unsigned int nodes_size = static_cast<const unsigned int>( nodes.size() );
                for (unsigned int ni=0; ni<nodes_size; ni++)
                {
                    const ExtendedSpacePointNodeP &node = nodes.at(ni);
                    unsigned int node_nx = static_cast<unsigned int>(node.nx);
                    unsigned int node_ny = static_cast<unsigned int>(node.ny);
                    _p05[i] += p05[node_ny][node_nx] * (node.w * (hx*hy));
                }
            }
#endif

            for (unsigned int j=0; j<No; j++)
            {
                _w05[j] = 0.0;
                for (unsigned int i=0; i<Nc; i++)
                {
                    _w05[j] += mOptParameter.k[i][j] * ( _p05[i] + 2.0*r1*gpi(i, 2*l+1, u_info, mOptParameter)*sgn(g0i(i, 2*l+1, u_info, mOptParameter)));
                }
            }

            delete [] _p05;

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n] = 0.0;
#ifdef OH1
                    if (m>0 && m<M)  dx[n] += p_aa_ht05__hyhy*(p00[m-1][n] - 2.0*p00[m][n]   + p00[m+1][n]);
                    else if (m == 0) dx[n] += p_aa_ht05__hyhy*(p00[0][n]   - 2.0*p00[1][n]   + p00[2][n]);
                    else if (m == M) dx[n] += p_aa_ht05__hyhy*(p00[M-2][n] - 2.0*p00[M-1][n] + p00[M][n]);
#else
                    if (m>0 && m<M)  dx[n] += p_aa_ht05__hyhy*(p00[m-1][n] - 2.0*p00[m][n]   + p00[m+1][n]);
                    else if (m == 0) dx[n] += p_aa_ht05__hyhy*(2.0*p00[0][n] - 5.0*p00[1][n]   + 4.0*p00[2][n]   - p00[3][n]);
                    else if (m == M) dx[n] += p_aa_ht05__hyhy*(2.0*p00[M][n] - 5.0*p00[M-1][n] + 4.0*p00[M-2][n] - p00[M-3][n]);
#endif

                    dx[n] += p00[m][n];

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int j=0; j<No; j++)
                    {
                        const ExtendedSpacePointP &extendedSpacePoint = msnExtSpacePoints.at(j);
                        if (extendedSpacePoint.contains(sn))
                        {
                            const std::vector<ExtendedSpacePointNodeP> &nodes = extendedSpacePoint.nodes;
                            unsigned int nodes_size = static_cast<const unsigned int>( nodes.size() );
                            for (unsigned int nj=0; nj<nodes_size; nj++)
                            {
                                const ExtendedSpacePointNodeP &node = nodes.at(nj);
                                if (node.equals(sn)) dx[n] += ht05 * _w05[j] * node.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }




                tomasAlgorithm(ax, bx, cx, dx, rx, N+1);
                for (unsigned int n=0; n<=N; n++) p05[m][n] = rx[n];
            }

            delete [] _w05;
        }

        if (rows1.size() != 0 && rows2.size() != 0)
        {
            //throw std::exception();
            //throw grid_exception("backward x2");

            for (unsigned int m=0; m<rows1_size; m++) for (unsigned int n=0; n<rows1_size; n++) w1[m][n] = 0.0;

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    const unsigned int index = offset+n;
                    d1[index] = 0.0;
#ifdef OH1
                    if (m>0 && m<M)  d1[index] += p_aa_ht05__hyhy*(p00[m-1][n] - 2.0*p00[m][n]   + p00[m+1][n]);
                    else if (m == 0) d1[index] += p_aa_ht05__hyhy*(p00[0][n]   - 2.0*p00[1][n]   + p00[2][n]);
                    else if (m == M) d1[index] += p_aa_ht05__hyhy*(p00[M-2][n] - 2.0*p00[M-1][n] + p00[M][n]);
#else
                    if (m>0 && m<M)  d1[index] += p_aa_ht05__hyhy*(p00[m-1][n] - 2.0*p00[m][n]   + p00[m+1][n]);
                    else if (m == 0) d1[index] += p_aa_ht05__hyhy*(2.0*p00[0][n] - 5.0*p00[1][n]   + 4.0*p00[2][n]   - p00[3][n]);
                    else if (m == M) d1[index] += p_aa_ht05__hyhy*(2.0*p00[M][n] - 5.0*p00[M-1][n] + 4.0*p00[M-2][n] - p00[M-3][n]);
#endif

                    d1[index] += p00[m][n];

                    a1[index] = m_aa_ht05__hxhx;
                    b1[index] = p_aa_ht__hxhx___alpha_ht05;
                    c1[index] = m_aa_ht05__hxhx;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int j=0; j<No; j++)
                    {
                        const ExtendedSpacePointP &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                        if (mExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNodeP> &nodes1 = mExtendedSpacePoint.nodes;
                            for (unsigned int nj=0; nj<nodes1.size(); nj++)
                            {
                                const ExtendedSpacePointNodeP &node1 = nodes1.at(nj);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                const ExtendedSpacePointP &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                                const std::vector<ExtendedSpacePointNodeP> &nodes2 = cExtendedSpacePoint.nodes;
                                for (unsigned int ni=0; ni<nodes2.size(); ni++)
                                {
                                    const ExtendedSpacePointNodeP &node2 = nodes2.at(ni);
                                    unsigned int node2_nx = static_cast<unsigned int>(node2.nx);
                                    unsigned int node2_ny = static_cast<unsigned int>(node2.ny);

                                    bool found = false;
                                    for (unsigned int rs=0; rs<rows1.size(); rs++)
                                    {
                                        if (static_cast<unsigned int>(node2.ny) == rows1[rs])
                                        {
                                            found = true;
                                            w1[index][rs*(N+1)+node2_nx] -= ht05 * mOptParameter.k[i][j] * (node2.w * (hx*hy)) * w;
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d1[index] += ht05 * mOptParameter.k[i][j] * p05[node2_ny][node2_nx] * (node2.w * (hx*hy)) * w;
                                    }
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1[index] += 2.0 * r1 * ht05 *  mOptParameter.k[i][j] * gpi(i, 2*l+1, u_info, mOptParameter)*sgn(g0i(i, 2*l+1, u_info, mOptParameter)) * w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1[offset+0] = 0.0; b1[offset+0] += aa_lambda_ht__hx; c1[offset+0] += m_aa_ht05__hxhx;
                c1[offset+N] = 0.0; b1[offset+N] += aa_lambda_ht__hx; a1[offset+N] += m_aa_ht05__hxhx;

                offset += N+1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w1, x1, rows1_size);

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=0; n<=N; n++)
                {
                    p05[m][n] = x1[offset+n];
                }
                offset += N+1;
            }
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        if (cols0.size() != 0)
        {
            for (unsigned int col=0; col<cols0.size(); col++)
            {
                unsigned int n = cols0.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m] = 0.0;
#ifdef OH1
                    if (n>0 && n<N) dy[m] += p_aa_ht05__hxhx*(p05[m][n-1] - 2.0*p05[m][n]   + p05[m][n+1]);
                    else if (n==0)  dy[m] += p_aa_ht05__hxhx*(p05[m][0]   - 2.0*p05[m][1]   + p05[m][2]);
                    else if (n==N)  dy[m] += p_aa_ht05__hxhx*(p05[m][N-2] - 2.0*p05[m][N-1] + p05[m][N]);
#else
                    if (n>0 && n<N)  dy[m] += p_aa_ht05__hxhx*(p05[m][n-1] - 2.0*p05[m][n]   + p05[m][n+1]);
                    else if (n == 0) dy[m] += p_aa_ht05__hxhx*(2.0*p05[m][0] - 5.0*p05[m][1]   + 4.0*p05[m][2]   - p05[m][3]);
                    else if (n == N) dy[m] += p_aa_ht05__hxhx*(2.0*p05[m][N] - 5.0*p05[m][N-1] + 4.0*p05[m][N-2] - p05[m][N-3]);
#endif

                    dy[m] += p05[m][n];
                }




                tomasAlgorithm(ay, by, cy, dy, ry, M+1);
                for (unsigned int m=0; m<=M; m++) p10[m][n] = ry[m];
            }
        }

        if (cols1.size() != 0 && cols2.size() == 0)
        {
            //throw std::exception();
            //throw grid_exception("backward y1");

            double *_w10 = new double[No];

            double* _p10 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p10[i] = 0.0;

#ifdef TIME_DISCRETE
            if (l == current_time)
            {
                unsigned int taus = discrete_times[current_indx+0];
                unsigned int taup = discrete_times[current_indx+1];
                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfoP &pi = p_info[i];
                    double pi_in = 0.0;
                    pi_in += 0.5 * ht * pi.vl[2*taus];
                    for (unsigned int ln=taus+1; ln<=taup-1; ln++) pi_in += ht * pi.vl[2*ln];
                    pi_in += 0.5 * ht * pi.vl[2*taup];

                    _p10[i] = pi_in * (1.0/ht);
                }


                //                unsigned int start = discrete_times[current_indx];
                //                unsigned int end   = 0;
                //                if (current_indx == discrete_times.size()-1)
                //                {
                //                    end = L;
                //                }
                //                else
                //                {
                //                    end = discrete_times[current_indx+1];
                //                }
                //                for (unsigned int i=0; i<Nc; i++)
                //                {
                //                    for (unsigned int ln=start; ln<end; ln++) _p10[i] += ht*p_info[i].vl[2*ln];
                //                    _p10[i] *= (1.0/ht);
                //                }
            }
#else
            for (unsigned int i=0; i<Nc; i++)
            {
                const ExtendedSpacePointP &extendedSpacePoint = cntExtSpacePoints.at(i);
                const std::vector<ExtendedSpacePointNodeP> &nodes = extendedSpacePoint.nodes;
                const unsigned int nodes_size = static_cast<const unsigned int>( nodes.size() );
                for (unsigned int ni=0; ni<nodes_size; ni++)
                {
                    const ExtendedSpacePointNodeP &node = nodes.at(ni);
                    unsigned int node_nx = static_cast<unsigned int>(node.nx);
                    unsigned int node_ny = static_cast<unsigned int>(node.ny);
                    _p10[i] += p10[node_ny][node_nx] * (node.w * (hx*hy));
                }
            }
#endif

            for (unsigned int j=0; j<No; j++)
            {
                _w10[j] = 0.0;
                for (unsigned int i=0; i<Nc; i++)
                {
                    _w10[j] += mOptParameter.k[i][j] * (_p10[i] + 2.0*r1*gpi(i, 2*l+0, u_info, mOptParameter)*sgn(g0i(i, 2*l+0, u_info, mOptParameter)));
                }
            }
            delete [] _p10;

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m] = 0.0;
#ifdef OH1
                    if (n>0 && n<N) dy[m] += p_aa_ht05__hxhx*(p05[m][n-1] - 2.0*p05[m][n]   + p05[m][n+1]);
                    else if (n==0)  dy[m] += p_aa_ht05__hxhx*(p05[m][0]   - 2.0*p05[m][1]   + p05[m][2]);
                    else if (n==N)  dy[m] += p_aa_ht05__hxhx*(p05[m][N-2] - 2.0*p05[m][N-1] + p05[m][N]);
#else
                    if (n>0 && n<N)  dy[m] += p_aa_ht05__hxhx*(p05[m][n-1] - 2.0*p05[m][n]   + p05[m][n+1]);
                    else if (n == 0) dy[m] += p_aa_ht05__hxhx*(2.0*p05[m][0] - 5.0*p05[m][1]   + 4.0*p05[m][2]   - p05[m][3]);
                    else if (n == N) dy[m] += p_aa_ht05__hxhx*(2.0*p05[m][N] - 5.0*p05[m][N-1] + 4.0*p05[m][N-2] - p05[m][N-3]);
#endif


                    dy[m] += p05[m][n];

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int j=0; j<No; j++)
                    {
                        const ExtendedSpacePointP &extendedSpacePoint = msnExtSpacePoints.at(j);
                        if (extendedSpacePoint.contains(sn))
                        {
                            const std::vector<ExtendedSpacePointNodeP> &nodes = extendedSpacePoint.nodes;
                            const unsigned int nodes_size = static_cast<const unsigned int>( nodes.size() );
                            for (unsigned int nj=0; nj<nodes_size; nj++)
                            {
                                const ExtendedSpacePointNodeP &node = nodes.at(nj);
                                if (node.equals(sn)) dy[m] += ht05 * _w10[j] * node.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }




                tomasAlgorithm(ay, by, cy, dy, ry, M+1);
                for (unsigned int m=0; m<=M; m++) p10[m][n] = ry[m];
            }

            delete [] _w10;
        }

        if (cols1.size() != 0 && cols2.size() != 0)
        {
            //throw std::exception();
            //throw grid_exception("backward y2");

            for (unsigned int m=0; m<cols1_size; m++) for (unsigned int n=0; n<cols1_size; n++) w2[m][n] = 0.0;

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    const unsigned int index = offset+m;
                    d2[index] = 0.0;
#ifdef OH1
                    if (n>0 && n<N) d2[index] = p_aa_ht05__hxhx*(p05[m][n-1] - 2.0*p05[m][n]   + p05[m][n+1]);
                    else if (n==0)  d2[index] = p_aa_ht05__hxhx*(p05[m][0]   - 2.0*p05[m][1]   + p05[m][2]);
                    else if (n==N)  d2[index] = p_aa_ht05__hxhx*(p05[m][N-2] - 2.0*p05[m][N-1] + p05[m][N]);
#else
                    if (n>0 && n<N)  d2[index] += p_aa_ht05__hxhx*(p05[m][n-1] - 2.0*p05[m][n]   + p05[m][n+1]);
                    else if (n == 0) d2[index] += p_aa_ht05__hxhx*(2.0*p05[m][0] - 5.0*p05[m][1]   + 4.0*p05[m][2]   - p05[m][3]);
                    else if (n == N) d2[index] += p_aa_ht05__hxhx*(2.0*p05[m][N] - 5.0*p05[m][N-1] + 4.0*p05[m][N-2] - p05[m][N-3]);
#endif

                    d2[index] += p05[m][n];

                    a2[index] = m_aa_ht05__hyhy;
                    b2[index] = p_aa_ht__hyhy___alpha_ht05;
                    c2[index] = m_aa_ht05__hyhy;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int j=0; j<No; j++)
                    {
                        const ExtendedSpacePointP &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                        if (mExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNodeP> &nodes1 = mExtendedSpacePoint.nodes;
                            for (unsigned int nj=0; nj<nodes1.size(); nj++)
                            {
                                const ExtendedSpacePointNodeP &node1 = nodes1.at(nj);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                const ExtendedSpacePointP &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                                const std::vector<ExtendedSpacePointNodeP> &nodes2 = cExtendedSpacePoint.nodes;
                                for (unsigned int ni=0; ni<nodes2.size(); ni++)
                                {
                                    const ExtendedSpacePointNodeP &node2 = nodes2.at(ni);
                                    unsigned int node2_nx = static_cast<unsigned int>(node2.nx);
                                    unsigned int node2_ny = static_cast<unsigned int>(node2.ny);

                                    bool found = false;
                                    for (unsigned int cs=0; cs<cols1.size(); cs++)
                                    {
                                        if (static_cast<unsigned int>(node2.nx) == cols1[cs])
                                        {
                                            found = true;
                                            w2[index][cs*(M+1)+node2_ny] -= ht05 * mOptParameter.k[i][i] * (node2.w * (hx*hy)) * w;
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d2[index] += ht05 * mOptParameter.k[i][i] * p10[node2_ny][node2_nx] * (node2.w * (hx*hy)) * w;
                                    }
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d2[index] += 2.0 * r1 * ht05 *  mOptParameter.k[i][j] * gpi(i, 2*l, u_info, mOptParameter)*sgn(g0i(i, 2*l, u_info, mOptParameter)) * w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0] = 0.0; b2[offset+0] += aa_lambda_ht__hy; c2[offset+0] += m_aa_ht05__hyhy;
                c2[offset+M] = 0.0; b2[offset+M] += aa_lambda_ht__hy; a2[offset+M] += m_aa_ht05__hyhy;

                offset += M+1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2, x2, cols1_size);

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=0; m<=M; m++)
                {
                    p10[m][n] = x2[offset+m];
                }
                offset += M+1;
            }
        }

        /**************************************************** y direction apprx ***************************************************/

        if (use == true) b_add2Info(p05, p_info, 2*l+1, hx, hy, cntExtSpacePoints); b_layerInfo(p05, 2*l+1);
        if (use == true) b_add2Info(p10, p_info, 2*l+0, hx, hy, cntExtSpacePoints); b_layerInfo(p10, 2*l+0);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
            }
        }
    }

    if (rows1.size() != 0 && rows2.size() != 0)
    {
        for (unsigned int row=0; row < rows1_size; row++) free(w1[row]); free(w1);
        free(x1);
        free(d1);
        free(c1);
        free(b1);
        free(a1);
    }

    if (cols1.size() != 0 && cols2.size() != 0)
    {
        for (unsigned int col=0; col < cols1_size; col++) free(w2[col]); free(w2);
        free(x2);
        free(d2);
        free(c2);
        free(b2);
        free(a2);
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

    msnExtSpacePoints.clear();
    cntExtSpacePoints.clear();

    p00.clear();
    p05.clear();
    p10.clear();
}

auto Problem2PNeumann::f_initialLayers(DoubleMatrix &u00, spif_vector &u_info, bool use, unsigned int N, unsigned int M,
                                       double hx, double hy, const std::vector<ExtendedSpacePointP> &msnExtSpacePoints) const -> void
{
    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u00[m][n] = f_initial(sn);
        }
    }

    if (use == true) f_add2Info(u00, u_info, 0, hx, hy, msnExtSpacePoints);
    f_layerInfo(u00, 0);
}

auto Problem2PNeumann::b_initialLayers(DoubleMatrix &p00, spif_vector &p_info, bool use, unsigned int N, unsigned int M,
                                       double hx, double hy, const std::vector<ExtendedSpacePointP> &cntExtSpacePoints,
                                       const DoubleMatrix &u) const -> void
{
    const unsigned int L = static_cast<const unsigned int>( mtimeDimension.sizeN() );

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p00[m][n] = b_initial(sn, u);
        }
    }

    if (use == true) b_add2Info(p00, p_info, 2*L, hx, hy, cntExtSpacePoints);
    b_layerInfo(p00, 2*L);
}

auto Problem2PNeumann::f_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2,
                                      uint_vector &cols0, uint_vector &cols1, uint_vector &cols2,
                                      unsigned int N, unsigned int M,
                                      const std::vector<ExtendedSpacePointP> &cntExtSpacePoints,
                                      const std::vector<ExtendedSpacePointP> &msnExtSpacePoints) const -> void
{
    for (unsigned int m=0; m<=M; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for(std::vector<ExtendedSpacePointP>::const_iterator csp_it=cntExtSpacePoints.begin(); csp_it != cntExtSpacePoints.end(); csp_it++)
        {
            const ExtendedSpacePointP &cxsp = *csp_it;
            const std::vector<ExtendedSpacePointNodeP> &c_nodes = cxsp.nodes;
            for (std::vector<ExtendedSpacePointNodeP>::const_iterator cnode_it=c_nodes.begin(); cnode_it != c_nodes.end(); cnode_it++)
            {
                const ExtendedSpacePointNodeP &cnode = *cnode_it;
                if (static_cast<unsigned int>(cnode.ny) == m)
                {
                    found1 = true;
                    for(std::vector<ExtendedSpacePointP>::const_iterator msp_it=msnExtSpacePoints.begin(); msp_it != msnExtSpacePoints.end(); msp_it++)
                    {
                        const ExtendedSpacePointP &mxsp = *msp_it;
                        const std::vector<ExtendedSpacePointNodeP> &mnodes = mxsp.nodes;
                        for (std::vector<ExtendedSpacePointNodeP>::const_iterator mnode_it=mnodes.begin(); mnode_it != mnodes.end(); mnode_it++)
                        {
                            const ExtendedSpacePointNodeP &mnode = *mnode_it;
                            if (static_cast<unsigned int>(mnode.ny) == m)
                            {
                                found2 = true;
                                break;
                            }
                        }
                        if (found2) break;
                    }
                    break;
                }
            }
            if (found1) break;
        }
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
    }

    for (unsigned int n=0; n<=N; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for(std::vector<ExtendedSpacePointP>::const_iterator csp_it=cntExtSpacePoints.begin(); csp_it != cntExtSpacePoints.end(); csp_it++)
        {
            const ExtendedSpacePointP &cxsp = *csp_it;
            const std::vector<ExtendedSpacePointNodeP> &c_nodes = cxsp.nodes;
            for (std::vector<ExtendedSpacePointNodeP>::const_iterator cnode_it=c_nodes.begin(); cnode_it != c_nodes.end(); cnode_it++)
            {
                const ExtendedSpacePointNodeP &cnode = *cnode_it;
                if (static_cast<unsigned int>(cnode.nx) == n)
                {
                    found1 = true;
                    for(std::vector<ExtendedSpacePointP>::const_iterator msp_it=msnExtSpacePoints.begin(); msp_it != msnExtSpacePoints.end(); msp_it++)
                    {
                        const ExtendedSpacePointP &mxsp = *msp_it;
                        const std::vector<ExtendedSpacePointNodeP> &mnodes = mxsp.nodes;
                        for (std::vector<ExtendedSpacePointNodeP>::const_iterator mnode_it=mnodes.begin(); mnode_it != mnodes.end(); mnode_it++)
                        {
                            const ExtendedSpacePointNodeP &mnode = *mnode_it;
                            if (static_cast<unsigned int>(mnode.nx) == n)
                            {
                                found2 = true;
                                break;
                            }
                        }
                        if (found2) break;
                    }
                    break;
                }
            }
            if (found1) break;
        }
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
    }
}

auto Problem2PNeumann::b_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2,
                                      uint_vector &cols0, uint_vector &cols1, uint_vector &cols2,
                                      unsigned int N, unsigned int M,
                                      const std::vector<ExtendedSpacePointP> &msnExtSpacePoints,
                                      const std::vector<ExtendedSpacePointP> &cntExtSpacePoints) const -> void
{
    for (unsigned int m=0; m<=M; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for(std::vector<ExtendedSpacePointP>::const_iterator csp_it=msnExtSpacePoints.begin(); csp_it != msnExtSpacePoints.end(); csp_it++)
        {
            const ExtendedSpacePointP &cxsp = *csp_it;
            const std::vector<ExtendedSpacePointNodeP> &c_nodes = cxsp.nodes;
            for (std::vector<ExtendedSpacePointNodeP>::const_iterator cnode_it=c_nodes.begin(); cnode_it != c_nodes.end(); cnode_it++)
            {
                const ExtendedSpacePointNodeP &cnode = *cnode_it;
                if (static_cast<unsigned int>(cnode.ny) == m)
                {
                    found1 = true;
                    for(std::vector<ExtendedSpacePointP>::const_iterator msp_it=cntExtSpacePoints.begin(); msp_it != cntExtSpacePoints.end(); msp_it++)
                    {
                        const ExtendedSpacePointP &mxsp = *msp_it;
                        const std::vector<ExtendedSpacePointNodeP> &mnodes = mxsp.nodes;
                        for (std::vector<ExtendedSpacePointNodeP>::const_iterator mnode_it=mnodes.begin(); mnode_it != mnodes.end(); mnode_it++)
                        {
                            const ExtendedSpacePointNodeP &mnode = *mnode_it;
                            if (static_cast<unsigned int>(mnode.ny) == m)
                            {
                                found2 = true;
                                break;
                            }
                        }
                        if (found2) break;
                    }
                    break;
                }
            }
            if (found1) break;
        }
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
    }

    for (unsigned int n=0; n<=N; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for(std::vector<ExtendedSpacePointP>::const_iterator csp_it=msnExtSpacePoints.begin(); csp_it != msnExtSpacePoints.end(); csp_it++)
        {
            const ExtendedSpacePointP &cxsp = *csp_it;
            const std::vector<ExtendedSpacePointNodeP> &c_nodes = cxsp.nodes;
            for (std::vector<ExtendedSpacePointNodeP>::const_iterator cnode_it=c_nodes.begin(); cnode_it != c_nodes.end(); cnode_it++)
            {
                const ExtendedSpacePointNodeP &cnode = *cnode_it;
                if (static_cast<unsigned int>(cnode.nx) == n)
                {
                    found1 = true;
                    for(std::vector<ExtendedSpacePointP>::const_iterator msp_it=cntExtSpacePoints.begin(); msp_it != cntExtSpacePoints.end(); msp_it++)
                    {
                        const ExtendedSpacePointP &mxsp = *msp_it;
                        const std::vector<ExtendedSpacePointNodeP> &mnodes = mxsp.nodes;
                        for (std::vector<ExtendedSpacePointNodeP>::const_iterator mnode_it=mnodes.begin(); mnode_it != mnodes.end(); mnode_it++)
                        {
                            const ExtendedSpacePointNodeP &mnode = *mnode_it;
                            if (static_cast<unsigned int>(mnode.nx) == n)
                            {
                                found2 = true;
                                break;
                            }
                        }
                        if (found2) break;
                    }
                    break;
                }
            }
            if (found1) break;
        }
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
    }
}

auto Problem2PNeumann::f_prepareInfo(unsigned int No, const std::vector<SpacePoint> &points, spif_vector &u_info, unsigned int L) const -> void
{
    u_info.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        SpacePointInfoP &inf = u_info[j];
        const SpacePoint &sp = points[j];
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(2*L+1);
    }
}

auto Problem2PNeumann::b_prepareInfo(unsigned int Nc, const std::vector<SpacePoint> &points, spif_vector &p_info, unsigned int L) const -> void
{
    p_info.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        SpacePointInfoP &inf = p_info[i];
        const SpacePoint &sp = points[i];
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(2*L+1);
    }
}

auto Problem2PNeumann::f_add2Info(const DoubleMatrix &u, spif_vector &u_info, unsigned int ln, double hx, double hy, const std::vector<ExtendedSpacePointP> &extMsmnts, int method) const -> void
{
    if (method == 1 || method == 2 || method == 4)
    {
        unsigned int No = static_cast<unsigned int>(extMsmnts.size());
        for (unsigned int j=0; j<No; j++)
        {
            const ExtendedSpacePointP &xsp = extMsmnts.at(j);
            SpacePointInfoP &ui = u_info[j];
            const unsigned int nodes_size = static_cast<const unsigned int>( xsp.nodes.size() );
            for (unsigned int i=0; i<nodes_size; i++)
            {
                const ExtendedSpacePointNodeP &node = xsp.nodes.at(i);
                const unsigned int nx = static_cast<const unsigned int>(node.nx);
                const unsigned int ny = static_cast<const unsigned int>(node.ny);
                ui.vl[ln] += u[ny][nx] * (node.w * (hx*hy));
                if (node.isCenter)
                {
                    const unsigned int rx = static_cast<const unsigned int>(xsp.rx);
                    const unsigned int ry = static_cast<const unsigned int>(xsp.ry);

                    ui.dx[ln] = (u[ry][rx+1] - u[ry][rx-1])/(2.0*hx);
                    ui.dy[ln] = (u[ry+1][rx] - u[ry-1][rx])/(2.0*hy);

                    ui.dx[ln] += ((xsp.x-rx*hx)/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
                    ui.dy[ln] += ((xsp.y-ry*hy)/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);

                    ui.dxx[ln] = (1.0/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
                    ui.dyy[ln] = (1.0/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);
                }
            }
        }
    }
}

auto Problem2PNeumann::b_add2Info(const DoubleMatrix &p, spif_vector &p_info, unsigned int ln, double hx, double hy, const std::vector<ExtendedSpacePointP> &extCntrls, int method) const -> void
{
    if (method == 1 || method == 2 || method == 4)
    {
        unsigned int Nc = static_cast<unsigned int>(extCntrls.size());
        for (unsigned int i=0; i<Nc; i++)
        {
            const ExtendedSpacePointP &xsp = extCntrls.at(i);
            SpacePointInfoP &pi = p_info[i];
            const unsigned int nodes_size = static_cast<const unsigned int>( xsp.nodes.size() );
            for (unsigned int i=0; i<nodes_size; i++)
            {
                const ExtendedSpacePointNodeP &node = xsp.nodes.at(i);
                const unsigned int nx = static_cast<const unsigned int>(node.nx);
                const unsigned int ny = static_cast<const unsigned int>(node.ny);
                pi.vl[ln] += p[ny][nx] * (node.w * (hx*hy));
                if (node.isCenter)
                {
                    const unsigned int rx = static_cast<const unsigned int>(xsp.rx);
                    const unsigned int ry = static_cast<const unsigned int>(xsp.ry);

                    pi.dx[ln] = (p[ry][rx+1] - p[ry][rx-1])/(2.0*hx);
                    pi.dy[ln] = (p[ry+1][rx] - p[ry-1][rx])/(2.0*hy);

                    pi.dx[ln] += ((xsp.x-rx*hx)/(hx*hx))*(p[ry][rx+1] - 2.0*p[ry][rx] + p[ry][rx-1]);
                    pi.dy[ln] += ((xsp.y-ry*hy)/(hy*hy))*(p[ry+1][rx] - 2.0*p[ry][rx] + p[ry-1][rx]);
                }
            }
        }
    }
}

auto Problem2PNeumann::f_initial(const SpaceNodePDE &) const -> double
{
    return mEquParameter.phi;
}

auto Problem2PNeumann::b_initial(const SpaceNodePDE &sn, const DoubleMatrix &u) const -> double
{
    return -2.0*mu(sn.x, sn.y)*(u[sn.j][sn.i]-U[sn.j][sn.i]);
}

auto Problem2PNeumann::f_layerInfo(const DoubleMatrix &u, unsigned int ln) const -> void
{
    C_UNUSED(u);
    C_UNUSED(ln);
    //    QPixmap pic;
    //    visualizeMatrixHeat(u, u.min(), u.max(), pic);
    //    pic.save("images/problem2P/f/100/pic_"+QString("%1").arg(ln)+".png", "PNG");
}

auto Problem2PNeumann::b_layerInfo(const DoubleMatrix &p, unsigned int ln) const -> void
{
    C_UNUSED(p);
    C_UNUSED(ln);
    //    QPixmap pic;
    //    visualizeMatrixHeat(p, p.min(), p.max(), pic);
    //    pic.save("images/problem2P/b/100/pic_"+QString("%1").arg(ln)+".png", "PNG");
}

auto Problem2PNeumann::newDistributeDeltaGaussCntrl(const std::vector<SpacePoint> &cntrls,
                                                    std::vector<ExtendedSpacePointP> &extCntrls,
                                                    const Dimension &dimX, const Dimension &dimY) const -> void
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.sizeN() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.sizeN() );
    unsigned int Nc = static_cast<unsigned int> ( cntrls.size() );

    extCntrls.clear();
    extCntrls.resize(Nc);

    int k = 4;
    double sigmaX = hx;
    double sigmaY = hy;

    for (unsigned int c=0; c<Nc; c++)
    {
        const SpacePoint &cntrl = cntrls.at(c);
        ExtendedSpacePointP &extCntrl = extCntrls.at(c);

        extCntrl.x = cntrl.x;
        extCntrl.y = cntrl.y;
        extCntrl.rx = static_cast<int> ( round(extCntrl.x*Nx) );
        extCntrl.ry = static_cast<int> ( round(extCntrl.y*Ny) );
        extCntrl.k = k;
        extCntrl.minX = extCntrl.rx - extCntrl.k;
        extCntrl.maxX = extCntrl.rx + extCntrl.k;
        extCntrl.minY = extCntrl.ry - extCntrl.k;
        extCntrl.maxY = extCntrl.ry + extCntrl.k;

        double sumX = 0.0;
        for (int n=extCntrl.minX; n<=extCntrl.maxX; n++) sumX += exp(-((n*hx-cntrl.x)*(n*hx-cntrl.x))/(2.0*sigmaX*sigmaX));
        sumX *= hx;

        double sumY = 0.0;
        for (int m=extCntrl.minY; m<=extCntrl.maxY; m++) sumY += exp(-((m*hy-cntrl.y)*(m*hy-cntrl.y))/(2.0*sigmaY*sigmaY));
        sumY *= hy;

        double sigma = (sumX*sumY) / (2.0*M_PI);
        double factor = 1.0/((2.0*M_PI)*sigma);

        for (int m=extCntrl.minY; m<=extCntrl.maxY; m++)
        {
            for (int n=extCntrl.minX; n<=extCntrl.maxX; n++)
            {
                ExtendedSpacePointNodeP node;
                node.nx = n; node.x = n*hx;
                node.ny = m; node.y = m*hy;
                node.w = factor*exp(-0.5*(((node.x-cntrl.x)*(node.x-cntrl.x))/(sigmaX*sigmaX)+((node.y-cntrl.y)*(node.y-cntrl.y))/(sigmaY*sigmaY)));
                node.isCenter = ( m==extCntrl.ry && n==extCntrl.rx );
                extCntrl.nodes.push_back(node);
            }
        }
    }
}

auto Problem2PNeumann::newDistributeDeltaGaussMsmnt(const std::vector<SpacePoint> &msmnts,
                                                    std::vector<ExtendedSpacePointP> &extMsmnts,
                                                    const Dimension &dimX, const Dimension &dimY) const -> void
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.sizeN() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.sizeN() );
    unsigned int Nc = static_cast<unsigned int> ( msmnts.size() );

    extMsmnts.clear();
    extMsmnts.resize(Nc);

    int k = 4;
    double sigmaX = hx;
    double sigmaY = hy;

    for (unsigned int c=0; c<Nc; c++)
    {
        const SpacePoint &msmnt = msmnts.at(c);
        ExtendedSpacePointP &extMsmnt = extMsmnts.at(c);

        extMsmnt.x = msmnt.x;
        extMsmnt.y = msmnt.y;
        extMsmnt.rx = static_cast<int> ( round(extMsmnt.x*Nx) );
        extMsmnt.ry = static_cast<int> ( round(extMsmnt.y*Ny) );
        extMsmnt.k = k;
        extMsmnt.minX = extMsmnt.rx - extMsmnt.k;
        extMsmnt.maxX = extMsmnt.rx + extMsmnt.k;
        extMsmnt.minY = extMsmnt.ry - extMsmnt.k;
        extMsmnt.maxY = extMsmnt.ry + extMsmnt.k;

        double sumX = 0.0;
        for (int n=extMsmnt.minX; n<=extMsmnt.maxX; n++) sumX += exp(-((n*hx-msmnt.x)*(n*hx-msmnt.x))/(2.0*sigmaX*sigmaX));
        sumX *= hx;

        double sumY = 0.0;
        for (int m=extMsmnt.minY; m<=extMsmnt.maxY; m++) sumY += exp(-((m*hy-msmnt.y)*(m*hy-msmnt.y))/(2.0*sigmaY*sigmaY));
        sumY *= hy;

        double sigma = (sumX*sumY) / (2.0*M_PI);
        double factor = 1.0/((2.0*M_PI)*sigma);

        for (int m=extMsmnt.minY; m<=extMsmnt.maxY; m++)
        {
            for (int n=extMsmnt.minX; n<=extMsmnt.maxX; n++)
            {
                ExtendedSpacePointNodeP node;
                node.nx = n; node.x = n*hx;
                node.ny = m; node.y = m*hy;
                node.w = factor*exp(-0.5*(((node.x-msmnt.x)*(node.x-msmnt.x))/(sigmaX*sigmaX)+((node.y-msmnt.y)*(node.y-msmnt.y))/(sigmaY*sigmaY)));
                node.isCenter = ( m==extMsmnt.ry && n==extMsmnt.rx );
                extMsmnt.nodes.push_back(node);
            }
        }
    }
}

auto Problem2PNeumann::VectorToPrm(const DoubleVector &pv, OptimizeParameterP &prm) const -> void
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;
#ifdef TIME_DISCRETE
    unsigned int Nt = mEquParameter.Nt;
#endif

    unsigned int index = 0;

    prm.k.clear();
    prm.k.resize(Nc, No);

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            prm.k[i][j] = pv[index]; index++;
        }
    }

    prm.z.clear();
    prm.z.resize(Nc, No);

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            prm.z[i][j] = pv[index]; index++;
        }
    }

    prm.xi.clear();
    prm.xi.resize(No);

    for (unsigned int j=0; j<No; j++)
    {
        prm.xi[j].x = pv[index]; index++;
        prm.xi[j].y = pv[index]; index++;
    }

    prm.eta.clear();
    prm.eta.resize(Nc);

    for (unsigned int i=0; i<Nc; i++)
    {
        prm.eta[i].x = pv[index]; index++;
        prm.eta[i].y = pv[index]; index++;
    }

#ifdef TIME_DISCRETE
    for (unsigned int s=0; s<Nt; s++)
    {
        const double tau = pv[index]; index++;
        prm.tau.push_back(tau);
    }
#endif
}

auto Problem2PNeumann::PrmToVector(const OptimizeParameterP &prm, DoubleVector &pv) const -> void
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;
#ifdef TIME_DISCRETE
    unsigned int Nt = mEquParameter.Nt;
#endif

    pv.clear();
#ifdef TIME_DISCRETE
    pv.resize(2*Nc*No+2*No+2*Nc+Nt);
#else
    pv.resize(2*Nc*No+2*No+2*Nc);
#endif

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            pv[i*No + j] = prm.k[i][j];
        }
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            pv[i*No + j + Nc*No] = prm.z[i][j];
        }
    }

    for (unsigned int j=0; j<No; j++)
    {
        pv[2*j + 0 + 2*Nc*No] = prm.xi[j].x;
        pv[2*j + 1 + 2*Nc*No] = prm.xi[j].y;
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        pv[2*i + 0 + 2*No + 2*Nc*No] = prm.eta[i].x;
        pv[2*i + 1 + 2*No + 2*Nc*No] = prm.eta[i].y;
    }

#ifdef TIME_DISCRETE
    for (unsigned int s=0; s<Nt; s++)
    {
        pv[s + 2*Nc + 2*No + 2*Nc*No] = prm.tau[s];
    }
#endif
}
