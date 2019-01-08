#include "problem2h_solver1.h"

auto Problem2HDirichlet1::f_layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{}

auto Problem2HDirichlet1::b_layerInfo(const DoubleMatrix &p UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{}

auto Problem2HDirichlet1::checkGradient1(const Problem2HDirichlet1 &prob) -> void
{
    EquationParameterH e_prm = prob.mEquParameter;
    OptimizeParameterH o_prm = prob.mOptParameter;
    OptimizeParameterH r_prm = prob.mRegParameter;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prob.PrmToVector(o_prm, pv);
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();
    DoubleVector r_pv;
    prob.PrmToVector(r_prm, r_pv);
    IPrinter::print(r_pv, r_pv.length(), 6, 4);
    IPrinter::printSeperatorLine();
    DoubleVector ag(pv.length());
    double functional = prob.fx(pv);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    prob.gradient(pv, ag);
    puts("Gradients are calculated.");

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    puts("Calculating numerical gradients.... dh=0.01");
    puts("*** Calculating numerical gradients for k...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for z...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for xi..... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    puts("*** Calculating numerical gradients for eta.... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    puts("Numerical gradients are calculated.");

    puts("Calculating numerical gradients.... hx=0.001");
    puts("*** Calculating numerical gradients for k...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for z...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for xi..... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    puts("*** Calculating numerical gradients for eta.... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    puts("Numerical gradients are calculated.");

    //k------------------------------------------------------//
    IPrinter::printSeperatorLine("k");
    DoubleVector pk0 = pv.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector ak0 = ag.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nk1 = ng1.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nk2 = ng2.mid(0, e_prm.Nc*e_prm.No-1);

    IPrinter::print(pk0,pk0.length(),14,4);
    IPrinter::print(ak0,ak0.length(),14,4); ak0.L2Normalize();
    IPrinter::print(nk1,nk1.length(),14,4); nk1.L2Normalize();
    IPrinter::print(nk2,nk2.length(),14,4); nk2.L2Normalize();
    IPrinter::print(ak0,ak0.length(),14,4);
    IPrinter::print(nk1,nk1.length(),14,4);
    IPrinter::print(nk2,nk2.length(),14,4);

    //z------------------------------------------------------//
    IPrinter::printSeperatorLine("z");
    DoubleVector pz0 = pv.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector az0 = ag.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nz1 = ng1.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nz2 = ng2.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);

    IPrinter::print(pz0,pz0.length(),14,4);
    IPrinter::print(az0,az0.length(),14,4); az0.L2Normalize();
    IPrinter::print(nz1,nz1.length(),14,4); nz1.L2Normalize();
    IPrinter::print(nz2,nz2.length(),14,4); nz2.L2Normalize();
    IPrinter::print(az0,az0.length(),14,4);
    IPrinter::print(nz1,nz1.length(),14,4);
    IPrinter::print(nz2,nz2.length(),14,4);

    //xi------------------------------------------------------//
    IPrinter::printSeperatorLine("xi");
    DoubleVector pe0 = pv.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector ae0 = ag.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector ne1 = ng1.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector ne2 = ng2.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);

    IPrinter::print(pe0,pe0.length(),14,4);
    IPrinter::print(ae0,ae0.length(),14,4); ae0.L2Normalize();
    IPrinter::print(ne1,ne1.length(),14,4); ne1.L2Normalize();
    IPrinter::print(ne2,ne2.length(),14,4); ne2.L2Normalize();
    IPrinter::print(ae0,ae0.length(),14,4);
    IPrinter::print(ne1,ne1.length(),14,4);
    IPrinter::print(ne2,ne2.length(),14,4);

    //eta------------------------------------------------------//
    IPrinter::printSeperatorLine("eta");
    DoubleVector px0 = pv.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector ax0 = ag.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nx1 = ng1.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nx2 = ng2.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);

    IPrinter::print(px0,px0.length(),14,4);
    IPrinter::print(ax0,ax0.length(),14,4); ax0.L2Normalize();
    IPrinter::print(nx1,nx1.length(),14,4); nx1.L2Normalize();
    IPrinter::print(nx2,nx2.length(),14,4); nx2.L2Normalize();
    IPrinter::print(ax0,ax0.length(),14,4);
    IPrinter::print(nx1,nx1.length(),14,4);
    IPrinter::print(nx2,nx2.length(),14,4);
    IPrinter::printSeperatorLine();
}

auto Problem2HDirichlet1::checkGradient2(const Problem2HDirichlet1 &prob) -> void
{
    EquationParameterH e_prm = prob.mEquParameter;
    OptimizeParameterH o_prm = prob.mOptParameter;
    OptimizeParameterH r_prm = prob.mRegParameter;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prob.PrmToVector(o_prm, pv);
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();
    DoubleVector r_pv;
    prob.PrmToVector(r_prm, r_pv);
    IPrinter::print(r_pv, r_pv.length(), 6, 4);
    IPrinter::printSeperatorLine();
    DoubleVector ag(pv.length());
    double functional = prob.fx(pv);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    prob.gradient(pv, ag);
    puts("Gradients are calculated.");

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    puts("Calculating numerical gradients.... dh=0.01");
    puts("*** Calculating numerical gradients for k...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for z...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for xi..... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    puts("*** Calculating numerical gradients for eta.... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    puts("Numerical gradients are calculated.");

    puts("Calculating numerical gradients.... hx=0.001");
    puts("*** Calculating numerical gradients for k...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for z...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for xi..... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    puts("*** Calculating numerical gradients for eta.... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    puts("Numerical gradients are calculated.");

    DoubleVector nag  = ag;  nag.EuclideanNormalize();
    DoubleVector nng1 = ng1; nng1.EuclideanNormalize();
    DoubleVector nng2 = ng2; nng2.EuclideanNormalize();

    //k------------------------------------------------------//
    IPrinter::printSeperatorLine("k");
    DoubleVector pk0 = pv.mid(0, e_prm.Nc*e_prm.No-1);
    IPrinter::print(pk0,pk0.length(),14,4);

    DoubleVector ak0 = ag.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nk1 = ng1.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nk2 = ng2.mid(0, e_prm.Nc*e_prm.No-1);
    IPrinter::print(ak0,ak0.length(),14,4);
    IPrinter::print(nk1,nk1.length(),14,4);
    IPrinter::print(nk2,nk2.length(),14,4);

    DoubleVector nak0 = nag.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nnk1 = nng1.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nnk2 = nng2.mid(0, e_prm.Nc*e_prm.No-1);

    IPrinter::print(nak0,nak0.length(),14,4);
    IPrinter::print(nnk1,nnk1.length(),14,4);
    IPrinter::print(nnk2,nnk2.length(),14,4);

    //z------------------------------------------------------//
    IPrinter::printSeperatorLine("z");
    DoubleVector pz0 = pv.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    IPrinter::print(pz0,pz0.length(),14,4);

    DoubleVector az0 = ag.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nz1 = ng1.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nz2 = ng2.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    IPrinter::print(az0,az0.length(),14,4);
    IPrinter::print(nz1,nz1.length(),14,4);
    IPrinter::print(nz2,nz2.length(),14,4);

    DoubleVector naz0 = nag.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nnz1 = nng1.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nnz2 = nng2.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    IPrinter::print(naz0,naz0.length(),14,4);
    IPrinter::print(nnz1,nnz1.length(),14,4);
    IPrinter::print(nnz2,nnz2.length(),14,4);

    //xi------------------------------------------------------//
    IPrinter::printSeperatorLine("xi");
    DoubleVector pe0 = pv.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    IPrinter::print(pe0,pe0.length(),14,4);

    DoubleVector ae0 = ag.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector ne1 = ng1.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector ne2 = ng2.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    IPrinter::print(ae0,ae0.length(),14,4);
    IPrinter::print(ne1,ne1.length(),14,4);
    IPrinter::print(ne2,ne2.length(),14,4);

    DoubleVector nae0 = nag.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector nne1 = nng1.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector nne2 = nng2.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    IPrinter::print(nae0,nae0.length(),14,4);
    IPrinter::print(nne1,nne1.length(),14,4);
    IPrinter::print(nne2,nne2.length(),14,4);

    //eta------------------------------------------------------//
    IPrinter::printSeperatorLine("eta");
    DoubleVector px0 = pv.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    IPrinter::print(px0,px0.length(),14,4);

    DoubleVector ax0 = ag.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nx1 = ng1.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nx2 = ng2.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    IPrinter::print(ax0,ax0.length(),14,4);
    IPrinter::print(nx1,nx1.length(),14,4);
    IPrinter::print(nx2,nx2.length(),14,4);

    DoubleVector nax0 = nag.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nnx1 = nng1.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nnx2 = nng2.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    IPrinter::print(nax0,nax0.length(),14,4);
    IPrinter::print(nnx1,nnx1.length(),14,4);
    IPrinter::print(nnx2,nnx2.length(),14,4);
    IPrinter::printSeperatorLine();
}

Problem2HDirichlet1::Problem2HDirichlet1()
{
    r = 0.0;
    regEpsilon = 0.0;
}

Problem2HDirichlet1::~Problem2HDirichlet1()
{}

auto Problem2HDirichlet1::fx(const DoubleVector &pv) const -> double
{
    OptimizeParameterH o_prm;

    VectorToPrm(pv, o_prm);

    Problem2HDirichlet1* prob = const_cast<Problem2HDirichlet1*>(this);
    prob->mOptParameter = o_prm;

    std::vector<DoubleMatrix> u;
    spif_vectorH u_info;
    prob->solveForwardIBVP(u, u_info, true);

    double intgrl = integral(u);
    double nrm = norm(mEquParameter, o_prm, mRegParameter);
    double pnt = penalty(u_info, o_prm);
    double sum = intgrl + regEpsilon*nrm + r*pnt;

    for (unsigned int i=0; i<u.size(); i++)      u[i].clear();      u.clear();
    for (unsigned int j=0; j<u_info.size(); j++) u_info[j].clear(); u_info.clear();

    return sum;
}

double Problem2HDirichlet1::integral(const std::vector<DoubleMatrix> &vu) const
{
    const double ht = timeDimension().step();
    double sum = 0.0;
    const DoubleMatrix &u0 = vu[0]; sum += 0.5*integralU(u0);
    for (unsigned int l=1; l<=LD-1; l++)
    {
        const DoubleMatrix &u1 = vu[2*l]; sum += integralU(u1);
    }
    const DoubleMatrix &u2 = vu[2*LD]; sum += 0.5*integralU(u2);
    return sum*ht;
}

double Problem2HDirichlet1::integralU(const DoubleMatrix &u) const
{
    const double hx = spaceDimension(Dimension::DimensionX).step();
    const double hy = spaceDimension(Dimension::DimensionY).step();
    const unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );
    const unsigned int M = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionY).size() );

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u[0][0]; usum += 0.25 * udiff * udiff * mu(0, 0);
    udiff = u[0][N]; usum += 0.25 * udiff * udiff * mu(N, 0);
    udiff = u[M][0]; usum += 0.25 * udiff * udiff * mu(0, M);
    udiff = u[M][N]; usum += 0.25 * udiff * udiff * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u[0][n]; usum += 0.5 * udiff * udiff * mu(n, 0);
        udiff = u[M][n]; usum += 0.5 * udiff * udiff * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u[m][0]; usum += 0.5 * udiff * udiff * mu(0, m);
        udiff = u[m][N]; usum += 0.5 * udiff * udiff * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u[m][n]; usum += udiff * udiff * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

double Problem2HDirichlet1::norm(const EquationParameterH& e_prm, const OptimizeParameterH &o_prm, const OptimizeParameterH &r_prm) const
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

double Problem2HDirichlet1::penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const
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

double Problem2HDirichlet1::gpi(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const
{
    double gpi_ln = fabs( g0i(i, layer, u_info, o_prm) ) - ( vmax.at(i) - vmin.at(i) )/2.0;
    return gpi_ln > 0.0 ? gpi_ln : 0.0;
}

double Problem2HDirichlet1::g0i(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const
{
    double vi = 0.0;
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        const SpacePointInfoH &u_xij = u_info[j];
        vi += o_prm.k[i][j] * ( u_xij.vl[layer] - o_prm.z[i][j] );
    }
    return ( vmax.at(i) + vmin.at(i) )/2.0 - vi;
}

double Problem2HDirichlet1::sign(double x) const
{
    if (x < 0.0)       return -1.0;
    else if (x > 0.0)  return +1.0;
    else               return  0.0;
}

auto Problem2HDirichlet1::gradient(const DoubleVector & pv, DoubleVector &g) const -> void
{
    const unsigned int L   = static_cast<const unsigned int>(mtimeDimension.size());
    const double ht        = mtimeDimension.step();
    const unsigned int Nc  = mEquParameter.Nc;
    const unsigned int No  = mEquParameter.No;
    const unsigned int LLD = L + LD;

    OptimizeParameterH o_prm;
    VectorToPrm(pv, o_prm);

    Problem2HDirichlet1* prob = const_cast<Problem2HDirichlet1*>(this);
    prob->mOptParameter = o_prm;

    std::vector<DoubleMatrix> u;

    spif_vectorH u_info;
    solveForwardIBVP(u, u_info, true);
    spif_vectorH p_info;
    solveBackwardIBVP(u, p_info, true, u_info);

    g.clear();
    g.resize(pv.length(), 0.0);
    unsigned int gi = 0;

    // k
    if (optimizeK)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfoH &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointInfoH &uj = u_info[j];

                double grad_Kij = 0.0;
                double zij = o_prm.z[i][j];

                grad_Kij += 0.5 * (pi.vl[0] + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * (uj.vl[0] - zij); //printf("%f %f\n", pi.vl[0], uj.vl[0]);
                for (unsigned int ln=1; ln<=LLD-1; ln++)
                {
                    grad_Kij += (pi.vl[2*ln] + 2.0*r*gpi(i,2*ln,u_info,o_prm)*sgn(g0i(i,2*ln,u_info,o_prm))) * (uj.vl[2*ln] - zij); //printf("%f %f\n", pi.vl[2*ln], uj.vl[2*ln]);
                }
                grad_Kij += 0.5 * (pi.vl[2*LLD] + 2.0*r*gpi(i,2*LLD,u_info,o_prm)*sgn(g0i(i,2*LLD,u_info,o_prm))) * (uj.vl[2*LLD] - zij); //printf("%f %f\n", pi.vl[2*LLD], uj.vl[2*LLD]);

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
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            g[gi++] = 0.0;
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

auto Problem2HDirichlet1::norm(const DoubleVector &v) const -> double
{
    return EuclideanNorm(v);
}

auto Problem2HDirichlet1::normalize(DoubleVector &v) const -> void
{
    if (optimizeK) { DoubleVector kv = v.mid(0, 3);   IVectorNormalizer::EuclideanNormalize(kv); v[0]  = kv[0]; v[1]  = kv[1];  v[2]  = kv[2]; v[3]  = kv[3]; kv.clear(); }
    if (optimizeZ) { DoubleVector zv = v.mid(4, 7);   IVectorNormalizer::EuclideanNormalize(zv); v[4]  = zv[0]; v[5]  = zv[1];  v[6]  = zv[2]; v[7]  = zv[3]; zv.clear(); }
    if (optimizeO) { DoubleVector ov = v.mid(8, 11);  IVectorNormalizer::EuclideanNormalize(ov); v[8]  = ov[0]; v[9]  = ov[1];  v[10] = ov[2]; v[11] = ov[3]; ov.clear(); }
    if (optimizeZ) { DoubleVector cv = v.mid(12, 15); IVectorNormalizer::EuclideanNormalize(cv); v[12] = cv[0]; v[13] = cv[1];  v[14] = cv[2]; v[15] = cv[3]; cv.clear(); }
}

auto Problem2HDirichlet1::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const -> void
{
    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(f); C_UNUSED(alpha); C_UNUSED(result);
    const char* msg = nullptr; C_UNUSED(msg);
    if (result == GradientMethod::MethodResult::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION   ";
    if (result == GradientMethod::MethodResult::FIRST_ITERATION)          msg = "FIRST_ITERATION         ";
    if (result == GradientMethod::MethodResult::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    if (result == GradientMethod::MethodResult::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS     ";
    if (result == GradientMethod::MethodResult::NEXT_ITERATION)           msg = "NEXT_ITERATION          ";

    Problem2HDirichlet1* prob = const_cast<Problem2HDirichlet1*>(this);
    //prob->gm->setStepTolerance(10.0*alpha);
    OptimizeParameterH o_prm;
    VectorToPrm(x, o_prm);
    prob->mOptParameter = o_prm;
    std::vector<DoubleMatrix> u;
    spif_vectorH u_info;
    solveForwardIBVP(u, u_info, true);
    double ing = integral(u);
    double pnt = penalty(u_info, o_prm);
    double nrm = norm(prob->mEquParameter, prob->mOptParameter, prob->mRegParameter);

    const unsigned int v_length = static_cast<const unsigned int>(timeDimension().size()) + LD;
    DoubleVector v1(v_length+1);
    DoubleVector v2(v_length+1);

    for (unsigned int ln=0; ln<=v_length; ln++)
    {
        v1[ln] = v(0, o_prm, mEquParameter, u_info, 2*ln);
        v2[ln] = v(1, o_prm, mEquParameter, u_info, 2*ln);
    }

    IPrinter::printVector(v1, "v1", 10);
    IPrinter::printVector(v2, "v2", 10);

    printf("I[%3d]: F:%.6f I:%.6f P:%.6f N:%.5f R:%.3f e:%.3f a:%.6f ", i, f, ing, pnt, nrm, r, regEpsilon, alpha);
    printf("min:%.6f max:%.6f min:%.6f max:%.6f U0:%.8f UT:%.8f", u.at(0).min(), u.at(0).max(), u.at(2*LD).min(), u.at(2*LD).max(), integralU(u[0]), integralU(u[LD]));
    printf("\n");
    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
    DoubleVector n = g; n.L2Normalize();
    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15]);

    u.clear();
    u_info.clear();

    C_UNUSED(prob);
    IPrinter::printSeperatorLine();

    //prob->optimizeK = i%4 == 3;
    //prob->optimizeZ = i%4 == 0;
    //prob->optimizeO = i%4 == 1;
    //prob->optimizeC = i%4 == 2;
}

auto Problem2HDirichlet1::project(DoubleVector &, unsigned int) -> void {}

auto Problem2HDirichlet1::project(DoubleVector &pv) const -> void
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

    //    for (unsigned int index = start; index <= end; index++)
    //    {
    //        if (index ==  8) { if (pv[ 8] < 0.15) pv[ 8] = 0.15; if (pv[ 8] > 0.45) pv[ 8] = 0.45; }
    //        if (index ==  9) { if (pv[ 9] < 0.55) pv[ 9] = 0.55; if (pv[ 9] > 0.95) pv[ 9] = 0.95; }
    //        if (index == 10) { if (pv[10] < 0.55) pv[10] = 0.55; if (pv[10] > 0.85) pv[10] = 0.85; }
    //        if (index == 11) { if (pv[11] < 0.15) pv[11] = 0.15; if (pv[11] > 0.55) pv[11] = 0.55; }

    //        if (index == 12) { if (pv[12] < 0.05) pv[12] = 0.05; if (pv[12] > 0.35) pv[12] = 0.35; }
    //        if (index == 13) { if (pv[13] < 0.05) pv[13] = 0.05; if (pv[13] > 0.45) pv[13] = 0.45; }
    //        if (index == 14) { if (pv[14] < 0.55) pv[14] = 0.55; if (pv[14] > 0.85) pv[14] = 0.85; }
    //        if (index == 15) { if (pv[15] < 0.65) pv[15] = 0.65; if (pv[15] > 0.95) pv[15] = 0.95; }
    //    }

    //IPrinter::print(pv.mid(start, end));
    for (unsigned int index = start; index <=end; index++)
    {
        //projectControlPoints(pv, index);
        projectMeasurePoints(pv, index);
    }
}

auto Problem2HDirichlet1::projectControlPoints(DoubleVector &pv, unsigned int index) const -> void
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

auto Problem2HDirichlet1::projectMeasurePoints(DoubleVector &pv, unsigned int index) const -> void
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

auto Problem2HDirichlet1::solveForwardIBVP1(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use) const -> void
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
    const double lambda = 0.25;

    const double a        = mEquParameter.a;
    const double alpha    = mEquParameter.lambda;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;

#ifdef TIME_DISCRETE_H
    const unsigned int Nt = mEquParameter.Nt;
#endif

    const double ht_ht_025 = ht*ht*0.25;
    const double alpha_ht_025 = alpha*ht*0.25;

    const double m_aa_htht__hxhx_025_lambda = -(0.25*a*a)*((ht*ht)/(hx*hx))*lambda;
    const double b_aa_htht__hxhx = +(1.0 + 0.5*(a*a)*((ht*ht)/(hx*hx))*lambda + alpha_ht_025);
    const double p_aa_htht__hyhy_025 = +(0.25*a*a)*((ht*ht)/(hy*hy));
    const double p_aa_htht__hyhy_025_lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*lambda;
    const double p_aa_htht__hyhy_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*(1.0-2.0*lambda);

    const double m_aa_htht__hyhy_025_lambda = -(0.25*a*a)*((ht*ht)/(hy*hy))*lambda;
    const double b_aa_htht__hyhy = +(1.0 + 0.5*(a*a)*((ht*ht)/(hy*hy))*lambda + alpha_ht_025);
    const double p_aa_htht__hxhx_025 = +(0.25*a*a)*((ht*ht)/(hx*hx));
    const double p_aa_htht__hxhx_025_lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*lambda;
    const double p_aa_htht__hxhx_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*(1.0-2.0*lambda);

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

    //------------------------------------- initial conditions -------------------------------------//
    initiatePulseGrid();
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
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;
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
        currentLayerFGrid(u10, cntrlDeltaGridList, measuremntGirdList);
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
        currentLayerFGrid(u15, cntrlDeltaGridList, measuremntGirdList);
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

auto Problem2HDirichlet1::f_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichlet1::initiatePulseGrid() const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    const_cast<Problem2HDirichlet1*>(this)->mPulseWeightMatrix.clear();
    const_cast<Problem2HDirichlet1*>(this)->mPulseWeightMatrix.resize(M+1, N+1, 0.0);

    const unsigned int Ns = mEquParameter.Ns;
    std::vector<DeltaGrid2D> deltaGrids(Ns);
    for (unsigned int s=0; s<Ns; s++)
    {
        deltaGrids[s].initGrid(N, hx, M, hy);
        deltaGrids[s].distributeGauss(mEquParameter.theta[s], 1, 1);
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            const_cast<Problem2HDirichlet1*>(this)->mPulseWeightMatrix[m][n] = 0.0;
            for (unsigned int s=0; s<Ns; s++)
            {
                const_cast<Problem2HDirichlet1*>(this)->mPulseWeightMatrix[m][n] += mEquParameter.q[s]*deltaGrids[s].weight(n, m);
            }
        }
    }

    for (unsigned int s=0; s<Ns; s++)
    {
        deltaGrids[s].cleanGrid();
    }
    deltaGrids.clear();
}

auto Problem2HDirichlet1::currentLayerFGrid(const DoubleMatrix &u,
                                            const std::vector<DeltaGrid2D> &controlDeltaGrids,
                                            const std::vector<DeltaGrid2D> &measurementDeltaGrids) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    const_cast<Problem2HDirichlet1*>(this)->mCrFfxWeightMatrix.clear();
    const_cast<Problem2HDirichlet1*>(this)->mCrFfxWeightMatrix.resize(M+1, N+1, 0.0);

    const unsigned int Nc = mEquParameter.Nc;
    const unsigned int No = mEquParameter.No;

    double* _u = new double[No];
    for (unsigned int j=0; j<No; j++)
    {
        _u[j] = 0.0;
        const DeltaGrid2D &dg = measurementDeltaGrids[j];
        for (unsigned int m=dg.minY(); m<=dg.maxY(); m++)
        {
            for (unsigned int n=dg.minX(); n<=dg.maxX(); n++)
            {
                _u[j] += u[m][n] * (dg.weight(n,m) * (hx*hy));
            }
        }
        _u[j] *= (1.0 + noise);
    }

    double *_v = new double[Nc];
    for (unsigned int i=0; i<Nc; i++)
    {
        _v[i] = 0.0;
        for (unsigned int j=0; j<No; j++)
        {
            _v[i] += mOptParameter.k[i][j] * (_u[j] - mOptParameter.z[i][j]);
        }
    }
    delete [] _u;

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            const_cast<Problem2HDirichlet1*>(this)->mCrFfxWeightMatrix[m][n] = 0.0;
            for (unsigned int i=0; i<Nc; i++)
            {
                const DeltaGrid2D &dg = controlDeltaGrids[i];
                const_cast<Problem2HDirichlet1*>(this)->mCrFfxWeightMatrix[m][n] += _v[i] * dg.weight(n,m);
            }
        }
    }

    delete [] _v;
}

auto Problem2HDirichlet1::currentLayerBGrid(const DoubleMatrix &p, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids,
                                            double ln, const spif_vectorH &u_info) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    const_cast<Problem2HDirichlet1*>(this)->mCrBfxWeightMatrix.clear();
    const_cast<Problem2HDirichlet1*>(this)->mCrBfxWeightMatrix.resize(M+1, N+1, 0.0);

    const unsigned int Nc = mEquParameter.Nc;
    const unsigned int No = mEquParameter.No;

    double* _p = new double[Nc];
    for (unsigned int i=0; i<Nc; i++)
    {
        _p[i] = 0.0;
        const DeltaGrid2D &dg = controlDeltaGrids[i];
        for (unsigned int m=dg.minY(); m<=dg.maxY(); m++)
        {
            for (unsigned int n=dg.minX(); n<=dg.maxX(); n++)
            {
                _p[i] += p[m][n] * (dg.weight(n,m) * (hx*hy));
            }
        }
    }

    double *_w = new double[No];
    for (unsigned int j=0; j<No; j++)
    {
        _w[j] = 0.0;
        for (unsigned int i=0; i<Nc; i++)
        {
            _w[j] += mOptParameter.k[i][j] * (_p[i] + 2.0*r*gpi(i, ln, u_info, mOptParameter)*sgn(g0i(i, ln, u_info, mOptParameter)) );
        }
    }
    delete [] _p;

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            const_cast<Problem2HDirichlet1*>(this)->mCrBfxWeightMatrix[m][n] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                const DeltaGrid2D &dg = measurementDeltaGrids[j];
                const_cast<Problem2HDirichlet1*>(this)->mCrBfxWeightMatrix[m][n] += _w[j] * dg.weight(n,m);
            }
        }
    }

    delete [] _w;
}

auto Problem2HDirichlet1::f_initial2(const SpaceNodePDE &sn) const -> double
{
    //return 0.0;
    return mPulseWeightMatrix[sn.j][sn.i];
}

auto Problem2HDirichlet1::f_boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichlet1::prepareInfo(unsigned int N, const std::vector<SpacePoint> &points, spif_vectorH &info, unsigned int size) const -> void
{
    info.resize(N);
    for (unsigned int i=0; i<N; i++)
    {
        SpacePointInfoH &inf = info[i];
        const SpacePoint &sp = points[i];
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(size);
    }
}

auto Problem2HDirichlet1::add2Info(const DoubleMatrix &u, spif_vectorH &info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &deltaList) const -> void
{
    unsigned int N = static_cast<unsigned int>(deltaList.size());

    for (unsigned int i=0; i<N; i++)
    {
        const DeltaGrid2D &dg = deltaList.at(i);
        SpacePointInfoH &ui = info[i];

        for (unsigned int m=dg.minY(); m<=dg.maxY(); m++)
        {
            for (unsigned int n=dg.minX(); n<=dg.maxX(); n++)
            {
                ui.vl[ln] += u[m][n] * (dg.weight(n,m) * (hx*hy));
                if (dg.isCenter(n,m))
                {
                    const unsigned int rx = static_cast<const unsigned int>(dg.rx());
                    const unsigned int ry = static_cast<const unsigned int>(dg.ry());

                    ui.dx[ln] = (u[ry][rx+1] - u[ry][rx-1])/(2.0*hx);
                    ui.dy[ln] = (u[ry+1][rx] - u[ry-1][rx])/(2.0*hy);

                    ui.dx[ln] += ((dg.p().x-rx*hx)/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
                    ui.dy[ln] += ((dg.p().y-ry*hy)/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);

                    //ui.dxx[ln] = (1.0/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
                    //ui.dyy[ln] = (1.0/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);
                }
            }
        }
    }
}

auto Problem2HDirichlet1::solveBackwardIBVP1(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info) const -> void
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
    const double lambda = 0.25;

    const double a        = mEquParameter.a;
    const double alpha    = mEquParameter.lambda;
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
    TimeNodePDE tn05; tn05.i = LLD-1; tn05.t = LLD*ht + 0.5*ht;
    TimeNodePDE tn10; tn10.i = LLD-1; tn10.t = LLD*ht;
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
            _w[j] += mOptParameter.k[i][j] * (_p00[i] + 2.0*r*gpi(i, 2*LLD, u_info, mOptParameter)*sgn(g0i(i,2*LLD, u_info, mOptParameter)));
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
            sum += lambda*b_initial2(sn10);

            for (unsigned int j=0; j<No; j++)
            {
                const DeltaGrid2D &dg = measuremntGirdList[j];
                sum += _w[j] * dg.weight(n,m);
            }

            p05[m][n] = p00[m][n] - (ht*0.5) * b_initial2(sn10) + 0.125*ht*ht*sum;
            p10[m][n] = p00[m][n] - (ht*1.0) * b_initial2(sn10) + 0.500*ht*ht*sum;
        }
    }
    if (use == true) add2Info(p05, p_info, 2*(LLD-1)+1, hx, hy, cntrlDeltaGridList); b_layerInfo(p05, 2*LLD-1);
    if (use == true) add2Info(p10, p_info, 2*(LLD-1)+0, hx, hy, cntrlDeltaGridList); b_layerInfo(p10, 2*LLD-2);
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

auto Problem2HDirichlet1::b_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichlet1::b_initial2(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichlet1::b_boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double
{
    return 0.0;
}

auto Problem2HDirichlet1::b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const -> double
{
    return -2.0*mu(n,m)*(u[m][n]);
}

auto Problem2HDirichlet1::PrmToVector(const OptimizeParameterH &prm, DoubleVector &pv) const -> void
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    pv.clear();
    pv.resize(2*Nc*No+2*No+2*Nc);

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
}

auto Problem2HDirichlet1::VectorToPrm(const DoubleVector &pv, OptimizeParameterH &prm) const -> void
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

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
}

auto Problem2HDirichlet1::v(unsigned int i, OptimizeParameterH o_prm, EquationParameterH e_prm, const spif_vectorH &u_info, unsigned int ln) const -> double
{
    const unsigned int No = static_cast<const unsigned int>(e_prm.No);
    double v = 0.0;
    for (unsigned int j=0; j<No; j++)
    {
        v += o_prm.k[i][j] * (u_info[j].vl[ln]-o_prm.z[i][j]);
    }
    return v;
}
