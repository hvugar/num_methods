#include "problem2h_solver_base.h"

Problem2hDirichletBase::Problem2hDirichletBase()
{
    r = 0.0;
    regEpsilon = 0.0;
}

Problem2hDirichletBase::~Problem2hDirichletBase()
{
    mPulseWeightMatrix.clear();
    mCrFfxWeightMatrix.clear();
    mCrBfxWeightMatrix.clear();
}

auto Problem2hDirichletBase::initPulseWeightMatrix(const std::vector<SpacePoint> &theta,
                                                   const std::vector<double> q) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    const unsigned int Ns = static_cast<unsigned int>(theta.size());
    DoubleMatrix &pulseWeightMatrix = const_cast<Problem2hDirichletBase*>(this)->mPulseWeightMatrix;

    std::vector<DeltaGrid2D> deltaGrids(Ns);
    for (unsigned int s=0; s<Ns; s++)
    {
        deltaGrids[s].initGrid(N, hx, M, hy);
        deltaGrids[s].distributeGauss(theta[s], 8, 8);
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            pulseWeightMatrix[m][n] = 0.0;
            for (unsigned int s=0; s<Ns; s++)
            {
                pulseWeightMatrix[m][n] += q[s]*deltaGrids[s].weight(n, m);
            }
        }
    }

    for (unsigned int s=0; s<Ns; s++) deltaGrids[s].cleanGrid();
    deltaGrids.clear();
}

double Problem2hDirichletBase::boundary(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return NAN;
}

auto Problem2hDirichletBase::f_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2hDirichletBase::f_initial2(const SpaceNodePDE &sn) const -> double
{
    unsigned int m = static_cast<unsigned int>(sn.j);
    unsigned int n = static_cast<unsigned int>(sn.i);
    return mPulseWeightMatrix[m][n];
}

auto Problem2hDirichletBase::f_boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2hDirichletBase::b_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2hDirichletBase::b_initial2(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2hDirichletBase::b_boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2hDirichletBase::f_layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{
    //if (ln == 1000) IPrinter::printMatrix(u);
    //printf("%4d %.8f %.8f\n", ln, u.min(), u.max());
    return;

    {
        Problem2hDirichletBase* tmp = const_cast<Problem2hDirichletBase*>(this);
        std::vector<DoubleMatrix> &rvu = tmp->vu;

        rvu.push_back(u);
        //if (rvu.size() > 3) rvu.erase(rvu.begin());
        if (rvu.size() > 2*LD+1) rvu.erase(rvu.begin());

        if (ln == 1701)
        {
            tmp->mOptParameter.k[0][0] = 0.0;
            tmp->mOptParameter.k[0][1] = 0.0;
            tmp->mOptParameter.k[1][0] = 0.0;
            tmp->mOptParameter.k[1][1] = 0.0;
        }

        if (rvu.size() == 2*LD+1)
//        if (rvu.size() == 3)
        {
            //std::cout << ln << " " << rvu.size() << std::endl;
            if (ln%2==0)
            {
                double fx = integral(rvu);
                //printf("%d %d %.10f\n", ln, ln/2-50, fx);
                printf("%.10f\n", fx);
            }
            //printf("%.10f\n", fx);
        }
    }
}

auto Problem2hDirichletBase::b_layerInfo(const DoubleMatrix &p UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{
    //    if (ln == 0) IPrinter::printMatrix(p);
}

auto Problem2hDirichletBase::checkGradient1(const Problem2hDirichletBase &prob) -> void
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
    IPrinter::printSeperatorLine();printf("%f\n", o_prm.eta[1].x);
}

auto Problem2hDirichletBase::checkGradient2(const Problem2hDirichletBase &prob) -> void
{
    EquationParameterH e_prm = prob.mEquParameter;
    OptimizeParameterH o_prm = prob.mOptParameter;
    OptimizeParameterH r_prm = prob.mRegParameter;

    printf("%f\n", prob.mOptParameter.eta[1].x);

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

auto Problem2hDirichletBase::PrmToVector(const OptimizeParameterH &prm, DoubleVector &pv) const -> void
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

auto Problem2hDirichletBase::VectorToPrm(const DoubleVector &pv, OptimizeParameterH &prm) const -> void
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

auto Problem2hDirichletBase::project(DoubleVector &, unsigned int) -> void {}

auto Problem2hDirichletBase::project(DoubleVector &pv) const -> void
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    unsigned int start = 2*Nc*No;
    unsigned int end  =  2*Nc*No + 2*No + 2*Nc - 1;

    for (unsigned int index = start; index <= end; index++)
    {
        if (pv[index] <= 0.05) pv[index] = 0.05;
        if (pv[index] >= 0.95) pv[index] = 0.95;
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
    //for (unsigned int index = start; index <=end; index++)
    //{
    //    //projectControlPoints(pv, index);
    //    projectMeasurePoints(pv, index);
    //}
}

auto Problem2hDirichletBase::projectControlPoints(DoubleVector &pv, unsigned int index) const -> void
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

auto Problem2hDirichletBase::projectMeasurePoints(DoubleVector &pv, unsigned int index) const -> void
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

auto Problem2hDirichletBase::sign(double x) const -> double
{
    if (x < 0.0)       return -1.0;
    else if (x > 0.0)  return +1.0;
    else               return  0.0;
}

auto Problem2hDirichletBase::v(unsigned int i, OptimizeParameterH o_prm, EquationParameterH e_prm, const spif_vectorH &u_info, unsigned int ln) const -> double
{
    const unsigned int No = static_cast<const unsigned int>(e_prm.No);
    double v = 0.0;
    for (unsigned int j=0; j<No; j++)
    {
        v += o_prm.k[i][j] * (u_info[j].vl[ln]-o_prm.z[i][j]);
    }
    return v;
}

auto Problem2hDirichletBase::gpi(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const -> double
{
    double gpi_ln = fabs( g0i(i, layer, u_info, o_prm) ) - ( vmax.at(i) - vmin.at(i) )/2.0;
    return gpi_ln > 0.0 ? gpi_ln : 0.0;
}

auto Problem2hDirichletBase::g0i(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const -> double
{
    double vi = 0.0;
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        const SpacePointInfoH &u_xij = u_info[j];
        vi += o_prm.k[i][j] * ( u_xij.vl[layer] - o_prm.z[i][j] );
    }
    return ( vmax.at(i) + vmin.at(i) )/2.0 - vi;
}

auto Problem2hDirichletBase::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const -> void
{
    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(f); C_UNUSED(alpha); C_UNUSED(result);
    const char* msg = nullptr; C_UNUSED(msg);
    if (result == GradientMethod::MethodResult::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION   ";
    if (result == GradientMethod::MethodResult::FIRST_ITERATION)          msg = "FIRST_ITERATION         ";
    if (result == GradientMethod::MethodResult::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    if (result == GradientMethod::MethodResult::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS     ";
    if (result == GradientMethod::MethodResult::NEXT_ITERATION)           msg = "NEXT_ITERATION          ";

    Problem2hDirichletBase* prob = const_cast<Problem2hDirichletBase*>(this);
    //prob->gm->setStepTolerance(10.0*alpha);
    OptimizeParameterH o_prm;
    VectorToPrm(x, o_prm);
    prob->mOptParameter = o_prm;
    std::vector<DoubleMatrix> u;
    spif_vectorH u_info;
    solveForwardIBVP(u, u_info, true, x);
    double ing = integral(u);
    double pnt = penalty(u_info, o_prm);
    double nrm = norm(prob->mEquParameter, prob->mOptParameter, prob->mRegParameter);

//    DoubleVector uf, um, ux;
//    for (unsigned int i=0; i<=50; i+=5)
//    {
//        uf << integralU(u[i]);
//        um << u[i].min();
//        ux << u[i].max();
//    }

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
//    IPrinter::printVector(uf, "uf", 10);
//    IPrinter::printVector(um, "um", 10);
//    IPrinter::printVector(ux, "ux", 10);

    printf("I[%3d]: F:%.6f I:%.6f P:%.6f N:%.5f R:%.3f e:%.3f a:%.6f ", i, f, ing, pnt, nrm, r, regEpsilon, alpha);
    printf("min:%.6f max:%.6f U:%.8f", u.at(1).min(), u.at(1).max(), integralU(u[0]));
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

auto Problem2hDirichletBase::fx(const DoubleVector &pv) const -> double
{
    OptimizeParameterH o_prm;

    VectorToPrm(pv, o_prm);

    Problem2hDirichletBase* prob = const_cast<Problem2hDirichletBase*>(this);
    prob->mOptParameter = o_prm;

    std::vector<DoubleMatrix> u;
    spif_vectorH u_info;
    prob->solveForwardIBVP(u, u_info, true, pv);

    double intgrl = integral(u);

    double nrm = norm(mEquParameter, o_prm, mRegParameter);
    double pnt = penalty(u_info, o_prm);
    double sum = intgrl + regEpsilon*nrm + r*pnt;

    for (unsigned int i=0; i<u.size(); i++)      u[i].clear();      u.clear();
    for (unsigned int j=0; j<u_info.size(); j++) u_info[j].clear(); u_info.clear();

    return sum;
}

double Problem2hDirichletBase::norm(const EquationParameterH& e_prm, const OptimizeParameterH &o_prm, const OptimizeParameterH &r_prm) const
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

auto Problem2hDirichletBase::norm(const DoubleVector &v) const -> double
{
    return EuclideanNorm(v);
}

auto Problem2hDirichletBase::normalize(DoubleVector &v) const -> void
{
    if (optimizeK) { DoubleVector kv = v.mid(0, 3);   IVectorNormalizer::EuclideanNormalize(kv); v[0]  = kv[0]; v[1]  = kv[1];  v[2]  = kv[2]; v[3]  = kv[3]; kv.clear(); }
    if (optimizeZ) { DoubleVector zv = v.mid(4, 7);   IVectorNormalizer::EuclideanNormalize(zv); v[4]  = zv[0]; v[5]  = zv[1];  v[6]  = zv[2]; v[7]  = zv[3]; zv.clear(); }
    if (optimizeO) { DoubleVector ov = v.mid(8, 11);  IVectorNormalizer::EuclideanNormalize(ov); v[8]  = ov[0]; v[9]  = ov[1];  v[10] = ov[2]; v[11] = ov[3]; ov.clear(); }
    if (optimizeZ) { DoubleVector cv = v.mid(12, 15); IVectorNormalizer::EuclideanNormalize(cv); v[12] = cv[0]; v[13] = cv[1];  v[14] = cv[2]; v[15] = cv[3]; cv.clear(); }
}

auto Problem2hDirichletBase::b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const -> double
{
    return -2.0*mu(n,m)*(u[m][n]);
}

double Problem2hDirichletBase::integralU(const DoubleMatrix &u) const
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

auto Problem2hDirichletBase::prepareInfo(unsigned int N, const std::vector<SpacePoint> &points, spif_vectorH &info, unsigned int size) const -> void
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

auto Problem2hDirichletBase::add2Info(const DoubleMatrix &u, spif_vectorH &info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &deltaList) const -> void
{
    const unsigned int N = static_cast<unsigned int>(deltaList.size());

    for (unsigned int i=0; i<N; i++)
    {
        const DeltaGrid2D &deltagrid = deltaList[i];
        SpacePointInfoH &ui = info[i];
        ui.vl[ln] = deltagrid.consentrateInPoint(u, ui.dx[ln], ui.dy[ln]);

//        for (unsigned int m=deltagrid.minY(); m<=deltagrid.maxY(); m++)
//        {
//            for (unsigned int n=deltagrid.minX(); n<=deltagrid.maxX(); n++)
//            {
//                ui.vl[ln] += u[m][n] * (deltagrid.weight(n,m) * (hx*hy));
//            }
//        }

//        const unsigned int rx = static_cast<unsigned int>(deltagrid.rx());
//        const unsigned int ry = static_cast<unsigned int>(deltagrid.ry());
//        const double px = deltagrid.p().x;
//        const double py = deltagrid.p().y;

//        ui.dx[ln] = (u[ry][rx+1] - u[ry][rx-1])/(2.0*hx);
//        ui.dy[ln] = (u[ry+1][rx] - u[ry-1][rx])/(2.0*hy);

//        ui.dx[ln] += ((px-rx*hx)/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
//        ui.dy[ln] += ((py-ry*hy)/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);
    }
}

auto Problem2hDirichletBase::currentLayerFGrid(const DoubleMatrix &u,
                                            const std::vector<DeltaGrid2D> &controlDeltaGrids,
                                            const std::vector<DeltaGrid2D> &measurementDeltaGrids) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    Problem2hDirichletBase* const_this = const_cast<Problem2hDirichletBase*>(this);

//    const_this->mCrFfxWeightMatrix.clear();
//    const_this->mCrFfxWeightMatrix.resize(M+1, N+1, 0.0);

    const unsigned int Nc = mEquParameter.Nc;
    const unsigned int No = mEquParameter.No;

    double* _u = new double[No];
    for (unsigned int j=0; j<No; j++)
    {
        _u[j] = measurementDeltaGrids[j].consentrateInPoint(u);
//        const DeltaGrid2D &dg = measurementDeltaGrids[j];
//        for (unsigned int m=dg.minY(); m<=dg.maxY(); m++)
//        {
//            for (unsigned int n=dg.minX(); n<=dg.maxX(); n++)
//            {
//                _u[j] += u[m][n] * dg.weight(n,m) * (hx*hy);
//            }
//        }
        //_u[j] *= (1.0 + noise);
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
            const_this->mCrFfxWeightMatrix[m][n] = 0.0;
            for (unsigned int i=0; i<Nc; i++)
            {
                const DeltaGrid2D &dg = controlDeltaGrids[i];
                const_this->mCrFfxWeightMatrix[m][n] += _v[i] * dg.weight(n,m);
            }
        }
    }

    delete [] _v;
}

auto Problem2hDirichletBase::currentLayerBGrid(const DoubleMatrix &p, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids,
                                            double ln, const spif_vectorH &u_info) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    Problem2hDirichletBase* const_this = const_cast<Problem2hDirichletBase*>(this);
//    const_this->mCrBfxWeightMatrix.clear();
//    const_this->mCrBfxWeightMatrix.resize(M+1, N+1, 0.0);

    const unsigned int Nc = mEquParameter.Nc;
    const unsigned int No = mEquParameter.No;

    double* _p = new double[Nc];
    for (unsigned int i=0; i<Nc; i++)
    {
        _p[i] = controlDeltaGrids[i].consentrateInPoint(p);
//        const DeltaGrid2D &dg = controlDeltaGrids[i];
//        for (unsigned int m=dg.minY(); m<=dg.maxY(); m++)
//        {
//            for (unsigned int n=dg.minX(); n<=dg.maxX(); n++)
//            {
//                _p[i] += p[m][n] * (dg.weight(n,m) * (hx*hy));
//            }
//        }
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
            const_this->mCrBfxWeightMatrix[m][n] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                const DeltaGrid2D &dg = measurementDeltaGrids[j];
                const_this->mCrBfxWeightMatrix[m][n] += _w[j] * dg.weight(n,m);
            }
        }
    }

    delete [] _w;
}

void Problem2hDirichletBase::setGridDimensions(const Dimension &time, const Dimension &dimX, const Dimension &dimY)
{
    setTimeDimension(time);
    addSpaceDimension(dimX);
    addSpaceDimension(dimY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    mPulseWeightMatrix.resize(M+1, N+1);
    mCrFfxWeightMatrix.resize(M+1, N+1);
    mCrBfxWeightMatrix.resize(M+1, N+1);
}


