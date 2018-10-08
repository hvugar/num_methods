#include "problem2h_solver.h"
#include "problem2h_example.h"
#include <map>
#include <utility>

double MIN = +100000.0;
double MAX = -100000.0;
void Problem2HNDirichlet::f_layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const
{
    return;
    if (ln == 0 || ln == 2 || ln == 2*timeDimension().sizeN())
    {
        IPrinter::printSeperatorLine();
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
        printf("%f\n", sqrt(integralU(u)));
    }
    return;

    if (ln%2 == 0 && false)
    {
        char filename1[40];
        int size1 = sprintf(filename1, "e:/data/img/image%d.png", ln);
        filename1[size1] = 0;
        char filename2[40];
        int size2 = sprintf(filename2, "e:/data/txt/image%d.txt", ln);
        filename2[size2] = 0;

        double min = u.min();
        double max = u.max();
        if (MIN>min) MIN = min;
        if (MAX<max) MAX = max;

        puts("Generating image...");
        QPixmap pxm;
        visualGrayScale(u, min, max, pxm, 0, 0);
        pxm.save(QString(filename1), "PNG");
        printf("Image generated. ln: %d min: %f max: %f MIN: %f MAX: %f\n", ln/2, min, max, MIN, MAX);
        FILE* file = fopen(filename2, "w");
        IPrinter::print(u, u.rows(), u.cols(), 10, 8, file);
        fclose(file);
    }

    //    if (ln == 2*timeDimension().sizeN())
    //    {
    //        puts("Generating image...");
    //        QPixmap pxm;
    //        visualGrayScale(u, u.min(), u.max(), pxm, 0, 0);
    //        pxm.save("E:/image1000.png", "PNG");
    //        printf("Image generated. ln: %d min: %f max: %f\n", 1000, u.min(), u.max());
    //        FILE* file = fopen("E:/image1000.txt", "w");
    //        IPrinter::print(u, u.rows(), u.cols(), 10, 8, file);
    //        fclose(file);
    //    }

#ifdef SAVE_TO_IMG
    //if (ln != 1 && ln != timeDimension().sizeN()+LD) return;
    //if (ln < timeDimension().sizeN()) return;

    double min = u.min();
    double max = u.max();
    if (MIN>min) MIN = min;
    if (MAX<max) MAX = max;

    //    double norm = 0.0;
    //    for (unsigned int m=0; m<u.rows(); m++)
    //    {
    //        for (unsigned int n=0; n<u.cols(); n++)
    //        {
    //            norm += u[m][n]*u[m][n];
    //        }
    //    }

    QPixmap pic;
    visualizeMatrixHeat(u, min, max, pic);
    pic.save("images/f/100/pic_"+QString("%1").arg(ln)+".png", "PNG");
    printf("Layer: %4d min: %8.4f max: %8.4f min: %8.4f max: %8.4f diff: %8.4f norm: %10.8f\n", ln, min, max, MIN, MAX, fabs(max-min), sqrt(integralU(u)));
#endif
}

void Problem2HNDirichlet::checkGradient1(const Problem2HNDirichlet &prob)
{
    EquationParameter e_prm = prob.mEquParameter;
    OptimizeParameter o_prm = prob.mOptParameter;
    OptimizeParameter r_prm = prob.mRegParameter;

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

void Problem2HNDirichlet::checkGradient2(const Problem2HNDirichlet &prob)
{
    EquationParameter e_prm = prob.mEquParameter;
    OptimizeParameter o_prm = prob.mOptParameter;
    OptimizeParameter r_prm = prob.mRegParameter;

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

Problem2HNDirichlet::Problem2HNDirichlet()
{
    r = 0.0;
    regEpsilon = 0.0;
}

Problem2HNDirichlet::~Problem2HNDirichlet()
{}

double Problem2HNDirichlet::fx(const DoubleVector &pv) const
{
    OptimizeParameter o_prm;

    VectorToPrm(pv, o_prm);

    Problem2HNDirichlet* prob = const_cast<Problem2HNDirichlet*>(this);
    prob->mOptParameter = o_prm;

    std::vector<DoubleMatrix> u;
    spif_vector u_info;
    prob->solveForwardIBVP(u, u_info, true);

    double intgrl = integral(u);

    for (unsigned int i=0; i<=LD; i++) u[i].clear();
    u.clear();

    double nrm = norm(mEquParameter, o_prm, mRegParameter);
    double pnt = penalty(u_info, o_prm);

    double sum = intgrl + regEpsilon*nrm + r*pnt;

    for (unsigned int j=0; j<u_info.size(); j++)
    {
        u_info[j].clear();
    }
    u_info.clear();
    //printf("%f %f %f\n", intgrl, pnt, nrm);

    return sum;
}

double Problem2HNDirichlet::integral(const std::vector<DoubleMatrix> &vu) const
{
    const double ht = timeDimension().step();

    const unsigned int L = static_cast<const unsigned int>( timeDimension().sizeN() );
    const unsigned int LLD = L + LD;

    double sum = 0.0;

    const DoubleMatrix &uL = vu.at(0);
    sum += 0.5*integralU(uL);
    for (unsigned int l=L+1; l<=LLD-1; l++)
    {
        const DoubleMatrix &u = vu.at(2*(l-L));
        sum += integralU(u);
    }
    const DoubleMatrix &uLLD = vu.at(2*LD);
    sum += 0.5*integralU(uLLD);

    return sum*ht;
}

double Problem2HNDirichlet::integralU(const DoubleMatrix &u) const
{
    const double hx = spaceDimension(Dimension::DimensionX).step();
    const double hy = spaceDimension(Dimension::DimensionY).step();
    const unsigned int N = static_cast<const unsigned int> ( spaceDimension(Dimension::DimensionX).sizeN() );
    const unsigned int M = static_cast<const unsigned int> ( spaceDimension(Dimension::DimensionY).sizeN() );

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u[0][0]; usum += 0.25 * udiff * udiff;
    udiff = u[0][N]; usum += 0.25 * udiff * udiff;
    udiff = u[M][0]; usum += 0.25 * udiff * udiff;
    udiff = u[M][N]; usum += 0.25 * udiff * udiff;

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u[0][n]; usum += 0.5 * udiff * udiff;
        udiff = u[M][n]; usum += 0.5 * udiff * udiff;
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u[m][0]; usum += 0.5 * udiff * udiff;
        udiff = u[m][N]; usum += 0.5 * udiff * udiff;
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u[m][n]; usum += udiff * udiff;
        }
    }

    return usum*(hx*hy);
}

double Problem2HNDirichlet::norm(const EquationParameter& e_prm, const OptimizeParameter &o_prm, const OptimizeParameter &o_prm0) const
{
    double norm = 0.0;

    for (unsigned int i=0; i<e_prm.Nc; i++)
    {
        norm += (o_prm.eta[i].x - o_prm0.eta[i].x)*(o_prm.eta[i].x - o_prm0.eta[i].x);
        norm += (o_prm.eta[i].y - o_prm0.eta[i].y)*(o_prm.eta[i].y - o_prm0.eta[i].y);

        for (unsigned int j=0; j<e_prm.No; j++)
        {
            norm += (o_prm.k[i][j] - o_prm0.k[i][j])*(o_prm.k[i][j] - o_prm0.k[i][j]);
            norm += (o_prm.z[i][j] - o_prm0.z[i][j])*(o_prm.z[i][j] - o_prm0.z[i][j]);

            norm += (o_prm.xi[j].x - o_prm0.xi[j].x)*(o_prm.xi[j].x - o_prm0.xi[j].x);
            norm += (o_prm.xi[j].y - o_prm0.xi[j].y)*(o_prm.xi[j].y - o_prm0.xi[j].y);
        }
    }

    return norm;
}

double Problem2HNDirichlet::penalty(const spif_vector &info, const OptimizeParameter &o_prm) const
{
    const double ht = mtimeDimension.step();
    const unsigned int L = static_cast<const unsigned int> ( mtimeDimension.sizeN() );

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

double Problem2HNDirichlet::gpi(unsigned int i, unsigned int layer, const spif_vector &u_info, const OptimizeParameter &o_prm) const
{
    double gpi_ln = fabs( g0i(i, layer, u_info, o_prm) ) - ( vmax.at(i) - vmin.at(i) )/2.0;
    return gpi_ln > 0.0 ? gpi_ln : 0.0;
}

double Problem2HNDirichlet::g0i(unsigned int i, unsigned int layer, const spif_vector &u_info, const OptimizeParameter &o_prm) const
{
    double vi = 0.0;
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        const SpacePointInfo &u_xij = u_info[j];
        vi += o_prm.k[i][j] * ( u_xij.vl[layer] - o_prm.z[i][j] );
    }
    return ( vmax.at(i) + vmin.at(i) )/2.0 - vi;
}

double Problem2HNDirichlet::sign(double x) const
{
    if (x < 0.0)       return -1.0;
    else if (x > 0.0)  return +1.0;
    else               return  0.0;
}

void Problem2HNDirichlet::gradient(const DoubleVector & pv, DoubleVector &g) const
{
    const unsigned int L   = mtimeDimension.sizeN();
    const double ht        = mtimeDimension.step();
    const unsigned int Nc  = mEquParameter.Nc;
    const unsigned int No  = mEquParameter.No;
    const unsigned int LLD = L + LD;

    OptimizeParameter o_prm;
    VectorToPrm(pv, o_prm);

    Problem2HNDirichlet* prob = const_cast<Problem2HNDirichlet*>(this);
    prob->mOptParameter = o_prm;

    std::vector<DoubleMatrix> u;

    spif_vector u_info;
    solveForwardIBVP(u, u_info, true);
    spif_vector p_info;
    solveBackwardIBVP(u, p_info, true, u_info);

    g.clear();
    g.resize(pv.length(), 0.0);
    unsigned int gi = 0;

    // k
    if (optimizeK)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointInfo &uj = u_info[j];

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
            const SpacePointInfo &pi = p_info[i];

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
            const SpacePointInfo &uj = u_info[j];

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
            const SpacePointInfo &pi = p_info[i];

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

auto Problem2HNDirichlet::norm(const DoubleVector &v) const -> double
{
    return EuclideanNorm(v);
}

auto Problem2HNDirichlet::normalize(DoubleVector &v) const -> void
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

void Problem2HNDirichlet::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const
{
    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(f); C_UNUSED(alpha); C_UNUSED(result);
    const char* msg = 0; C_UNUSED(msg);
    if (result == GradientMethod::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION   ";
    if (result == GradientMethod::FIRST_ITERATION)          msg = "FIRST_ITERATION         ";
    if (result == GradientMethod::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    if (result == GradientMethod::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS     ";
    if (result == GradientMethod::NEXT_ITERATION)           msg = "NEXT_ITERATION          ";

    Problem2HNDirichlet* prob = const_cast<Problem2HNDirichlet*>(this);
    OptimizeParameter o_prm;
    VectorToPrm(x, o_prm);
    prob->mOptParameter = o_prm;
    std::vector<DoubleMatrix> u;
    spif_vector u_info;
    solveForwardIBVP(u, u_info, true);
    double ing = integral(u);
    double pnt = penalty(u_info, o_prm);
    double nrm = norm(prob->mEquParameter, prob->mOptParameter, prob->mRegParameter);

    unsigned int v_length = timeDimension().sizeN() + LD;
    DoubleVector v1(v_length+1);
    DoubleVector v2(v_length+1);

    for (unsigned int ln=0; ln<=v_length; ln++)
    {
        v1[ln] = o_prm.k[0][0] * (u_info.at(0).vl[2*ln] - o_prm.z[0][0]) + o_prm.k[0][1] * (u_info.at(1).vl[2*ln] - o_prm.z[0][1]);
        v2[ln] = o_prm.k[1][0] * (u_info.at(0).vl[2*ln] - o_prm.z[1][0]) + o_prm.k[1][1] * (u_info.at(1).vl[2*ln] - o_prm.z[1][1]);
    }

    IPrinter::printVector(v1, "v1", 10);
    IPrinter::printVector(v2, "v2", 10);

    printf("I[%3d]: F:%10.6f I:%10.6f P:%12.6f N:%10.6f R:%7.3f e:%5.3f a:%10.6f min:%10.6f max:%10.6f \n",
           i, f, ing, pnt, nrm, r, regEpsilon, alpha, u.at(0).min(), u.at(0).max());
    //if (result == GradientMethod::BREAK_GRADIENT_NORM_LESS || result == GradientMethod::BREAK_DISTANCE_LESS)
    printf("k:%7.4f %7.4f %7.4f %7.4f z:%7.4f %7.4f %7.4f %7.4f o: %6.4f %6.4f %6.4f %6.4f c: %6.4f %6.4f %6.4f %6.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    //printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o:%8.4f %8.4f %8.4f %8.4f c:%8.4f %8.4f %8.4f %8.4f\n", g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
    //DoubleVector n = g;
    //n.L2Normalize();
    //printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o:%8.4f %8.4f %8.4f %8.4f c:%8.4f %8.4f %8.4f %8.4f\n", n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15]);

    u.clear();
    u_info.clear();

    C_UNUSED(prob);
    IPrinter::printSeperatorLine();

    //    prob->optimizeK = i%4 == 3;
    //    prob->optimizeZ = i%4 == 0;
    //    prob->optimizeC = i%4 == 1;
    //    prob->optimizeO = i%4 == 2;
    //    if (alpha > 0.00001) prob->gm->setR1MinimizeEpsilon(alpha, 0.0001);
    //    exit(-1);
}

auto Problem2HNDirichlet::project(DoubleVector &pv, unsigned int index) -> void
{
    C_UNUSED(pv);
    C_UNUSED(index);
    return;

    //    unsigned int Nc = mEquParameter.Nc;
    //    unsigned int No = mEquParameter.No;

    //    unsigned int offset = 2*Nc*No;

    //    // xi
    //    if ( offset <= index && index <= offset + 2*No - 1 )
    //    {
    //        if (pv[index] < 0.05) pv[index] = 0.05;
    //        if (pv[index] > 0.95) pv[index] = 0.95;
    //    }

    //    // eta
    //    if ( offset + 2*No <= index && index <= offset + 2*No + 2*Nc - 1 )
    //    {
    //        if (pv[index] < 0.05) pv[index] = 0.05;
    //        if (pv[index] > 0.95) pv[index] = 0.95;
    //    }

    //return;
    //projectControlPoints(pv, index);
    //projectMeasurePoints(pv, index);
}

auto Problem2HNDirichlet::project(DoubleVector &pv) const -> void
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

auto Problem2HNDirichlet::projectControlPoints(DoubleVector &pv, unsigned int index) const -> void
{
    double dist = 0.10;

    if (index == 12)
    {
        if (fabs(pv[12] - pv[ 8])<dist) pv[12] = pv[ 8] + sign(pv[12] - pv[ 8])*dist;
        if (fabs(pv[12] - pv[10])<dist) pv[12] = pv[10] + sign(pv[12] - pv[10])*dist;
        if (pv[12] < 0.05) pv[12] = 0.05;
        if (pv[12] > 0.95) pv[12] = 0.95;

        //        if (fabs(pv[12] - pv[8])<dx)
        //        {
        //            pv[12] = pv[8] + dx;
        //            if (pv[12] > 0.95) pv[12] = pv[8] - dx;
        //        }
        //        if (fabs(pv[12] - pv[10])<dx)
        //        {
        //            pv[12] = pv[10] + dx;
        //            if (pv[12] > 0.95) pv[12] = pv[10] - dx;
        //        }
    }

    if (index == 13)
    {
        if (fabs(pv[13] - pv[ 9])<dist) pv[13] = pv[ 9] + sign(pv[13] - pv[ 9])*dist;
        if (fabs(pv[13] - pv[11])<dist) pv[13] = pv[11] + sign(pv[13] - pv[11])*dist;
        if (pv[13] < 0.05) pv[13] = 0.05;
        if (pv[13] > 0.95) pv[13] = 0.95;

        //        if (fabs(pv[13] - pv[9])<dx)
        //        {
        //            pv[13] = pv[9] + dx;
        //            if (pv[13] > 0.95) pv[13] = pv[9] - dx;
        //        }
        //        if (fabs(pv[13] - pv[11])<dx)
        //        {
        //            pv[13] = pv[11] + dx;
        //            if (pv[13] > 0.95) pv[13] = pv[11] - dx;
        //        }
    }

    if (index == 14)
    {
        if (fabs(pv[14] - pv[ 8])<dist) pv[14] = pv[ 8] + sign(pv[14] - pv[ 8])*dist;
        if (fabs(pv[14] - pv[10])<dist) pv[14] = pv[10] + sign(pv[14] - pv[10])*dist;
        if (pv[14] < 0.05) pv[14] = 0.05;
        if (pv[14] > 0.95) pv[14] = 0.95;

        //        if (fabs(pv[14] - pv[8])<dx)
        //        {
        //            pv[14] = pv[8] + dx;
        //            if (pv[14] > 0.95) pv[14] = pv[8] - dx;
        //        }
        //        if (fabs(pv[14] - pv[10])<dx)
        //        {
        //            pv[14] = pv[10] + dx;
        //            if (pv[14] > 0.95) pv[14] = pv[10] - dx;
        //        }
    }

    if (index == 15)
    {
        if (fabs(pv[15] - pv[ 9])<dist) pv[15] = pv[ 9] + sign(pv[15] - pv[ 9])*dist;
        if (fabs(pv[15] - pv[11])<dist) pv[15] = pv[11] + sign(pv[15] - pv[11])*dist;
        if (pv[15] < 0.05) pv[15] = 0.05;
        if (pv[15] > 0.95) pv[15] = 0.95;

        //        if (fabs(pv[15] - pv[9])<dx)
        //        {
        //            pv[15] = pv[9] + dx;
        //            if (pv[15] > 0.95) pv[15] = pv[9] - dx;
        //        }
        //        if (fabs(pv[15] - pv[11])<dx)
        //        {
        //            pv[15] = pv[11] + dx;
        //            if (pv[15] > 0.95) pv[15] = pv[11] - dx;
        //        }
    }
}

auto Problem2HNDirichlet::projectMeasurePoints(DoubleVector &pv, unsigned int index) const -> void
{
    double dist = 0.10;

    if (index == 8)
    {
        if (fabs(pv[ 8] - pv[12])<dist) pv[8] = pv[12] + sign(pv[ 8] - pv[12])*dist;
        if (fabs(pv[ 8] - pv[14])<dist) pv[8] = pv[14] + sign(pv[ 8] - pv[14])*dist;
        if (pv[8] < 0.05) pv[8] = 0.05;
        if (pv[8] > 0.95) pv[8] = 0.95;

        //        if ( fabs(pv[8] - pv[12]) < dist)
        //        {
        //            pv[8] = pv[12] + dist;
        //            if (pv[8] > 0.95) pv[8] = pv[12] - dist;
        //        }
        //        if (fabs(pv[8] - pv[14])<dist)
        //        {
        //            pv[8] = pv[14] + dist;
        //            if (pv[8] > 0.95) pv[8] = pv[14] - dist;
        //        }
    }

    if (index == 9)
    {
        if (fabs(pv[ 9] - pv[13])<dist) pv[9] = pv[13] + sign(pv[ 9] - pv[13])*dist;
        if (fabs(pv[ 9] - pv[15])<dist) pv[9] = pv[15] + sign(pv[ 9] - pv[15])*dist;
        if (pv[9] < 0.05) pv[9] = 0.05;
        if (pv[9] > 0.95) pv[9] = 0.95;

        //        if (fabs(pv[9] - pv[13])<dist)
        //        {
        //            pv[9] = pv[13] + dist;
        //            if (pv[9] > 0.95) pv[9] = pv[13] - dist;
        //        }
        //        if (fabs(pv[9] - pv[15])<dist)
        //        {
        //            pv[9] = pv[15] + dist;
        //            if (pv[9] > 0.95) pv[9] = pv[15] - dist;
        //        }
    }

    if (index == 10)
    {
        if (fabs(pv[10] - pv[12])<dist) pv[10] = pv[12] + sign(pv[10] - pv[12])*dist;
        if (fabs(pv[10] - pv[14])<dist) pv[10] = pv[14] + sign(pv[10] - pv[14])*dist;
        if (pv[10] < 0.05) pv[10] = 0.05;
        if (pv[10] > 0.95) pv[10] = 0.95;

        //        if (fabs(pv[10] - pv[12])<dist)
        //        {
        //            pv[10] = pv[12] + dist;
        //            if (pv[10] > 0.95) pv[10] = pv[12] - dist;
        //        }
        //        if (fabs(pv[10] - pv[14])<dist)
        //        {
        //            pv[10] = pv[14] + dist;
        //            if (pv[10] > 0.95) pv[10] = pv[12] - dist;
        //        }
    }

    if (index == 11)
    {
        if (fabs(pv[11] - pv[13])<dist) pv[11] = pv[13] + sign(pv[11] - pv[13])*dist;
        if (fabs(pv[11] - pv[15])<dist) pv[11] = pv[15] + sign(pv[11] - pv[15])*dist;
        if (pv[11] < 0.05) pv[11] = 0.05;
        if (pv[11] > 0.95) pv[11] = 0.95;

        //        if (fabs(pv[11] - pv[13])<dist)
        //        {
        //            pv[11] = pv[13] + dist;
        //            if (pv[11] > 0.95) pv[11] = pv[13] - dist;
        //        }
        //        if (fabs(pv[11] - pv[15])<dist)
        //        {
        //            pv[11] = pv[15] + dist;
        //            if (pv[11] > 0.95) pv[11] = pv[15] - dist;
        //        }
    }
}

auto Problem2HNDirichlet::solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vector &u_info, bool use) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<const unsigned int> ( dimX.sizeN() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.sizeN() );
    const unsigned int L = static_cast<const unsigned int> ( time.sizeN() );
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a        = mEquParameter.a;
    const double lambda   = mEquParameter.lambda;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;
    const unsigned int Ns = mEquParameter.Ns;

    const double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    const double p_aa_htht__hxhx___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hx*hx) + (lambda*ht);
    const double p_aa_htht__hyhy = +(a*a*ht*ht)/(hy*hy);

    const double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    const double p_aa_htht__hyhy___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hy*hy) + (lambda*ht);
    const double p_aa_htht__hxhx = +(a*a*ht*ht)/(hx*hx);

    const double htht = ht*ht;
    const double lambda_ht = lambda*ht;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    for (unsigned int l=0; l<u.size(); l++) u[l].clear(); u.clear();
    unsigned int u_size = 2*LD + 1;
    u.resize(u_size); for (unsigned int l=0; l<u_size; l++) u[l].resize(M+1, N+1);

    //----------------------------------------------------------------------------------------------//
    espn_vector obsPointNodes, cntDeltaNodes, qPointNodes;
    std::vector<ExtendedSpacePoint> msnExtSpacePoints, cntExtSpacePoints, qExtSpacePoints;
#if defined(OLD_VERSION) || defined(NEW_VERSION)
    for (unsigned int j=0; j<No; j++) distributeDelta0(mOptParameter.xi[j], j, obsPointNodes, dimX, dimY, 4);
    for (unsigned int i=0; i<Nc; i++) distributeDelta0(mOptParameter.eta[i], i, cntDeltaNodes, dimX, dimY, 4);
    for (unsigned int s=0; s<Ns; s++) distributeDeltaGaussPulse(mEquParameter.theta[s], s, qPointNodes, dimX, dimY);
#endif
#ifdef NEW_VERSION
    newDistributeDeltaGaussPulse(mEquParameter.theta, qExtSpacePoints, dimX, dimY);
    newDistributeDeltaGauseCntrl(mOptParameter.eta, cntExtSpacePoints, dimX, dimY);
    newDistributeDeltaGauseMsmnt(mOptParameter.xi,  msnExtSpacePoints, dimX, dimY);
#endif

    //----------------------------------------------------------------------------------------------//
    uint_vector rows0, rows1, rows2, cols0, cols1, cols2;
    f_findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, obsPointNodes, cntDeltaNodes, N, M);
    //cntDeltaNodes.clear();
    //obsPointNodes.clear();;
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) f_prepareInfo(No, mOptParameter.xi, u_info, LLD, dimX, dimY);
    //----------------------------------------------------------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    f_initialLayers(u00, u10, u_info, use, obsPointNodes, cntDeltaNodes, qPointNodes, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, lambda,
                    qExtSpacePoints, msnExtSpacePoints);
    //------------------------------------- initial conditions -------------------------------------//

    double *ax = static_cast<double*>( malloc(sizeof(double)*(N-1)) ); for (unsigned int n=1; n<=N-1; n++) ax[n-1] = m_aa_htht__hxhx;
    double *bx = static_cast<double*>( malloc(sizeof(double)*(N-1)) ); for (unsigned int n=1; n<=N-1; n++) bx[n-1] = p_aa_htht__hxhx___lambda_ht;
    double *cx = static_cast<double*>( malloc(sizeof(double)*(N-1)) ); for (unsigned int n=1; n<=N-1; n++) cx[n-1] = m_aa_htht__hxhx;
    double *dx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *rx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    ax[0] = cx[N-2] = 0.0;

    double *ay = static_cast<double*>( malloc(sizeof(double)*(M-1)) ); for (unsigned int m=1; m<=M-1; m++) ay[m-1] = m_aa_htht__hyhy;
    double *by = static_cast<double*>( malloc(sizeof(double)*(M-1)) ); for (unsigned int m=1; m<=M-1; m++) by[m-1] = p_aa_htht__hyhy___lambda_ht;
    double *cy = static_cast<double*>( malloc(sizeof(double)*(M-1)) ); for (unsigned int m=1; m<=M-1; m++) cy[m-1] = m_aa_htht__hyhy;
    double *dy = static_cast<double*>( malloc(sizeof(double)*(M-1)) );
    double *ry = static_cast<double*>( malloc(sizeof(double)*(M-1)) );
    ay[0] = cy[M-2] = 0.0;

    size_t rows1_size = rows1.size()*(N-1);
    double *a1, *b1, *c1, *d1, *x1, **w1;
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

    size_t cols1_size = cols1.size()*(M-1);
    double *a2, *b2, *c2, *d2, *x2, **w2;
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

    for (unsigned int l=2; l<=LLD; l++)
    {
        TimeNodePDE tn20; tn20.i = l; tn20.t = l*ht;
        TimeNodePDE tn15; tn15.i = l; tn15.t = l*ht-0.5*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0; sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy; u15[m][0] = f_boundary(sn0, tn15); u20[m][0] = f_boundary(sn0, tn20);
            sn1.j = m; sn1.y = m*hy; u15[m][N] = f_boundary(sn1, tn15); u20[m][N] = f_boundary(sn1, tn20);
        }

        sn0.j = 0; sn0.y = 0.0; sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx; u15[0][n] = f_boundary(sn0, tn15); u20[0][n] = f_boundary(sn0, tn20);
            sn1.i = n; sn1.x = n*hx; u15[M][n] = f_boundary(sn1, tn15); u20[M][n] = f_boundary(sn1, tn20);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = rows0.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    dx[n-1] += lambda_ht*(u10[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    dx[n-1] += 2.0*u10[m][n] + (u10[m][n]-u00[m][n]);

                    //------------------------------------- Adding time delta part --------------------------------//
                    dx[n-1] += htht*distributeTimeDelta(tn15.t, ht, 2*l-1, qPointNodes, sn, qExtSpacePoints);
                    //------------------------------------- Adding time delta part --------------------------------//
                }

                dx[0]   -= m_aa_htht__hxhx * u15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * u15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0 && rows2.size() == 0)
        {
            //throw std::exception();
            double *_v15 = new double[Nc];

            double* _u15 = new double[No]; for (unsigned int j=0; j<No; j++) _u15[j] = 0.0;
#ifdef OLD_VERSION
            for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[opj];
                _u15[opn.id] += u15[opn.j][opn.i] * (opn.w * (hx*hy));
            }
#endif
#ifdef NEW_VERSION1
            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointExt &spx = obsPointNodeExts.at(j);
                const SpacePointExt::GridNodeMap &map = spx.distPoints();
                for (SpacePointExt::GridNodeMap::const_iterator it=map.begin(); it != map.end(); it++)
                {
                    const SpacePointExt::GridNodeWeight &gnw = it->second;
                    _u15[gnw.id] += u15[gnw.j][gnw.i] * (gnw.w * (hx*hy));
                }
            }
#endif
#ifdef NEW_VERSION
            for (unsigned int j=0; j<No; j++)
            {
                const ExtendedSpacePoint &extendedSpacePoint = msnExtSpacePoints.at(j);
                const std::vector<ExtendedSpacePointNode1> &nodes = extendedSpacePoint.nodes;
                unsigned int nodes_size = nodes.size();
                for (unsigned int ni=0; ni<nodes_size; ni++)
                {
                    const ExtendedSpacePointNode1 &node = nodes.at(ni);
                    _u15[j] += u15[node.ny][node.nx] * (node.w * (hx*hy));
                }
            }
#endif
            for (unsigned int i=0; i<Nc; i++)
            {
                _v15[i] = 0.0;
                for (unsigned int j=0; j<No; j++)
                {
                    _v15[i] += mOptParameter.k[i][j] * (_u15[j] - mOptParameter.z[i][j]);
                }
            }
            delete [] _u15;

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    dx[n-1] += lambda_ht*(u10[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    dx[n-1] += 2.0*u10[m][n] + (u10[m][n]-u00[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
#ifdef OLD_VERSION
                    for (unsigned int cdi=0; cdi<cntDeltaNodes.size(); cdi++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cdi);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            dx[n-1] += htht * _v15[cdn.id] * cdn.w;
                        }
                    }
#endif
#ifdef NEW_VERSION1
                    for (unsigned int c=0; c<Nc; c++)
                    {
                        const SpacePointExt &spx = cntDeltaNodeExts.at(c);
                        SpacePointExt::GridNodePair node(sn.j, sn.i);
                        SpacePointExt::GridNodeMap::const_iterator it = spx.distPoints().find( node);
                        dx[n-1] += htht * _v15[c] * it->second.w;
                    }
#endif
#ifdef NEW_VERSION
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePoint &extendedSpacePoint = cntExtSpacePoints.at(i);
                        if (extendedSpacePoint.contains(sn))
                        {
                            const std::vector<ExtendedSpacePointNode1> &nodes = extendedSpacePoint.nodes;
                            unsigned int nodes_size = nodes.size();
                            for (unsigned int ni=0; ni<nodes_size; ni++)
                            {
                                const ExtendedSpacePointNode1 &node = nodes.at(ni);
                                if (node.equals(sn)) dx[n-1] += htht * _v15[i] * node.w;
                            }
                        }
                    }
#endif
                    //------------------------------------- Adding delta part -------------------------------------//

                    //------------------------------------- Adding time delta part --------------------------------//
                    dx[n-1] += htht*distributeTimeDelta(tn15.t, ht, 2*l-1, qPointNodes, sn, qExtSpacePoints);
                    //------------------------------------- Adding time delta part --------------------------------//
                }

                dx[0]   -= m_aa_htht__hxhx * u15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * u15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }

            delete [] _v15;
        }

        if (rows1.size() != 0 && rows2.size() != 0)
        {
            //throw std::exception();

            for (unsigned int i=0; i < rows1_size; i++) for (unsigned int j=0; j < rows1_size; j++) w1[i][j] = 0.0;

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    const unsigned int index = offset+(n-1);
                    d1[index] = 0.0;
                    if (m>0 && m<M) d1[index] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    else if (m==0)  d1[index] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    else if (m==M)  d1[index] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1[index] += lambda_ht*(u10[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    d1[index] += 2.0*u10[m][n] + (u10[m][n]-u00[m][n]);

                    a1[index] = m_aa_htht__hxhx;
                    b1[index] = p_aa_htht__hxhx___lambda_ht;
                    c1[index] = m_aa_htht__hxhx;

                    //------------------------------------- Adding delta part -------------------------------------//
#ifdef OLD_VERSION
                    for (unsigned int cdi=0; cdi<cntDeltaNodes.size(); cdi++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cdi);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(opj);

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (opn.j == rows1[rs])
                                    {
                                        found = true;
                                        w1[index][rs*(N-1)+(opn.i-1)] -= htht * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d1[index] += htht * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[index] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
#endif
#ifdef NEW_VERSION
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePoint &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                        if (cExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNode1> &nodes1 = cExtendedSpacePoint.nodes;
                            for (unsigned int ni=0; ni<nodes1.size(); ni++)
                            {
                                const ExtendedSpacePointNode1 &node1 = nodes1.at(ni);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                const ExtendedSpacePoint &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                                const std::vector<ExtendedSpacePointNode1> &nodes2 = mExtendedSpacePoint.nodes;
                                for (unsigned int nj=0; nj<nodes2.size(); nj++)
                                {
                                    const ExtendedSpacePointNode1 &node2 = nodes2.at(nj);

                                    bool found = false;
                                    for (unsigned int rs=0; rs<rows1.size(); rs++)
                                    {
                                        if (node2.ny == rows1[rs])
                                        {
                                            found = true;
                                            w1[index][rs*(N-1)+(node2.nx-1)] -= htht * mOptParameter.k[i][j] * (node2.w * (hx*hy)) * w;
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d1[index] += htht * mOptParameter.k[i][j] * u15[node2.ny][node2.nx] * (node2.w * (hx*hy)) * w;
                                    }
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[index] -= htht * mOptParameter.k[i][j] * mOptParameter.z[i][j] * w;
                            }
                        }
                    }
#endif
                    //------------------------------------- Adding delta part -------------------------------------//

                    //------------------------------------- Adding time delta part --------------------------------//
                    d1[index] += htht*distributeTimeDelta(tn15.t, ht, 2*l-1, qPointNodes, sn, qExtSpacePoints);
                    //------------------------------------- Adding time delta part --------------------------------//
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx * u15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx * u15[m][N];

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w1, x1, rows1_size);

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
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        if (cols0.size() != 0)
        {
            for (unsigned int col=0; col<cols0.size(); col++)
            {
                unsigned int n = cols0.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    dy[m-1] += lambda_ht*(u15[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    dy[m-1] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]);

                    //------------------------------------- Adding time delta part --------------------------------//
                    dy[m-1] += htht*distributeTimeDelta(tn20.t, ht, 2*l, qPointNodes, sn, qExtSpacePoints);
                    //------------------------------------- Adding time delta part --------------------------------//
                }

                dy[0]   -= m_aa_htht__hyhy * u20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * u20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0 && cols2.size() == 0)
        {
            //throw std::exception();
            double *_v20 = new double[Nc];

            double* _u20 = new double[No]; for (unsigned int j=0; j<No; j++) _u20[j] = 0.0;
#ifdef OLD_VERSION
            for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[opj];
                _u20[opn.id] += u20[opn.j][opn.i] * (opn.w * (hx*hy));
            }
#endif
#ifdef NEW_VERSION1
            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointExt &spx = obsPointNodeExts.at(j);
                const SpacePointExt::GridNodeMap &map = spx.distPoints();
                for (SpacePointExt::GridNodeMap::const_iterator it=map.begin(); it != map.end(); it++)
                {
                    const SpacePointExt::GridNodeWeight &gnw = it->second;
                    _u20[gnw.id] += u20[gnw.j][gnw.i] * (gnw.w * (hx*hy));
                }
            }
#endif
#ifdef NEW_VERSION
            for (unsigned int j=0; j<No; j++)
            {
                const ExtendedSpacePoint &extendedSpacePoint = msnExtSpacePoints.at(j);
                const std::vector<ExtendedSpacePointNode1> &nodes = extendedSpacePoint.nodes;
                unsigned int nodes_size = nodes.size();
                for (unsigned int ni=0; ni<nodes_size; ni++)
                {
                    const ExtendedSpacePointNode1 &node = nodes.at(ni);
                    _u20[j] += u20[node.ny][node.nx] * (node.w * (hx*hy));
                }
            }
#endif
            for (unsigned int i=0; i<Nc; i++)
            {
                _v20[i] = 0.0;
                for (unsigned int j=0; j<No; j++)
                {
                    _v20[i] += mOptParameter.k[i][j] * ( _u20[j] - mOptParameter.z[i][j] );
                }
            }
            delete [] _u20;

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    dy[m-1] += lambda_ht*(u15[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    dy[m-1] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
#ifdef OLD_VERSION
                    for (unsigned int cdi=0; cdi<cntDeltaNodes.size(); cdi++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cdi);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            dy[m-1] += htht * _v20[cdn.id] * cdn.w;
                        }
                    }
#endif
#ifdef NEW_VERSION1
                    for (unsigned int c=0; c<Nc; c++)
                    {
                        const SpacePointExt &spx = cntDeltaNodeExts.at(c);
                        SpacePointExt::GridNodePair node(sn.j, sn.i);
                        SpacePointExt::GridNodeMap::const_iterator it = spx.distPoints().find( node);
                        dy[m-1] += htht * _v20[c] * it->second.w;
                    }
#endif
#ifdef NEW_VERSION
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePoint &extendedSpacePoint = cntExtSpacePoints.at(i);
                        if (extendedSpacePoint.contains(sn))
                        {
                            const std::vector<ExtendedSpacePointNode1> &nodes = extendedSpacePoint.nodes;
                            unsigned int nodes_size = nodes.size();
                            for (unsigned int ni=0; ni<nodes_size; ni++)
                            {
                                const ExtendedSpacePointNode1 &node = nodes.at(ni);
                                if (node.equals(sn)) dy[m-1] += htht * _v20[i] * node.w;
                            }
                        }
                    }
#endif
                    //------------------------------------- Adding delta part -------------------------------------//

                    //------------------------------------- Adding time delta part --------------------------------//
                    dy[m-1] += htht*distributeTimeDelta(tn20.t, ht, 2*l, qPointNodes, sn, qExtSpacePoints);
                    //------------------------------------- Adding time delta part --------------------------------//
                }

                dy[0]   -= m_aa_htht__hyhy * u20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * u20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
            }

            delete [] _v20;
        }

        if (cols1.size() != 0 && cols2.size() != 0)
        {
            //throw std::exception();

            for (unsigned int i=0; i < cols1_size; i++) for (unsigned int j=0; j < cols1_size; j++) w2[i][j] = 0.0;

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    const unsigned int index = offset+(m-1);
                    d2[index] = 0.0;
                    if (n>0 && n<N) d2[index] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    else if (n==0)  d2[index] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    else if (n==N)  d2[index] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d2[index] += lambda_ht*(u15[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    d2[index] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]);

                    a2[index] = m_aa_htht__hyhy;
                    b2[index] = p_aa_htht__hyhy___lambda_ht;
                    c2[index] = m_aa_htht__hyhy;

                    //------------------------------------- Adding delta part -------------------------------------//
#ifdef OLD_VERSION
                    for (unsigned int cdi=0; cdi<cntDeltaNodes.size(); cdi++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cdi);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(opj);

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (opn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[index][cs*(M-1)+(opn.j-1)] -= htht * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d2[index] += htht * mOptParameter.k[cdn.id][opn.id] * u20[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[index] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
#endif
#ifdef NEW_VERSION1
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const SpacePointExt &cspx = cntDeltaNodeExts.at(i);
                        const SpacePointExt::GridNodeMap &cmap = cspx.distPoints();
                        SpacePointExt::GridNodePair cnode(sn.j, sn.i);
                        SpacePointExt::GridNodeMap::const_iterator cit = cmap.find( cnode);
                        if ( cit != cmap.end() )
                        {
                            double cw = cit->second.w;

                            for (unsigned int j=0; j<No; j++)
                            {
                                const SpacePointExt &ospx = obsPointNodeExts.at(j);
                                const SpacePointExt::GridNodeMap &omap = ospx.distPoints();
                                for (SpacePointExt::GridNodeMap::const_iterator oit=omap.begin(); oit!=omap.end(); oit++)
                                {

                                    double ow = oit->second.w;

                                    bool found = false;
                                    for (unsigned int cs=0; cs<cols1.size(); cs++)
                                    {
                                        if (oit->second.i == cols1[cs])
                                        {
                                            found = true;
                                            w2[index][cs*(M-1)+(oit->second.j-1)] -= htht * mOptParameter.k[i][j] * (ow * (hx*hy)) * cw;
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d2[index] += htht * mOptParameter.k[i][j] * u20[oit->second.j][oit->second.i] * (ow * (hx*hy)) * cw;
                                    }
                                }
                            }

                            //                            for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
                            //                            {
                            //                                const ExtendedSpacePointNode &opn = obsPointNodes.at(opj);

                            //                                bool found = false;
                            //                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                            //                                {
                            //                                    if (opn.i == cols1[cs])
                            //                                    {
                            //                                        found = true;
                            //                                        w2[index][cs*(M-1)+(opn.j-1)] -= htht * mOptParameter.k[i][opn.id] * (opn.w * (hx*hy)) * cw;
                            //                                        break;
                            //                                    }
                            //                                }

                            //                                if (!found)
                            //                                {
                            //                                    d2[index] += htht * mOptParameter.k[i][opn.id] * u20[opn.j][opn.i] * (opn.w * (hx*hy)) * cw;
                            //                                }
                            //                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[index] -= htht * mOptParameter.k[i][j] * mOptParameter.z[i][j] * cw;
                            }
                        }
                    }
#endif
#ifdef NEW_VERSION
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePoint &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                        if (cExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNode1> &nodes1 = cExtendedSpacePoint.nodes;
                            for (unsigned int ni=0; ni<nodes1.size(); ni++)
                            {
                                const ExtendedSpacePointNode1 &node1 = nodes1.at(ni);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                const ExtendedSpacePoint &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                                const std::vector<ExtendedSpacePointNode1> &nodes2 = mExtendedSpacePoint.nodes;
                                for (unsigned int nj=0; nj<nodes2.size(); nj++)
                                {
                                    const ExtendedSpacePointNode1 &node2 = nodes2.at(nj);

                                    bool found = false;
                                    for (unsigned int cs=0; cs<cols1.size(); cs++)
                                    {
                                        if (node2.nx == cols1[cs])
                                        {
                                            found = true;
                                            w2[index][cs*(M-1)+(node2.ny-1)] -= htht * mOptParameter.k[i][j] * (node2.w * (hx*hy)) * w;
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d2[index] += htht * mOptParameter.k[i][j] * u20[node2.ny][node2.nx] * (node2.w * (hx*hy)) * w;
                                    }
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[index] -= htht * mOptParameter.k[i][j] * mOptParameter.z[i][j] * w;
                            }
                        }
                    }
#endif
                    //------------------------------------- Adding delta part -------------------------------------//

                    //------------------------------------- Adding time delta part --------------------------------//
                    d2[index] += htht*distributeTimeDelta(tn20.t, ht, 2*l, qPointNodes, sn, qExtSpacePoints);
                    //------------------------------------- Adding time delta part --------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy * u20[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy * u20[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2, x2, cols1_size);

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    u20[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }
        }

        /**************************************************** y direction apprx ***************************************************/

        if (use == true) f_add2Info(u15, u_info, obsPointNodes, 2*l-1, hx, hy, msnExtSpacePoints); f_layerInfo(u15, 2*l-1);
        if (use == true) f_add2Info(u20, u_info, obsPointNodes, 2*l+0, hx, hy, msnExtSpacePoints); f_layerInfo(u20, 2*l+0);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }

        /**************************************************** saving last LD layers ***********************************************/

        if (L == l)
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    u[l-L][m][n] = u20[m][n];
                }
            }
        }

        if ( L+1 <= l && l <= LLD )
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    u[2*(l-L)-1][m][n] = u15[m][n];
                    u[2*(l-L)+0][m][n] = u20[m][n];
                }
            }
        }

        /**************************************************** saving last LD layers ***********************************************/
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

    qPointNodes.clear();
    obsPointNodes.clear();
    cntDeltaNodes.clear();

    u00.clear();
    u10.clear();
    u15.clear();
    u20.clear();
}

void Problem2HNDirichlet::f_initialLayers(DoubleMatrix &u00, DoubleMatrix &u10, spif_vector &u_info, bool use,
                                          espn_vector &obsPointNodes, espn_vector &cntDeltaNodes UNUSED_PARAM,
                                          espn_vector &qPointNodes UNUSED_PARAM, unsigned int N, unsigned int M,
                                          double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda,
                                          const std::vector<ExtendedSpacePoint> &qExtSpacePoints,
                                          std::vector<ExtendedSpacePoint> &msnExtSpacePoints) const
{
    DoubleMatrix u05 = u00;



    /************************************************************************/
    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u00[m][n] = f_initial1(sn);
        }
    }

    /************************************************************************/
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*ht;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; u10[m][0] = f_boundary(sn0, tn10); u05[m][0] = f_boundary(sn0, tn05);
        sn1.j = m; sn1.y = m*hy; u10[m][N] = f_boundary(sn1, tn10); u05[m][N] = f_boundary(sn1, tn05);
    }

    sn0.j = 0; sn0.y = 0.0;
    sn1.j = M; sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; u10[0][n] = f_boundary(sn0, tn10); u05[0][n] = f_boundary(sn0, tn05);
        sn1.i = n; sn1.x = n*hx; u10[M][n] = f_boundary(sn1, tn10); u05[M][n] = f_boundary(sn1, tn05);
    }

    /************************************************************************/

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= lambda*f_initial2(sn);

            u05[m][n] = u00[m][n] + (0.5*ht) * f_initial2(sn) + (0.125*ht*ht) * sum;
            u10[m][n] = u00[m][n] + (1.0*ht) * f_initial2(sn) + (0.500*ht*ht) * sum;

            //            double sum1 = 0.0;
            //            for (unsigned int cdi=0; cdi<cntDeltaNodes.size(); cdi++)
            //            {
            //                const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cdi);
            //                if (cdn.i == n && cdn.j == m)
            //                {
            //                    for (unsigned int j=0; j<mEquParameter.No; j++)
            //                        sum1 -= mOptParameter.z[cdn.id][j] * cdn.w;
            //                }
            //            }

            //            u05[m][n] += (0.125*ht*ht) * sum1;
            //            u10[m][n] += (0.500*ht*ht) * sum1;

            u05[m][n] += (0.125*ht*ht)*distributeTimeDelta(0.5*ht, ht, 1, qPointNodes, sn, qExtSpacePoints);
            u10[m][n] += (0.500*ht*ht)*distributeTimeDelta(1.0*ht, ht, 2, qPointNodes, sn, qExtSpacePoints);
        }
    }

    /************************************************************************/

    if (use == true) f_add2Info(u00, u_info, obsPointNodes, 0, hx, hy, msnExtSpacePoints);
    f_layerInfo(u00, 0);

    if (use == true) f_add2Info(u05, u_info, obsPointNodes, 1, hx, hy, msnExtSpacePoints);
    f_layerInfo(u05, 1);

    if (use == true) f_add2Info(u10, u_info, obsPointNodes, 2, hx, hy, msnExtSpacePoints);
    f_layerInfo(u10, 2);
}

double Problem2HNDirichlet::f_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return 0.0;
}

double Problem2HNDirichlet::f_initial2(const SpaceNodePDE &) const
{
    return 0.0;
}

double Problem2HNDirichlet::f_boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType) const
{
    return 0.0;
}

void Problem2HNDirichlet::f_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2, uint_vector &cols0, uint_vector &cols1, uint_vector &cols2,
                                         espn_vector &obsPointNodes, espn_vector &cntDeltaNodes, unsigned int N, unsigned int M) const
{
#ifdef OLD_VERSION
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
#endif
#ifdef NEW_VERSION
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
#endif
}

void Problem2HNDirichlet::f_borderLayer(DoubleMatrix &u, DoubleMatrix &um5, unsigned int ln) const
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = dimX.sizeN();
    const unsigned int M = dimY.sizeN();

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    TimeNodePDE tn00; tn00.i = ln; tn00.t = ln*ht;
    TimeNodePDE tnm5; tnm5.i = ln; tnm5.t = ln*ht-0.5*ht;

    SpaceNodePDE sn0;
    SpaceNodePDE sn1;

    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = hx*N;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; um5[m][0] = f_boundary(sn0, tnm5); u[m][0] = f_boundary(sn0, tn00);
        sn1.j = m; sn1.y = m*hy; um5[m][N] = f_boundary(sn1, tnm5); u[m][N] = f_boundary(sn1, tn00);
    }

    sn0.j = 0; sn0.y = 0.0;
    sn1.j = M; sn1.y = hy*M;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; um5[0][n] = f_boundary(sn0, tnm5); u[0][n] = f_boundary(sn0, tn00);
        sn1.i = n; sn1.x = n*hx; um5[M][n] = f_boundary(sn1, tnm5); u[M][n] = f_boundary(sn1, tn00);
    }
}

void Problem2HNDirichlet::f_prepareInfo(unsigned int No, const std::vector<SpacePoint> &points, spif_vector &u_info,
                                        unsigned int LLD, const Dimension &dimX, const Dimension &dimY) const
{
    u_info.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        SpacePointInfo &inf = u_info[j];
        const SpacePoint &sp = points[j];
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(2*LLD+1);
    }
}

void Problem2HNDirichlet::f_add2Info(const DoubleMatrix &u, spif_vector &u_info, const espn_vector &obsPointNodes, unsigned int ln,
                                     double hx, double hy, std::vector<ExtendedSpacePoint> &extMsmnts, int method) const
{
    if (method == 1)
    {
#ifdef OLD_VERSION
        for (unsigned int j=0; j<u_info.size(); j++) u_info[j]._vl[ln] = u_info[j]._dx[ln] = u_info[j]._dy[ln] = 0.0;

        for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
        {
            const ExtendedSpacePointNode &opn = obsPointNodes[opj];
            const SpacePointInfo &ui = u_info[opn.id];

            ui._vl[ln] += u[opn.j][opn.i] * (opn.w * (hx*hy));
            if (opn.isCenter)
            {
                ui._dx[ln] = (u[ui.j][ui.i+1] - u[ui.j][ui.i-1])/(2.0*hx);
                ui._dy[ln] = (u[ui.j+1][ui.i] - u[ui.j-1][ui.i])/(2.0*hy);
            }
        }
#endif
    }

    if (method == 2)
    {
#ifdef OLD_VERSION
        for (unsigned int j=0; j<u_info.size(); j++) u_info[j]._vl[ln] = u_info[j]._dx[ln] = u_info[j]._dy[ln] = 0.0;

        for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
        {
            const ExtendedSpacePointNode &opn = obsPointNodes[opj];
            const SpacePointInfo &ui = u_info[opn.id];

            ui._vl[ln] += u[opn.j][opn.i] * (opn.w * (hx*hy));
            if (opn.isCenter)
            {
                ui._dx[ln] = (u[ui.j][ui.i+1] - u[ui.j][ui.i-1])/(2.0*hx);
                ui._dy[ln] = (u[ui.j+1][ui.i] - u[ui.j-1][ui.i])/(2.0*hy);
            }
        }
#endif
    }

    if (method == 4)
    {
#ifdef OLD_VERSION
        for (unsigned int j=0; j<u_info.size(); j++) u_info[j]._vl[ln] = u_info[j]._dx[ln] = u_info[j]._dy[ln] = 0.0;

        for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
        {
            const ExtendedSpacePointNode &opn = obsPointNodes[opj];
            const SpacePointInfo &ui = u_info[opn.id];

            ui._vl[ln] += u[opn.j][opn.i] * (opn.w * (hx*hy));
            if (opn.isCenter)
            {
                ui._dx[ln] = (u[ui.j][ui.i+1] - u[ui.j][ui.i-1])/(2.0*hx);
                ui._dy[ln] = (u[ui.j+1][ui.i] - u[ui.j-1][ui.i])/(2.0*hy);

                ui._dx[ln] += ((ui.x-ui.i*hx)/(hx*hx))*(u[ui.j][ui.i+1] - 2.0*u[ui.j][ui.i] + u[ui.j][ui.i-1]);
                ui._dy[ln] += ((ui.y-ui.j*hy)/(hy*hy))*(u[ui.j+1][ui.i] - 2.0*u[ui.j][ui.i] + u[ui.j-1][ui.i]);
            }
        }
#endif
#ifdef NEW_VERSION
        unsigned int No = static_cast<unsigned int>(extMsmnts.size());
        for (unsigned int j=0; j<No; j++)
        {
            ExtendedSpacePoint &xsp = extMsmnts.at(j);
            SpacePointInfo &ui = u_info[j];
            unsigned int nodes_size = xsp.nodes.size();
            for (unsigned int i=0; i<nodes_size; i++)
            {
                const ExtendedSpacePointNode1 &node = xsp.nodes.at(i);
                ui.vl[ln] += u[node.ny][node.nx] * (node.w * (hx*hy));
                if (node.isCenter)
                {
                    ui.dx[ln] = (u[xsp.ry][xsp.rx+1] - u[xsp.ry][xsp.rx-1])/(2.0*hx);
                    ui.dy[ln] = (u[xsp.ry+1][xsp.rx] - u[xsp.ry-1][xsp.rx])/(2.0*hy);

                    ui.dx[ln] += ((xsp.x-xsp.rx*hx)/(hx*hx))*(u[xsp.ry][xsp.rx+1] - 2.0*u[xsp.ry][xsp.rx] + u[xsp.ry][xsp.rx-1]);
                    ui.dy[ln] += ((xsp.y-xsp.ry*hy)/(hy*hy))*(u[xsp.ry+1][xsp.rx] - 2.0*u[xsp.ry][xsp.rx] + u[xsp.ry-1][xsp.rx]);
                }
            }
        }
#endif
    }
}

void Problem2HNDirichlet::solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const int unsigned N = dimX.sizeN();
    const int unsigned M = dimY.sizeN();
    const int unsigned L = time.sizeN();
    const int unsigned LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a = mEquParameter.a;
    const double lambda = mEquParameter.lambda;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;

    const double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    const double p_aa_htht__hxhx___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hx*hx) + (lambda*ht);
    const double p_aa_htht__hxhx = +(a*a*ht*ht)/(hx*hx);

    const double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    const double p_aa_htht__hyhy___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hy*hy) + (lambda*ht);
    const double p_aa_htht__hyhy = +(a*a*ht*ht)/(hy*hy);

    const double htht = ht*ht;
    const double lambda_ht = lambda*ht;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    espn_vector obsDeltaNodes, cntPointNodes;
    std::vector<ExtendedSpacePoint> msnExtSpacePoints, cntExtSpacePoints, qExtSpacePoints;
#if defined(OLD_VERSION) || defined(NEW_VERSION)
    for (unsigned int j=0; j<No; j++) distributeDelta0(mOptParameter.xi[j], j, obsDeltaNodes, dimX, dimY, 4);
    for (unsigned int i=0; i<Nc; i++) distributeDelta0(mOptParameter.eta[i], i, cntPointNodes, dimX, dimY, 4);
#endif
#ifdef NEW_VERSION
    newDistributeDeltaGauseCntrl(mOptParameter.eta, cntExtSpacePoints, dimX, dimY);
    newDistributeDeltaGauseMsmnt(mOptParameter.xi,  msnExtSpacePoints, dimX, dimY);
#endif

    //----------------------------------------------------------------------------------------------//
    uint_vector rows0, rows1, rows2, cols0, cols1, cols2;
    b_findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, cntPointNodes, obsDeltaNodes, N, M);

    //-------------------------------------------- info --------------------------------------------//
    if (use == true) b_prepareInfo(Nc, mOptParameter.eta, p_info, LLD, dimX, dimY);
    //-------------------------------------------- info --------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    b_initialLayers(p00, p10, p_info, use, u_info, cntPointNodes, obsDeltaNodes, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, lambda, cntExtSpacePoints);
    //------------------------------------- initial conditions -------------------------------------//

    double *ax = static_cast<double*>( malloc(sizeof(double)*(N-1)) ); for (unsigned int n=1; n<=N-1; n++) ax[n-1] = m_aa_htht__hxhx;
    double *bx = static_cast<double*>( malloc(sizeof(double)*(N-1)) ); for (unsigned int n=1; n<=N-1; n++) bx[n-1] = p_aa_htht__hxhx___lambda_ht;
    double *cx = static_cast<double*>( malloc(sizeof(double)*(N-1)) ); for (unsigned int n=1; n<=N-1; n++) cx[n-1] = m_aa_htht__hxhx;
    double *dx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    double *rx = static_cast<double*>( malloc(sizeof(double)*(N-1)) );
    ax[0] = cx[N-2] = 0.0;

    double *ay = static_cast<double*>( malloc(sizeof(double)*(M-1)) ); for (unsigned int m=1; m<=M-1; m++) ay[m-1] = m_aa_htht__hyhy;
    double *by = static_cast<double*>( malloc(sizeof(double)*(M-1)) ); for (unsigned int m=1; m<=M-1; m++) by[m-1] = p_aa_htht__hyhy___lambda_ht;
    double *cy = static_cast<double*>( malloc(sizeof(double)*(M-1)) ); for (unsigned int m=1; m<=M-1; m++) cy[m-1] = m_aa_htht__hyhy;
    double *dy = static_cast<double*>( malloc(sizeof(double)*(M-1)) );
    double *ry = static_cast<double*>( malloc(sizeof(double)*(M-1)) );
    ay[0] = cy[M-2] = 0.0;

    unsigned int row1_size = rows1.size()*(N-1);
    double* a1, *b1, *c1, *d1, *x1, **w1;
    if (rows1.size() != 0 && rows2.size() != 0)
    {
        a1 = static_cast<double*>( malloc(sizeof(double)*row1_size) );
        b1 = static_cast<double*>( malloc(sizeof(double)*row1_size) );
        c1 = static_cast<double*>( malloc(sizeof(double)*row1_size) );
        d1 = static_cast<double*>( malloc(sizeof(double)*row1_size) );
        x1 = static_cast<double*>( malloc(sizeof(double)*row1_size) );
        w1 = static_cast<double**>( malloc(sizeof(double*)*row1_size) );
        for (unsigned int row=0; row < row1_size; row++) w1[row] = static_cast<double*>( malloc(sizeof(double)*row1_size) );
    }

    unsigned int cols1_size = cols1.size()*(M-1);
    double *a2, *b2, *c2, *d2, *x2, **w2;
    if (cols1.size() != 0 && cols2.size() != 0)
    {
        a2 = static_cast<double*>( malloc(sizeof(double)*cols1_size) );
        b2 = static_cast<double*>( malloc(sizeof(double)*cols1_size) );
        c2 = static_cast<double*>( malloc(sizeof(double)*cols1_size) );
        d2 = static_cast<double*>( malloc(sizeof(double)*cols1_size) );
        x2 = static_cast<double*>( malloc(sizeof(double)*cols1_size) );
        w2 = static_cast<double**>( malloc(sizeof(double*)*cols1_size) );
        for (unsigned int col=0; col < cols1_size; col++) w2[col] = static_cast<double*>( malloc(sizeof(double)*cols1_size) );
    }

    SpaceNodePDE sn;

    for (unsigned int l=LLD-2; l != (unsigned)0-1; l--)
    {
        TimeNodePDE tn20; tn20.i = l; tn20.t = l*ht;
        TimeNodePDE tn15; tn15.i = l; tn15.t = l*ht+0.5*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy; p15[m][0] = b_boundary(sn0, tn15); p20[m][0] = b_boundary(sn0, tn20);
            sn1.j = m; sn1.y = m*hy; p15[m][N] = b_boundary(sn1, tn15); p20[m][N] = b_boundary(sn1, tn20);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx; p15[0][n] = b_boundary(sn0, tn15); p20[0][n] = b_boundary(sn0, tn20);
            sn1.i = n; sn1.x = n*hx; p15[M][n] = b_boundary(sn1, tn15); p20[M][n] = b_boundary(sn1, tn20);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = rows0.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    dx[n-1] += lambda_ht*(p10[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dx[n-1] += 2.0*p10[m][n] + (p10[m][n]-p00[m][n]);

                    //------------------------------------- Adding functional part --------------------------------//
                    if (L <= l && l <= LLD) dx[n-1] += -2.0*(u.at(2*(l-L)+1)[m][n]) * htht;
                    //------------------------------------- Adding functional part --------------------------------//
                }

                dx[0]   -= m_aa_htht__hxhx * p15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * p15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0 && rows2.size() == 0)
        {
            //throw std::exception();
            double *_w15 = new double[No];

            double* _p15 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p15[i] = 0.0;
            for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
            {
                const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];
                _p15[cpn.id] += p15[cpn.j][cpn.i] * (cpn.w * (hx*hy));
            }

            for (unsigned int j=0; j<No; j++)
            {
                _w15[j] = 0.0;
                for (unsigned int i=0; i<Nc; i++)
                {
                    _w15[j] += mOptParameter.k[i][j] * ( _p15[i] + 2.0*r*gpi(i, 2*l+1, u_info, mOptParameter)*sgn(g0i(i, 2*l+1, u_info, mOptParameter)));
                }
            }
            delete [] _p15;

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    dx[n-1] += lambda_ht*(p10[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dx[n-1] += 2.0*p10[m][n] + (p10[m][n]-p00[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int odj=0; odj<obsDeltaNodes.size(); odj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(odj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            dx[n-1] += htht * _w15[odn.id] * odn.w;
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//

                    //------------------------------------- Adding functional part --------------------------------//

                    if (L <= l && l <= LLD) dx[n-1] += -2.0*(u.at(2*(l-L))[m][n]) * htht;

                    //------------------------------------- Adding functional part --------------------------------//
                }

                dx[0]   -= m_aa_htht__hxhx * p15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * p15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
            }

            delete [] _w15;
        }

        if (rows1.size() != 0 && rows2.size() != 0)
        {
            //throw std::exception();

            for (unsigned int i=0; i < row1_size; i++) for (unsigned int j=0; j < row1_size; j++) w1[i][j] = 0.0;

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    const unsigned int index = offset+(n-1);
                    d1[index] = 0.0;
                    if (m>0 && m<M)  d1[index] = p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    else if (m == 0) d1[index] = p_aa_htht__hyhy*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    else if (m == M) d1[index] = p_aa_htht__hyhy*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    d1[index] += lambda_ht*(p10[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    d1[index] += 2.0*p10[m][n] + (p10[m][n]-p00[m][n]);

                    a1[index] = m_aa_htht__hxhx;
                    b1[index] = p_aa_htht__hxhx___lambda_ht;
                    c1[index] = m_aa_htht__hxhx;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int odj=0; odj<obsDeltaNodes.size(); odj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(odj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cdi=0; cdi<cntPointNodes.size(); cdi++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cdi];

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (cpn.j == rows1[rs])
                                    {
                                        found = true;
                                        w1[index][rs*(N-1)+(cpn.i-1)] -= htht * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d1[index] += htht * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1[index] += 2.0 * r * htht *  mOptParameter.k[i][odn.id] * gpi(i, 2*l+1, u_info, mOptParameter)*sgn(g0i(i, 2*l+1, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//

                    //------------------------------------- Adding functional part --------------------------------//

                    if (L <= l && l <= LLD) d1[offset+(n-1)] += -2.0*(u.at(2*(l-L)+1)[m][n]) * htht;

                    //------------------------------------- Adding functional part --------------------------------//
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx * p15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx * p15[m][N];

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w1, x1, row1_size);

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
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        if (cols0.size() != 0)
        {
            for (unsigned int col=0; col<cols0.size(); col++)
            {
                unsigned int n = cols0.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    dy[m-1] += lambda_ht*(p15[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dy[m-1] += 2.0*p15[m][n] + (p10[m][n]-p00[m][n]);

                    //------------------------------------- Adding functional part --------------------------------//
                    if (L <= l && l <= LLD) dy[m-1] += -2.0*(u[2*(l-L)][m][n]) * htht;
                    //------------------------------------- Adding functional part --------------------------------//
                }

                dy[0]   -= m_aa_htht__hyhy * p20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * p20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0 && cols2.size() == 0)
        {
            //throw std::exception();
            double *_w20 = new double[No];

            double* _p20 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p20[i] = 0.0;
            for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
            {
                const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];
                _p20[cpn.id] += p20[cpn.j][cpn.i] * (cpn.w * (hx*hy));
            }

            for (unsigned int j=0; j<No; j++)
            {
                _w20[j] = 0.0;
                for (unsigned int i=0; i<Nc; i++)
                {
                    _w20[j] += mOptParameter.k[i][j] * (_p20[i] + 2.0*r*gpi(i, 2*l+0, u_info, mOptParameter)*sgn(g0i(i, 2*l+0, u_info, mOptParameter)));
                }
            }
            delete [] _p20;

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    dy[m-1] += lambda_ht*(p15[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dy[m-1] += 2.0*p15[m][n] + (p10[m][n]-p00[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int odj=0; odj<obsDeltaNodes.size(); odj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(odj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            dy[m-1] += htht * _w20[odn.id] * odn.w;
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//

                    //------------------------------------- Adding functional part --------------------------------//

                    if (L <= l && l <= LLD) dy[m-1] += -2.0*(u[2*(l-L)][m][n]) * htht;

                    //------------------------------------- Adding functional part --------------------------------//
                }

                dy[0]   -= m_aa_htht__hyhy * p20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * p20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
            }

            delete [] _w20;
        }

        if (cols1.size() != 0 && cols2.size() != 0)
        {
            //throw std::exception();

            for (unsigned int i=0; i < cols1_size; i++) for (unsigned int j=0; j < cols1_size; j++) w2[i][j] = 0.0;

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    const unsigned int index = offset+(m-1);
                    d2[index] = 0.0;
                    if (n>0 && n<N) d2[index] = p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    else if (n==0)  d2[index] = p_aa_htht__hxhx*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    else if (n==N)  d2[index] = p_aa_htht__hxhx*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    d2[index] += lambda_ht*(p15[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    d2[index] += 2.0*p15[m][n] + (p10[m][n]-p00[m][n]);

                    a2[index] = m_aa_htht__hyhy;
                    b2[index] = p_aa_htht__hyhy___lambda_ht;
                    c2[index] = m_aa_htht__hyhy;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int odj=0; odj<obsDeltaNodes.size(); odj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(odj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (cpn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[index][cs*(M-1)+(cpn.j-1)] -= htht * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d2[index] += htht * mOptParameter.k[cpn.id][odn.id] * p20[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d2[index] += 2.0 * r * htht *  mOptParameter.k[i][odn.id] * gpi(i, 2*l, u_info, mOptParameter)*sgn(g0i(i, 2*l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//

                    //------------------------------------- Adding functional part --------------------------------//

                    if (L <= l && l <= LLD) d2[index] += -2.0*(u[2*(l-L)][m][n]) * htht;

                    //------------------------------------- Adding functional part --------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy * p20[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy * p20[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2, x2, cols1_size);

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    p20[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }
        }

        /**************************************************** y direction apprx ***************************************************/

        if (use == true) b_add2Info(p15, p_info, cntPointNodes, 2*l+1, hx, hy); b_layerInfo(p15, 2*l+1);
        if (use == true) b_add2Info(p20, p_info, cntPointNodes, 2*l+0, hx, hy); b_layerInfo(p20, 2*l+0);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p20[m][n];
            }
        }
    }

    if (rows1.size() != 0 && rows2.size() != 0)
    {
        for (unsigned int row=0; row < row1_size; row++) free(w1[row]); free(w1);
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

    obsDeltaNodes.clear();
    cntPointNodes.clear();

    p00.clear();
    p10.clear();
    p15.clear();
    p20.clear();
}

void Problem2HNDirichlet::b_initialLayers(DoubleMatrix &p00, DoubleMatrix &p10, spif_vector &p_info, bool use, const spif_vector &u_info,
                                          espn_vector &cntPointNodes, espn_vector &obsDeltaNodes, unsigned int N, unsigned int M, double
                                          hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda,
                                          std::vector<ExtendedSpacePoint> &cntExtSpacePoints) const
{
    DoubleMatrix p05 = p00;
    unsigned int L = mtimeDimension.sizeN();
    unsigned int LLD = L+LD;

    /************************************************************************/
    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p00[m][n] = b_initial1(sn);
        }
    }

    /************************************************************************/
    TimeNodePDE tn05; tn05.i = LLD-1; tn05.t = LLD*ht + 0.5*ht;
    TimeNodePDE tn10; tn10.i = LLD-1; tn10.t = LLD*ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; p10[m][0] = b_boundary(sn0, tn10); p05[m][0] = b_boundary(sn0, tn05);
        sn1.j = m; sn1.y = m*hy; p10[m][N] = b_boundary(sn1, tn10); p05[m][N] = b_boundary(sn1, tn05);
    }

    sn0.j = 0;  sn0.y = 0.0;
    sn1.j = M; sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; p10[0][n] = b_boundary(sn0, tn10); p05[0][n] = b_boundary(sn0, tn05);
        sn1.i = n; sn1.x = n*hx; p10[M][n] = b_boundary(sn1, tn10); p05[M][n] = b_boundary(sn1, tn05);
    }

    /************************************************************************/

    unsigned int No = mEquParameter.No;
    unsigned int Nc = mEquParameter.Nc;
    double *_w = new double[No];

    double* _p = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p[i] = 0.0;

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        const ExtendedSpacePoint &esp = cntExtSpacePoints.at(i);
//        const std::vector<ExtendedSpacePointNode1> &nodes = esp.nodes;
//        for (unsigned int ni=0; ni<nodes.size(); ni++)
//        {
//            const ExtendedSpacePointNode1 &node = nodes.at(ni);
//            _p[cpn.id] += p00[node.ny][node.nx] * (node.w * (hx*hy));
//        }
//    }
    for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
    {
        const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];
        _p[cpn.id] += p00[cpn.j][cpn.i] * (cpn.w * (hx*hy));
    }

    for (unsigned int j=0; j<No; j++)
    {
        _w[j] = 0.0;
        for (unsigned int i=0; i<Nc; i++)
        {
            _w[j] += mOptParameter.k[i][j] * (_p[i] + 2.0*r*gpi(i, 2*LLD, u_info, mOptParameter)*sgn(g0i(i,2*LLD, u_info, mOptParameter)));
        }
    }
    delete [] _p;

    /************************************************************************/

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);
            sum += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);
            sum += lambda*b_initial2(sn);

            p05[m][n] = p00[m][n] - (ht*0.5) * b_initial2(sn) + 0.125*ht*ht*sum;
            p10[m][n] = p00[m][n] - (ht*1.0) * b_initial2(sn) + 0.500*ht*ht*sum;

            double sum1 = 0.0;
            for (unsigned int odj=0; odj<obsDeltaNodes.size(); odj++)
            {
                const ExtendedSpacePointNode &odn = obsDeltaNodes.at(odj);
                if (odn.i == n && odn.j == m)
                {
                    sum1 += _w[odn.id] * odn.w;
                }
            }

            p05[m][n] += (0.125*ht*ht)*sum1;
            p10[m][n] += (0.500*ht*ht)*sum1;
        }
    }

    /************************************************************************/

    if (use == true) b_add2Info(p00, p_info, cntPointNodes, 2*LLD, hx, hy);
    b_layerInfo(p00, 2*LLD);

    if (use == true) b_add2Info(p05, p_info, cntPointNodes, 2*(LLD-1)+1, hx, hy);
    b_layerInfo(p05, 2*(LLD-1)+1);

    if (use == true) b_add2Info(p10, p_info, cntPointNodes, 2*(LLD-1)+0, hx, hy);
    b_layerInfo(p10, 2*(LLD-1)+0);

    delete [] _w;
}

double Problem2HNDirichlet::b_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return 0.0;
}

double Problem2HNDirichlet::b_initial2(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return 0.0;
}

double Problem2HNDirichlet::b_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const
{
    return 0.0;
}

double Problem2HNDirichlet::b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const
{
    return -2.0*(u[m][n]);
}

void Problem2HNDirichlet::b_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2, uint_vector &cols0, uint_vector &cols1, uint_vector &cols2,
                                         espn_vector &cntPointNodes, espn_vector &obsDeltaNodes, unsigned int N, unsigned int M) const
{
#ifdef OLD_VERSION
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
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
    }

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
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
    }
#endif
#ifdef NEW_VERSION
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
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
    }

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
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
    }
#endif
}

void Problem2HNDirichlet::b_prepareInfo(unsigned int Nc, const std::vector<SpacePoint> &points, spif_vector &p_info,
                                        unsigned int LLD, const Dimension &dimX, const Dimension &dimY) const
{
    p_info.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        SpacePointInfo &inf = p_info[i];
        const SpacePoint &sp = points[i];
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(2*LLD+1);
    }
}

void Problem2HNDirichlet::b_add2Info(const DoubleMatrix &p, spif_vector &p_info, const espn_vector &cntPointNodes, unsigned int ln, double hx, double hy, int method) const
{
    if (method == 1)
    {
        for (unsigned int j=0; j<p_info.size(); j++) p_info[j].vl[ln] = p_info[j].dx[ln] = p_info[j].dy[ln] = 0.0;

        for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
        {
            const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];
            SpacePointInfo &pi = p_info[cpn.id];

            pi.vl[ln] += p[cpn.j][cpn.i] * (cpn.w * (hx*hy));
            if (cpn.isCenter)
            {
                pi.dx[ln] = (p[cpn.j][cpn.i+1] - p[cpn.j][cpn.i-1])/(2.0*hx);
                pi.dy[ln] = (p[cpn.j+1][cpn.i] - p[cpn.j-1][cpn.i])/(2.0*hy);
            }
        }
    }

    if (method == 2)
    {
        for (unsigned int j=0; j<p_info.size(); j++) p_info[j].vl[ln] = p_info[j].dx[ln] = p_info[j].dy[ln] = 0.0;

        for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
        {
            const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];
            SpacePointInfo &pi = p_info[cpn.id];

            pi.vl[ln] += p[cpn.j][cpn.i] * (cpn.w * (hx*hy));
            if (cpn.isCenter)
            {
                pi.dx[ln] = (p[cpn.j][cpn.i+1] - p[cpn.j][cpn.i-1])/(2.0*hx);
                pi.dy[ln] = (p[cpn.j+1][cpn.i] - p[cpn.j-1][cpn.i])/(2.0*hy);
            }
        }
    }

    if (method == 4)
    {
        for (unsigned int j=0; j<p_info.size(); j++) p_info[j].vl[ln] = p_info[j].dx[ln] = p_info[j].dy[ln] = 0.0;

        for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
        {
            const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];
            SpacePointInfo &pi = p_info[cpn.id];

            pi.vl[ln] += p[cpn.j][cpn.i] * (cpn.w * (hx*hy));
            if (cpn.isCenter)
            {
                pi.dx[ln] = (p[cpn.j][cpn.i+1] - p[cpn.j][cpn.i-1])/(2.0*hx);
                pi.dy[ln] = (p[cpn.j+1][cpn.i] - p[cpn.j-1][cpn.i])/(2.0*hy);

                pi.dx[ln] += ((pi.x-cpn.i*hx)/(hx*hx))*(p[cpn.j][cpn.i+1] - 2.0*p[cpn.j][cpn.i] + p[cpn.j][cpn.i-1]);
                pi.dy[ln] += ((pi.y-cpn.j*hy)/(hy*hy))*(p[cpn.j+1][cpn.i] - 2.0*p[cpn.j][cpn.i] + p[cpn.j-1][cpn.i]);
            }
        }
    }
}

void Problem2HNDirichlet::b_layerInfo(const DoubleMatrix &p UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const
{
    //double min = p.min();
    //double max = p.max();

    //QPixmap pic;
    //visualizeMatrixHeat(p, min, max, pic);
    //pic.save("images/b/pic"+QString("%1").arg(ln)+".png", "PNG");
    //printf("Layer: %d min: %f max: %f min: %f max: %f norm: %f\n", ln, min, max, min, max, fabs(max-min));

    if (ln == 0)
    {
        IPrinter::printSeperatorLine();
        IPrinter::printMatrix(p);
        IPrinter::printSeperatorLine();
        printf("Backward: %f\n", sqrt(integralU(p)));
    }
}

// backward -----------------------------------

void Problem2HNDirichlet::distributeDelta0(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, int method, unsigned int k) const
{
    if (method == 1) distributeDeltaP(pt, id, nodes, dimX, dimY);

    if (method == 2) distributeDeltaR(pt, id, nodes, dimX, dimY);

    if (method == 4) distributeDeltaG(pt, id, nodes, dimX, dimY, k);
}

void Problem2HNDirichlet::distributeDeltaP(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY) const
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.sizeN() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.sizeN() );

    unsigned int rx = static_cast<unsigned int> ( round( pt.x * Nx ) );
    unsigned int ry = static_cast<unsigned int> ( round( pt.y * Ny ) );

    ExtendedSpacePointNode node; node.id = id; node.pt = pt; node.i = rx; node.j = ry; node.w = 1.0/(hx*hy); nodes.push_back(node);
}

void Problem2HNDirichlet::distributeDeltaR(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY) const
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.sizeN() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.sizeN() );

    unsigned int rx = static_cast<unsigned int> ( floor( pt.x * Nx ) );
    unsigned int ry = static_cast<unsigned int> ( floor( pt.y * Ny ) );

    double h1x = fabs(pt.x - rx*hx);
    double h1y = fabs(pt.y - ry*hy);
    double h2x = hx-fabs(pt.x - rx*hx);
    double h2y = hx-fabs(pt.y - ry*hy);

    ExtendedSpacePointNode node00; node00.id = id; node00.pt = pt; node00.i = rx+0; node00.j = ry+0; node00.w = ((h2x/hx)*(h2y/hy))/(hx*hy); nodes.push_back(node00);
    ExtendedSpacePointNode node01; node01.id = id; node01.pt = pt; node01.i = rx+0; node01.j = ry+1; node01.w = ((h2x/hx)*(h1y/hy))/(hx*hy); nodes.push_back(node01);
    ExtendedSpacePointNode node11; node11.id = id; node11.pt = pt; node11.i = rx+1; node11.j = ry+1; node11.w = ((h1x/hx)*(h1y/hy))/(hx*hy); nodes.push_back(node11);
    ExtendedSpacePointNode node10; node10.id = id; node10.pt = pt; node10.i = rx+1; node10.j = ry+0; node10.w = ((h1x/hx)*(h2y/hy))/(hx*hy); nodes.push_back(node10);
}

void Problem2HNDirichlet::distributeDeltaG(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k) const
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.sizeN() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.sizeN() );

    unsigned int rx = static_cast<unsigned int> ( round(pt.x*Nx) );
    unsigned int ry = static_cast<unsigned int> ( round(pt.y*Ny) );

    double sigmaX = 1.0*hx;
    double sigmaY = 1.0*hy;

    double sumX = 0.0;
    for (unsigned int n=rx-k; n<=rx+k; n++) sumX += exp(-((n*hx-pt.x)*(n*hx-pt.x))/(2.0*sigmaX*sigmaX));
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=ry-k; m<=ry+k; m++) sumY += exp(-((m*hy-pt.y)*(m*hy-pt.y))/(2.0*sigmaY*sigmaY));
    sumY *= hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/((2.0*M_PI)*sigma);

    for (unsigned int m=ry-k; m<=ry+k; m++)
    {
        for (unsigned int n=rx-k; n<=rx+k; n++)
        {
            ExtendedSpacePointNode node;
            node.i = n; node.x = n*hx;
            node.j = m; node.y = m*hy;
            node.pt = pt; node.id = id;
            node.w = factor*exp(-0.5*(((node.x-pt.x)*(node.x-pt.x))/(sigmaX*sigmaX)+((node.y-pt.y)*(node.y-pt.y))/(sigmaY*sigmaY)));
            node.isCenter = ( m==ry && n==rx );
            nodes.push_back(node);
        }
    }
}

void Problem2HNDirichlet::distributeDeltaGaussPulse(const SpacePoint &pt, unsigned id, std::vector<ExtendedSpacePointNode> &qPointNodes,
                                                    const Dimension &dimX, const Dimension &dimY, unsigned int k) const
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.sizeN() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.sizeN() );

    unsigned int rx = static_cast<unsigned int> ( round(pt.x*Nx) );
    unsigned int ry = static_cast<unsigned int> ( round(pt.y*Ny) );

    k *= 5;
    double sigmaX = hx*5.0;
    double sigmaY = hy*5.0;

    double sumX = 0.0;
    for (unsigned int n=rx-k; n<=rx+k; n++) sumX += exp(-((n*hx-pt.x)*(n*hx-pt.x))/(2.0*sigmaX*sigmaX));
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=ry-k; m<=ry+k; m++) sumY += exp(-((m*hy-pt.y)*(m*hy-pt.y))/(2.0*sigmaY*sigmaY));
    sumY *= hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/((2.0*M_PI)*sigma);

    for (unsigned int m=ry-k; m<=ry+k; m++)
    {
        for (unsigned int n=rx-k; n<=rx+k; n++)
        {
            ExtendedSpacePointNode node;
            node.i = n; node.x = n*hx;
            node.j = m; node.y = m*hy;
            node.pt = pt; node.id = id;
            node.w = factor*exp(-0.5*(((node.x-pt.x)*(node.x-pt.x))/(sigmaX*sigmaX)+((node.y-pt.y)*(node.y-pt.y))/(sigmaY*sigmaY)));
            node.isCenter = ( m==ry && n==rx );
            qPointNodes.push_back(node);
        }
    }
}

void Problem2HNDirichlet::newDistributeDeltaGaussPulse(const std::vector<SpacePoint> &thetas, std::vector<ExtendedSpacePoint> &extThetas, const Dimension &dimX, const Dimension &dimY) const
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.sizeN() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.sizeN() );
    unsigned int Ns = static_cast<unsigned int> ( thetas.size() );

    extThetas.clear();
    extThetas.resize(Ns);

    int k = 4;
    double sigmaX = hx*5.0;
    double sigmaY = hy*5.0;

    for (unsigned int s=0; s<Ns; s++)
    {
        const SpacePoint &theta = thetas.at(s);
        ExtendedSpacePoint &extTheta = extThetas.at(s);

        extTheta.x = theta.x;
        extTheta.y = theta.y;
        extTheta.rx = static_cast<unsigned int> ( round(extTheta.x*Nx) );
        extTheta.ry = static_cast<unsigned int> ( round(extTheta.y*Ny) );
        extTheta.k = k*5;
        extTheta.minX = extTheta.rx - extTheta.k;
        extTheta.maxX = extTheta.rx + extTheta.k;
        extTheta.minY = extTheta.ry - extTheta.k;
        extTheta.maxY = extTheta.ry + extTheta.k;

        double sumX = 0.0;
        for (int n=extTheta.minX; n<=extTheta.maxX; n++) sumX += exp(-((n*hx-theta.x)*(n*hx-theta.x))/(2.0*sigmaX*sigmaX));
        sumX *= hx;

        double sumY = 0.0;
        for (int m=extTheta.minY; m<=extTheta.maxY; m++) sumY += exp(-((m*hy-theta.y)*(m*hy-theta.y))/(2.0*sigmaY*sigmaY));
        sumY *= hy;

        double sigma = (sumX*sumY) / (2.0*M_PI);
        double factor = 1.0/((2.0*M_PI)*sigma);

        for (int m=extTheta.minY; m<=extTheta.maxY; m++)
        {
            for (int n=extTheta.minX; n<=extTheta.maxX; n++)
            {
                ExtendedSpacePointNode1 node;
                node.nx = n; node.x = n*hx;
                node.ny = m; node.y = m*hy;
                node.w = factor*exp(-0.5*(((node.x-theta.x)*(node.x-theta.x))/(sigmaX*sigmaX)+((node.y-theta.y)*(node.y-theta.y))/(sigmaY*sigmaY)));
                node.isCenter = ( m==extTheta.ry && n==extTheta.rx );
                extTheta.nodes.push_back(node);
            }
        }
    }
}

void Problem2HNDirichlet::newDistributeDeltaGauseCntrl(const std::vector<SpacePoint> &cntrls, std::vector<ExtendedSpacePoint> &extCntrls, const Dimension &dimX, const Dimension &dimY) const
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
        ExtendedSpacePoint &extCntrl = extCntrls.at(c);

        extCntrl.x = cntrl.x;
        extCntrl.y = cntrl.y;
        extCntrl.rx = static_cast<unsigned int> ( round(extCntrl.x*Nx) );
        extCntrl.ry = static_cast<unsigned int> ( round(extCntrl.y*Ny) );
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
                ExtendedSpacePointNode1 node;
                node.nx = n; node.x = n*hx;
                node.ny = m; node.y = m*hy;
                node.w = factor*exp(-0.5*(((node.x-cntrl.x)*(node.x-cntrl.x))/(sigmaX*sigmaX)+((node.y-cntrl.y)*(node.y-cntrl.y))/(sigmaY*sigmaY)));
                node.isCenter = ( m==extCntrl.ry && n==extCntrl.rx );
                extCntrl.nodes.push_back(node);
            }
        }
    }
}

void Problem2HNDirichlet::newDistributeDeltaGauseMsmnt(const std::vector<SpacePoint> &msmnts, std::vector<ExtendedSpacePoint> &extMsmnts, const Dimension &dimX, const Dimension &dimY) const
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
        ExtendedSpacePoint &extMsmnt = extMsmnts.at(c);

        extMsmnt.x = msmnt.x;
        extMsmnt.y = msmnt.y;
        extMsmnt.rx = static_cast<unsigned int> ( round(extMsmnt.x*Nx) );
        extMsmnt.ry = static_cast<unsigned int> ( round(extMsmnt.y*Ny) );
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
                ExtendedSpacePointNode1 node;
                node.nx = n; node.x = n*hx;
                node.ny = m; node.y = m*hy;
                node.w = factor*exp(-0.5*(((node.x-msmnt.x)*(node.x-msmnt.x))/(sigmaX*sigmaX)+((node.y-msmnt.y)*(node.y-msmnt.y))/(sigmaY*sigmaY)));
                node.isCenter = ( m==extMsmnt.ry && n==extMsmnt.rx );
                extMsmnt.nodes.push_back(node);
            }
        }
    }
}

double Problem2HNDirichlet::distributeTimeDelta(double t, double ht, unsigned int ln, const espn_vector &qPointNodes, const SpaceNodePDE &sn, const std::vector<ExtendedSpacePoint> &xsps) const
{
    if ( ln >= 40 ) return 0.0;

    double Q = 0.0;
#ifdef OLD_VERSION
    for (unsigned int si=0; si<qPointNodes.size(); si++)
    {
        const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
        if (qNode.i == sn.i && qNode.j == sn.j)
        {
            Q += mEquParameter.q[qNode.id] * qNode.w;
        }
    }
#endif
#ifdef NEW_VERSION
    unsigned int Ns = xsps.size();
    for (unsigned int s=0; s<Ns; s++)
    {
        double q = mEquParameter.q[s];
        const ExtendedSpacePoint &extendedSpacePoint = xsps.at(s);
        if (extendedSpacePoint.contains(sn.i, sn.j))
        {
            const std::vector<ExtendedSpacePointNode1> &nodes = extendedSpacePoint.nodes;
            unsigned int nodes_size = nodes.size();
            for (unsigned int i=0; i<nodes_size; i++)
            {
                const ExtendedSpacePointNode1 &node = nodes.at(i);
                if (node.equals(sn))
                {
                    Q += q * node.w;
                }
            }
        }
    }
#endif
    const double sigma = 5.0*ht;
    const double mu = 20.0*ht;
    const double factor = 1.0 / (sqrt(2*M_PI)*sigma);
    return factor * exp( -0.5*((t - mu)*(t - mu))/(sigma*sigma) ) * Q;
}

void Problem2HNDirichlet::PrmToVector(const OptimizeParameter &prm, DoubleVector &pv) const
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

void Problem2HNDirichlet::VectorToPrm(const DoubleVector &pv, OptimizeParameter &prm) const
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
