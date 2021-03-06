﻿#include "problem2h_solver.h"

void Problem2HDirichlet::f_layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const
{
    //IPrinter::printMatrix(u);
    //IPrinter::printSeperatorLine();
    return;

    //if (ln%2==0) printf("ln: %d fx: %8.6f min: %8.6f max: %8.6f\n", (ln/2), integralU(u), u.min(), u.max());
    return;

    static double MIN = +100000.0;
    static double MAX = -100000.0;

    if (ln%2 == 0 || ln == 3)
    {
        double min = u.min();
        double max = u.max();
        if (MIN>min) MIN = min;
        if (MAX<max) MAX = max;
        printf("ln: %d min: %f max: %f MIN: %f MAX: %f %.10f\n", (ln/2), min, max, MIN, MAX, integralU(u));
    }
    //    return;

    if (ln == 10 || ln%20 == 0)
    {
        char filename1[40];
        int size1 = sprintf(filename1, "txt/h_layer%d.txt", (ln/2));
        filename1[size1] = 0;

        FILE* file = fopen(filename1, "w");
        //IPrinter::printSeperatorLine();
        IPrinter::printMatrix(u, u.rows(), u.cols(), nullptr, file);
        //IPrinter::printSeperatorLine();
        //printf("%f\n", sqrt(integralU(u)));
        fclose(file);
    }
    return;

    if (ln%2 == 0)
    {
        //#ifdef USE_IMAGING
        char filename1[40];
        int size1 = sprintf(filename1, "e:/data/img/image%d.png", ln);
        filename1[size1] = 0;
        //char filename2[40];
        //int size2 = sprintf(filename2, "e:/data/txt/image%d.txt", ln);
        //filename2[size2] = 0;

        double min = u.min();
        double max = u.max();
        if (MIN>min) MIN = min;
        if (MAX<max) MAX = max;

        //puts("Generating image...");
        //QPixmap pxm;
        //visualGrayScale(u, min, max, pxm, 0, 0);
        //visualizeMatrixHeat(u, min, max, pxm, 0, 0);
        //pxm.save(QString(filename1), "PNG");
        printf("Image generated. ln: %d min: %f max: %f MIN: %f MAX: %f\n", ln/2, min, max, MIN, MAX);
        //FILE* file = fopen(filename2, "w");
        //IPrinter::print(u, u.rows(), u.cols(), 10, 8, file);
        //fclose(file);
        //#endif
    }

    //    //    if (ln == 2*timeDimension().sizeN())
    //    //    {
    //    //        puts("Generating image...");
    //    //        QPixmap pxm;
    //    //        visualGrayScale(u, u.min(), u.max(), pxm, 0, 0);
    //    //        pxm.save("E:/image1000.png", "PNG");
    //    //        printf("Image generated. ln: %d min: %f max: %f\n", 1000, u.min(), u.max());
    //    //        FILE* file = fopen("E:/image1000.txt", "w");
    //    //        IPrinter::print(u, u.rows(), u.cols(), 10, 8, file);
    //    //        fclose(file);
    //    //    }

    //#ifdef SAVE_TO_IMG
    //    //if (ln != 1 && ln != timeDimension().sizeN()+LD) return;
    //    //if (ln < timeDimension().sizeN()) return;

    //    double min = u.min();
    //    double max = u.max();
    //    if (MIN>min) MIN = min;
    //    if (MAX<max) MAX = max;

    //    //    double norm = 0.0;
    //    //    for (unsigned int m=0; m<u.rows(); m++)
    //    //    {
    //    //        for (unsigned int n=0; n<u.cols(); n++)
    //    //        {
    //    //            norm += u[m][n]*u[m][n];
    //    //        }
    //    //    }

    //    QPixmap pic;
    //    visualizeMatrixHeat(u, min, max, pic);
    //    pic.save("images/f/100/pic_"+QString("%1").arg(ln)+".png", "PNG");
    //    printf("Layer: %4d min: %8.4f max: %8.4f min: %8.4f max: %8.4f diff: %8.4f norm: %10.8f\n", ln, min, max, MIN, MAX, fabs(max-min), sqrt(integralU(u)));
    //#endif
}

void Problem2HDirichlet::checkGradient1(const Problem2HDirichlet &prob)
{
//    EquationParameterH e_prm = prob.mEquParameter;
//    OptimizeParameterH o_prm = prob.mOptParameter;
//    OptimizeParameterH r_prm = prob.mRegParameter;

//    IPrinter::printSeperatorLine();
//    DoubleVector pv;
//    prob.PrmToVector(o_prm, pv);
//    IPrinter::print(pv, pv.length(), 6, 4);
//    IPrinter::printSeperatorLine();
//    DoubleVector r_pv;
//    prob.PrmToVector(r_prm, r_pv);
//    IPrinter::print(r_pv, r_pv.length(), 6, 4);
//    IPrinter::printSeperatorLine();
//    DoubleVector ag(pv.length());
//    double functional = prob.fx(pv);
//    printf("Functional: %f\n", functional);
//    puts("Calculating gradients....");
//    prob.gradient(pv, ag);
//    puts("Gradients are calculated.");

//    DoubleVector ng1(pv.length(), 0.0);
//    DoubleVector ng2(pv.length(), 0.0);

//    puts("Calculating numerical gradients.... dh=0.01");
//    puts("*** Calculating numerical gradients for k...... dh=0.01");
//    IGradient::Gradient(&prob, 0.01, pv, ng1, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
//    puts("*** Calculating numerical gradients for z...... dh=0.01");
//    IGradient::Gradient(&prob, 0.01, pv, ng1, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
//    puts("*** Calculating numerical gradients for xi..... dh=0.01");
//    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    puts("*** Calculating numerical gradients for eta.... dh=0.01");
//    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    puts("Numerical gradients are calculated.");

//    puts("Calculating numerical gradients.... hx=0.001");
//    puts("*** Calculating numerical gradients for k...... dh=0.001");
//    IGradient::Gradient(&prob, 0.001, pv, ng2, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
//    puts("*** Calculating numerical gradients for z...... dh=0.001");
//    IGradient::Gradient(&prob, 0.001, pv, ng2, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
//    puts("*** Calculating numerical gradients for xi..... dh=0.001");
//    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    puts("*** Calculating numerical gradients for eta.... dh=0.001");
//    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    puts("Numerical gradients are calculated.");

//    //k------------------------------------------------------//
//    IPrinter::printSeperatorLine("k");
//    DoubleVector pk0 = pv.mid(0, e_prm.Nc*e_prm.No-1);
//    DoubleVector ak0 = ag.mid(0, e_prm.Nc*e_prm.No-1);
//    DoubleVector nk1 = ng1.mid(0, e_prm.Nc*e_prm.No-1);
//    DoubleVector nk2 = ng2.mid(0, e_prm.Nc*e_prm.No-1);

//    IPrinter::print(pk0,pk0.length(),14,4);
//    IPrinter::print(ak0,ak0.length(),14,4); ak0.L2Normalize();
//    IPrinter::print(nk1,nk1.length(),14,4); nk1.L2Normalize();
//    IPrinter::print(nk2,nk2.length(),14,4); nk2.L2Normalize();
//    IPrinter::print(ak0,ak0.length(),14,4);
//    IPrinter::print(nk1,nk1.length(),14,4);
//    IPrinter::print(nk2,nk2.length(),14,4);

//    //z------------------------------------------------------//
//    IPrinter::printSeperatorLine("z");
//    DoubleVector pz0 = pv.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    DoubleVector az0 = ag.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    DoubleVector nz1 = ng1.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    DoubleVector nz2 = ng2.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);

//    IPrinter::print(pz0,pz0.length(),14,4);
//    IPrinter::print(az0,az0.length(),14,4); az0.L2Normalize();
//    IPrinter::print(nz1,nz1.length(),14,4); nz1.L2Normalize();
//    IPrinter::print(nz2,nz2.length(),14,4); nz2.L2Normalize();
//    IPrinter::print(az0,az0.length(),14,4);
//    IPrinter::print(nz1,nz1.length(),14,4);
//    IPrinter::print(nz2,nz2.length(),14,4);

//    //xi------------------------------------------------------//
//    IPrinter::printSeperatorLine("xi");
//    DoubleVector pe0 = pv.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    DoubleVector ae0 = ag.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    DoubleVector ne1 = ng1.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    DoubleVector ne2 = ng2.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);

//    IPrinter::print(pe0,pe0.length(),14,4);
//    IPrinter::print(ae0,ae0.length(),14,4); ae0.L2Normalize();
//    IPrinter::print(ne1,ne1.length(),14,4); ne1.L2Normalize();
//    IPrinter::print(ne2,ne2.length(),14,4); ne2.L2Normalize();
//    IPrinter::print(ae0,ae0.length(),14,4);
//    IPrinter::print(ne1,ne1.length(),14,4);
//    IPrinter::print(ne2,ne2.length(),14,4);

//    //eta------------------------------------------------------//
//    IPrinter::printSeperatorLine("eta");
//    DoubleVector px0 = pv.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    DoubleVector ax0 = ag.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    DoubleVector nx1 = ng1.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    DoubleVector nx2 = ng2.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);

//    IPrinter::print(px0,px0.length(),14,4);
//    IPrinter::print(ax0,ax0.length(),14,4); ax0.L2Normalize();
//    IPrinter::print(nx1,nx1.length(),14,4); nx1.L2Normalize();
//    IPrinter::print(nx2,nx2.length(),14,4); nx2.L2Normalize();
//    IPrinter::print(ax0,ax0.length(),14,4);
//    IPrinter::print(nx1,nx1.length(),14,4);
//    IPrinter::print(nx2,nx2.length(),14,4);
//    IPrinter::printSeperatorLine();
}

void Problem2HDirichlet::checkGradient2(const Problem2HDirichlet &prob)
{
//    EquationParameterH e_prm = prob.mEquParameter;
//    OptimizeParameterH o_prm = prob.mOptParameter;
//    OptimizeParameterH r_prm = prob.mRegParameter;

//    IPrinter::printSeperatorLine();
//    DoubleVector pv;
//    prob.PrmToVector(o_prm, pv);
//    IPrinter::print(pv, pv.length(), 6, 4);
//    IPrinter::printSeperatorLine();
//    DoubleVector r_pv;
//    prob.PrmToVector(r_prm, r_pv);
//    IPrinter::print(r_pv, r_pv.length(), 6, 4);
//    IPrinter::printSeperatorLine();
//    DoubleVector ag(pv.length());
//    double functional = prob.fx(pv);
//    printf("Functional: %f\n", functional);
//    puts("Calculating gradients....");
//    prob.gradient(pv, ag);
//    puts("Gradients are calculated.");

//    DoubleVector ng1(pv.length(), 0.0);
//    DoubleVector ng2(pv.length(), 0.0);

//    puts("Calculating numerical gradients.... dh=0.01");
//    puts("*** Calculating numerical gradients for k...... dh=0.01");
//    IGradient::Gradient(&prob, 0.01, pv, ng1, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
//    puts("*** Calculating numerical gradients for z...... dh=0.01");
//    IGradient::Gradient(&prob, 0.01, pv, ng1, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
//    puts("*** Calculating numerical gradients for xi..... dh=0.01");
//    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    puts("*** Calculating numerical gradients for eta.... dh=0.01");
//    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    puts("Numerical gradients are calculated.");

//    puts("Calculating numerical gradients.... hx=0.001");
//    puts("*** Calculating numerical gradients for k...... dh=0.001");
//    IGradient::Gradient(&prob, 0.001, pv, ng2, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
//    puts("*** Calculating numerical gradients for z...... dh=0.001");
//    IGradient::Gradient(&prob, 0.001, pv, ng2, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
//    puts("*** Calculating numerical gradients for xi..... dh=0.001");
//    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    puts("*** Calculating numerical gradients for eta.... dh=0.001");
//    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    puts("Numerical gradients are calculated.");

//    DoubleVector nag  = ag;  nag.EuclideanNormalize();
//    DoubleVector nng1 = ng1; nng1.EuclideanNormalize();
//    DoubleVector nng2 = ng2; nng2.EuclideanNormalize();

//    //k------------------------------------------------------//
//    IPrinter::printSeperatorLine("k");
//    DoubleVector pk0 = pv.mid(0, e_prm.Nc*e_prm.No-1);
//    IPrinter::print(pk0,pk0.length(),14,4);

//    DoubleVector ak0 = ag.mid(0, e_prm.Nc*e_prm.No-1);
//    DoubleVector nk1 = ng1.mid(0, e_prm.Nc*e_prm.No-1);
//    DoubleVector nk2 = ng2.mid(0, e_prm.Nc*e_prm.No-1);
//    IPrinter::print(ak0,ak0.length(),14,4);
//    IPrinter::print(nk1,nk1.length(),14,4);
//    IPrinter::print(nk2,nk2.length(),14,4);

//    DoubleVector nak0 = nag.mid(0, e_prm.Nc*e_prm.No-1);
//    DoubleVector nnk1 = nng1.mid(0, e_prm.Nc*e_prm.No-1);
//    DoubleVector nnk2 = nng2.mid(0, e_prm.Nc*e_prm.No-1);

//    IPrinter::print(nak0,nak0.length(),14,4);
//    IPrinter::print(nnk1,nnk1.length(),14,4);
//    IPrinter::print(nnk2,nnk2.length(),14,4);

//    //z------------------------------------------------------//
//    IPrinter::printSeperatorLine("z");
//    DoubleVector pz0 = pv.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    IPrinter::print(pz0,pz0.length(),14,4);

//    DoubleVector az0 = ag.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    DoubleVector nz1 = ng1.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    DoubleVector nz2 = ng2.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    IPrinter::print(az0,az0.length(),14,4);
//    IPrinter::print(nz1,nz1.length(),14,4);
//    IPrinter::print(nz2,nz2.length(),14,4);

//    DoubleVector naz0 = nag.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    DoubleVector nnz1 = nng1.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    DoubleVector nnz2 = nng2.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
//    IPrinter::print(naz0,naz0.length(),14,4);
//    IPrinter::print(nnz1,nnz1.length(),14,4);
//    IPrinter::print(nnz2,nnz2.length(),14,4);

//    //xi------------------------------------------------------//
//    IPrinter::printSeperatorLine("xi");
//    DoubleVector pe0 = pv.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    IPrinter::print(pe0,pe0.length(),14,4);

//    DoubleVector ae0 = ag.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    DoubleVector ne1 = ng1.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    DoubleVector ne2 = ng2.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    IPrinter::print(ae0,ae0.length(),14,4);
//    IPrinter::print(ne1,ne1.length(),14,4);
//    IPrinter::print(ne2,ne2.length(),14,4);

//    DoubleVector nae0 = nag.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    DoubleVector nne1 = nng1.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    DoubleVector nne2 = nng2.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
//    IPrinter::print(nae0,nae0.length(),14,4);
//    IPrinter::print(nne1,nne1.length(),14,4);
//    IPrinter::print(nne2,nne2.length(),14,4);

//    //eta------------------------------------------------------//
//    IPrinter::printSeperatorLine("eta");
//    DoubleVector px0 = pv.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    IPrinter::print(px0,px0.length(),14,4);

//    DoubleVector ax0 = ag.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    DoubleVector nx1 = ng1.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    DoubleVector nx2 = ng2.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    IPrinter::print(ax0,ax0.length(),14,4);
//    IPrinter::print(nx1,nx1.length(),14,4);
//    IPrinter::print(nx2,nx2.length(),14,4);

//    DoubleVector nax0 = nag.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    DoubleVector nnx1 = nng1.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    DoubleVector nnx2 = nng2.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
//    IPrinter::print(nax0,nax0.length(),14,4);
//    IPrinter::print(nnx1,nnx1.length(),14,4);
//    IPrinter::print(nnx2,nnx2.length(),14,4);
//    IPrinter::printSeperatorLine();
}

Problem2HDirichlet::Problem2HDirichlet()
{
//    r = 0.0;
//    regEpsilon = 0.0;
}

Problem2HDirichlet::~Problem2HDirichlet()
{}

double Problem2HDirichlet::fx(const DoubleVector &pv) const
{
    OptimizeParameterH o_prm;

    VectorToPrm(pv, o_prm);

    Problem2HDirichlet* prob = const_cast<Problem2HDirichlet*>(this);
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

double Problem2HDirichlet::integral(const std::vector<DoubleMatrix> &vu) const
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

double Problem2HDirichlet::integralU(const DoubleMatrix &u) const
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

double Problem2HDirichlet::norm(const EquationParameterH& e_prm, const OptimizeParameterH &o_prm, const OptimizeParameterH &r_prm) const
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

double Problem2HDirichlet::penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const
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

double Problem2HDirichlet::gpi(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const
{
    double gpi_ln = fabs( g0i(i, layer, u_info, o_prm) ) - ( vmax.at(i) - vmin.at(i) )/2.0;
    return gpi_ln > 0.0 ? gpi_ln : 0.0;
}

double Problem2HDirichlet::g0i(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const
{
    double vi = 0.0;
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        const SpacePointInfoH &u_xij = u_info[j];
        vi += o_prm.k[i][j] * ( u_xij.vl[layer] - o_prm.z[i][j] );
    }
    return ( vmax.at(i) + vmin.at(i) )/2.0 - vi;
}

double Problem2HDirichlet::sign(double x) const
{
    if (x < 0.0)       return -1.0;
    else if (x > 0.0)  return +1.0;
    else               return  0.0;
}

void Problem2HDirichlet::gradient(const DoubleVector & pv, DoubleVector &g) const
{
    const unsigned int L   = static_cast<const unsigned int>(mtimeDimension.size());
    const double ht        = mtimeDimension.step();
    const unsigned int Nc  = mEquParameter.Nc;
    const unsigned int No  = mEquParameter.No;
    const unsigned int LLD = L + LD;

    OptimizeParameterH o_prm;
    VectorToPrm(pv, o_prm);

    Problem2HDirichlet* prob = const_cast<Problem2HDirichlet*>(this);
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

auto Problem2HDirichlet::norm(const DoubleVector &v) const -> double
{
    return EuclideanNorm(v);
}

auto Problem2HDirichlet::normalize(DoubleVector &v) const -> void
{
    if (optimizeK) { DoubleVector kv = v.mid(0, 3);   IVectorNormalizer::EuclideanNormalize(kv); v[0]  = kv[0]; v[1]  = kv[1];  v[2]  = kv[2]; v[3]  = kv[3]; kv.clear(); }
    if (optimizeZ) { DoubleVector zv = v.mid(4, 7);   IVectorNormalizer::EuclideanNormalize(zv); v[4]  = zv[0]; v[5]  = zv[1];  v[6]  = zv[2]; v[7]  = zv[3]; zv.clear(); }
    if (optimizeO) { DoubleVector ov = v.mid(8, 11);  IVectorNormalizer::EuclideanNormalize(ov); v[8]  = ov[0]; v[9]  = ov[1];  v[10] = ov[2]; v[11] = ov[3]; ov.clear(); }
    if (optimizeZ) { DoubleVector cv = v.mid(12, 15); IVectorNormalizer::EuclideanNormalize(cv); v[12] = cv[0]; v[13] = cv[1];  v[14] = cv[2]; v[15] = cv[3]; cv.clear(); }
}

void Problem2HDirichlet::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const
{
    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(f); C_UNUSED(alpha); C_UNUSED(result);
    const char* msg = nullptr; C_UNUSED(msg);
    if (result == GradientMethod::MethodResult::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION   ";
    if (result == GradientMethod::MethodResult::FIRST_ITERATION)          msg = "FIRST_ITERATION         ";
    if (result == GradientMethod::MethodResult::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    if (result == GradientMethod::MethodResult::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS     ";
    if (result == GradientMethod::MethodResult::NEXT_ITERATION)           msg = "NEXT_ITERATION          ";

    Problem2HDirichlet* prob = const_cast<Problem2HDirichlet*>(this);
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

    //IPrinter::printVector(v1, "v1", 10);
    //IPrinter::printVector(v2, "v2", 10);

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

auto Problem2HDirichlet::project(DoubleVector &, unsigned int) -> void
{}

auto Problem2HDirichlet::project(DoubleVector &pv) const -> void
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

auto Problem2HDirichlet::projectControlPoints(DoubleVector &pv, unsigned int index) const -> void
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

auto Problem2HDirichlet::projectMeasurePoints(DoubleVector &pv, unsigned int index) const -> void
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

auto Problem2HDirichlet::solveForwardIBVP1(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use) const -> void
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

    const unsigned int Ns = mEquParameter.Ns;
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
//    std::vector<DeltaGrid2D> pulseDeltaGridList; pulseDeltaGridList.resize(Ns);
//    std::vector<DeltaGrid2D> measuremntGirdList; measuremntGirdList.resize(No);
//    std::vector<DeltaGrid2D> cntrlDeltaGridList; cntrlDeltaGridList.resize(Nc);

//    for (unsigned int s=0; s<Ns; s++)
//    {
//        pulseDeltaGridList[s].initGrid(N, hx, M, hy);
//        pulseDeltaGridList[s].distributeGauss(mEquParameter.theta[s], 8, 8);
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        measuremntGirdList[j].initGrid(N, hx, M, hy);
//        measuremntGirdList[j].distributeGauss(mOptParameter.xi[j], 1, 1);
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        cntrlDeltaGridList[i].initGrid(N, hx, M, hy);
//        cntrlDeltaGridList[i].distributeGauss(mOptParameter.eta[i], 1, 1);
//    }

//    std::vector<unsigned int> rows0, rows1, rows2, cols0, cols1, cols2;

//    for (unsigned int m=1; m<=M-1; m++)
//    {
//        bool foundC = false;
//        bool foundM = false;
//        for (unsigned int i=0; i<Nc; i++) { foundC |= cntrlDeltaGridList[i]._rows[m]; }
//        for (unsigned int j=0; j<No; j++) { foundM |= measuremntGirdList[j]._rows[m]; }
//        if ( foundC==false && foundM==false && std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
//        if ( foundC==true                   && std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
//        if ( foundC==true  && foundM==true  && std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
//    }

//    for (unsigned int n=1; n<=N-1; n++)
//    {
//        bool foundC = false;
//        bool foundM = false;
//        for (unsigned int i=0; i<Nc; i++) { foundC |= cntrlDeltaGridList[i]._cols[n]; }
//        for (unsigned int j=0; j<No; j++) { foundM |= measuremntGirdList[j]._cols[n]; }
//        if ( foundC==false && foundM==false && std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
//        if ( foundC==true                   && std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
//        if ( foundC==true  && foundM==true  && std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
//    }

    std::vector<ExtendedSpacePointH> msnExtSpacePoints, cntExtSpacePoints, qExtSpacePoints;
    newDistributeDeltaGaussPulse(mEquParameter.theta, qExtSpacePoints, dimX, dimY);
    newDistributeDeltaGaussCntrl(mOptParameter.eta, cntExtSpacePoints, dimX, dimY);
    newDistributeDeltaGaussMsmnt(mOptParameter.xi,  msnExtSpacePoints, dimX, dimY);

    //----------------------------------------------------------------------------------------------//
    GridH grid;
    uint_vectorH rows0, rows1, rows2, cols0, cols1, cols2;
    f_findRowsCols(grid, rows0, rows1, rows2, cols0, cols1, cols2, N, M, cntExtSpacePoints, msnExtSpacePoints);
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) f_prepareInfo(No, mOptParameter.xi, u_info, LLD);
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

    const unsigned int rows1_size = static_cast<const unsigned int>( rows1.size()*(N-1) );
    double *a1=nullptr, *b1=nullptr, *c1=nullptr, *d1=nullptr, *x1=nullptr, **w1=nullptr;
    if (rows1.size() != 0 && rows2.size() != 0)
    {
        a1 = static_cast<double*>(malloc(sizeof(double)*rows1_size));
        b1 = static_cast<double*>(malloc(sizeof(double)*rows1_size));
        c1 = static_cast<double*>(malloc(sizeof(double)*rows1_size));
        d1 = static_cast<double*>(malloc(sizeof(double)*rows1_size));
        x1 = static_cast<double*>(malloc(sizeof(double)*rows1_size));
        w1 = static_cast<double**>(malloc(sizeof(double*)*rows1_size));
        for (unsigned int row=0; row < rows1_size; row++) w1[row] = static_cast<double*> ( malloc(sizeof(double)*rows1_size) );
    }

    const unsigned int cols1_size = static_cast<const unsigned int>( cols1.size()*(M-1) );
    double *a2=nullptr, *b2=nullptr, *c2=nullptr, *d2=nullptr, *x2=nullptr, **w2=nullptr;
    if (cols1.size() != 0 && cols2.size() != 0)
    {
        a2 = static_cast<double*> (malloc(sizeof(double)*cols1_size) );
        b2 = static_cast<double*> (malloc(sizeof(double)*cols1_size) );
        c2 = static_cast<double*> (malloc(sizeof(double)*cols1_size) );
        d2 = static_cast<double*> (malloc(sizeof(double)*cols1_size) );
        x2 = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
        w2 = static_cast<double**> ( malloc(sizeof(double*)*cols1_size) );
        for (unsigned int col=0; col < cols1_size; col++) w2[col] = static_cast<double*> ( malloc(sizeof(double)*cols1_size) );
    }

    //------------------------------------- initial conditions -------------------------------------//
    //f_initialLayers(u00, u05, u10, u_info, use, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, alpha, qExtSpacePoints, msnExtSpacePoints, cntExtSpacePoints);
    /************************************************************************/
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
    /************************************************************************/
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*ht;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;
    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0; sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u10[m][0] = f_boundary(sn0, tn10); u05[m][0] = f_boundary(sn0, tn05);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u10[m][N] = f_boundary(sn1, tn10); u05[m][N] = f_boundary(sn1, tn05);
    }

    sn0.j = 0; sn0.y = 0.0; sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u10[0][n] = f_boundary(sn0, tn10); u05[0][n] = f_boundary(sn0, tn05);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u10[M][n] = f_boundary(sn1, tn10); u05[M][n] = f_boundary(sn1, tn05);
    }
    SpaceNodePDE sn10;
    for (unsigned int m=1; m<=M-1; m++)
    {
        sn10.j = static_cast<int>(m); sn10.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn10.i = static_cast<int>(n); sn10.x = n*hx;

            double Q = 0.0;
//            for (unsigned int s=0; s<Ns; s++)
//            {
//                Q += mEquParameter.q[s]*pulseDeltaGridList[s].weight(sn10);
//            }

            double sum = 0.0;
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= alpha*(f_initial2(sn10)+Q);

            u05[m][n] = u00[m][n] + (0.5*ht) * (f_initial2(sn10)+Q) + (0.125*ht*ht) * sum;
            u10[m][n] = u00[m][n] + (1.0*ht) * (f_initial2(sn10)+Q) + (0.500*ht*ht) * sum;
        }
    }
    /************************************************************************/
    if (use == true) f_add2Info(u00, u_info, 0, hx, hy, msnExtSpacePoints);
    f_layerInfo(u00, 0);

    if (use == true) f_add2Info(u05, u_info, 1, hx, hy, msnExtSpacePoints);
    f_layerInfo(u05, 1);

    if (use == true) f_add2Info(u10, u_info, 2, hx, hy, msnExtSpacePoints);
    f_layerInfo(u10, 2);
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
        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = rows0.at(row);
                sn.j = static_cast<int>(m); sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = static_cast<int>(n); sn.x = n*hx;
                    dx[n-1] = 0.0;
                    dx[n-1] = (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025;
                    dx[n-1] += 2.0*u10[m][n] - u05[m][n] + alpha_ht_025*u05[m][n];
                    dx[n-1] += (u10[m][n-1] - 2.0*u10[m][n] + u10[m][n+1])*p_aa_htht__hxhx_025_1m2lambda;
                    dx[n-1] += (u05[m][n-1] - 2.0*u05[m][n] + u05[m][n+1])*p_aa_htht__hxhx_025_lambda;
                }
                dx[0]   -= u15[m][0]*m_aa_htht__hxhx_025_lambda;
                dx[N-2] -= u15[m][N]*m_aa_htht__hxhx_025_lambda;
                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0/* && rows2.size() == 0*/)
        {
            //throw grid_exception("forward x1");
            currentLayerFGrid(u10);

//            double* _u15 = new double[No];
//            for (unsigned int j=0; j<No; j++)
//            {
//                const ExtendedSpacePointH &extendedSpacePoint = msnExtSpacePoints.at(j);
//                const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
//                const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
//                _u15[j] = 0.0;
//                for (unsigned int nj=0; nj<nodes_size; nj++)
//                {
//                    const ExtendedSpacePointNodeH &node = nodes.at(nj);
//                    const unsigned int node_nx = static_cast<unsigned int>(node.nx);
//                    const unsigned int node_ny = static_cast<unsigned int>(node.ny);
//                    _u15[j] += u10[node_ny][node_nx] * (node.w * (hx*hy));
//                }
//                _u15[j] *= (1.0 + noise);
//            }

//            double *_v15 = new double[Nc];
//            for (unsigned int i=0; i<Nc; i++)
//            {
//                _v15[i] = 0.0;
//                for (unsigned int j=0; j<No; j++)
//                {
//                    _v15[i] += mOptParameter.k[i][j] * (_u15[j] - mOptParameter.z[i][j]);
//                }
//            }
//            delete [] _u15;

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = static_cast<int>(m); sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = static_cast<int>(n); sn.x = n*hx;
                    dx[n-1] = 0.0;
                    dx[n-1] += (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025;
                    dx[n-1] += 2.0*u10[m][n] - u05[m][n] + alpha_ht_025*u05[m][n];
                    dx[n-1] += (u10[m][n-1] - 2.0*u10[m][n] + u10[m][n+1])*p_aa_htht__hxhx_025_1m2lambda;
                    dx[n-1] += (u05[m][n-1] - 2.0*u05[m][n] + u05[m][n+1])*p_aa_htht__hxhx_025_lambda;
                    //------------------------------------- Adding delta part -------------------------------------//
                    double fx15 = mCrFfxWeightMatrix[m][n];
//                    for (unsigned int i=0; i<Nc; i++)
//                    {
//                        const ExtendedSpacePointH &extendedSpacePoint = cntExtSpacePoints.at(i);
//                        if (extendedSpacePoint.contains(sn))
//                        {
//                            const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
//                            const unsigned int nodes_size = static_cast<unsigned int>( nodes.size() );
//                            for (unsigned int ni=0; ni<nodes_size; ni++)
//                            {
//                                const ExtendedSpacePointNodeH &node = nodes.at(ni);
//                                if (node.equals(sn)) fx15 += _v15[i] * node.w;
//                            }
//                        }
//                    }
                    dx[n-1] += ht_ht_025 * fx15;
                    //------------------------------------- Adding delta part -------------------------------------//
                }
                dx[0]   -= u15[m][0]*m_aa_htht__hxhx_025_lambda;
                dx[N-2] -= u15[m][N]*m_aa_htht__hxhx_025_lambda;
                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }
//            delete [] _v15;
        }

        if (rows1.size() != 0 && rows2.size() != 0 && false)
        {
            //throw grid_exception("forward x2");
            throw exception();

            for (unsigned int m=0; m<rows1_size; m++) for (unsigned int n=0; n<rows1_size; n++) w1[m][n] = 0.0;

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = static_cast<int>(m); sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = static_cast<int>(n); sn.x = n*hx;
                    const unsigned int index = offset+(n-1);
                    d1[index] = 0.0;
                    d1[index] += (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025;
                    d1[index] += 2.0*u10[m][n] - u05[m][n] + alpha_ht_025*u05[m][n];
                    d1[index] += (u10[m][n-1] - 2.0*u10[m][n] + u10[m][n+1])*p_aa_htht__hxhx_025_1m2lambda;
                    d1[index] += (u05[m][n-1] - 2.0*u05[m][n] + u05[m][n+1])*p_aa_htht__hxhx_025_lambda;

                    a1[index] = m_aa_htht__hxhx_025_lambda;
                    b1[index] = b_aa_htht__hxhx;
                    c1[index] = m_aa_htht__hxhx_025_lambda;
                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePointH &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                        if (cExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNodeH> &nodes1 = cExtendedSpacePoint.nodes;
                            for (unsigned int ni=0; ni<nodes1.size(); ni++)
                            {
                                const ExtendedSpacePointNodeH &node1 = nodes1.at(ni);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                const ExtendedSpacePointH &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                                const std::vector<ExtendedSpacePointNodeH> &nodes2 = mExtendedSpacePoint.nodes;
                                for (unsigned int nj=0; nj<nodes2.size(); nj++)
                                {
                                    const ExtendedSpacePointNodeH &node2 = nodes2.at(nj);
                                    const unsigned int node2_nx = static_cast<unsigned int>(node2.nx);
                                    const unsigned int node2_ny = static_cast<unsigned int>(node2.ny);

                                    bool found = false;
                                    for (unsigned int rs=0; rs<rows1.size(); rs++)
                                    {
                                        if (node2_ny == rows1[rs])
                                        {
                                            found = true;
                                            w1[index][rs*(N-1)+(node2_nx-1)] -= ht_ht_025 * mOptParameter.k[i][j] * (node2.w * (hx*hy)) * w * (1.0+noise);
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d1[index] += ht_ht_025 * mOptParameter.k[i][j] * u15[node2_ny][node2_nx] * (node2.w * (hx*hy)) * w * (1.0+noise);
                                    }
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[index] -= ht_ht_025 * mOptParameter.k[i][j] * mOptParameter.z[i][j] * w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }
                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= u15[m][0]*m_aa_htht__hxhx_025_lambda;
                d1[offset+N-2] -= u15[m][N]*m_aa_htht__hxhx_025_lambda;

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
                sn.i = static_cast<int>(n); sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = static_cast<int>(m); sn.y = m*hy;
                    dy[m-1] = 0.0;
                    dy[m-1] += (u15[m][n-1] - 2.0*u15[m][n] + u15[m][n+1])*p_aa_htht__hxhx_025;
                    dy[m-1] += 2.0*u15[m][n] - u10[m][n] + alpha_ht_025*u10[m][n];
                    dy[m-1] += (u15[m-1][n] - 2.0*u15[m][n] + u15[m+1][n])*p_aa_htht__hyhy_025_1m2lambda;
                    dy[m-1] += (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025_lambda;
                }
                dy[0]   -= u20[0][n]*m_aa_htht__hyhy_025_lambda;
                dy[M-2] -= u20[M][n]*m_aa_htht__hyhy_025_lambda;
                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0/* && cols2.size() == 0*/)
        {
            //throw grid_exception("forward y1");
            currentLayerFGrid(u15);

//            double* _u20 = new double[No];

//            for (unsigned int j=0; j<No; j++)
//            {
//                const ExtendedSpacePointH &extendedSpacePoint = msnExtSpacePoints.at(j);
//                const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
//                const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
//                _u20[j] = 0.0;
//                for (unsigned int nj=0; nj<nodes_size; nj++)
//                {
//                    const ExtendedSpacePointNodeH &node = nodes.at(nj);
//                    const unsigned int node_nx = static_cast<unsigned int>(node.nx);
//                    const unsigned int node_ny = static_cast<unsigned int>(node.ny);
//                    _u20[j] += u15[node_ny][node_nx] * (node.w * (hx*hy));
//                }
//                _u20[j] *= (1.0+noise);
//            }

//            double *_v20 = new double[Nc];
//            for (unsigned int i=0; i<Nc; i++)
//            {
//                _v20[i] = 0.0;
//                for (unsigned int j=0; j<No; j++)
//                {
//                    _v20[i] += mOptParameter.k[i][j] * ( _u20[j] - mOptParameter.z[i][j] );
//                }
//            }

//            delete [] _u20;

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = static_cast<int>(n); sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = static_cast<int>(m); sn.y = m*hy;

                    dy[m-1] = 0.0;
                    dy[m-1] += (u15[m][n-1] - 2.0*u15[m][n] + u15[m][n+1])*p_aa_htht__hxhx_025;
                    dy[m-1] += 2.0*u15[m][n] - u10[m][n] + alpha_ht_025*u10[m][n];
                    dy[m-1] += (u15[m-1][n] - 2.0*u15[m][n] + u15[m+1][n])*p_aa_htht__hyhy_025_1m2lambda;
                    dy[m-1] += (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025_lambda;
                    //------------------------------------- Adding delta part -------------------------------------//
                    double fx20 =mCrFfxWeightMatrix[m][n];
//                    for (unsigned int i=0; i<Nc; i++)
//                    {
//                        const ExtendedSpacePointH &extendedSpacePoint = cntExtSpacePoints.at(i);
//                        if (extendedSpacePoint.contains(sn))
//                        {
//                            const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
//                            const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
//                            for (unsigned int ni=0; ni<nodes_size; ni++)
//                            {
//                                const ExtendedSpacePointNodeH &node = nodes.at(ni);
//                                if (node.equals(sn)) fx20 += _v20[i] * node.w;
//                            }
//                        }
//                    }
                    dy[m-1] += ht_ht_025 * fx20;
                    //------------------------------------- Adding delta part -------------------------------------//
                }
                dy[0]   -= u20[0][n]*m_aa_htht__hyhy_025_lambda;
                dy[M-2] -= u20[M][n]*m_aa_htht__hyhy_025_lambda;
                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
            }

//            delete [] _v20;
        }

        if (cols1.size() != 0 && cols2.size() != 0 && false)
        {
            //throw grid_exception("forward y2");
            throw exception();

            for (unsigned int m=0; m<cols1_size; m++) for (unsigned int n=0; n<cols1_size; n++) w2[m][n] = 0.0;

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = static_cast<int>(n); sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = static_cast<int>(m); sn.y = m*hy;

                    const unsigned int index = offset+(m-1);
                    d2[index] = 0.0;
                    d2[index] += (u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1])*p_aa_htht__hxhx_025;
                    d2[index] += 2.0*u15[m][n] - u10[m][n] + alpha_ht_025*u10[m][n];
                    d2[index] += (u15[m-1][n] - 2.0*u15[m][n] + u15[m+1][n])*p_aa_htht__hyhy_025_1m2lambda;
                    d2[index] += (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025_lambda;

                    a2[index] = m_aa_htht__hyhy_025_lambda;
                    b2[index] = b_aa_htht__hyhy;
                    c2[index] = m_aa_htht__hyhy_025_lambda;
                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const ExtendedSpacePointH &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                        if (cExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNodeH> &nodes1 = cExtendedSpacePoint.nodes;
                            for (unsigned int ni=0; ni<nodes1.size(); ni++)
                            {
                                const ExtendedSpacePointNodeH &node1 = nodes1.at(ni);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                const ExtendedSpacePointH &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                                const std::vector<ExtendedSpacePointNodeH> &nodes2 = mExtendedSpacePoint.nodes;
                                for (unsigned int nj=0; nj<nodes2.size(); nj++)
                                {
                                    const ExtendedSpacePointNodeH &node2 = nodes2.at(nj);
                                    const unsigned int node2_nx = static_cast<unsigned int>(node2.nx);
                                    const unsigned int node2_ny = static_cast<unsigned int>(node2.ny);

                                    bool found = false;
                                    for (unsigned int cs=0; cs<cols1.size(); cs++)
                                    {
                                        if (node2_nx == cols1[cs])
                                        {
                                            found = true;
                                            w2[index][cs*(M-1)+(node2_ny-1)] -= ht_ht_025 * mOptParameter.k[i][j] * (node2.w * (hx*hy)) * w * (1.0+noise);
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d2[index] += ht_ht_025 * mOptParameter.k[i][j] * u20[node2_ny][node2_nx] * (node2.w * (hx*hy)) * w * (1.0+noise);
                                    }
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[index] -= ht_ht_025 * mOptParameter.k[i][j] * mOptParameter.z[i][j] * w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= u20[0][n]*m_aa_htht__hyhy_025_lambda;
                d2[offset+M-2] -= u20[M][n]*m_aa_htht__hyhy_025_lambda;

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

        if (use == true) f_add2Info(u15, u_info, 2*ln-1, hx, hy, msnExtSpacePoints); f_layerInfo(u15, 2*ln-1);
        if (use == true) f_add2Info(u20, u_info, 2*ln+0, hx, hy, msnExtSpacePoints); f_layerInfo(u20, 2*ln+0);

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

    qExtSpacePoints.clear();
    msnExtSpacePoints.clear();
    cntExtSpacePoints.clear();

    u00.clear();
    u10.clear();
    u15.clear();
    u20.clear();
}

auto Problem2HDirichlet::f_initialLayers(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10,
                                          spif_vectorH &u_info, bool use, unsigned int N, unsigned int M,
                                          double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda,
                                          const std::vector<ExtendedSpacePointH> &qExtSpacePoints,
                                          const std::vector<ExtendedSpacePointH> &msnExtSpacePoints,
                                          const std::vector<ExtendedSpacePointH> &/*cntExtSpacePoints*/) const -> void
{
    /************************************************************************/
    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            u00[m][n] = f_initial1(sn);
        }
    }

    /************************************************************************/
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*ht;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0; sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u10[m][0] = f_boundary(sn0, tn10); u05[m][0] = f_boundary(sn0, tn05);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u10[m][N] = f_boundary(sn1, tn10); u05[m][N] = f_boundary(sn1, tn05);
    }

    sn0.j = 0; sn0.y = 0.0; sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u10[0][n] = f_boundary(sn0, tn10); u05[0][n] = f_boundary(sn0, tn05);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u10[M][n] = f_boundary(sn1, tn10); u05[M][n] = f_boundary(sn1, tn05);
    }

    /************************************************************************/

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;

            double Q = 0.0;
            const unsigned int Ns = static_cast<const unsigned int>(qExtSpacePoints.size());
            for (unsigned int s=0; s<Ns; s++)
            {
                double q = mEquParameter.q[s];
                const ExtendedSpacePointH &extendedSpacePoint = qExtSpacePoints.at(s);
                if (extendedSpacePoint.contains(sn))
                {
                    const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
                    const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
                    for (unsigned int i=0; i<nodes_size; i++)
                    {
                        const ExtendedSpacePointNodeH &node = nodes.at(i);
                        if (node.equals(sn))
                        {
                            Q += q * node.w;
                        }
                    }
                }
            }

            double sum = 0.0;
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= lambda*(f_initial2(sn)+Q);

            u05[m][n] = u00[m][n] + (0.5*ht) * (f_initial2(sn)+Q) + (0.125*ht*ht) * sum;
            u10[m][n] = u00[m][n] + (1.0*ht) * (f_initial2(sn)+Q) + (0.500*ht*ht) * sum;
        }
    }

    /************************************************************************/
    if (use == true) f_add2Info(u00, u_info, 0, hx, hy, msnExtSpacePoints);
    f_layerInfo(u00, 0);

    if (use == true) f_add2Info(u05, u_info, 1, hx, hy, msnExtSpacePoints);
    f_layerInfo(u05, 1);

    if (use == true) f_add2Info(u10, u_info, 2, hx, hy, msnExtSpacePoints);
    f_layerInfo(u10, 2);
}

auto Problem2HDirichlet::f_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichlet::initiatePulseGrid() const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    const_cast<Problem2HDirichlet*>(this)->mPulseWeightMatrix.clear();
    const_cast<Problem2HDirichlet*>(this)->mPulseWeightMatrix.resize(M+1, N+1, 0.0);

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
            const_cast<Problem2HDirichlet*>(this)->mPulseWeightMatrix[m][n] = 0.0;
            for (unsigned int s=0; s<Ns; s++)
            {
                const_cast<Problem2HDirichlet*>(this)->mPulseWeightMatrix[m][n] += mEquParameter.q[s]*deltaGrids[s].weight(n, m);
            }
        }
    }

    for (unsigned int s=0; s<Ns; s++)
    {
        deltaGrids[s].cleanGrid();
    }
    deltaGrids.clear();
}

auto Problem2HDirichlet::currentLayerFGrid(const DoubleMatrix &u) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    const_cast<Problem2HDirichlet*>(this)->mCrFfxWeightMatrix.clear();
    const_cast<Problem2HDirichlet*>(this)->mCrFfxWeightMatrix.resize(M+1, N+1, 0.0);

    const unsigned int Nc = mEquParameter.Nc;
    std::vector<DeltaGrid2D> controlDeltaGrids(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        controlDeltaGrids[i].initGrid(N, hx, M, hy);
        controlDeltaGrids[i].distributeGauss(mOptParameter.eta[i], 1, 1);
    }

    const unsigned int No = mEquParameter.No;
    std::vector<DeltaGrid2D> measurementDeltaGrids(No);
    for (unsigned int j=0; j<No; j++)
    {
        measurementDeltaGrids[j].initGrid(N, hx, M, hy);
        measurementDeltaGrids[j].distributeGauss(mOptParameter.xi[j], 1, 1);
    }

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
            const_cast<Problem2HDirichlet*>(this)->mCrFfxWeightMatrix[m][n] = 0.0;
            for (unsigned int i=0; i<Nc; i++)
            {
                const DeltaGrid2D &dg = controlDeltaGrids[i];
                const_cast<Problem2HDirichlet*>(this)->mCrFfxWeightMatrix[m][n] += _v[i] * dg.weight(n,m);
            }
        }
    }

    delete [] _v;

    for (unsigned int j=0; j<No; j++) measurementDeltaGrids[j].cleanGrid(); measurementDeltaGrids.clear();
    for (unsigned int i=0; i<Nc; i++) controlDeltaGrids[i].cleanGrid(); controlDeltaGrids.clear();
}

auto Problem2HDirichlet::currentLayerBGrid(const DoubleMatrix &p, double ln, const spif_vectorH &u_info) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    const_cast<Problem2HDirichlet*>(this)->mCrBfxWeightMatrix.clear();
    const_cast<Problem2HDirichlet*>(this)->mCrBfxWeightMatrix.resize(M+1, N+1, 0.0);

    const unsigned int Nc = mEquParameter.Nc;
    std::vector<DeltaGrid2D> controlDeltaGrids(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        controlDeltaGrids[i].initGrid(N, hx, M, hy);
        controlDeltaGrids[i].distributeGauss(mOptParameter.eta[i], 1, 1);
    }

    const unsigned int No = mEquParameter.No;
    std::vector<DeltaGrid2D> measurementDeltaGrids(No);
    for (unsigned int j=0; j<No; j++)
    {
        measurementDeltaGrids[j].initGrid(N, hx, M, hy);
        measurementDeltaGrids[j].distributeGauss(mOptParameter.xi[j], 1, 1);
    }

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
            const_cast<Problem2HDirichlet*>(this)->mCrBfxWeightMatrix[m][n] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                const DeltaGrid2D &dg = measurementDeltaGrids[j];
                const_cast<Problem2HDirichlet*>(this)->mCrBfxWeightMatrix[m][n] += _w[j] * dg.weight(n,m);
            }
        }
    }

    delete [] _w;

    for (unsigned int j=0; j<No; j++) measurementDeltaGrids[j].cleanGrid(); measurementDeltaGrids.clear();
    for (unsigned int i=0; i<Nc; i++) controlDeltaGrids[i].cleanGrid(); controlDeltaGrids.clear();
}

auto Problem2HDirichlet::f_initial2(const SpaceNodePDE &sn) const -> double
{
    //return 0.0;
    return mPulseWeightMatrix[sn.j][sn.i];
}

double Problem2HDirichlet::f_boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 0.0;
}

void Problem2HDirichlet::f_findRowsCols(GridH &grid, uint_vectorH &rows0, uint_vectorH &rows1, uint_vectorH &rows2, uint_vectorH &cols0, uint_vectorH &cols1, uint_vectorH &cols2, unsigned int N, unsigned int M,
                                         const std::vector<ExtendedSpacePointH> &cntExtSpacePoints,
                                         const std::vector<ExtendedSpacePointH> &msnExtSpacePoints) const
{
    for (unsigned int m=1; m<=M-1; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for(std::vector<ExtendedSpacePointH>::const_iterator csp_it=cntExtSpacePoints.begin(); csp_it != cntExtSpacePoints.end(); csp_it++)
        {
            const ExtendedSpacePointH &cxsp = *csp_it;
            const std::vector<ExtendedSpacePointNodeH> &c_nodes = cxsp.nodes;
            for (std::vector<ExtendedSpacePointNodeH>::const_iterator cnode_it=c_nodes.begin(); cnode_it != c_nodes.end(); cnode_it++)
            {
                const ExtendedSpacePointNodeH &cnode = *cnode_it;
                if (static_cast<unsigned int>(cnode.ny) == m)
                {
                    found1 = true;
                    for(std::vector<ExtendedSpacePointH>::const_iterator msp_it=msnExtSpacePoints.begin(); msp_it != msnExtSpacePoints.end(); msp_it++)
                    {
                        const ExtendedSpacePointH &mxsp = *msp_it;
                        const std::vector<ExtendedSpacePointNodeH> &mnodes = mxsp.nodes;
                        for (std::vector<ExtendedSpacePointNodeH>::const_iterator mnode_it=mnodes.begin(); mnode_it != mnodes.end(); mnode_it++)
                        {
                            const ExtendedSpacePointNodeH &mnode = *mnode_it;
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

        if (found1 == false && found2 == false) if(std::find(grid.rows0.begin(), grid.rows0.end(), m) == grid.rows0.end()) grid.rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(grid.rows2.begin(), grid.rows2.end(), m) == grid.rows2.end()) grid.rows2.push_back(m);
        if (found1 == true)                     if(std::find(grid.rows1.begin(), grid.rows1.end(), m) == grid.rows1.end()) grid.rows1.push_back(m);
    }

    for (unsigned int n=1; n<=N-1; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for(std::vector<ExtendedSpacePointH>::const_iterator csp_it=cntExtSpacePoints.begin(); csp_it != cntExtSpacePoints.end(); csp_it++)
        {
            const ExtendedSpacePointH &cxsp = *csp_it;
            const std::vector<ExtendedSpacePointNodeH> &c_nodes = cxsp.nodes;
            for (std::vector<ExtendedSpacePointNodeH>::const_iterator cnode_it=c_nodes.begin(); cnode_it != c_nodes.end(); cnode_it++)
            {
                const ExtendedSpacePointNodeH &cnode = *cnode_it;
                if (static_cast<unsigned int>(cnode.nx) == n)
                {
                    found1 = true;
                    for(std::vector<ExtendedSpacePointH>::const_iterator msp_it=msnExtSpacePoints.begin(); msp_it != msnExtSpacePoints.end(); msp_it++)
                    {
                        const ExtendedSpacePointH &mxsp = *msp_it;
                        const std::vector<ExtendedSpacePointNodeH> &mnodes = mxsp.nodes;
                        for (std::vector<ExtendedSpacePointNodeH>::const_iterator mnode_it=mnodes.begin(); mnode_it != mnodes.end(); mnode_it++)
                        {
                            const ExtendedSpacePointNodeH &mnode = *mnode_it;
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

        if (found1 == false && found2 == false) if(std::find(grid.cols0.begin(), grid.cols0.end(), n) == grid.cols0.end()) grid.cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(grid.cols2.begin(), grid.cols2.end(), n) == grid.cols2.end()) grid.cols2.push_back(n);
        if (found1 == true)                     if(std::find(grid.cols1.begin(), grid.cols1.end(), n) == grid.cols1.end()) grid.cols1.push_back(n);
    }
}

void Problem2HDirichlet::b_findRowsCols(GridH &grid, uint_vectorH &rows0, uint_vectorH &rows1, uint_vectorH &rows2, uint_vectorH &cols0, uint_vectorH &cols1, uint_vectorH &cols2, unsigned int N, unsigned int M,
                                         const std::vector<ExtendedSpacePointH> &msnExtSpacePoints,
                                         const std::vector<ExtendedSpacePointH> &cntExtSpacePoints) const
{
    for (unsigned int m=1; m<=M-1; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for(std::vector<ExtendedSpacePointH>::const_iterator csp_it=msnExtSpacePoints.begin(); csp_it != msnExtSpacePoints.end(); csp_it++)
        {
            const ExtendedSpacePointH &cxsp = *csp_it;
            const std::vector<ExtendedSpacePointNodeH> &c_nodes = cxsp.nodes;
            for (std::vector<ExtendedSpacePointNodeH>::const_iterator cnode_it=c_nodes.begin(); cnode_it != c_nodes.end(); cnode_it++)
            {
                const ExtendedSpacePointNodeH &cnode = *cnode_it;
                if (static_cast<unsigned int>(cnode.ny) == m)
                {
                    found1 = true;
                    for(std::vector<ExtendedSpacePointH>::const_iterator msp_it=cntExtSpacePoints.begin(); msp_it != cntExtSpacePoints.end(); msp_it++)
                    {
                        const ExtendedSpacePointH &mxsp = *msp_it;
                        const std::vector<ExtendedSpacePointNodeH> &mnodes = mxsp.nodes;
                        for (std::vector<ExtendedSpacePointNodeH>::const_iterator mnode_it=mnodes.begin(); mnode_it != mnodes.end(); mnode_it++)
                        {
                            const ExtendedSpacePointNodeH &mnode = *mnode_it;
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

        if (found1 == false && found2 == false) if(std::find(grid.rows0.begin(), grid.rows0.end(), m) == grid.rows0.end()) grid.rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(grid.rows2.begin(), grid.rows2.end(), m) == grid.rows2.end()) grid.rows2.push_back(m);
        if (found1 == true)                     if(std::find(grid.rows1.begin(), grid.rows1.end(), m) == grid.rows1.end()) grid.rows1.push_back(m);
    }

    for (unsigned int n=1; n<=N-1; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for(std::vector<ExtendedSpacePointH>::const_iterator csp_it=msnExtSpacePoints.begin(); csp_it != msnExtSpacePoints.end(); csp_it++)
        {
            const ExtendedSpacePointH &cxsp = *csp_it;
            const std::vector<ExtendedSpacePointNodeH> &c_nodes = cxsp.nodes;
            for (std::vector<ExtendedSpacePointNodeH>::const_iterator cnode_it=c_nodes.begin(); cnode_it != c_nodes.end(); cnode_it++)
            {
                const ExtendedSpacePointNodeH &cnode = *cnode_it;
                if (static_cast<unsigned int>(cnode.nx) == n)
                {
                    found1 = true;
                    for(std::vector<ExtendedSpacePointH>::const_iterator msp_it=cntExtSpacePoints.begin(); msp_it != cntExtSpacePoints.end(); msp_it++)
                    {
                        const ExtendedSpacePointH &mxsp = *msp_it;
                        const std::vector<ExtendedSpacePointNodeH> &mnodes = mxsp.nodes;
                        for (std::vector<ExtendedSpacePointNodeH>::const_iterator mnode_it=mnodes.begin(); mnode_it != mnodes.end(); mnode_it++)
                        {
                            const ExtendedSpacePointNodeH &mnode = *mnode_it;
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

        if (found1 == false && found2 == false) if(std::find(grid.cols0.begin(), grid.cols0.end(), n) == grid.cols0.end()) grid.cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(grid.cols2.begin(), grid.cols2.end(), n) == grid.cols2.end()) grid.cols2.push_back(n);
        if (found1 == true)                     if(std::find(grid.cols1.begin(), grid.cols1.end(), n) == grid.cols1.end()) grid.cols1.push_back(n);
    }
}

void Problem2HDirichlet::f_borderLayer(DoubleMatrix &u, DoubleMatrix &um5, unsigned int ln) const
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<const unsigned int>(dimX.size());
    const unsigned int M = static_cast<const unsigned int>(dimY.size());

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

auto Problem2HDirichlet::f_prepareInfo(unsigned int No, const std::vector<SpacePoint> &points, spif_vectorH &u_info, unsigned int LLD) const -> void
{
    u_info.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        SpacePointInfoH &inf = u_info[j];
        const SpacePoint &sp = points[j];
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(2*LLD+1);
    }
}

auto Problem2HDirichlet::b_prepareInfo(unsigned int Nc, const std::vector<SpacePoint> &points, spif_vectorH &p_info, unsigned int LLD) const -> void
{
    p_info.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        SpacePointInfoH &inf = p_info[i];
        const SpacePoint &sp = points[i];
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(2*LLD+1);
    }
}

auto Problem2HDirichlet::f_add2Info(const DoubleMatrix &u, spif_vectorH &u_info, unsigned int ln, double hx, double hy, const std::vector<ExtendedSpacePointH> &extMsmnts, int method) const -> void
{
    if (method == 1 || method == 2 || method == 4)
    {
        unsigned int No = static_cast<unsigned int>(extMsmnts.size());
        for (unsigned int j=0; j<No; j++)
        {
            const ExtendedSpacePointH &xsp = extMsmnts.at(j);
            SpacePointInfoH &ui = u_info[j];
            const unsigned int nodes_size = static_cast<const unsigned int>( xsp.nodes.size() );
            for (unsigned int i=0; i<nodes_size; i++)
            {
                const ExtendedSpacePointNodeH &node = xsp.nodes.at(i);
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

                    //ui.dxx[ln] = (1.0/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
                    //ui.dyy[ln] = (1.0/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);
                }
            }
        }
    }
}

auto Problem2HDirichlet::b_add2Info(const DoubleMatrix &p, spif_vectorH &p_info, unsigned int ln, double hx, double hy, const std::vector<ExtendedSpacePointH> &extCntrls, int method) const -> void
{
    if (method == 1 || method == 2 || method == 4)
    {
        unsigned int Nc = static_cast<unsigned int>(extCntrls.size());
        for (unsigned int i=0; i<Nc; i++)
        {
            const ExtendedSpacePointH &xsp = extCntrls.at(i);
            SpacePointInfoH &pi = p_info[i];
            const unsigned int nodes_size = static_cast<const unsigned int>( xsp.nodes.size() );
            for (unsigned int i=0; i<nodes_size; i++)
            {
                const ExtendedSpacePointNodeH &node = xsp.nodes.at(i);
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

auto Problem2HDirichlet::solveBackwardIBVP1(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info) const -> void
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
    const double alpha   = mEquParameter.lambda;
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
    std::vector<ExtendedSpacePointH> cntExtSpacePoints, msnExtSpacePoints;
    newDistributeDeltaGaussCntrl(mOptParameter.eta, cntExtSpacePoints, dimX, dimY);
    newDistributeDeltaGaussMsmnt(mOptParameter.xi,  msnExtSpacePoints, dimX, dimY);
    //----------------------------------------------------------------------------------------------//
    GridH grid;
    uint_vectorH rows0, rows1, rows2, cols0, cols1, cols2;
    b_findRowsCols(grid, rows0, rows1, rows2, cols0, cols1, cols2, N, M, msnExtSpacePoints, cntExtSpacePoints);
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) b_prepareInfo(Nc, mOptParameter.eta, p_info, LLD);
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

    const unsigned int row1_size = static_cast<const unsigned int>(grid.rows1.size()*(N-1));
    double* a1=nullptr, *b1=nullptr, *c1=nullptr, *d1=nullptr, *x1=nullptr, **w1=nullptr;
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

    const unsigned int cols1_size = static_cast<const unsigned int>(cols1.size()*(M-1));
    double *a2=nullptr, *b2=nullptr, *c2=nullptr, *d2=nullptr, *x2=nullptr, **w2=nullptr;
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

    //------------------------------------- initial conditions -------------------------------------//
    b_initialLayers(p00, p05, p10, p_info, use, u_info, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, alpha, cntExtSpacePoints, msnExtSpacePoints);
    //------------------------------------- initial conditions -------------------------------------//
    SpaceNodePDE sn00;
    for (unsigned int m=0; m<=M; m++)
    {
        sn00.j = static_cast<int>(m); sn00.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn00.i = static_cast<int>(n); sn00.x = n*hx;
            p00[m][n] = f_initial1(sn00);
        }
    }
    /************************************************************************/
    TimeNodePDE tn05; tn05.i = LLD-1; tn05.t = LLD*ht + 0.5*ht;
    TimeNodePDE tn10; tn10.i = LLD-1; tn10.t = LLD*ht;
    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0; sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; p10[m][0] = b_boundary(sn0, tn10); p05[m][0] = b_boundary(sn0, tn05);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; p10[m][N] = b_boundary(sn1, tn10); p05[m][N] = b_boundary(sn1, tn05);
    }
    sn0.j = 0;  sn0.y = 0.0; sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; p10[0][n] = b_boundary(sn0, tn10); p05[0][n] = b_boundary(sn0, tn05);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; p10[M][n] = b_boundary(sn1, tn10); p05[M][n] = b_boundary(sn1, tn05);
    }
    /************************************************************************/
    double *_w = new double[No];

    double* _p00 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p00[i] = 0.0;

    for (unsigned int i=0; i<Nc; i++)
    {
        const ExtendedSpacePointH &extendedSpacePoint = cntExtSpacePoints.at(i);
        const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
        const unsigned int nodes_size = static_cast<unsigned int>( nodes.size() );
        for (unsigned int ni=0; ni<nodes_size; ni++)
        {
            const ExtendedSpacePointNodeH &node = nodes.at(ni);
            const unsigned int node_nx = static_cast<unsigned int>(node.nx);
            const unsigned int node_ny = static_cast<unsigned int>(node.ny);
            _p00[i] += p00[node_ny][node_nx] * (node.w * (hx*hy));
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

    SpaceNodePDE sn2;
    for (unsigned int m=1; m<=M-1; m++)
    {
        sn2.j = static_cast<int>(m); sn2.y = static_cast<int>(m)*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn2.i = static_cast<int>(n); sn2.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);
            sum += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);
            sum += alpha*b_initial2(sn2);

            p05[m][n] = p00[m][n] - (ht*0.5) * b_initial2(sn2) + 0.125*ht*ht*sum;
            p10[m][n] = p00[m][n] - (ht*1.0) * b_initial2(sn2) + 0.500*ht*ht*sum;

            double sum1 = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                const ExtendedSpacePointH &extendedSpacePoint = msnExtSpacePoints.at(j);
                if (extendedSpacePoint.contains(sn2))
                {
                    const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
                    const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
                    for (unsigned int nj=0; nj<nodes_size; nj++)
                    {
                        const ExtendedSpacePointNodeH &node = nodes.at(nj);
                        if (node.equals(sn2)) sum1 += _w[j] * node.w;
                    }
                }
            }

            p05[m][n] += (0.125*ht*ht)*sum1;
            p10[m][n] += (0.500*ht*ht)*sum1;
        }
    }
    /************************************************************************/

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
        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = grid.rows0.at(row);
                sn.j = static_cast<int>(m); sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = static_cast<int>(n); sn.x = n*hx;
                    dx[n-1] = 0.0;
                    dx[n-1] = (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n]) * p_aa_htht__hyhy_025;
                    dx[n-1] += 2.0*p10[m][n] - p05[m][n] + alpha_ht_025*p05[m][n];
                    dx[n-1] += (p10[m][n-1] - 2.0*p10[m][n] + p10[m][n+1]) * p_aa_htht__hxhx_025_1m2lambda;
                    dx[n-1] += (p05[m][n-1] - 2.0*p05[m][n] + p05[m][n+1]) * p_aa_htht__hxhx_025_lambda;
                    //------------------------------------- Adding functional part --------------------------------//
                    if (L <= ln && ln <= LLD) dx[n-1] += -2.0*mu(sn.i,sn.j)*(u.at(2*(ln-L+1))[m][n]) * ht_ht_025;
                    //------------------------------------- Adding functional part --------------------------------//
                }
                dx[0]   -= p15[m][0] * m_aa_htht__hxhx_025_lambda;
                dx[N-2] -= p15[m][N] * m_aa_htht__hxhx_025_lambda;
                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0 /*&& rows2.size() == 0*/)
        {
            //throw grid_exception("forward x1");
            currentLayerBGrid(p10, 2*ln+2, u_info);

//            double* _p15 = new double[Nc];
//            for (unsigned int i=0; i<Nc; i++)
//            {
//                const ExtendedSpacePointH &extendedSpacePoint = cntExtSpacePoints.at(i);
//                const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
//                const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
//                _p15[i] = 0.0;
//                for (unsigned int ni=0; ni<nodes_size; ni++)
//                {
//                    const ExtendedSpacePointNodeH &node = nodes.at(ni);
//                    const unsigned int node_nx = static_cast<unsigned int>(node.nx);
//                    const unsigned int node_ny = static_cast<unsigned int>(node.ny);
//                    _p15[i] += p10[node_ny][node_nx] * (node.w * (hx*hy));
//                }
//            }

//            double *_w15 = new double[No];
//            for (unsigned int j=0; j<No; j++)
//            {
//                _w15[j] = 0.0;
//                for (unsigned int i=0; i<Nc; i++)
//                {
//                    _w15[j] += mOptParameter.k[i][j] * ( _p15[i] + 2.0*r*gpi(i, 2*ln+1, u_info, mOptParameter)*sgn(g0i(i, 2*ln+1, u_info, mOptParameter)));
//                }
//            }
//            delete [] _p15;

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = static_cast<int>(m); sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = static_cast<int>(n); sn.x = n*hx;
                    dx[n-1] = 0.0;
                    dx[n-1] += (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n]) * p_aa_htht__hyhy_025;
                    dx[n-1] += 2.0*p10[m][n] - p05[m][n] + alpha_ht_025*p05[m][n];
                    dx[n-1] += (p10[m][n-1] - 2.0*p10[m][n] + p10[m][n+1]) * p_aa_htht__hxhx_025_1m2lambda;
                    dx[n-1] += (p05[m][n-1] - 2.0*p05[m][n] + p05[m][n+1]) * p_aa_htht__hxhx_025_lambda;
                    //------------------------------------- Adding delta part -------------------------------------//
                    double fx15 = mCrBfxWeightMatrix[m][n];
//                    for (unsigned int j=0; j<No; j++)
//                    {
//                        const ExtendedSpacePointH &extendedSpacePoint = msnExtSpacePoints.at(j);
//                        if (extendedSpacePoint.contains(sn))
//                        {
//                            const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
//                            const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
//                            for (unsigned int nj=0; nj<nodes_size; nj++)
//                            {
//                                const ExtendedSpacePointNodeH &node = nodes.at(nj);
//                                if (node.equals(sn)) fx15 += _w15[j] * node.w;
//                            }
//                        }
//                    }
                    dx[n-1] += ht_ht_025 * fx15;
                    //------------------------------------- Adding delta part -------------------------------------//
                    //------------------------------------- Adding functional part --------------------------------//
                    if (L <= ln && ln <= LLD) dx[n-1] += -2.0*mu(sn.i,sn.j)*(u.at(2*(ln-L)+2)[m][n]) * ht_ht_025;
                    //------------------------------------- Adding functional part --------------------------------//
                }
                dx[0]   -= p15[m][0] * m_aa_htht__hxhx_025_lambda;
                dx[N-2] -= p15[m][N] * m_aa_htht__hxhx_025_lambda;
                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
            }
//            delete [] _w15;
        }

        if (rows1.size() != 0 && rows2.size() != 0 && false)
        {
            //throw grid_exception("forward x2");
            throw exception();

            for (unsigned int m=0; m<row1_size; m++) for (unsigned int n=0; n<row1_size; n++) w1[m][n] = 0.0;

            unsigned int offset = 0;
            for (unsigned int row=0; row<grid.rows1.size(); row++)
            {
                unsigned int m = grid.rows1.at(row);
                sn.j = static_cast<int>(m); sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = static_cast<int>(n); sn.x = n*hx;

                    const unsigned int index = offset+(n-1);
                    d1[index] = 0.0;
                    d1[index] += (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n])*p_aa_htht__hyhy_025;
                    d1[index] += 2.0*p10[m][n] - p05[m][n] + alpha_ht_025*p05[m][n];
                    d1[index] += (p10[m][n-1] - 2.0*p10[m][n] + p10[m][n+1])*p_aa_htht__hxhx_025_1m2lambda;
                    d1[index] += (p05[m][n-1] - 2.0*p05[m][n] + p05[m][n+1])*p_aa_htht__hxhx_025_lambda;

                    a1[index] = m_aa_htht__hxhx_025_lambda;
                    b1[index] = b_aa_htht__hxhx;
                    c1[index] = m_aa_htht__hxhx_025_lambda;
                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int j=0; j<No; j++)
                    {
                        const ExtendedSpacePointH &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                        if (mExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNodeH> &nodes1 = mExtendedSpacePoint.nodes;
                            for (unsigned int nj=0; nj<nodes1.size(); nj++)
                            {
                                const ExtendedSpacePointNodeH &node1 = nodes1.at(nj);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                const ExtendedSpacePointH &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                                const std::vector<ExtendedSpacePointNodeH> &nodes2 = cExtendedSpacePoint.nodes;
                                for (unsigned int ni=0; ni<nodes2.size(); ni++)
                                {
                                    const ExtendedSpacePointNodeH &node2 = nodes2.at(ni);
                                    const unsigned int node2_nx = static_cast<unsigned int>(node2.nx);
                                    const unsigned int node2_ny = static_cast<unsigned int>(node2.ny);

                                    bool found = false;
                                    for (unsigned int rs=0; rs<grid.rows1.size(); rs++)
                                    {
                                        if (node2_ny == grid.rows1[rs])
                                        {
                                            found = true;
                                            w1[index][rs*(N-1)+(node2_nx-1)] -= ht_ht_025 * mOptParameter.k[i][j] * (node2.w * (hx*hy)) * w;
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d1[index] += ht_ht_025 * mOptParameter.k[i][j] * p15[node2_ny][node2_nx] * (node2.w * (hx*hy)) * w;
                                    }
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1[index] += 2.0 * r * ht_ht_025 *  mOptParameter.k[i][j] * gpi(i, 2*ln+1, u_info, mOptParameter)*sgn(g0i(i, 2*ln+1, u_info, mOptParameter)) * w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                    //------------------------------------- Adding functional part --------------------------------//
                    if (L <= ln && ln <= LLD) d1[index] += -2.0*mu(sn.i,sn.j)*(u.at(2*(ln-L)+1)[m][n]) * ht_ht_025;
                    //------------------------------------- Adding functional part --------------------------------//
                }
                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= p15[m][0] * m_aa_htht__hxhx_025_lambda;
                d1[offset+N-2] -= p15[m][N] * m_aa_htht__hxhx_025_lambda;

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w1, x1, row1_size);

            offset = 0;
            for (unsigned int row=0; row<grid.rows1.size(); row++)
            {
                unsigned int m=grid.rows1.at(row);
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
            for (unsigned int col=0; col<grid.cols0.size(); col++)
            {
                unsigned int n = grid.cols0.at(col);
                sn.i = sn.i = static_cast<int>(n); sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = sn.i = static_cast<int>(m); sn.y = m*hy;

                    dy[m-1] = 0.0;
                    dy[m-1] += (p15[m][n-1] - 2.0*p15[m][n] + p15[m][n+1]) * p_aa_htht__hxhx_025;
                    dy[m-1] += 2.0*p15[m][n] - p10[m][n] + alpha_ht_025*p10[m][n];
                    dy[m-1] += (p15[m-1][n] - 2.0*p15[m][n] + p15[m+1][n]) * p_aa_htht__hyhy_025_1m2lambda;
                    dy[m-1] += (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n]) * p_aa_htht__hyhy_025_lambda;
                    //------------------------------------- Adding functional part --------------------------------//
                    if (L <= ln && ln <= LLD) dy[m-1] += -2.0*mu(sn.i,sn.j)*(u.at(2*(ln-L+1)-1)[m][n]) * ht_ht_025;
                    //------------------------------------- Adding functional part --------------------------------//
                }
                dy[0]   -= p20[0][n] * m_aa_htht__hyhy_025_lambda;
                dy[M-2] -= p20[M][n] * m_aa_htht__hyhy_025_lambda;
                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0 /*&& cols2.size() == 0*/)
        {
            //throw grid_exception("forward y1");
            currentLayerBGrid(p15, 2*ln+1, u_info);

//            double* _p20 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p20[i] = 0.0;

//            for (unsigned int i=0; i<Nc; i++)
//            {
//                const ExtendedSpacePointH &extendedSpacePoint = cntExtSpacePoints.at(i);
//                const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
//                const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
//                _p20[i] = 0.0;
//                for (unsigned int ni=0; ni<nodes_size; ni++)
//                {
//                    const ExtendedSpacePointNodeH &node = nodes.at(ni);
//                    const unsigned int node_nx = static_cast<const unsigned int>(node.nx);
//                    const unsigned int node_ny = static_cast<const unsigned int>(node.ny);
//                    _p20[i] += p15[node_ny][node_nx] * (node.w * (hx*hy));
//                }
//            }

//            double *_w20 = new double[No];
//            for (unsigned int j=0; j<No; j++)
//            {
//                _w20[j] = 0.0;
//                for (unsigned int i=0; i<Nc; i++)
//                {
//                    _w20[j] += mOptParameter.k[i][j] * (_p20[i] + 2.0*r*gpi(i, 2*ln+0, u_info, mOptParameter)*sgn(g0i(i, 2*ln+0, u_info, mOptParameter)));
//                }
//            }

//            delete [] _p20;

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = grid.cols1.at(col);
                sn.i = static_cast<int>(n); sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = static_cast<int>(m); sn.y = m*hy;

                    dy[m-1] = 0.0;
                    dy[m-1] += (p15[m][n-1] - 2.0*p15[m][n] + p15[m][n+1]) * p_aa_htht__hxhx_025;
                    dy[m-1] += 2.0*p15[m][n] - p10[m][n] + alpha_ht_025*p10[m][n];
                    dy[m-1] += (p15[m-1][n] - 2.0*p15[m][n] + p15[m+1][n]) * p_aa_htht__hyhy_025_1m2lambda;
                    dy[m-1] += (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n]) * p_aa_htht__hyhy_025_lambda;
                    //------------------------------------- Adding delta part -------------------------------------//
                    double fx20 = mCrBfxWeightMatrix[m][n];
//                    for (unsigned int j=0; j<No; j++)
//                    {
//                        const ExtendedSpacePointH &extendedSpacePoint = msnExtSpacePoints.at(j);
//                        if (extendedSpacePoint.contains(sn))
//                        {
//                            const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
//                            const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
//                            for (unsigned int nj=0; nj<nodes_size; nj++)
//                            {
//                                const ExtendedSpacePointNodeH &node = nodes.at(nj);
//                                if (node.equals(sn)) fx20 += _w20[j] * node.w;
//                            }
//                        }
//                    }
                    dy[m-1] += ht_ht_025 * fx20;
                    //------------------------------------- Adding delta part -------------------------------------//
                    //------------------------------------- Adding functional part --------------------------------//
                    if (L <= ln && ln <= LLD) dy[m-1] += -2.0*mu(sn.i,sn.j)*(u.at(2*(ln-L)+1)[m][n]) * ht_ht_025;
                    //------------------------------------- Adding functional part --------------------------------//
                }
                dy[0]   -= p20[0][n] * m_aa_htht__hyhy_025_lambda;
                dy[M-2] -= p20[M][n] * m_aa_htht__hyhy_025_lambda;
                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
            }

//            delete [] _w20;
        }

        if (cols1.size() != 0 && cols2.size() != 0 && false)
        {
            //throw grid_exception("forward y2");
            throw exception();

            for (unsigned int m=0; m<cols1_size; m++) for (unsigned int n=0; n<cols1_size; n++) w2[m][n] = 0.0;

            unsigned int offset = 0;
            for (unsigned int col=0; col<grid.cols1.size(); col++)
            {
                unsigned int n = grid.cols1.at(col);
                sn.i = static_cast<int>(n); sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = static_cast<int>(m); sn.y = m*hy;

                    const unsigned int index = offset+(m-1);
                    d2[index] = 0.0;
                    d2[index] += (p15[m][n-1] - 2.0*p15[m][n] + p15[m][n+1]) * p_aa_htht__hxhx_025;
                    d2[index] += 2.0*p15[m][n] - p10[m][n] + alpha_ht_025*p10[m][n];
                    d2[index] += (p15[m-1][n] - 2.0*p15[m][n] + p15[m+1][n]) * p_aa_htht__hyhy_025_1m2lambda;
                    d2[index] += (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n]) * p_aa_htht__hyhy_025_lambda;

                    a2[index] = p_aa_htht__hyhy_025_lambda;
                    b2[index] = b_aa_htht__hyhy;
                    c2[index] = m_aa_htht__hyhy_025_lambda;
                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int j=0; j<No; j++)
                    {
                        const ExtendedSpacePointH &mExtendedSpacePoint = msnExtSpacePoints.at(j);
                        if (mExtendedSpacePoint.contains(sn))
                        {
                            double w = 0.0;
                            const std::vector<ExtendedSpacePointNodeH> &nodes1 = mExtendedSpacePoint.nodes;
                            for (unsigned int nj=0; nj<nodes1.size(); nj++)
                            {
                                const ExtendedSpacePointNodeH &node1 = nodes1.at(nj);
                                if (node1.equals(sn))
                                {
                                    w = node1.w;
                                    break;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                const ExtendedSpacePointH &cExtendedSpacePoint = cntExtSpacePoints.at(i);
                                const std::vector<ExtendedSpacePointNodeH> &nodes2 = cExtendedSpacePoint.nodes;
                                for (unsigned int ni=0; ni<nodes2.size(); ni++)
                                {
                                    const ExtendedSpacePointNodeH &node2 = nodes2.at(ni);
                                    const unsigned int node2_nx = static_cast<const unsigned int>(node2.nx);
                                    const unsigned int node2_ny = static_cast<const unsigned int>(node2.ny);

                                    bool found = false;
                                    for (unsigned int cs=0; cs<grid.cols1.size(); cs++)
                                    {
                                        if (node2_nx == grid.cols1[cs])
                                        {
                                            found = true;
                                            w2[index][cs*(M-1)+(node2_ny-1)] -= ht_ht_025 * mOptParameter.k[i][i] * (node2.w * (hx*hy)) * w;
                                            break;
                                        }
                                    }

                                    if (!found)
                                    {
                                        d2[index] += ht_ht_025 * mOptParameter.k[i][i] * p20[node2_ny][node2_nx] * (node2.w * (hx*hy)) * w;
                                    }
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d2[index] += 2.0 * r * ht_ht_025 *  mOptParameter.k[i][j] * gpi(i, 2*ln, u_info, mOptParameter)*sgn(g0i(i, 2*ln, u_info, mOptParameter)) * w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                    //------------------------------------- Adding functional part --------------------------------//
                    if (L <= ln && ln <= LLD) d2[index] += -2.0*mu(sn.i,sn.j)*(u[2*(ln-L)][m][n]) * ht_ht_025;
                    //------------------------------------- Adding functional part --------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= p20[0][n] * p_aa_htht__hyhy_025_lambda;
                d2[offset+M-2] -= p20[M][n] * p_aa_htht__hyhy_025_lambda;

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2, x2, cols1_size);

            offset = 0;
            for (unsigned int col=0; col<grid.cols1.size(); col++)
            {
                unsigned int n=grid.cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    p20[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }
        }

        /**************************************************** y direction apprx ***************************************************/

        if (use == true) b_add2Info(p15, p_info, 2*ln+1, hx, hy, cntExtSpacePoints); b_layerInfo(p15, 2*ln+1);
        if (use == true) b_add2Info(p20, p_info, 2*ln+0, hx, hy, cntExtSpacePoints); b_layerInfo(p20, 2*ln+0);

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

    if (grid.rows1.size() != 0 && grid.rows2.size() != 0)
    {
        for (unsigned int row=0; row < row1_size; row++) free(w1[row]); free(w1);
        free(x1);
        free(d1);
        free(c1);
        free(b1);
        free(a1);
    }

    if (grid.cols1.size() != 0 && grid.cols2.size() != 0)
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

    grid.rows0.clear();
    grid.rows1.clear();
    grid.rows2.clear();

    grid.cols0.clear();
    grid.cols1.clear();
    grid.cols2.clear();

    msnExtSpacePoints.clear();
    cntExtSpacePoints.clear();

    p00.clear();
    p10.clear();
    p15.clear();
    p20.clear();
}

auto Problem2HDirichlet::b_initialLayers(DoubleMatrix &p00, DoubleMatrix &p05, DoubleMatrix &p10,
                                          spif_vectorH &p_info, bool use, const spif_vectorH &u_info,
                                          unsigned int N, unsigned int M, double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda,
                                          const std::vector<ExtendedSpacePointH> &cntExtSpacePoints,
                                          const std::vector<ExtendedSpacePointH> &msnExtSpacePoints) const -> void
{
    /************************************************************************/
    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            p00[m][n] = b_initial1(sn);
        }
    }

    /************************************************************************/
    const unsigned int L = static_cast<unsigned int>(mtimeDimension.size());
    const unsigned int LLD = L+LD;

    TimeNodePDE tn05; tn05.i = LLD-1; tn05.t = LLD*ht + 0.5*ht;
    TimeNodePDE tn10; tn10.i = LLD-1; tn10.t = LLD*ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0; sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; p10[m][0] = b_boundary(sn0, tn10); p05[m][0] = b_boundary(sn0, tn05);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; p10[m][N] = b_boundary(sn1, tn10); p05[m][N] = b_boundary(sn1, tn05);
    }

    sn0.j = 0;  sn0.y = 0.0; sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; p10[0][n] = b_boundary(sn0, tn10); p05[0][n] = b_boundary(sn0, tn05);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; p10[M][n] = b_boundary(sn1, tn10); p05[M][n] = b_boundary(sn1, tn05);
    }
    /************************************************************************/

    unsigned int No = mEquParameter.No;
    unsigned int Nc = mEquParameter.Nc;

    double *_w = new double[No];

    double* _p00 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p00[i] = 0.0;

    for (unsigned int i=0; i<Nc; i++)
    {
        const ExtendedSpacePointH &extendedSpacePoint = cntExtSpacePoints.at(i);
        const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
        const unsigned int nodes_size = static_cast<unsigned int>( nodes.size() );
        for (unsigned int ni=0; ni<nodes_size; ni++)
        {
            const ExtendedSpacePointNodeH &node = nodes.at(ni);
            const unsigned int node_nx = static_cast<unsigned int>(node.nx);
            const unsigned int node_ny = static_cast<unsigned int>(node.ny);
            _p00[i] += p00[node_ny][node_nx] * (node.w * (hx*hy));
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
        sn.j = static_cast<int>(m); sn.y = static_cast<int>(m)*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);
            sum += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);
            sum += lambda*b_initial2(sn);

            p05[m][n] = p00[m][n] - (ht*0.5) * b_initial2(sn) + 0.125*ht*ht*sum;
            p10[m][n] = p00[m][n] - (ht*1.0) * b_initial2(sn) + 0.500*ht*ht*sum;

            double sum1 = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                const ExtendedSpacePointH &extendedSpacePoint = msnExtSpacePoints.at(j);
                if (extendedSpacePoint.contains(sn))
                {
                    const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
                    const unsigned int nodes_size = static_cast<unsigned int>(nodes.size());
                    for (unsigned int nj=0; nj<nodes_size; nj++)
                    {
                        const ExtendedSpacePointNodeH &node = nodes.at(nj);
                        if (node.equals(sn)) sum1 += _w[j] * node.w;
                    }
                }
            }

            p05[m][n] += (0.125*ht*ht)*sum1;
            p10[m][n] += (0.500*ht*ht)*sum1;
        }
    }

    /************************************************************************/

    if (use == true) b_add2Info(p00, p_info, 2*LLD, hx, hy, cntExtSpacePoints);
    b_layerInfo(p00, 2*LLD);

    if (use == true) b_add2Info(p05, p_info, 2*(LLD-1)+1, hx, hy, cntExtSpacePoints);
    b_layerInfo(p05, 2*(LLD-1)+1);

    if (use == true) b_add2Info(p10, p_info, 2*(LLD-1)+0, hx, hy, cntExtSpacePoints);
    b_layerInfo(p10, 2*(LLD-1)+0);

    delete [] _w;
}

auto Problem2HDirichlet::b_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichlet::b_initial2(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichlet::b_boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double
{
    return 0.0;
}

auto Problem2HDirichlet::b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const -> double
{
    return -2.0*mu(n,m)*(u[m][n]);
}

auto Problem2HDirichlet::b_layerInfo(const DoubleMatrix &p UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{
//    IPrinter::printMatrix(p);
//    IPrinter::printSeperatorLine();
//    return;

//    double min = p.min();
//    double max = p.max();

//    QPixmap pic;
//    visualizeMatrixHeat(p, min, max, pic);
//    pic.save("images/b/pic"+QString("%1").arg(ln)+".png", "PNG");
//    printf("Layer: %d min: %f max: %f min: %f max: %f norm: %f\n", ln, min, max, min, max, fabs(max-min));

//    if (ln == 0)
//    {
//        IPrinter::printSeperatorLine();
//        IPrinter::printMatrix(p);
//        IPrinter::printSeperatorLine();
//        printf("Backward: %f\n", sqrt(integralU(p)));
//    }
}

auto Problem2HDirichlet::newDistributeDeltaGaussPulse(const std::vector<SpacePoint> &thetas, std::vector<ExtendedSpacePointH> &extThetas, const Dimension &dimX, const Dimension &dimY) const -> void
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.size() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.size() );
    unsigned int Ns = static_cast<unsigned int> ( thetas.size() );

    extThetas.clear();
    extThetas.resize(Ns);

    const int k = 4;
    const double sigmaX = 1.0*hx;
    const double sigmaY = 1.0*hy;

    for (unsigned int s=0; s<Ns; s++)
    {
        const SpacePoint &theta = thetas.at(s);
        ExtendedSpacePointH &extTheta = extThetas.at(s);

        extTheta.x = theta.x;
        extTheta.y = theta.y;
        extTheta.rx = static_cast<int> ( round(extTheta.x*Nx) );
        extTheta.ry = static_cast<int> ( round(extTheta.y*Ny) );
        extTheta.k = k;
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
        //double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);

        for (int m=extTheta.minY; m<=extTheta.maxY; m++)
        {
            for (int n=extTheta.minX; n<=extTheta.maxX; n++)
            {
                ExtendedSpacePointNodeH node;
                node.nx = n; node.x = n*hx;
                node.ny = m; node.y = m*hy;
                node.w = factor*exp(-0.5*(((node.x-theta.x)*(node.x-theta.x))/(sigmaX*sigmaX)+((node.y-theta.y)*(node.y-theta.y))/(sigmaY*sigmaY)));
                node.isCenter = ( m==extTheta.ry && n==extTheta.rx );
                extTheta.nodes.push_back(node);
            }
        }
    }
}

auto Problem2HDirichlet::newDistributeDeltaGaussCntrl(const std::vector<SpacePoint> &cntrls,std::vector<ExtendedSpacePointH> &extCntrls,
                                                       const Dimension &dimX, const Dimension &dimY) const -> void
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.size() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.size() );
    unsigned int Nc = static_cast<unsigned int> ( cntrls.size() );

    extCntrls.clear();
    extCntrls.resize(Nc);

    int k = 4;
    double sigmaX = hx;
    double sigmaY = hy;

    for (unsigned int c=0; c<Nc; c++)
    {
        const SpacePoint &cntrl = cntrls.at(c);
        ExtendedSpacePointH &extCntrl = extCntrls.at(c);

        extCntrl.x = cntrl.x;
        extCntrl.y = cntrl.y;
        extCntrl.rx = static_cast<int>(round(extCntrl.x*Nx));
        extCntrl.ry = static_cast<int>(round(extCntrl.y*Ny));
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
                ExtendedSpacePointNodeH node;
                node.nx = n; node.x = n*hx;
                node.ny = m; node.y = m*hy;
                node.w = factor*exp(-0.5*(((node.x-cntrl.x)*(node.x-cntrl.x))/(sigmaX*sigmaX)+((node.y-cntrl.y)*(node.y-cntrl.y))/(sigmaY*sigmaY)));
                node.isCenter = ( m==extCntrl.ry && n==extCntrl.rx );
                extCntrl.nodes.push_back(node);
            }
        }
    }

    //    EquationParameterHE &prm = const_cast<EquationParameterHE&>(this->mParameter);
    //    for (unsigned int i=0; i<Nc; i++)
    //    {
    //        SpacePoint &cntrl = prm.eta[i];
    //        EquationParameterHE::SpacePointExt &spe = prm.eta_ext[i];

    //        spe.rx = static_cast<unsigned int>( round(cntrl.x*Nx));
    //        spe.ry = static_cast<unsigned int>( round(cntrl.y*Ny));
    //        spe.minX = spe.rx - k;
    //        spe.maxX = spe.rx + k;
    //        spe.minY = spe.ry - k;
    //        spe.maxY = spe.ry + k;

    //        double sumX = 0.0;
    //        for (unsigned int n=spe.minX; n<=spe.maxX; n++) sumX += exp(-((n*hx-cntrl.x)*(n*hx-cntrl.x))/(2.0*sigmaX*sigmaX));
    //        sumX *= hx;

    //        double sumY = 0.0;
    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++) sumY += exp(-((m*hy-cntrl.y)*(m*hy-cntrl.y))/(2.0*sigmaY*sigmaY));
    //        sumY *= hy;

    //        double sigma = (sumX*sumY) / (2.0*M_PI);
    //        double factor = 1.0/((2.0*M_PI)*sigma);

    //        spe.nodes.clear();
    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++)
    //        {
    //            for (unsigned int n=spe.minX; n<=spe.maxX; n++)
    //            {
    //                EquationParameterHE::SpacePointExt::Node node;
    //                node.nx = n; node.x = n*hx;
    //                node.ny = m; node.y = m*hy;
    //                node.w = factor*exp(-0.5*(((node.x-cntrl.x)*(node.x-cntrl.x))/(sigmaX*sigmaX)+((node.y-cntrl.y)*(node.y-cntrl.y))/(sigmaY*sigmaY)));
    //                node.isCenter = ( m==spe.ry && n==spe.rx );
    //                spe.nodes.push_back(node);
    //            }
    //        }
    //    }
}

auto Problem2HDirichlet::newDistributeDeltaGaussMsmnt(const std::vector<SpacePoint> &msmnts, std::vector<ExtendedSpacePointH> &extMsmnts, const Dimension &dimX, const Dimension &dimY) const -> void
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = static_cast<unsigned int> ( dimX.size() );
    unsigned int Ny = static_cast<unsigned int> ( dimY.size() );
    unsigned int No = static_cast<unsigned int> ( msmnts.size() );

    extMsmnts.clear();
    extMsmnts.resize(No);

    int k = 4;
    double sigmaX = hx;
    double sigmaY = hy;

    for (unsigned int c=0; c<No; c++)
    {
        const SpacePoint &msmnt = msmnts.at(c);
        ExtendedSpacePointH &extMsmnt = extMsmnts.at(c);

        extMsmnt.x = msmnt.x;
        extMsmnt.y = msmnt.y;
        extMsmnt.rx = static_cast<int>(round(extMsmnt.x*Nx));
        extMsmnt.ry = static_cast<int>(round(extMsmnt.y*Ny));
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
                ExtendedSpacePointNodeH node;
                node.nx = n; node.x = n*hx;
                node.ny = m; node.y = m*hy;
                node.w = factor*exp(-0.5*(((node.x-msmnt.x)*(node.x-msmnt.x))/(sigmaX*sigmaX)+((node.y-msmnt.y)*(node.y-msmnt.y))/(sigmaY*sigmaY)));
                node.isCenter = ( m==extMsmnt.ry && n==extMsmnt.rx );
                extMsmnt.nodes.push_back(node);
            }
        }
    }

    //    EquationParameterHE &prm = const_cast<EquationParameterHE&>(this->mParameter);
    //    for (unsigned int j=0; j<No; j++)
    //    {
    //        SpacePoint &msmnt = prm.xi[j];
    //        EquationParameterHE::SpacePointExt &spe = prm.xi_ext[j];

    //        spe.rx = static_cast<unsigned int>(round(msmnt.x*Nx));
    //        spe.ry = static_cast<unsigned int>(round(msmnt.y*Ny));
    //        spe.minX = spe.rx - k;
    //        spe.maxX = spe.rx + k;
    //        spe.minY = spe.ry - k;
    //        spe.maxY = spe.ry + k;

    //        double sumX = 0.0;
    //        for (unsigned int n=spe.minX; n<=spe.maxX; n++) sumX += exp(-((n*hx-msmnt.x)*(n*hx-msmnt.x))/(2.0*sigmaX*sigmaX));
    //        sumX *= hx;

    //        double sumY = 0.0;
    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++) sumY += exp(-((m*hy-msmnt.y)*(m*hy-msmnt.y))/(2.0*sigmaY*sigmaY));
    //        sumY *= hy;

    //        double sigma = (sumX*sumY) / (2.0*M_PI);
    //        double factor = 1.0/((2.0*M_PI)*sigma);

    //        spe.nodes.clear();
    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++)
    //        {
    //            for (unsigned int n=spe.minX; n<=spe.maxX; n++)
    //            {
    //                EquationParameterHE::SpacePointExt::Node node;
    //                node.nx = n; node.x = n*hx;
    //                node.ny = m; node.y = m*hy;
    //                node.w = factor*exp(-0.5*(((node.x-msmnt.x)*(node.x-msmnt.x))/(sigmaX*sigmaX)+((node.y-msmnt.y)*(node.y-msmnt.y))/(sigmaY*sigmaY)));
    //                node.isCenter = ( m==spe.ry && n==spe.rx );
    //                spe.nodes.push_back(node);
    //            }
    //        }
    //    }
}

double Problem2HDirichlet::distributeTimeDelta(double t, double ht, unsigned int ln, const SpaceNodePDE &sn, const std::vector<ExtendedSpacePointH> &xsps) const
{
    return 0.0;

    if ( ln >= 40 ) return 0.0;

    double Q = 0.0;

    const unsigned int Ns = static_cast<const unsigned int>(xsps.size());

    for (unsigned int s=0; s<Ns; s++)
    {
        double q = mEquParameter.q[s];
        const ExtendedSpacePointH &extendedSpacePoint = xsps.at(s);
        if (extendedSpacePoint.contains(sn))
        {
            const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
            const unsigned int nodes_size = static_cast<const unsigned int>(nodes.size());
            //printf_s("%d\n", nodes_size);
            for (unsigned int i=0; i<nodes_size; i++)
            {
                const ExtendedSpacePointNodeH &node = nodes.at(i);
                if (node.equals(sn))
                {
                    Q += q * node.w;
                }
            }
        }
    }

    const double sigma = 5.0*ht;
    const double mu = 20.0*ht;
    const double factor = 1.0 / (sqrt(2*M_PI)*sigma);
    return factor * exp( -0.5*((t - mu)*(t - mu))/(sigma*sigma) ) * Q;
}

void Problem2HDirichlet::PrmToVector(const OptimizeParameterH &prm, DoubleVector &pv) const
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

void Problem2HDirichlet::VectorToPrm(const DoubleVector &pv, OptimizeParameterH &prm) const
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

auto Problem2HDirichlet::v(unsigned int i, OptimizeParameterH o_prm, EquationParameterH e_prm, const spif_vectorH &u_info, unsigned int ln) const -> double
{
    const unsigned int No = static_cast<const unsigned int>(e_prm.No);
    double v = 0.0;
    for (unsigned int j=0; j<No; j++)
    {
        v += o_prm.k[i][j] * (u_info[j].vl[ln]-o_prm.z[i][j]);
    }
    return v;
}
