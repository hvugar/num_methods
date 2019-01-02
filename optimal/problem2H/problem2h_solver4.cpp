#include "problem2h_solver4.h"

auto Problem2HNDirichlet4::f_layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{
    return;
    if (ln==500)
    {
        IPrinter::printSeperatorLine("+layerInfo: ln == 2");
        IPrinter::printMatrix(u);
    }

    return;
    {
        Problem2HNDirichlet4* tmp = const_cast<Problem2HNDirichlet4*>(this);
        std::vector<DoubleMatrix> &rvu = tmp->vu;
        rvu.push_back(u);
        if (ln > LD) rvu.erase(rvu.begin());
        if (ln == 501)
        {
            tmp->mOptParameter.k[0][0] = 0.0;
            tmp->mOptParameter.k[0][1] = 0.0;
            tmp->mOptParameter.k[1][0] = 0.0;
            tmp->mOptParameter.k[1][1] = 0.0;
            //tmp->mEquParameter.lambda = 0.01;
        }
        //std::cout << ln << " " << rvu.size() << std::endl;
        if (rvu.size() == 51)
        {
            double fx = integral(rvu);
            //printf("%d,%.10f\n", ln-50, fx);
            printf("%.10f\n", fx);
        }
    }


    //printf("ln: %d fx: %8.6f min: %8.6f max: %8.6f\n", ln, integralU(u), u.min(), u.max());
    return;

    static double MIN = +100000.0;
    static double MAX = -100000.0;

    //if (ln%2 == 0 || ln == 3)
    //{
    //    double min = u.min();
    //    double max = u.max();
    //    if (MIN>min) MIN = min;
    //    if (MAX<max) MAX = max;
    //    printf("ln: %d min: %f max: %f MIN: %f MAX: %f %.10f\n", (ln/2), min, max, MIN, MAX, integralU(u));
    //}
    //return;

    if (ln == 0 || ln == 1 || ln <= 50 || ln%100 == 0)
    {
        char filename1[40];
        int size1 = sprintf(filename1, "txt/h_layer%d.txt", ln);
        filename1[size1] = 0;

        FILE* file = fopen(filename1, "w");
        //IPrinter::printSeperatorLine();
        IPrinter::printMatrix(u, u.rows(), u.cols(), nullptr, file);
        //IPrinter::printSeperatorLine();
        printf("%d %f %f %f %f %f %f\n", ln, sqrt(integralU(u)), integralU(u), integralU(u)*integralU(u), u.min(), u.max(), integralU(u));
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

auto Problem2HNDirichlet4::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const -> void
{
    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(f); C_UNUSED(alpha); C_UNUSED(result);
    const char* msg = nullptr; C_UNUSED(msg);
    if (result == GradientMethod::MethodResult::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION   ";
    if (result == GradientMethod::MethodResult::FIRST_ITERATION)          msg = "FIRST_ITERATION         ";
    if (result == GradientMethod::MethodResult::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    if (result == GradientMethod::MethodResult::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS     ";
    if (result == GradientMethod::MethodResult::NEXT_ITERATION)           msg = "NEXT_ITERATION          ";

    Problem2HNDirichlet4* prob = const_cast<Problem2HNDirichlet4*>(this);
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
        v1[ln] = v(0, o_prm, mEquParameter, u_info, ln);
        v2[ln] = v(1, o_prm, mEquParameter, u_info, ln);
    }

    IPrinter::printVector(v1, "v1", 10);
    IPrinter::printVector(v2, "v2", 10);

    //    printf("I[%3d]: F:%.5f I:%.5f P:%.5f N:%.5f R:%.3f e:%.3f a:%.6f ", i, f, ing, pnt, nrm, r, regEpsilon, alpha);
    //    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %6.4f %6.4f %6.4f %6.4f c: %6.4f %6.4f %6.4f %6.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);

    printf("I[%3d]: F:%.5f I:%.5f P:%.5f N:%.5f R:%.3f e:%.3f a:%.6f ", i, f, ing, pnt, nrm, r, regEpsilon, alpha);
    printf("min:%.6f max:%.6f min:%.6f max:%.6f U0:%.8f UT:%.8f", u.at(0).min(), u.at(0).max(), u.at(LD).min(), u.at(LD).max(), integralU(u[0]), integralU(u[LD]));
    printf("\n");
    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
    DoubleVector n = g; n.L2Normalize();
    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15]);

    u.clear();
    u_info.clear();

    C_UNUSED(prob);
    IPrinter::printSeperatorLine();

    //    prob->optimizeK = i%4 == 3;
    //    prob->optimizeZ = i%4 == 0;
    //    prob->optimizeO = i%4 == 1;
    //    prob->optimizeC = i%4 == 2;
}

auto Problem2HNDirichlet4::checkGradient1(const Problem2HNDirichlet4 &prob) -> void
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

auto Problem2HNDirichlet4::checkGradient2(const Problem2HNDirichlet4 &prob) -> void
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

Problem2HNDirichlet4::Problem2HNDirichlet4()
{
    r = 0.0;
    regEpsilon = 0.0;
}

Problem2HNDirichlet4::~Problem2HNDirichlet4() {}

double Problem2HNDirichlet4::fx(const DoubleVector &pv) const
{
    OptimizeParameterH o_prm;

    VectorToPrm(pv, o_prm);

    //    Problem2HNDirichletForward1 forward;
    //    forward.setTimeDimension(timeDimension());
    //    forward.addSpaceDimension()
    //    forward.setParameters(mEquParameter, o_prm, LD);

    Problem2HNDirichlet4* prob = const_cast<Problem2HNDirichlet4*>(this);
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

double Problem2HNDirichlet4::integral(const std::vector<DoubleMatrix> &vu) const
{
    unsigned int size = static_cast<unsigned int>(vu.size());
    const double ht = timeDimension().step();
    double sum = 0.0;
    sum += 0.5*integralU(vu[0]);
    for (unsigned int l=1; l<=size-1; l++) { sum += integralU(vu[l]); }
    sum += 0.5*integralU(vu[size-1]);
    return sum*ht;
}

double Problem2HNDirichlet4::integralU(const DoubleMatrix &u) const
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

double Problem2HNDirichlet4::norm(const EquationParameterH& e_prm, const OptimizeParameterH &o_prm, const OptimizeParameterH &r_prm) const
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

double Problem2HNDirichlet4::penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const
{
    const double ht = mtimeDimension.step();
    const unsigned int L = static_cast<const unsigned int> ( mtimeDimension.size() );

    double pnlt = 0.0;
    for (unsigned int i=0; i<mEquParameter.Nc; i++)
    {
        double pnlt_i = 0.0;
        double _gpi_0 = gpi(i, 0, info, o_prm);
        pnlt_i += 0.5*_gpi_0*_gpi_0;
        for (unsigned int ln=1; ln<=L+LD-1; ln++)
        {
            double _gpi_l = gpi(i, ln, info, o_prm);
            pnlt_i += _gpi_l*_gpi_l;
        }
        double _gpi_L = gpi(i, (L+LD), info, o_prm);
        pnlt_i += 0.5*_gpi_L*_gpi_L;

        pnlt += pnlt_i*ht;
    }

    return pnlt;
}

double Problem2HNDirichlet4::gpi(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const
{
    double gpi_ln = fabs( g0i(i, layer, u_info, o_prm) ) - ( vmax.at(i) - vmin.at(i) )/2.0;
    return gpi_ln > 0.0 ? gpi_ln : 0.0;
}

double Problem2HNDirichlet4::g0i(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const
{
    double vi = 0.0;
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        const SpacePointInfoH &u_xij = u_info[j];
        vi += o_prm.k[i][j] * ( u_xij.vl[layer] - o_prm.z[i][j] );
    }
    return ( vmax.at(i) + vmin.at(i) )/2.0 - vi;
}

double Problem2HNDirichlet4::sign(double x) const
{
    if (x < 0.0)       return -1.0;
    else if (x > 0.0)  return +1.0;
    else               return  0.0;
}

void Problem2HNDirichlet4::gradient(const DoubleVector & pv, DoubleVector &g) const
{
    const unsigned int L   = static_cast<const unsigned int>(mtimeDimension.size());
    const double ht        = mtimeDimension.step();
    const unsigned int Nc  = mEquParameter.Nc;
    const unsigned int No  = mEquParameter.No;
    const unsigned int LLD = L + LD;

    OptimizeParameterH o_prm;
    VectorToPrm(pv, o_prm);

    Problem2HNDirichlet4* prob = const_cast<Problem2HNDirichlet4*>(this);
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

                grad_Kij += 0.5 * (pi.vl[0] + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * (uj.vl[0] - zij);
                for (unsigned int ln=1; ln<=LLD-1; ln++)
                {
                    grad_Kij += (pi.vl[ln] + 2.0*r*gpi(i,ln,u_info,o_prm)*sgn(g0i(i,ln,u_info,o_prm))) * (uj.vl[ln] - zij);
                }
                grad_Kij += 0.5 * (pi.vl[LLD] + 2.0*r*gpi(i,LLD,u_info,o_prm)*sgn(g0i(i,LLD,u_info,o_prm))) * (uj.vl[LLD] - zij);

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
                for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].vl[ln] + 2.0*r*gpi(i,ln,u_info,o_prm)*sgn(g0i(i,ln,u_info,o_prm)));
                gradXijX += uj.dx[ln] * vi;
                gradXijY += uj.dy[ln] * vi;
            }

            vi = 0.0;
            for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].vl[LLD] + 2.0*r*gpi(i,LLD,u_info,o_prm)*sgn(g0i(i,LLD,u_info,o_prm)));
            gradXijX += 0.5 * uj.dx[LLD] * vi;
            gradXijY += 0.5 * uj.dy[LLD] * vi;

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
                for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].vl[ln] - o_prm.z[i][j]);
                gradEtaiX += pi.dx[ln] * vi;
                gradEtaiY += pi.dy[ln] * vi;
            }

            vi = 0.0;
            for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].vl[LLD] - o_prm.z[i][j]);
            gradEtaiX += 0.5 * pi.dx[LLD] * vi;
            gradEtaiY += 0.5 * pi.dy[LLD] * vi;

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

auto Problem2HNDirichlet4::norm(const DoubleVector &v) const -> double
{
    return EuclideanNorm(v);
}

auto Problem2HNDirichlet4::normalize(DoubleVector &v) const -> void
{
    if (optimizeK) { DoubleVector kv = v.mid(0, 3);   IVectorNormalizer::EuclideanNormalize(kv); v[0]  = kv[0]; v[1]  = kv[1];  v[2]  = kv[2]; v[3]  = kv[3]; kv.clear(); }
    if (optimizeZ) { DoubleVector zv = v.mid(4, 7);   IVectorNormalizer::EuclideanNormalize(zv); v[4]  = zv[0]; v[5]  = zv[1];  v[6]  = zv[2]; v[7]  = zv[3]; zv.clear(); }
    if (optimizeO) { DoubleVector ov = v.mid(8, 11);  IVectorNormalizer::EuclideanNormalize(ov); v[8]  = ov[0]; v[9]  = ov[1];  v[10] = ov[2]; v[11] = ov[3]; ov.clear(); }
    if (optimizeZ) { DoubleVector cv = v.mid(12, 15); IVectorNormalizer::EuclideanNormalize(cv); v[12] = cv[0]; v[13] = cv[1];  v[14] = cv[2]; v[15] = cv[3]; cv.clear(); }
}

auto Problem2HNDirichlet4::project(DoubleVector &, unsigned int) -> void {}

auto Problem2HNDirichlet4::project(DoubleVector &pv) const -> void
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

auto Problem2HNDirichlet4::projectControlPoints(DoubleVector &pv, unsigned int index) const -> void
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

auto Problem2HNDirichlet4::projectMeasurePoints(DoubleVector &pv, unsigned int index) const -> void
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

void Problem2HNDirichlet4::initDeltaGrids(std::vector<DeltaGrid2D> &pulseDeltaGridList, std::vector<DeltaGrid2D> &msrntDeltaGridList, std::vector<DeltaGrid2D> &cntrlDeltaGridList,
                                          const EquationParameterH &equationParameter, const OptimizeParameterH &optimizeParameter) const
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();

    const unsigned int No = equationParameter.No;
    const unsigned int Nc = equationParameter.Nc;
    const unsigned int Ns = equationParameter.Ns;

    pulseDeltaGridList.resize(Ns);
    for (unsigned int s=0; s<Ns; s++)
    {
        SpacePoint sp = equationParameter.theta[s];
        pulseDeltaGridList[s].initGrid(N,hx,M,hy);
        pulseDeltaGridList[s].distributeGauss(sp, 8, 8);
    }

    msrntDeltaGridList.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        SpacePoint sp = optimizeParameter.xi[j];
        msrntDeltaGridList[j].initGrid(N,hx,M,hy);
        msrntDeltaGridList[j].distributeGauss(sp);
    }

    cntrlDeltaGridList.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        SpacePoint sp = optimizeParameter.eta[i];
        cntrlDeltaGridList[i].initGrid(N,hx,M,hy);
        cntrlDeltaGridList[i].distributeGauss(sp);
    }
}

void Problem2HNDirichlet4::releaseDeltaGrids(std::vector<DeltaGrid2D> &pulseDeltaGridList, std::vector<DeltaGrid2D> &msrntDeltaGridList, std::vector<DeltaGrid2D> &cntrlDeltaGridList) const
{
    for (auto it = pulseDeltaGridList.begin(); it != pulseDeltaGridList.end(); it++) it->cleanGrid(); pulseDeltaGridList.clear();
    //    for (unsigned int s=0; s<pulseDeltaGridList.size(); s++) pulseDeltaGridList[s].cleanGrid(); pulseDeltaGridList.clear();
    for (unsigned int j=0; j<msrntDeltaGridList.size(); j++) msrntDeltaGridList[j].cleanGrid(); msrntDeltaGridList.clear();
    for (unsigned int i=0; i<cntrlDeltaGridList.size(); i++) cntrlDeltaGridList[i].cleanGrid(); cntrlDeltaGridList.clear();
}

void Problem2HNDirichlet4::f_currentLayer(const DoubleMatrix &u10, const EquationParameterH &equationParameter, const OptimizeParameterH &optimizeParameter,
                                          const std::vector<DeltaGrid2D> &msrntDeltaGridList, const std::vector<DeltaGrid2D> &cntrlDeltaGridList,
                                          unsigned int N, double hx, unsigned int M, double hy, DoubleMatrix &fxv, unsigned int) const
{
    unsigned int No = equationParameter.No;
    unsigned int Nc = equationParameter.Nc;

    double* _u = new double[No];
    double *_v = new double[Nc];

    for (unsigned int j=0; j<No; j++)
    {
        _u[j] = 0.0;
        const DeltaGrid2D &mdg = msrntDeltaGridList[j];
        for (unsigned int m=mdg.minY(); m<=mdg.maxY(); m++)
        {
            for (unsigned int n=mdg.minX(); n<=mdg.maxX(); n++)
            {
                _u[j] += u10[m][n] * mdg.weight(n,m) * (hx*hy);
            }
        }
        //_u[j] *= (1.0 + noise * (rand()%2==0 ? +1.0 : -1.0));
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        _v[i] = 0.0;
        for (unsigned int j=0; j<No; j++)
        {
            _v[i] += mOptParameter.k[i][j] * (_u[j] - mOptParameter.z[i][j]);
        }
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            fxv[m][n] = 0.0;
            for (unsigned int i=0; i<Nc; i++)
            {
                fxv[m][n] += _v[i] * cntrlDeltaGridList[i].weight(n,m);
            }
        }
    }

    delete [] _u;
    delete [] _v;
}

void Problem2HNDirichlet4::b_currentLayer(const DoubleMatrix &p10, const EquationParameterH &equationParameter, const OptimizeParameterH &optimizeParameter,
                                          const std::vector<DeltaGrid2D> &msrntDeltaGridList, const std::vector<DeltaGrid2D> &cntrlDeltaGridList,
                                          unsigned int N, double hx, unsigned int M, double hy, DoubleMatrix &fxv, unsigned int ln, const spif_vectorH &u_info) const
{
    unsigned int No = equationParameter.No;
    unsigned int Nc = equationParameter.Nc;

    double* _p = new double[Nc];
    double *_w = new double[No];

    for (unsigned int i=0; i<Nc; i++)
    {
        _p[i] = 0.0;
        const DeltaGrid2D &cdg = cntrlDeltaGridList[i];
        for (unsigned int m=cdg.minY(); m<=cdg.maxY(); m++)
        {
            for (unsigned int n=cdg.minX(); n<=cdg.maxX(); n++)
            {
                _p[i] += p10[m][n] * cdg.weight(n,m) * (hx*hy);
            }
        }

    }

    for (unsigned int j=0; j<No; j++)
    {
        _w[j] = 0.0;
        for (unsigned int i=0; i<Nc; i++)
        {
            double _gpi = gpi(i, ln+1, u_info, mOptParameter);
            double _sgn = sgn(g0i(i, ln+1, u_info, mOptParameter));
            _w[j] += mOptParameter.k[j][i] * (_p[i] + 2.0*r*_gpi*_sgn);
        }
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            fxv[m][n] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                fxv[m][n] += _w[j] * msrntDeltaGridList[j].weight(n,m);
            }
        }
    }

    delete [] _p;
    delete [] _w;
}

auto Problem2HNDirichlet4::solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimY.size() );
    const unsigned int L = static_cast<unsigned int> ( time.size() );
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a        = mEquParameter.a;
    const double lambda   = mEquParameter.lambda;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;
    const unsigned int Ns = mEquParameter.Ns;

    const double lambda_ht_05 = lambda*ht*0.5;
    const double inv__alpha_ht = 1.0/(1.0 + lambda_ht_05);
    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double aa_htht__hyhy = ((a*a*ht*ht)/(hy*hy));

    const double ht_ht = ht*ht;
    const double lambda_ht = lambda*ht;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    //----------------------------------------------------------------------------------------------//
    std::vector<DeltaGrid2D> pulseDeltaGridList, msrntDeltaGridList, cntrlDeltaGridList;
    initDeltaGrids(pulseDeltaGridList, msrntDeltaGridList, cntrlDeltaGridList, mEquParameter, mOptParameter);
    //----------------------------------------------------------------------------------------------//
    if (use == true) f_prepareInfo(No, mOptParameter.xi, u_info, LLD);
    //----------------------------------------------------------------------------------------------//
    //------------------------------------- initial conditions -------------------------------------//
    f_initialLayers(u00, u10, u_info, use, N, hx, M, hy, ht, aa__hxhx, aa__hyhy, lambda, pulseDeltaGridList, msrntDeltaGridList, cntrlDeltaGridList);
    //------------------------------------- initial conditions -------------------------------------//

    for (unsigned int l=0; l<u.size(); l++) u[l].clear(); u.clear();
    unsigned int u_size = LD + 1;
    u.resize(u_size); for (unsigned int l=0; l<u_size; l++) u[l].resize(M+1, N+1);

    DoubleMatrix fxv(M+1, N+1);

    for (unsigned int ln=2; ln<=LLD; ln++)
    {
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;

        /**************************************************** border conditions ***************************************************/
        f_borderCalculate(u20, N, hx, M, hy, tn20);
        /**************************************************** border conditions ***************************************************/
        f_currentLayer(u10, mEquParameter, mOptParameter, msrntDeltaGridList, cntrlDeltaGridList, N, hx, M, hy, fxv, ln);
        /**************************************************************************************************************************/
        SpaceNodePDE sn;
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                u20[m][n] = inv__alpha_ht * (aa_htht__hxhx*(u10[m][n-1]-2.0*u10[m][n]+u10[m][n+1]) + aa_htht__hyhy*(u10[m-1][n]-2.0*u10[m][n]+u10[m+1][n])
                        + lambda_ht_05*u00[m][n] + 2.0*u10[m][n] - u00[m][n] + ht_ht*fxv[m][n]);
            }
        }
        /**************************************************************************************************************************/
        f_layerInfo(u20, ln);
        if (use == true) f_add2Info(u20, u_info, ln, hx, hy, msrntDeltaGridList);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }

        /**************************************************** saving last LD layers ***********************************************/
        if (L <= ln)
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    u[ln-L][m][n] = u20[m][n];
                }
            }
        }
        /**************************************************** saving last LD layers ***********************************************/
    }

    fxv.clear();
    releaseDeltaGrids(pulseDeltaGridList, msrntDeltaGridList, cntrlDeltaGridList);

    u00.clear();
    u10.clear();
    u20.clear();
}

auto Problem2HNDirichlet4::f_initialLayers(DoubleMatrix &u00, DoubleMatrix &u10, spif_vectorH &u_info, bool use,
                                           unsigned int N, double hx, unsigned int M, double hy,
                                           double ht, double aa__hxhx, double aa__hyhy, double lambda,
                                           const std::vector<DeltaGrid2D> &pulseDeltaGridList,
                                           const std::vector<DeltaGrid2D> &msrntDeltaGridList,
                                           const std::vector<DeltaGrid2D> &cntrlDeltaGridList) const -> void
{
    unsigned int Ns = mEquParameter.Ns;
    /***********************************************************************************************/
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
    f_layerInfo(u00, 0);
    if (use == true) f_add2Info(u00, u_info, 0, hx, hy, msrntDeltaGridList);
    /***********************************************************************************************/
    TimeNodePDE tn10; tn10.i = 1; tn10.t = tn10.i*ht;
    f_borderCalculate(u10, N, hx, M, hy, tn10);
    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;

            double q_value = 0.0;
            for (unsigned int s=0; s<Ns; s++)
            {
                q_value += mEquParameter.q[s]*pulseDeltaGridList[s].weight(sn);
            }

            double sum = 0.0;
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= lambda*(f_initial2(sn)+q_value);

            u10[m][n] = u00[m][n] + ht*(f_initial2(sn)+q_value) + 0.5*ht*ht*sum;
        }
    }
    f_layerInfo(u10, 1);
    if (use == true) f_add2Info(u10, u_info, 1, hx, hy, msrntDeltaGridList);
    /***********************************************************************************************/
}

auto Problem2HNDirichlet4::b_initialLayers(DoubleMatrix &p00, DoubleMatrix &p10, spif_vectorH &p_info, bool use,
                                           unsigned int N, double hx, unsigned int M, double hy,
                                           double ht, double aa__hxhx, double aa__hyhy, double lambda,
                                           const std::vector<DeltaGrid2D> &pulseDeltaGridList,
                                           const std::vector<DeltaGrid2D> &msrntDeltaGridList,
                                           const std::vector<DeltaGrid2D> &cntrlDeltaGridList) const -> void
{
    unsigned int Ns = mEquParameter.Ns;
    const Dimension time = timeDimension();
    const unsigned int L = static_cast<unsigned int> ( time.size() );
    const unsigned int LLD = L+LD;
    /***********************************************************************************************/
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
    b_layerInfo(p00, LLD);
    if (use == true) b_add2Info(p00, p_info, LLD, hx, hy, cntrlDeltaGridList);
    /***********************************************************************************************/
    //    double *_w = new double[No];
    //    double* _p00 = new double[Nc]; for (unsigned int i=0; i<Nc; i++) _p00[i] = 0.0;
    //    for (unsigned int i=0; i<Nc; i++)
    //    {
    //        const ExtendedSpacePointH &extendedSpacePoint = cntExtSpacePoints.at(i);
    //        const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
    //        const unsigned int nodes_size = static_cast<unsigned int>( nodes.size() );
    //        for (unsigned int ni=0; ni<nodes_size; ni++)
    //        {
    //            const ExtendedSpacePointNodeH &node = nodes.at(ni);
    //            const unsigned int node_nx = static_cast<unsigned int>(node.nx);
    //            const unsigned int node_ny = static_cast<unsigned int>(node.ny);
    //            _p00[i] += p00[node_ny][node_nx] * (node.w * (hx*hy));
    //        }
    //    }

    //    for (unsigned int j=0; j<No; j++)
    //    {
    //        _w[j] = 0.0;
    //        for (unsigned int i=0; i<Nc; i++)
    //        {
    //            _w[j] += mOptParameter.k[i][j] * (_p00[i] + 2.0*r*gpi(i, LLD, u_info, mOptParameter)*sgn(g0i(i,LLD, u_info, mOptParameter)));
    //        }
    //    }

    //    delete [] _p00;

    TimeNodePDE tn10; tn10.i = LLD-1; tn10.t = tn10.i*ht;
    f_borderCalculate(p10, N, hx, M, hy, tn10);
    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);
            sum += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);
            sum += lambda*b_initial2(sn);

            p10[m][n] = p00[m][n] - ht*b_initial2(sn) + 0.5*ht*ht*sum;
        }
    }
    b_layerInfo(p10, LLD-1);
    if (use == true) b_add2Info(p10, p_info, LLD-1, hx, hy, cntrlDeltaGridList);
    /***********************************************************************************************/
}

void Problem2HNDirichlet4::f_borderCalculate(DoubleMatrix &u, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const
{
    SpaceNodePDE sn0, sn1;
    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u[m][0] = f_boundary(sn0, tn);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u[m][N] = f_boundary(sn1, tn);
    }
    sn0.j = static_cast<int>(0); sn0.y = 0*hy;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u[0][n] = f_boundary(sn0, tn);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u[M][n] = f_boundary(sn1, tn);
    }
}

void Problem2HNDirichlet4::b_borderCalculate(DoubleMatrix &p, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const
{
    SpaceNodePDE sn0, sn1;
    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; p[m][0] = b_boundary(sn0, tn);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; p[m][N] = b_boundary(sn1, tn);
    }
    sn0.j = static_cast<int>(0); sn0.y = 0*hy;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; p[0][n] = b_boundary(sn0, tn);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; p[M][n] = b_boundary(sn1, tn);
    }
}

auto Problem2HNDirichlet4::solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int M = static_cast<unsigned int>( dimY.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a        = mEquParameter.a;
    const double lambda   = mEquParameter.lambda;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;

    const double lambda_ht_05 = lambda*ht*0.5;
    const double inv__alpha_ht = 1.0/(1.0 + lambda_ht_05);
    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double aa_htht__hyhy = ((a*a*ht*ht)/(hy*hy));

    const double ht_ht = ht*ht;
    const double lambda_ht = lambda*ht;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<DeltaGrid2D> pulseDeltaGridList, msrntDeltaGridList, cntrlDeltaGridList;
    initDeltaGrids(pulseDeltaGridList, msrntDeltaGridList, cntrlDeltaGridList, mEquParameter, mOptParameter);
    //----------------------------------------------------------------------------------------------//
    if (use == true) b_prepareInfo(Nc, mOptParameter.eta, p_info, LLD);
    //----------------------------------------------------------------------------------------------//
    //------------------------------------- initial conditions -------------------------------------//
    b_initialLayers(p00, p10, p_info, use, N, hx, M, hy, ht, aa__hxhx, aa__hyhy, lambda, pulseDeltaGridList, msrntDeltaGridList, cntrlDeltaGridList);
    //------------------------------------- initial conditions -------------------------------------//

    DoubleMatrix fxv(M+1, N+1);

    const unsigned int size_ln = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=LLD-2; ln != size_ln; ln--)
    {
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = ln+1; tn10.t = tn10.i*ht;

        /**************************************************** border conditions ***************************************************/
        b_borderCalculate(p20, N, hx, M, hy, tn20);
        /**************************************************** border conditions ***************************************************/
        b_currentLayer(p10, mEquParameter, mOptParameter, msrntDeltaGridList, cntrlDeltaGridList, N, hx, M, hy, fxv, ln, u_info);
        /**************************************************************************************************************************/
        SpaceNodePDE sn;
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                p20[m][n] = inv__alpha_ht * (aa_htht__hxhx*(p10[m][n-1]-2.0*p10[m][n]+p10[m][n+1]) + aa_htht__hyhy*(p10[m-1][n]-2.0*p10[m][n]+p10[m+1][n])
                        + lambda_ht_05*p00[m][n] + 2.0*p10[m][n] - p00[m][n] + ht_ht*fxv[m][n]);
                //------------------------------------- Adding functional part --------------------------------//
                double _mu = mu(static_cast<unsigned int>(sn.i),static_cast<unsigned int>(sn.j));
                if (L <= ln && ln <= LLD) p20[m][n] += -2.0*_mu*(u[ln-L+1][m][n]) * ht_ht * inv__alpha_ht;
                //------------------------------------- Adding functional part --------------------------------//
            }
        }
        /**************************************************************************************************************************/
        b_layerInfo(p20, ln);
        if (use == true) b_add2Info(p20, p_info, ln, hx, hy, cntrlDeltaGridList);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p20[m][n];
            }
        }
    }

    fxv.clear();
    releaseDeltaGrids(pulseDeltaGridList, msrntDeltaGridList, cntrlDeltaGridList);

    p00.clear();
    p10.clear();
    p20.clear();
}

auto Problem2HNDirichlet4::f_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HNDirichlet4::f_initial2(const SpaceNodePDE &) const -> double
{
    return 0.0;
}

auto Problem2HNDirichlet4::f_boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HNDirichlet4::f_borderLayer(DoubleMatrix &u, DoubleMatrix &um5, unsigned int ln) const -> void
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

auto Problem2HNDirichlet4::f_prepareInfo(unsigned int No, const std::vector<SpacePoint> &points, spif_vectorH &u_info, unsigned int LLD) const -> void
{
    u_info.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        SpacePointInfoH &inf = u_info[j];
        const SpacePoint &sp = points[j];
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(LLD+1);
    }
}

auto Problem2HNDirichlet4::b_prepareInfo(unsigned int Nc, const std::vector<SpacePoint> &points, spif_vectorH &p_info, unsigned int LLD) const -> void
{
    p_info.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        SpacePointInfoH &inf = p_info[i];
        const SpacePoint &sp = points[i];
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(LLD+1);
    }
}

auto Problem2HNDirichlet4::f_add2Info(const DoubleMatrix &u, spif_vectorH &u_info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &msrntDeltaGridList) const -> void
{
    unsigned int No = static_cast<unsigned int>(msrntDeltaGridList.size());
    for (unsigned int j=0; j<No; j++)
    {
        SpacePointInfoH &ui =  u_info[j];
        const DeltaGrid2D &mdg = msrntDeltaGridList[j];

        ui.vl[ln] = 0.0;
        for (unsigned int m=mdg.minY(); m<=mdg.maxY(); m++)
        {
            for (unsigned int n=mdg.minX(); n<=mdg.maxX(); n++)
            {
                ui.vl[ln] += u[m][n] * mdg.weight(n,m) * (hx*hy);
            }
        }

        double px = mdg.p().x;
        double py = mdg.p().y;
        unsigned int rx = mdg.rx();
        unsigned int ry = mdg.ry();

        ui.dx[ln] = (u[ry][rx+1] - u[ry][rx-1])/(2.0*hx);
        ui.dy[ln] = (u[ry+1][rx] - u[ry-1][rx])/(2.0*hy);

        ui.dx[ln] += ((px-rx*hx)/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
        ui.dy[ln] += ((py-ry*hy)/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);

        //ui.dxx[ln] = (1.0/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
        //ui.dyy[ln] = (1.0/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);
    }
}

auto Problem2HNDirichlet4::b_add2Info(const DoubleMatrix &p, spif_vectorH &p_info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &cntrlDeltaGridList) const -> void
{
    unsigned int Nc = static_cast<unsigned int>(cntrlDeltaGridList.size());
    for (unsigned int i=0; i<Nc; i++)
    {
        SpacePointInfoH &pi =  p_info[i];
        const DeltaGrid2D &cdg = cntrlDeltaGridList[i];

        pi.vl[ln] = 0.0;
        for (unsigned int m=cdg.minY(); m<=cdg.maxY(); m++)
        {
            for (unsigned int n=cdg.minX(); n<=cdg.maxX(); n++)
            {
                pi.vl[ln] += p[m][n] * cdg.weight(n,m) * (hx*hy);
            }
        }

        double px = cdg.p().x;
        double py = cdg.p().y;
        unsigned int rx = cdg.rx();
        unsigned int ry = cdg.ry();

        pi.dx[ln] = (p[ry][rx+1] - p[ry][rx-1])/(2.0*hx);
        pi.dy[ln] = (p[ry+1][rx] - p[ry-1][rx])/(2.0*hy);

        pi.dx[ln] += ((px-rx*hx)/(hx*hx))*(p[ry][rx+1] - 2.0*p[ry][rx] + p[ry][rx-1]);
        pi.dy[ln] += ((py-ry*hy)/(hy*hy))*(p[ry+1][rx] - 2.0*p[ry][rx] + p[ry-1][rx]);

        //ui.dxx[ln] = (1.0/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
        //ui.dyy[ln] = (1.0/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);


    }
}

auto Problem2HNDirichlet4::b_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HNDirichlet4::b_initial2(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HNDirichlet4::b_boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double {
    return 0.0;
}

auto Problem2HNDirichlet4::b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const -> double
{
    return -2.0*mu(n,m)*(u[m][n]);
}

auto Problem2HNDirichlet4::b_layerInfo(const DoubleMatrix &p UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{}

//auto Problem2HNDirichlet4::newDistributeDeltaGaussPulse(const std::vector<SpacePoint> &thetas, std::vector<ExtendedSpacePointH> &extThetas, const Dimension &dimX, const Dimension &dimY) const -> void
//{
//    double hx = dimX.step();
//    double hy = dimY.step();

//    unsigned int Nx = static_cast<unsigned int> ( dimX.size() );
//    unsigned int Ny = static_cast<unsigned int> ( dimY.size() );
//    unsigned int Ns = static_cast<unsigned int> ( thetas.size() );

//    extThetas.clear();
//    extThetas.resize(Ns);

//    //    const SpacePoint &theta = thetas.at(0);
//    //    ExtendedSpacePointH &extTheta = extThetas.at(0);
//    //    extTheta.x = theta.x;
//    //    extTheta.y = theta.y;
//    //    extTheta.rx = static_cast<int> ( round(extTheta.x*Nx) );
//    //    extTheta.ry = static_cast<int> ( round(extTheta.y*Ny) );
//    //    extTheta.k = 0;
//    //    extTheta.minX = extTheta.rx - extTheta.k;
//    //    extTheta.maxX = extTheta.rx + extTheta.k;
//    //    extTheta.minY = extTheta.ry - extTheta.k;
//    //    extTheta.maxY = extTheta.ry + extTheta.k;
//    //    ExtendedSpacePointNodeH node;
//    //    node.nx = 50; node.x = 50*hx;
//    //    node.ny = 50; node.y = 50*hy;
//    //    node.w = 1.0/(hx*hy);
//    //    node.isCenter = true;
//    //    extTheta.nodes.push_back(node);
//    //    return;

//    const int k = 24;
//    const double sigmaX = 8.0*hx;
//    const double sigmaY = 8.0*hy;

//    for (unsigned int s=0; s<Ns; s++)
//    {
//        const SpacePoint &theta = thetas.at(s);
//        ExtendedSpacePointH &extTheta = extThetas.at(s);

//        extTheta.x = theta.x;
//        extTheta.y = theta.y;
//        extTheta.rx = static_cast<int>( round(extTheta.x*Nx) );
//        extTheta.ry = static_cast<int>( round(extTheta.y*Ny) );
//        extTheta.k = k;
//        extTheta.minX = extTheta.rx - extTheta.k;
//        extTheta.maxX = extTheta.rx + extTheta.k;
//        extTheta.minY = extTheta.ry - extTheta.k;
//        extTheta.maxY = extTheta.ry + extTheta.k;

//        double sumX = 0.0;
//        for (int n=extTheta.minX; n<=extTheta.maxX; n++) sumX += exp(-((n*hx-theta.x)*(n*hx-theta.x))/(2.0*sigmaX*sigmaX));
//        sumX *= hx;

//        double sumY = 0.0;
//        for (int m=extTheta.minY; m<=extTheta.maxY; m++) sumY += exp(-((m*hy-theta.y)*(m*hy-theta.y))/(2.0*sigmaY*sigmaY));
//        sumY *= hy;

//        double sigma = (sumX*sumY) / (2.0*M_PI);
//        double factor = 1.0/(2.0*M_PI*sigma);
//        //double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);

//        for (int m=extTheta.minY; m<=extTheta.maxY; m++)
//        {
//            for (int n=extTheta.minX; n<=extTheta.maxX; n++)
//            {
//                ExtendedSpacePointNodeH node;
//                node.nx = n; node.x = n*hx;
//                node.ny = m; node.y = m*hy;
//                node.w = factor*exp(-0.5*(((node.x-theta.x)*(node.x-theta.x))/(sigmaX*sigmaX)+((node.y-theta.y)*(node.y-theta.y))/(sigmaY*sigmaY)));
//                node.isCenter = ( m==extTheta.ry && n==extTheta.rx );
//                extTheta.nodes.push_back(node);
//            }
//        }
//    }

//    //    EquationParameterHE &prm = const_cast<EquationParameterHE&>(this->mParameter);
//    //    for (unsigned int s=0; s<Ns; s++)
//    //    {
//    //        SpacePoint &msmnt = prm.theta[s];
//    //        EquationParameterHE::SpacePointExt &spe = prm.theta_ext[s];

//    //        spe.rx = static_cast<unsigned int>(round(msmnt.x*Nx));
//    //        spe.ry = static_cast<unsigned int>(round(msmnt.y*Ny));
//    //        spe.minX = spe.rx - k;
//    //        spe.maxX = spe.rx + k;
//    //        spe.minY = spe.ry - k;
//    //        spe.maxY = spe.ry + k;

//    //        double sumX = 0.0;
//    //        for (unsigned int n=spe.minX; n<=spe.maxX; n++) sumX += exp(-((n*hx-msmnt.x)*(n*hx-msmnt.x))/(2.0*sigmaX*sigmaX));
//    //        sumX *= hx;

//    //        double sumY = 0.0;
//    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++) sumY += exp(-((m*hy-msmnt.y)*(m*hy-msmnt.y))/(2.0*sigmaY*sigmaY));
//    //        sumY *= hy;

//    //        double sigma = (sumX*sumY) / (2.0*M_PI);
//    //        double factor = 1.0/((2.0*M_PI)*sigma);

//    //        spe.nodes.clear();
//    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++)
//    //        {
//    //            for (unsigned int n=spe.minX; n<=spe.maxX; n++)
//    //            {
//    //                EquationParameterHE::SpacePointExt::Node node;
//    //                node.nx = n; node.x = n*hx;
//    //                node.ny = m; node.y = m*hy;
//    //                node.w = factor*exp(-0.5*(((node.x-msmnt.x)*(node.x-msmnt.x))/(sigmaX*sigmaX)+((node.y-msmnt.y)*(node.y-msmnt.y))/(sigmaY*sigmaY)));
//    //                node.isCenter = ( m==spe.ry && n==spe.rx );
//    //                spe.nodes.push_back(node);
//    //            }
//    //        }
//    //    }
//}

//auto Problem2HNDirichlet4::newDistributeDeltaGaussCntrl(const std::vector<SpacePoint> &cntrls, std::vector<ExtendedSpacePointH> &extCntrls, const Dimension &dimX, const Dimension &dimY) const -> void
//{
//    double hx = dimX.step();
//    double hy = dimY.step();

//    unsigned int Nx = static_cast<unsigned int> ( dimX.size() );
//    unsigned int Ny = static_cast<unsigned int> ( dimY.size() );
//    unsigned int Nc = static_cast<unsigned int> ( cntrls.size() );

//    extCntrls.clear();
//    extCntrls.resize(Nc);

//    int k = 3;
//    double sigmaX = hx;
//    double sigmaY = hy;

//    for (unsigned int c=0; c<Nc; c++)
//    {
//        const SpacePoint &cntrl = cntrls.at(c);
//        ExtendedSpacePointH &extCntrl = extCntrls.at(c);

//        extCntrl.x = cntrl.x;
//        extCntrl.y = cntrl.y;
//        extCntrl.rx = static_cast<int>(round(extCntrl.x*Nx));
//        extCntrl.ry = static_cast<int>(round(extCntrl.y*Ny));
//        extCntrl.k = k;
//        extCntrl.minX = extCntrl.rx - extCntrl.k;
//        extCntrl.maxX = extCntrl.rx + extCntrl.k;
//        extCntrl.minY = extCntrl.ry - extCntrl.k;
//        extCntrl.maxY = extCntrl.ry + extCntrl.k;

//        double sumX = 0.0;
//        for (int n=extCntrl.minX; n<=extCntrl.maxX; n++) sumX += exp(-((n*hx-cntrl.x)*(n*hx-cntrl.x))/(2.0*sigmaX*sigmaX));
//        sumX *= hx;

//        double sumY = 0.0;
//        for (int m=extCntrl.minY; m<=extCntrl.maxY; m++) sumY += exp(-((m*hy-cntrl.y)*(m*hy-cntrl.y))/(2.0*sigmaY*sigmaY));
//        sumY *= hy;

//        double sigma = (sumX*sumY) / (2.0*M_PI);
//        double factor = 1.0/((2.0*M_PI)*sigma);

//        for (int m=extCntrl.minY; m<=extCntrl.maxY; m++)
//        {
//            for (int n=extCntrl.minX; n<=extCntrl.maxX; n++)
//            {
//                ExtendedSpacePointNodeH node;
//                node.nx = n; node.x = n*hx;
//                node.ny = m; node.y = m*hy;
//                node.w = factor*exp(-0.5*(((node.x-cntrl.x)*(node.x-cntrl.x))/(sigmaX*sigmaX)+((node.y-cntrl.y)*(node.y-cntrl.y))/(sigmaY*sigmaY)));
//                node.isCenter = ( m==extCntrl.ry && n==extCntrl.rx );
//                extCntrl.nodes.push_back(node);
//            }
//        }
//    }

//    //    EquationParameterHE &prm = const_cast<EquationParameterHE&>(this->mParameter);
//    //    for (unsigned int i=0; i<Nc; i++)
//    //    {
//    //        SpacePoint &cntrl = prm.eta[i];
//    //        EquationParameterHE::SpacePointExt &spe = prm.eta_ext[i];

//    //        spe.rx = static_cast<unsigned int>( round(cntrl.x*Nx));
//    //        spe.ry = static_cast<unsigned int>( round(cntrl.y*Ny));
//    //        spe.minX = spe.rx - k;
//    //        spe.maxX = spe.rx + k;
//    //        spe.minY = spe.ry - k;
//    //        spe.maxY = spe.ry + k;

//    //        double sumX = 0.0;
//    //        for (unsigned int n=spe.minX; n<=spe.maxX; n++) sumX += exp(-((n*hx-cntrl.x)*(n*hx-cntrl.x))/(2.0*sigmaX*sigmaX));
//    //        sumX *= hx;

//    //        double sumY = 0.0;
//    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++) sumY += exp(-((m*hy-cntrl.y)*(m*hy-cntrl.y))/(2.0*sigmaY*sigmaY));
//    //        sumY *= hy;

//    //        double sigma = (sumX*sumY) / (2.0*M_PI);
//    //        double factor = 1.0/((2.0*M_PI)*sigma);

//    //        spe.nodes.clear();
//    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++)
//    //        {
//    //            for (unsigned int n=spe.minX; n<=spe.maxX; n++)
//    //            {
//    //                EquationParameterHE::SpacePointExt::Node node;
//    //                node.nx = n; node.x = n*hx;
//    //                node.ny = m; node.y = m*hy;
//    //                node.w = factor*exp(-0.5*(((node.x-cntrl.x)*(node.x-cntrl.x))/(sigmaX*sigmaX)+((node.y-cntrl.y)*(node.y-cntrl.y))/(sigmaY*sigmaY)));
//    //                node.isCenter = ( m==spe.ry && n==spe.rx );
//    //                spe.nodes.push_back(node);
//    //            }
//    //        }
//    //    }
//}

//auto Problem2HNDirichlet4::newDistributeDeltaGaussMsmnt(const std::vector<SpacePoint> &msmnts, std::vector<ExtendedSpacePointH> &extMsmnts, const Dimension &dimX, const Dimension &dimY) const -> void
//{
//    double hx = dimX.step();
//    double hy = dimY.step();

//    unsigned int Nx = static_cast<unsigned int> ( dimX.size() );
//    unsigned int Ny = static_cast<unsigned int> ( dimY.size() );
//    unsigned int No = static_cast<unsigned int> ( msmnts.size() );

//    extMsmnts.clear();
//    extMsmnts.resize(No);

//    int k = 3;
//    double sigmaX = hx;
//    double sigmaY = hy;

//    for (unsigned int c=0; c<No; c++)
//    {
//        const SpacePoint &msmnt = msmnts.at(c);
//        ExtendedSpacePointH &extMsmnt = extMsmnts.at(c);

//        extMsmnt.x = msmnt.x;
//        extMsmnt.y = msmnt.y;
//        extMsmnt.rx = static_cast<int>(round(extMsmnt.x*Nx));
//        extMsmnt.ry = static_cast<int>(round(extMsmnt.y*Ny));
//        extMsmnt.k = k;
//        extMsmnt.minX = extMsmnt.rx - extMsmnt.k;
//        extMsmnt.maxX = extMsmnt.rx + extMsmnt.k;
//        extMsmnt.minY = extMsmnt.ry - extMsmnt.k;
//        extMsmnt.maxY = extMsmnt.ry + extMsmnt.k;

//        double sumX = 0.0;
//        for (int n=extMsmnt.minX; n<=extMsmnt.maxX; n++) sumX += exp(-((n*hx-msmnt.x)*(n*hx-msmnt.x))/(2.0*sigmaX*sigmaX));
//        sumX *= hx;

//        double sumY = 0.0;
//        for (int m=extMsmnt.minY; m<=extMsmnt.maxY; m++) sumY += exp(-((m*hy-msmnt.y)*(m*hy-msmnt.y))/(2.0*sigmaY*sigmaY));
//        sumY *= hy;

//        double sigma = (sumX*sumY) / (2.0*M_PI);
//        double factor = 1.0/((2.0*M_PI)*sigma);

//        for (int m=extMsmnt.minY; m<=extMsmnt.maxY; m++)
//        {
//            for (int n=extMsmnt.minX; n<=extMsmnt.maxX; n++)
//            {
//                ExtendedSpacePointNodeH node;
//                node.nx = n; node.x = n*hx;
//                node.ny = m; node.y = m*hy;
//                node.w = factor*exp(-0.5*(((node.x-msmnt.x)*(node.x-msmnt.x))/(sigmaX*sigmaX)+((node.y-msmnt.y)*(node.y-msmnt.y))/(sigmaY*sigmaY)));
//                node.isCenter = ( m==extMsmnt.ry && n==extMsmnt.rx );
//                extMsmnt.nodes.push_back(node);
//            }
//        }
//    }

//    //    EquationParameterHE &prm = const_cast<EquationParameterHE&>(this->mParameter);
//    //    for (unsigned int j=0; j<No; j++)
//    //    {
//    //        SpacePoint &msmnt = prm.xi[j];
//    //        EquationParameterHE::SpacePointExt &spe = prm.xi_ext[j];

//    //        spe.rx = static_cast<unsigned int>(round(msmnt.x*Nx));
//    //        spe.ry = static_cast<unsigned int>(round(msmnt.y*Ny));
//    //        spe.minX = spe.rx - k;
//    //        spe.maxX = spe.rx + k;
//    //        spe.minY = spe.ry - k;
//    //        spe.maxY = spe.ry + k;

//    //        double sumX = 0.0;
//    //        for (unsigned int n=spe.minX; n<=spe.maxX; n++) sumX += exp(-((n*hx-msmnt.x)*(n*hx-msmnt.x))/(2.0*sigmaX*sigmaX));
//    //        sumX *= hx;

//    //        double sumY = 0.0;
//    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++) sumY += exp(-((m*hy-msmnt.y)*(m*hy-msmnt.y))/(2.0*sigmaY*sigmaY));
//    //        sumY *= hy;

//    //        double sigma = (sumX*sumY) / (2.0*M_PI);
//    //        double factor = 1.0/((2.0*M_PI)*sigma);

//    //        spe.nodes.clear();
//    //        for (unsigned int m=spe.minY; m<=spe.maxY; m++)
//    //        {
//    //            for (unsigned int n=spe.minX; n<=spe.maxX; n++)
//    //            {
//    //                EquationParameterHE::SpacePointExt::Node node;
//    //                node.nx = n; node.x = n*hx;
//    //                node.ny = m; node.y = m*hy;
//    //                node.w = factor*exp(-0.5*(((node.x-msmnt.x)*(node.x-msmnt.x))/(sigmaX*sigmaX)+((node.y-msmnt.y)*(node.y-msmnt.y))/(sigmaY*sigmaY)));
//    //                node.isCenter = ( m==spe.ry && n==spe.rx );
//    //                spe.nodes.push_back(node);
//    //            }
//    //        }
//    //    }
//}

//auto Problem2HNDirichlet4::distributeTimeDelta(double t, double ht, unsigned int ln, const SpaceNodePDE &sn, const std::vector<ExtendedSpacePointH> &xsps) const -> double
//{
//    return 0.0;

//    if ( ln >= 40 ) return 0.0;

//    double Q = 0.0;

//    const unsigned int Ns = static_cast<const unsigned int>(xsps.size());

//    for (unsigned int s=0; s<Ns; s++)
//    {
//        double q = mEquParameter.q[s];
//        const ExtendedSpacePointH &extendedSpacePoint = xsps.at(s);
//        if (extendedSpacePoint.contains(sn))
//        {
//            const std::vector<ExtendedSpacePointNodeH> &nodes = extendedSpacePoint.nodes;
//            const unsigned int nodes_size = static_cast<const unsigned int>(nodes.size());
//            //printf_s("%d\n", nodes_size);
//            for (unsigned int i=0; i<nodes_size; i++)
//            {
//                const ExtendedSpacePointNodeH &node = nodes.at(i);
//                if (node.equals(sn))
//                {
//                    Q += q * node.w;
//                }
//            }
//        }
//    }

//    const double sigma = 5.0*ht;
//    const double mu = 20.0*ht;
//    const double factor = 1.0 / (sqrt(2*M_PI)*sigma);
//    return factor * exp( -0.5*((t - mu)*(t - mu))/(sigma*sigma) ) * Q;
//}

auto Problem2HNDirichlet4::PrmToVector(const OptimizeParameterH &prm, DoubleVector &pv) const -> void
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

auto Problem2HNDirichlet4::VectorToPrm(const DoubleVector &pv, OptimizeParameterH &prm) const -> void
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

auto Problem2HNDirichlet4::v(unsigned int i, OptimizeParameterH o_prm, EquationParameterH e_prm, const spif_vectorH &u_info, unsigned int ln) const -> double
{
    const unsigned int No = static_cast<const unsigned int>(e_prm.No);
    double v = 0.0;
    for (unsigned int j=0; j<No; j++)
    {
        v += o_prm.k[i][j] * (u_info[j].vl[ln]-o_prm.z[i][j]);
    }
    return v;
}
