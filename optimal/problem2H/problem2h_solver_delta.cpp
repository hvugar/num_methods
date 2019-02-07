﻿#include "problem2h_solver_delta.h"

void Problem2HDirichletDelta::Main(int argc, char* argv[])
{
    QGuiApplication app(argc, argv);
    example1();
}

auto Problem2HDirichletDelta::example1() -> void
{
    // Equation parameters
    EquaParameter2H equaPrm;
    equaPrm.a = 1.0;
    equaPrm.alpha = +0.001;

    // Functional parameters
    FuncParameter2H funcPrm;
    funcPrm.Q1 << 0.2;// << 0.22 << 0.24;
    funcPrm.Q2 << 0.2;// << 0.27 << 0.29;

    equaPrm.Nt = 10;
    equaPrm.tm.resize(equaPrm.Nt);
    for (unsigned int s=0; s<equaPrm.Nt; s++)
    {
        equaPrm.tm[s].i = (s+1)*30;
        equaPrm.tm[s].t = (s+1)*0.3;
    }

    // Optimization parameters
    equaPrm.initParemeters(10, 2, 2);
//    for (unsigned int s=0; s<equaPrm.Nt; s++)
//    {
//        for (unsigned int r=0; r<equaPrm.Nc; r++)
//        {
//            for (unsigned int c=0; c<equaPrm.No; c++)
//            {
//                //o_prm.k[s][r][c] = -static_cast<double>((rand() % 1000))/1000.0;
//                //o_prm.k[s][r][c] = 1.0-static_cast<double>((rand() % 2000))/1000.0;
//                //o_prm.z[s][r][c] = +static_cast<double>((rand() % 1000))/100000.0;
//                //equaPrm.opt.k[s][r][c] = -0.1*fabs(sin((c+1)*10.0)*cos((r+1)*20.0)*sin((s+1)*0.1));
//                //equaPrm.opt.z[s][r][c] = 0.01*cos((c+1)*10.0)*sin((r+1)*20.0)*sin((s+1)*0.2);
//
//                //equaPrm.opt.k[s][r][c] = 0.500+(rand()%1000)*0.0001;
//                equaPrm.opt.k[s][r][c] = -0.010*sin(c+r+s+1.0);
//                equaPrm.opt.z[s][r][c] = +0.010*cos(c+r+s+1.0);
//                //printf("%d %d %d %f %f\n", s, r, c, -sin(c+r+s+1.0), -cos(c+r+s+1.0));
//            }
//        }
//    }
    equaPrm.initPulseParemeters(2);
    equaPrm.pulses[0] = InitialPulse2D(SpacePoint(0.2600, 0.2600), 0.2);
    equaPrm.pulses[1] = InitialPulse2D(SpacePoint(0.7300, 0.7200), 0.2);

    DoubleVector ox;
    ox << -0.008415 << -0.009093 << -0.009581 << -0.008411 << -0.002504 << -0.001411 << -0.002457 << -0.003568 << -0.002481 << -0.001568;
    ox << -0.007568 << -0.009259 << -0.006758 << -0.002589 << -0.005989 << -0.002794 << -0.002389 << -0.003974 << -0.004279 << -0.006570;
    ox <<  0.001794 << -0.004570 << -0.002476 << -0.006325 << -0.004125 << -0.002354 << -0.001286 << -0.002321 << -0.004694 << -0.001281;
    ox << -0.004281 <<  0.003440 << -0.001575 <<  0.005440 <<  0.001403 <<  0.000845 <<  0.006942 <<  0.004512 <<  0.003267 <<  0.004456;
    ox <<  0.005403 << -0.004161 << -0.001461 << -0.009245 << -0.008121 << -0.001997 << -0.001299 << -0.005663 << -0.008647 << -0.002357;
    ox << -0.004593 <<  0.002837 << -0.006536 <<  0.002975 <<  0.008237 <<  0.009602 <<  0.002237 <<  0.009512 <<  0.008612 <<  0.007539;
    ox <<  0.003602 <<  0.005529 <<  0.006529 << -0.001455 <<  0.002539 << -0.002425 << -0.003428 << -0.009751 << -0.004145 << -0.002851;
    ox << -0.005251 << -0.008391 << -0.009111 << -0.002839 << -0.002831 <<  0.005444 << -0.001839 <<  0.002364 <<  0.006544 <<  0.008439;
    ox <<  0.344500 <<  0.382500 <<  0.826700 <<  0.916400;
    ox <<  0.228200 <<  0.636700 <<  0.718400 <<  0.373600;

    DoubleVector rx;
    rx << -0.115452 << -0.111167 << -0.083476 << -0.075845 <<  0.015146 << -0.215321 << -0.017947 << -0.171416 <<  0.059772 <<  0.089042;
    rx <<  0.071801 <<  0.122163 <<  0.156964 <<  0.058127 <<  0.241945 <<  0.105073 << -0.073669 << -0.072246 << -0.028962 <<  0.025317;
    rx << -0.006582 <<  0.012226 << -0.054747 <<  0.089466 << -0.004514 <<  0.200153 << +0.008543 <<  0.230858 << -0.124258 << -0.033923;
    rx << -0.168729 << -0.065452 <<  0.127873 <<  0.025357 <<  0.063209 << -0.005781 << -0.044238 << -0.020863 << -0.037664 << -0.022886;
    rx << -0.008409 <<  0.016980 << -0.003193 << -0.010559 << -0.027209 <<  0.105051 << -0.000021 <<  0.103941 << -0.016936 << -0.021695;
    rx << -0.029625 << -0.048136 << -0.031988 << -0.009787 << -0.047384 << -0.033184 << -0.023529 << -0.038822 << -0.042849 << -0.042849;
    rx <<  0.005021 <<  0.004248 <<  0.009179 << -0.001834 <<  0.018598 << -0.157891 << -0.006210 << -0.147856 << -0.044797 << -0.012049;
    rx << -0.056247 << -0.025108 <<  0.018231 <<  0.002586 <<  0.040771 <<  0.013513 << -0.053478 << -0.028416 << -0.019315 << -0.008417;
    rx <<  0.350321 <<  0.361990 <<  0.798991 <<  0.870608;
    rx <<  0.344082 <<  0.627388 <<  0.687065 <<  0.348083;

    equaPrm.OptimalParameterFromVector(ox);
    equaPrm.RegularParameterFromVector(rx);

    // Penalty paramteres
    DoubleVector r; r << 1.0000 << 10.000 << 100.00;
    // Regularization coefficients
    DoubleVector e; e << 0.0000 << 0.0000 << 0.0000;

    DoubleVector e1; e1 << 0.1000 << 0.1000 << 0.1000;
    DoubleVector e2; e2 << 0.0100 << 0.0100 << 0.0100;

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem2HDirichletDelta prob;
        prob.setGridDimensions(Dimension(0.010, 0, 300), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
        prob.equaPrm = equaPrm;
        prob.funcPrm = funcPrm;
        prob.optimizeK = true;
        prob.optimizeZ = false;
        prob.optimizeO = false;
        prob.optimizeC = false;
        prob.funcPrm.vmin.resize(equaPrm.Nc, -0.05);
        prob.funcPrm.vmax.resize(equaPrm.Nc, +0.05);
        prob.LD = 30;
        prob.noise = 0.0;

        prob.funcPrm.regEpsilon = e[i];
        prob.funcPrm.r = r[i];

        if (i==0)
        {
            prob.equaPrm.OptimalParameterToVector(x);
            //prob.checkGradient1(prob);
            //return;

            //prob.checkGradient2(prob);
            //return;

            //prob.printLayers = false;
            //if (prob.printLayers)
            //{
            //    std::vector<DoubleVector> u;
            //    spif_vector1H u_info, p_info;
            //    prob.solveForwardIBVP(u, u_info, true);
            //    return;
            //}
            //IPrinter::printSeperatorLine();
            //prob.solveBackwardIBVP(u, p_info, true, u_info);
            //return;
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setProjection(new ProjectionEx1);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.00001);
        g.setFunctionTolerance(0.00001);
        g.setStepTolerance(0.00001);
        g.setR1MinimizeEpsilon(e1[i], e2[i]);
        g.setMaxIterations(50);
        g.setNormalize(true);
        g.showExitMessage(true);
        //prob.gm = &g;

        IPrinter::printSeperatorLine(nullptr, '*');
        printf("k : "); IPrinter::print(x.mid(00, 19), x.mid(00, 19).length(), 9, 6);
        printf("k : "); IPrinter::print(x.mid(20, 39), x.mid(20, 39).length(), 9, 6);
        printf("z : "); IPrinter::print(x.mid(40, 59), x.mid(40, 59).length(), 9, 6);
        printf("z : "); IPrinter::print(x.mid(60, 79), x.mid(60, 79).length(), 9, 6);
        printf("xy: "); IPrinter::print(x.mid(80, 87), x.mid(80, 87).length(), 9, 6);
        IPrinter::printSeperatorLine(nullptr, '*');

        g.calculate(x);

        IPrinter::printSeperatorLine(nullptr, '=');
        printf("k : "); IPrinter::print(x.mid(00, 19), x.mid(00, 19).length(), 9, 6);
        printf("k : "); IPrinter::print(x.mid(20, 39), x.mid(20, 39).length(), 9, 6);
        printf("z : "); IPrinter::print(x.mid(40, 59), x.mid(40, 59).length(), 9, 6);
        printf("z : "); IPrinter::print(x.mid(60, 79), x.mid(60, 79).length(), 9, 6);
        printf("xy: "); IPrinter::print(x.mid(80, 87), x.mid(80, 87).length(), 9, 6);
        IPrinter::printSeperatorLine(nullptr, '=');

        prob.equaPrm.OptimalParameterFromVector(x);
        prob.printLayers = true;
        if (prob.printLayers)
        {
            prob.setGridDimensions(Dimension(0.010, 0, 800), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
            std::vector<DoubleMatrix> u;
            spif_vectorH u_info, p_info;
            prob.solveForwardIBVP(u, u_info, true, x, 0.25);
        }
        puts("End");
    }
}

auto Problem2HDirichletDelta::checkGradient1(const Problem2HDirichletDelta &prob) -> void
{
    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prob.equaPrm.OptimalParameterToVector(pv);
    printf("ok: "); IPrinter::print(pv.mid(00, 19), pv.mid(00, 19).length(), 9, 6);
    printf("ok: "); IPrinter::print(pv.mid(20, 39), pv.mid(20, 39).length(), 9, 6);
    printf("oz: "); IPrinter::print(pv.mid(40, 59), pv.mid(40, 59).length(), 9, 6);
    printf("oz: "); IPrinter::print(pv.mid(60, 79), pv.mid(60, 79).length(), 9, 6);
    printf("xy: "); IPrinter::print(pv.mid(80, 87), pv.mid(80, 87).length(), 9, 6);
    IPrinter::printSeperatorLine();

    DoubleVector rv;
    prob.equaPrm.RegularParameterToVector(rv);
    printf("rk: "); IPrinter::print(rv.mid(00, 19), rv.mid(00, 19).length(), 9, 6);
    printf("rk: "); IPrinter::print(rv.mid(20, 39), rv.mid(20, 39).length(), 9, 6);
    printf("rz: "); IPrinter::print(rv.mid(40, 59), rv.mid(40, 59).length(), 9, 6);
    printf("rz: "); IPrinter::print(rv.mid(60, 79), rv.mid(60, 79).length(), 9, 6);
    printf("xy: "); IPrinter::print(rv.mid(80, 87), rv.mid(80, 87).length(), 9, 6);
    IPrinter::printSeperatorLine();

    DoubleVector ag(pv.length());
    double functional = prob.fx(pv);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    prob.gradient(pv, ag);
    puts("Gradients are calculated.");
    //return;

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    const unsigned int Nt = prob.equaPrm.Nt;
    const unsigned int Nc = prob.equaPrm.Nc;
    const unsigned int No = prob.equaPrm.No;
    const unsigned int offset = Nc*No*Nt;

    puts("Calculating numerical gradients.... dh=0.01");
    printf("*** Calculating numerical gradients for k...... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 0*offset, 1*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for z...... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 1*offset, 2*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for xi..... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*offset+0*No, 2*offset+2*No-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for eta.... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*offset+2*No, 2*offset+2*No+2*Nc-1);
    printf("Calculated.\n");

    puts("Calculating numerical gradients.... dh=0.001");
    printf("*** Calculating numerical gradients for k...... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 0*offset, 1*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for z...... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 1*offset, 2*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for xi..... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*offset+0*No, 2*offset+2*No-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for eta.... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*offset+2*No, 2*offset+2*No+2*Nc-1);
    printf("Calculated.\n");

    const unsigned int N = 20;
    const unsigned int W = 9;
    const unsigned int P = 6;
    //k------------------------------------------------------//
    IPrinter::printSeperatorLine("k");
    DoubleVector pk0 = pv.mid(0, offset-1);
    DoubleVector ak0 = ag.mid(0, offset-1);
    DoubleVector nk1 = ng1.mid(0, offset-1);
    DoubleVector nk2 = ng2.mid(0, offset-1);

    IPrinter::printVector(W,P,pk0,nullptr,N);
    IPrinter::printVector(W,P,ak0,nullptr,N); ak0.L2Normalize();
    IPrinter::printVector(W,P,nk1,nullptr,N); nk1.L2Normalize();
    IPrinter::printVector(W,P,nk2,nullptr,N); nk2.L2Normalize();
    IPrinter::printVector(W,P,ak0,nullptr,N);
    IPrinter::printVector(W,P,nk1,nullptr,N);
    IPrinter::printVector(W,P,nk2,nullptr,N);

    //z------------------------------------------------------//
    IPrinter::printSeperatorLine("z");
    DoubleVector pz0 = pv.mid(offset, 2*offset-1);
    DoubleVector az0 = ag.mid(offset, 2*offset-1);
    DoubleVector nz1 = ng1.mid(offset, 2*offset-1);
    DoubleVector nz2 = ng2.mid(offset, 2*offset-1);

    IPrinter::printVector(W,P,pz0,nullptr,N);
    IPrinter::printVector(W,P,az0,nullptr,N); az0.L2Normalize();
    IPrinter::printVector(W,P,nz1,nullptr,N); nz1.L2Normalize();
    IPrinter::printVector(W,P,nz2,nullptr,N); nz2.L2Normalize();
    IPrinter::printVector(W,P,az0,nullptr,N);
    IPrinter::printVector(W,P,nz1,nullptr,N);
    IPrinter::printVector(W,P,nz2,nullptr,N);

    //xi------------------------------------------------------//
    IPrinter::printSeperatorLine("xi");
    DoubleVector pe0 = pv.mid(2*offset, 2*offset+2*prob.equaPrm.No-1);
    DoubleVector ae0 = ag.mid(2*offset, 2*offset+2*prob.equaPrm.No-1);
    DoubleVector ne1 = ng1.mid(2*offset, 2*offset+2*prob.equaPrm.No-1);
    DoubleVector ne2 = ng2.mid(2*offset, 2*offset+2*prob.equaPrm.No-1);

    IPrinter::print(pe0,pe0.length(),14,4);
    IPrinter::print(ae0,ae0.length(),14,4); ae0.L2Normalize();
    IPrinter::print(ne1,ne1.length(),14,4); ne1.L2Normalize();
    IPrinter::print(ne2,ne2.length(),14,4); ne2.L2Normalize();
    IPrinter::print(ae0,ae0.length(),14,4);
    IPrinter::print(ne1,ne1.length(),14,4);
    IPrinter::print(ne2,ne2.length(),14,4);

    //eta------------------------------------------------------//
    IPrinter::printSeperatorLine("eta");
    DoubleVector px0 = pv.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
    DoubleVector ax0 = ag.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
    DoubleVector nx1 = ng1.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
    DoubleVector nx2 = ng2.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);

    IPrinter::print(px0,px0.length(),14,4);
    IPrinter::print(ax0,ax0.length(),14,4); ax0.L2Normalize();
    IPrinter::print(nx1,nx1.length(),14,4); nx1.L2Normalize();
    IPrinter::print(nx2,nx2.length(),14,4); nx2.L2Normalize();
    IPrinter::print(ax0,ax0.length(),14,4);
    IPrinter::print(nx1,nx1.length(),14,4);
    IPrinter::print(nx2,nx2.length(),14,4);
    IPrinter::printSeperatorLine();
}

auto Problem2HDirichletDelta::checkGradient2(const Problem2HDirichletDelta &prob) -> void
{
    /*
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
    */
}

Problem2HDirichletDelta::Problem2HDirichletDelta() {}

Problem2HDirichletDelta::~Problem2HDirichletDelta() {}

auto Problem2HDirichletDelta::fx(const DoubleVector &pv) const -> double
{
    Problem2HDirichletDelta* prob = const_cast<Problem2HDirichletDelta*>(this);
    prob->equaPrm.OptimalParameterFromVector(pv);

    const DoubleVector &Q1 = funcPrm.Q1;
    const DoubleVector &Q2 = funcPrm.Q2;
    std::vector<InitialPulse2D> &pulses = prob->equaPrm.pulses;
    const double regEpsilon = funcPrm.regEpsilon;
    const double r = funcPrm.r;

    double SUM = 0.0;
    for (unsigned int q1=0; q1<Q1.length(); q1++)
    {
        pulses[0].q = Q1[q1];
        for (unsigned int q2=0; q2<Q2.length(); q2++)
        {
            pulses[1].q = Q2[q2];

            std::vector<DoubleMatrix> u;
            spif_vectorH u_info;
            prob->solveForwardIBVP(u, u_info, true, pv);

            double intgrl = integral(u);
            double pnt = penalty(u_info, equaPrm);
            double nrm = norm(equaPrm);
            double sum = intgrl + regEpsilon*nrm + r*pnt;

            for (unsigned int i=0; i<u.size(); i++)      u[i].clear();      u.clear();
            for (unsigned int j=0; j<u_info.size(); j++) u_info[j].clear(); u_info.clear();

            SUM += sum * (1.0/(double(Q1.length())*double(Q2.length())));
        }
    }
    return SUM;
}

auto Problem2HDirichletDelta::integral(const std::vector<DoubleMatrix> &vu) const -> double
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

auto Problem2HDirichletDelta::integralU(const DoubleMatrix &u) const -> double
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

auto Problem2HDirichletDelta::norm(const EquaParameter2H &prm) const -> double
{
    double _norm = 0.0;
    const unsigned int Nc = prm.Nc;
    const unsigned int No = prm.No;
    const unsigned int Nt = prm.Nt;

    for (unsigned int s=0; s<Nt; s++)
    {
        const DoubleMatrix &ok = prm.opt.k[s];
        const DoubleMatrix &rk = prm.reg.k[s];
        const DoubleMatrix &oz = prm.opt.z[s];
        const DoubleMatrix &rz = prm.reg.z[s];

        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                _norm += (ok[i][j] - rk[i][j])*(ok[i][j] - rk[i][j]);
                _norm += (oz[i][j] - rz[i][j])*(oz[i][j] - rz[i][j]);
            }
        }
    }

    for (unsigned int j=0; j<No; j++)
    {
        const SpacePoint &oksi = prm.opt.ksi[j];
        const SpacePoint &rksi = prm.reg.ksi[j];

        _norm += (oksi.x - rksi.x)*(oksi.x - rksi.x);
        _norm += (oksi.y - rksi.y)*(oksi.y - rksi.y);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        const SpacePoint &oeta = prm.opt.eta[i];
        const SpacePoint &reta = prm.reg.eta[i];

        _norm += (oeta.x - reta.x)*(oeta.x - reta.x);
        _norm += (oeta.y - reta.y)*(oeta.y - reta.y);
    }

    return _norm;
}

auto Problem2HDirichletDelta::penalty(const spif_vectorH &info, const EquaParameter2H &prm) const -> double
{
    double pnlt = 0.0;

    const unsigned int Nt = prm.Nt;
    const unsigned int Nc = prm.Nc;
    //const unsigned int ln = mEquParameter.timeMoments[s].i;

    for (unsigned int s=0; s<Nt; s++)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            double _gpi_s = gpi(i, s, info, prm);
            pnlt += _gpi_s*_gpi_s;
        }
    }

    return pnlt;
}

auto Problem2HDirichletDelta::gpi(unsigned int i, unsigned int s, const spif_vectorH &u_info, const EquaParameter2H &prm) const -> double
{
    const Dimension time = timeDimension();
    const unsigned int L = static_cast<unsigned int>(time.size());
    const unsigned int LLD = L + LD;
    const unsigned int ln = prm.tm[s].i;

    if (ln>LLD) return 0.0;

    const DoubleVector &vmin = funcPrm.vmin;
    const DoubleVector &vmax = funcPrm.vmax;

    double gpi_ln = fabs( g0i(i, s, u_info, prm) ) - ( vmax.at(i) - vmin.at(i) )/2.0;
    return gpi_ln > 0.0 ? gpi_ln : 0.0;
}

auto Problem2HDirichletDelta::g0i(unsigned int i, unsigned int s, const spif_vectorH &u_info, const EquaParameter2H &prm) const -> double
{
    const unsigned int No = prm.No;
    const unsigned int ln = 2*prm.tm[s].i;
    const DoubleVector &vmin = funcPrm.vmin;
    const DoubleVector &vmax = funcPrm.vmax;

    double vi = 0.0;
    for (unsigned int j=0; j<No; j++)
    {
        const DoubleMatrix &ok = prm.opt.k[s];
        const DoubleMatrix &oz = prm.opt.z[s];

        const SpacePointInfoH &u_xij = u_info[j];
        vi += ok[i][j] * ( u_xij.vl[ln] - oz[i][j] );
    }

    return ( vmax.at(i) + vmin.at(i) )/2.0 - vi;
}

auto Problem2HDirichletDelta::gradient(const DoubleVector &pv, DoubleVector &g) const -> void
{
    g.clear();
    g.resize(pv.length(), 0.0);

    const unsigned int Nc = equaPrm.Nc;
    const unsigned int No = equaPrm.No;
    const unsigned int Nt = equaPrm.Nt;
    const double regEpsilon = funcPrm.regEpsilon;
    const double r = funcPrm.r;

    Problem2HDirichletDelta* prob = const_cast<Problem2HDirichletDelta*>(this);
    prob->equaPrm.OptimalParameterFromVector(pv);
    std::vector<InitialPulse2D> &pulses = prob->equaPrm.pulses;

    const DoubleVector &Q1 = funcPrm.Q1;
    const DoubleVector &Q2 = funcPrm.Q2;

    for (unsigned int q1=0; q1<Q1.length(); q1++)
    {
        pulses[0].q = Q1[q1];
        for (unsigned int q2=0; q2<Q2.length(); q2++)
        {
            pulses[1].q = Q2[q2];

            std::vector<DoubleMatrix> u;
            spif_vectorH u_info;
            solveForwardIBVP(u, u_info, true, pv);
            spif_vectorH p_info;
            solveBackwardIBVP(u, p_info, true, u_info, pv);

            unsigned int gi = 0;

            // k
            if (optimizeK)
            {
                //puts("Calculating k gradients...");
                for (unsigned int s=0; s<Nt; s++)
                {
                    const unsigned int ln = 2*equaPrm.tm[s].i;
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const SpacePointInfoH &pi = p_info[i];
                        for (unsigned int j=0; j<No; j++)
                        {
                            const SpacePointInfoH &uj = u_info[j];
                            double zij = equaPrm.opt.z[s][i][j];
                            double grad_Kij = 0.0;
                            grad_Kij += -(uj.vl[ln] - zij) * pi.vl[ln];
                            grad_Kij += -(uj.vl[ln] - zij) * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
                            grad_Kij += +2.0*regEpsilon*(equaPrm.opt.k[s][i][j] - equaPrm.reg.k[s][i][j]);
                            g[gi++] += grad_Kij * (1.0/(double(Q1.length())*double(Q2.length())));
                        }
                    }
                }
            }
            else
            {
                for (unsigned int s=0; s<Nt; s++)
                {
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        for (unsigned int j=0; j<No; j++)
                        {
                            g[gi++] = 0.0;
                        }
                    }
                }
            }

            // z
            if (optimizeZ)
            {
                //puts("Calculating z gradients...");
                for (unsigned int s=0; s<Nt; s++)
                {
                    const unsigned int ln = 2*equaPrm.tm[s].i;
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const SpacePointInfoH &pi = p_info[i];
                        for (unsigned int j=0; j<No; j++)
                        {
                            double kij = equaPrm.opt.k[s][i][j];
                            double grad_Zij = 0.0;
                            grad_Zij += kij * pi.vl[ln];
                            grad_Zij += kij * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
                            grad_Zij += +2.0*regEpsilon*(equaPrm.opt.z[s][i][j] - equaPrm.opt.z[s][i][j]);
                            g[gi++] += grad_Zij * (1.0/(double(Q1.length())*double(Q2.length())));
                        }
                    }
                }
            }
            else
            {
                for (unsigned int s=0; s<Nt; s++)
                {
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        for (unsigned int j=0; j<No; j++)
                        {
                            g[gi++] = 0.0;
                        }
                    }
                }
            }

            // xi
            if (optimizeO)
            {
                //puts("Calculating o gradients...");
                for (unsigned int j=0; j<No; j++)
                {
                    const SpacePointInfoH &uj = u_info[j];

                    double gradXijX = 0.0;
                    double gradXijY = 0.0;

                    for (unsigned int s=0; s<Nt; s++)
                    {
                        const unsigned int ln = 2*equaPrm.tm[s].i;
                        for (unsigned int i=0; i<Nc; i++)
                        {
                            gradXijX += -equaPrm.opt.k[s][i][j] * uj.dx[ln] * p_info[i].vl[ln];
                            gradXijY += -equaPrm.opt.k[s][i][j] * uj.dy[ln] * p_info[i].vl[ln];
                            gradXijX += -equaPrm.opt.k[s][i][j] * uj.dx[ln] * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
                            gradXijY += -equaPrm.opt.k[s][i][j] * uj.dy[ln] * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
                        }
                    }

                    gradXijX += 2.0*regEpsilon*(equaPrm.opt.ksi[j].x - equaPrm.reg.ksi[j].x);
                    gradXijY += 2.0*regEpsilon*(equaPrm.opt.ksi[j].y - equaPrm.reg.ksi[j].y);

                    g[gi++] += gradXijX * (1.0/(double(Q1.length())*double(Q2.length())));
                    g[gi++] += gradXijY * (1.0/(double(Q1.length())*double(Q2.length())));
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
                //puts("Calculating c gradients...");
                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfoH &pi = p_info[i];

                    double gradEtaiX = 0.0;
                    double gradEtaiY = 0.0;

                    for (unsigned int s=0; s<Nt; s++)
                    {
                        const unsigned int ln = 2*equaPrm.tm[s].i;
                        for (unsigned int j=0; j<No; j++)
                        {
                            gradEtaiX += -pi.dx[ln] * equaPrm.opt.k[s][i][j] * (u_info[j].vl[ln] - equaPrm.opt.z[s][i][j]);
                            gradEtaiY += -pi.dy[ln] * equaPrm.opt.k[s][i][j] * (u_info[j].vl[ln] - equaPrm.opt.z[s][i][j]);
                        }
                    }

                    gradEtaiX += 2.0*regEpsilon*(equaPrm.opt.eta[i].x - equaPrm.reg.eta[i].x);
                    gradEtaiY += 2.0*regEpsilon*(equaPrm.opt.eta[i].y - equaPrm.reg.eta[i].y);

                    g[gi++] += gradEtaiX * (1.0/(double(Q1.length())*double(Q2.length())));
                    g[gi++] += gradEtaiY * (1.0/(double(Q1.length())*double(Q2.length())));
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

            for (unsigned int i=0; i<u_info.size(); i++) u_info[i].clear(); u_info.clear();
            for (unsigned int i=0; i<p_info.size(); i++) p_info[i].clear(); p_info.clear();
            for (unsigned int i=0; i<u.size(); i++) u[i].clear(); u.clear();
        }
    }
}

auto Problem2HDirichletDelta::solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use, const DoubleVector &pv, double lambda) const -> void
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

    const double a        = equaPrm.a;
    const double alpha    = equaPrm.alpha;
    const unsigned int No = equaPrm.No;
    const unsigned int Nc = equaPrm.Nc;

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
        measuremntGirdList[j].distributeGauss(equaPrm.opt.ksi[j], MSRMT_SIGMA, MSRMT_SIGMA);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        cntrlDeltaGridList[i].initGrid(N, hx, M, hy);
        cntrlDeltaGridList[i].distributeGauss(equaPrm.opt.eta[i], CNTRL_SIGMA, CNTRL_SIGMA);
    }
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) prepareInfo(No,  equaPrm.opt.ksi, u_info, 2*LLD+1);
    //----------------------------------------------------------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    initPulseWeightMatrix(equaPrm.pulses);
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
        //currentLayerFGrid(u10, cntrlDeltaGridList, measuremntGirdList, 2*ln-1);
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
        //currentLayerFGrid(u15, cntrlDeltaGridList, measuremntGirdList, 2*ln-0);
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

auto Problem2HDirichletDelta::solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info, const DoubleVector &pv, double lambda) const -> void
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

    const double a        = equaPrm.a;
    const double alpha    = equaPrm.alpha;
    const unsigned int No = equaPrm.No;
    const unsigned int Nc = equaPrm.Nc;

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

    //--------------------------------------------------------------------------------------------//
    std::vector<DeltaGrid2D> measuremntGirdList(No);
    std::vector<DeltaGrid2D> cntrlDeltaGridList(Nc);
    for (unsigned int j=0; j<No; j++)
    {
        measuremntGirdList[j].initGrid(N, hx, M, hy);
        measuremntGirdList[j].distributeGauss(equaPrm.opt.ksi[j], MSRMT_SIGMA, MSRMT_SIGMA);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        cntrlDeltaGridList[i].initGrid(N, hx, M, hy);
        cntrlDeltaGridList[i].distributeGauss(equaPrm.opt.eta[i], CNTRL_SIGMA, CNTRL_SIGMA);
    }
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) prepareInfo(Nc, equaPrm.opt.eta, p_info, 2*LLD+1);
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
        //currentLayerBGrid(p10, cntrlDeltaGridList, measuremntGirdList, 2*ln+1, u_info);
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
        //currentLayerBGrid(p10, cntrlDeltaGridList, measuremntGirdList, 2*ln+0, u_info);
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

auto Problem2HDirichletDelta::initPulseWeightMatrix(const std::vector<InitialPulse2D> &pulses) const -> void
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<const unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<const unsigned int> ( dimY.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();


    DoubleMatrix &pulseWeightMatrix = const_cast<Problem2HDirichletDelta*>(this)->mPulseWeightMatrix;

    const unsigned int Ns = static_cast<unsigned int>(pulses.size());
    DeltaGrid2D *deltaGrids = new DeltaGrid2D[Ns];
    for (unsigned int s=0; s<Ns; s++)
    {
        deltaGrids[s].initGrid(N, hx, M, hy);
        deltaGrids[s].distributeGauss(pulses[s].theta, PULSE_SIGMA, PULSE_SIGMA);
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            pulseWeightMatrix[m][n] = 0.0;
            for (unsigned int s=0; s<Ns; s++)
            {
                pulseWeightMatrix[m][n] += pulses[s].q*deltaGrids[s].weight(n, m);
            }
        }
    }

    for (unsigned int s=0; s<Ns; s++) deltaGrids[s].cleanGrid();
    delete [] deltaGrids;
}

auto Problem2HDirichletDelta::f_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
    //unsigned int m = static_cast<unsigned int>(sn.j);
    //unsigned int n = static_cast<unsigned int>(sn.i);
    //return mPulseWeightMatrix[m][n];
}

auto Problem2HDirichletDelta::f_initial2(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    unsigned int m = static_cast<unsigned int>(sn.j);
    unsigned int n = static_cast<unsigned int>(sn.i);
    return mPulseWeightMatrix[m][n];
    //return 0.0;
}

auto Problem2HDirichletDelta::f_boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichletDelta::b_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichletDelta::b_initial2(const SpaceNodePDE &sn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichletDelta::b_boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const -> double
{
    return 0.0;
}

auto Problem2HDirichletDelta::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const -> void
{
    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(f); C_UNUSED(alpha); C_UNUSED(result);
    const char* msg = nullptr; C_UNUSED(msg);
    if (result == GradientMethod::MethodResult::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION   ";
    if (result == GradientMethod::MethodResult::FIRST_ITERATION)          msg = "FIRST_ITERATION         ";
    if (result == GradientMethod::MethodResult::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    if (result == GradientMethod::MethodResult::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS     ";
    if (result == GradientMethod::MethodResult::NEXT_ITERATION)           msg = "NEXT_ITERATION          ";

    Problem2HDirichletDelta* prob = const_cast<Problem2HDirichletDelta*>(this);
    prob->equaPrm.OptimalParameterFromVector(x);
    std::vector<DoubleMatrix> u;
    spif_vectorH u_info;
    solveForwardIBVP(u, u_info, true, x);
    double ing = integral(u);
    double pnt = penalty(u_info, prob->equaPrm);
    double nrm = norm(prob->equaPrm);

    //    DoubleVector uf, um, ux;
    //    for (unsigned int i=0; i<=50; i+=5)
    //    {
    //        uf << integralU(u[i]);
    //        um << u[i].min();
    //        ux << u[i].max();
    //    }

    //const unsigned int v_length = static_cast<const unsigned int>(timeDimension().size()) + LD;
    //DoubleVector v1(v_length+1);
    //DoubleVector v2(v_length+1);

    //const unsigned int Nt = static_cast<unsigned int>(equaPrm.Nt);
    //DoubleVector v1(Nt);
    //DoubleVector v2(Nt);
    //for (unsigned int s=0; s<Nt; s++)
    //{
    //    v1[s] = v(0, s, equaPrm, u_info);
    //    v2[s] = v(1, s, equaPrm, u_info);
    //}
    //IPrinter::printVector(v1, "v1", v1.length());
    //IPrinter::printVector(v2, "v2", v2.length());

    //    IPrinter::printVector(uf, "uf", 10);
    //    IPrinter::printVector(um, "um", 10);
    //    IPrinter::printVector(ux, "ux", 10);

    printf("I[%3d]: F:%.6f I:%.6f P:%.6f N:%.5f R:%.3f e:%.3f a:%10.6f ", i, f, ing, pnt, nrm, funcPrm.r, funcPrm.regEpsilon, alpha);
    printf("min:%10.6f max:%10.6f U:%.8f ", u.front().min(), u.front().max(), integralU(u.front()));
    printf("min:%10.6f max:%10.6f U:%.8f ", u.back().min(), u.back().max(), integralU(u.back()));
    //printf("\n");
    //printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    unsigned int s = 2*(equaPrm.Nc*equaPrm.No*equaPrm.Nt);
    printf("o: %6.4f %6.4f %6.4f %6.4f c: %6.4f %6.4f %6.4f %6.4f\n", x[s+0], x[s+1], x[s+2], x[s+3], x[s+4], x[s+5], x[s+6], x[s+7]);
    //printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
    //DoubleVector n = g; n.L2Normalize();
    //printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o: %8.4f %8.4f %8.4f %8.4f c: %8.4f %8.4f %8.4f %8.4f\n", n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15]);

    //IPrinter::printSeperatorLine("Vector", '-');
    //printf("k : "); IPrinter::print(x.mid(00, 19), x.mid(00, 19).length(), 9, 6);
    //printf("k : "); IPrinter::print(x.mid(20, 39), x.mid(20, 39).length(), 9, 6);
    //printf("z : "); IPrinter::print(x.mid(40, 59), x.mid(40, 59).length(), 9, 6);
    //printf("z : "); IPrinter::print(x.mid(60, 79), x.mid(60, 79).length(), 9, 6);
    //printf("xy: "); IPrinter::print(x.mid(80, 87), x.mid(80, 87).length(), 9, 6);
    //IPrinter::printSeperatorLine("Gradient", '-');
    //printf("k : "); IPrinter::print(g.mid(00, 19), g.mid(00, 19).length(), 9, 6);
    //printf("k : "); IPrinter::print(g.mid(20, 39), g.mid(20, 39).length(), 9, 6);
    //printf("z : "); IPrinter::print(g.mid(40, 59), g.mid(40, 59).length(), 9, 6);
    //printf("z : "); IPrinter::print(g.mid(60, 79), g.mid(60, 79).length(), 9, 6);
    //printf("xy: "); IPrinter::print(g.mid(80, 87), g.mid(80, 87).length(), 9, 6);
    //IPrinter::printSeperatorLine(nullptr, '-');

    //u.clear();
    //u_info.clear();

    //C_UNUSED(prob);
    //IPrinter::printSeperatorLine();

//    prob->optimizeK = (i%4 == 3);
//    prob->optimizeZ = (i%4 == 0);
//    prob->optimizeO = (i%4 == 1);
//    prob->optimizeC = (i%4 == 2);

    prob->optimizeK = (i%3 == 2);
    prob->optimizeZ = (i%3 == 0);
    prob->optimizeO = (i%3 == 1);
    prob->optimizeC = (i%3 == 1);
}

auto Problem2HDirichletDelta::currentLayerFGrid(const DoubleMatrix &u, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids, unsigned int ln) const -> void
{
    Problem2HDirichletDelta* const_this = const_cast<Problem2HDirichletDelta*>(this);
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();
    const double ht = time.step();

    const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimY.size() );
    const unsigned int Nc = equaPrm.Nc;
    const unsigned int No = equaPrm.No;
    const unsigned int Nt = equaPrm.Nt;

    for (unsigned int m=0; m<=M; m++) for (unsigned int n=0; n<=N; n++) const_this->mCrFfxWeightMatrix[m][n] = 0.0;

    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        const unsigned int sln = equaPrm.tm[s].i;
        //if ((ln == 2*sln+0) || (ln == 2*sln+1)) { wt = 1.0/ht; } else { continue; }
        //if ((ln == 2*sln-2) || (ln == 2*sln-1)) { wt = 1.0/ht; } else { continue; }
        //if ((ln == 2*sln+0) || (ln == 2*sln-1)) { wt = 1.0/ht; } else { continue; }
        if (ln == 2*sln) { wt = 2.0/ht; } else { continue; }
        //if (ln == 2*sln-1) { wt = 2.0/ht; } else { continue; }
        //printf("ln: %d %d\n", ln, sln);

        double* _u = new double[No];
        for (unsigned int j=0; j<No; j++)
        {
            _u[j] = measurementDeltaGrids[j].consentrateInPoint(u);
            //_u[j] *= (1.0 + noise);
        }

        double *_v = new double[Nc];
        for (unsigned int i=0; i<Nc; i++)
        {
            _v[i] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                _v[i] += equaPrm.opt.k[s][i][j] * (_u[j] - equaPrm.opt.z[s][i][j]);
            }
        }
        delete [] _u;

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                for (unsigned int i=0; i<Nc; i++)
                {
                    const DeltaGrid2D &dg = controlDeltaGrids[i];
                    const_this->mCrFfxWeightMatrix[m][n] += _v[i] * dg.weight(n,m) * wt;
                }
            }
        }
        delete [] _v;
    }
}

auto Problem2HDirichletDelta::currentLayerBGrid(const DoubleMatrix &p, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids, unsigned int ln, const spif_vectorH &u_info) const -> void
{
    Problem2HDirichletDelta* const_this = const_cast<Problem2HDirichletDelta*>(this);
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();
    const double ht = time.step();
    const double r = funcPrm.r;

    const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimY.size() );
    const unsigned int Nc = equaPrm.Nc;
    const unsigned int No = equaPrm.No;
    const unsigned int Nt = equaPrm.Nt;

    for (unsigned int m=0; m<=M; m++) for (unsigned int n=0; n<=N; n++) const_this->mCrBfxWeightMatrix[m][n] = 0.0;

    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        const unsigned int sln = equaPrm.tm[s].i;
        //if ((ln == 2*sln-0) || (ln == 2*sln-1)) { wt = 1.0/ht; } else { continue; }
        //if ((ln == 2*sln+2) || (ln == 2*sln+1)) { wt = 1.0/ht; } else { continue; }
        //if ((ln == 2*sln-0) || (ln == 2*sln+1)) { wt = 1.0/ht; } else { continue; }
        if (ln == 2*sln) { wt = 2.0/ht; } else { continue; }
        //if (ln == 2*sln+1) { wt = 2.0/ht; } else { continue; }
        //printf("ln: %d %d\n", ln, sln);

        double* _p = new double[Nc];
        for (unsigned int i=0; i<Nc; i++)
        {
            _p[i] = controlDeltaGrids[i].consentrateInPoint(p);
        }

        double *_w = new double[No];
        for (unsigned int j=0; j<No; j++)
        {
            _w[j] = 0.0;
            for (unsigned int i=0; i<Nc; i++)
            {
                _w[j] += equaPrm.opt.k[s][i][j] * _p[i];
                _w[j] += equaPrm.opt.k[s][i][j] * 2.0*r*gpi(i, s, u_info, equaPrm)*sgn(g0i(i, s, u_info, equaPrm));
            }
        }
        delete [] _p;

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                for (unsigned int j=0; j<No; j++)
                {
                    const DeltaGrid2D &dg = measurementDeltaGrids[j];
                    const_this->mCrBfxWeightMatrix[m][n] += _w[j] * dg.weight(n,m) * wt;
                }
            }
        }
        delete [] _w;
    }
}

auto Problem2HDirichletDelta::prepareInfo(unsigned int N, const std::vector<SpacePoint> &points, spif_vectorH &info, unsigned int size) const -> void
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

auto Problem2HDirichletDelta::add2Info(const DoubleMatrix &u, spif_vectorH &info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &deltaList) const -> void
{
    const unsigned int N = static_cast<unsigned int>(deltaList.size());

    for (unsigned int i=0; i<N; i++)
    {
        const DeltaGrid2D &deltagrid = deltaList[i];
        SpacePointInfoH &ui = info[i];
        ui.vl[ln] = deltagrid.consentrateInPoint(u, ui.dx[ln], ui.dy[ln]);
    }
}

auto Problem2HDirichletDelta::project(DoubleVector &, unsigned int) -> void {}

auto Problem2HDirichletDelta::project(DoubleVector &pv) const -> void
{
    unsigned int Nc = equaPrm.Nc;
    unsigned int No = equaPrm.No;
    unsigned int Nt = equaPrm.Nt;

    unsigned int start = 2*Nc*No*Nt;
    unsigned int end = 2*Nc*No*Nt + 2*No + 2*Nc - 1;

    for (unsigned int index = start; index <= end; index++)
    {
        if (pv[index] <= 0.05) pv[index] = 0.05;
        if (pv[index] >= 0.95) pv[index] = 0.95;
    }

    unsigned int s = start;
    if (s+0) { if (pv[s+0] < 0.05) pv[s+0] = 0.05; if (pv[s+0] > 0.45) pv[s+0] = 0.45; }
    if (s+1) { if (pv[s+1] < 0.05) pv[s+1] = 0.05; if (pv[s+1] > 0.45) pv[s+1] = 0.45; }
    if (s+2) { if (pv[s+2] < 0.55) pv[s+2] = 0.55; if (pv[s+2] > 0.95) pv[s+2] = 0.95; }
    if (s+3) { if (pv[s+3] < 0.55) pv[s+3] = 0.55; if (pv[s+3] > 0.95) pv[s+3] = 0.95; }

    if (s+4) { if (pv[s+4] < 0.05) pv[s+4] = 0.05; if (pv[s+4] > 0.45) pv[s+4] = 0.45; }
    if (s+5) { if (pv[s+5] < 0.55) pv[s+5] = 0.55; if (pv[s+5] > 0.95) pv[s+5] = 0.95; }
    if (s+6) { if (pv[s+6] < 0.55) pv[s+6] = 0.55; if (pv[s+6] > 0.95) pv[s+6] = 0.95; }
    if (s+7) { if (pv[s+7] < 0.05) pv[s+7] = 0.05; if (pv[s+7] > 0.45) pv[s+7] = 0.45; }

    //IPrinter::print(pv.mid(start, end));
    //for (unsigned int index = start; index <=end; index++)
    //{
    //    //projectControlPoints(pv, index);
    //    projectMeasurePoints(pv, index);
    //}
}

auto Problem2HDirichletDelta::projectControlPoints(DoubleVector &pv, unsigned int index) const -> void
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

auto Problem2HDirichletDelta::projectMeasurePoints(DoubleVector &pv, unsigned int index) const -> void
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

auto Problem2HDirichletDelta::f_layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{
    //return;
    if (ln%2==0 && printLayers)
    {
        QPixmap pxm;
        visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
        pxm.save(QString("data/images/all/%1_f.png").arg(ln/2, 4, 10, QChar('0')));
    }

    const Dimension time = timeDimension();
    const unsigned int L = static_cast<const unsigned int> ( time.size() );
    //IPrinter::printVector(u, nullptr, 100);
    //printf("%4d %.8f %.8f\n", ln, u.min(), u.max());
    //    return;

    if (printLayers)
    {
        FILE *file;
        if (ln == 0) file = fopen("data/data.txt", "w"); else file = fopen("data/data.txt", "a");

        Problem2HDirichletDelta* tmp = const_cast<Problem2HDirichletDelta*>(this);
        std::vector<DoubleMatrix> &rvu = tmp->vu;

        rvu.push_back(u);
        if (rvu.size() > 2*LD+1) rvu.erase(rvu.begin());

        if (ln%2==0)
        {
            std::string filename = "data/txt/f/" + std::to_string(ln/2) + ".txt";
            IPrinter::print(u, filename.data(), u.rows(), u.cols(), 14, 10);


            if (rvu.size() == 2*LD+1)
            {
                double fx = integral(rvu);
                //fprintf(file, "%.10f\n", fx);
                //printf("%d %d %.10f\n", ln, ln-LD, fx);
                fprintf(file, "%.10f %.10f %.10f\n", fx, u.min(), u.max());
            }
            else
            {
                fprintf(file, "%.10f %.10f %.10f\n", 0.0, u.min(), u.max());
            }
        }

        fclose(file);

        //printf("%d,%.10f,%.10f\n", ln, u.min(), u.max());
        //visualString1(u, -1.00, +1.00, 100, 100, Qt::white, Qt::blue, QString("d:/img/string/%1.png").arg(ln,5));
    }
}

auto Problem2HDirichletDelta::b_layerInfo(const DoubleMatrix &p UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{
    //return;
    if (ln%2==0 && printLayers)
    {
        QPixmap pxm;
        visualizeMatrixHeat(p, p.min(), p.max(), pxm, 101, 101);
        pxm.save(QString("data/images/all/%1_b.png").arg(ln/2, 4, 10, QChar('0')));
        return;
    }
}

void Problem2HDirichletDelta::setGridDimensions(const Dimension &time, const Dimension &dimX, const Dimension &dimY)
{
    setTimeDimension(time);
    addSpaceDimension(dimX);
    addSpaceDimension(dimY);

    const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimY.size() );

    mPulseWeightMatrix.resize(M+1, N+1);
    mCrFfxWeightMatrix.resize(M+1, N+1);
    mCrBfxWeightMatrix.resize(M+1, N+1);
}

auto Problem2HDirichletDelta::v(unsigned int i, unsigned int s, const EquaParameter2H &prm, const spif_vectorH &u_info) const -> double
{
    const unsigned int ln = 2*prm.tm[s].i;
    const unsigned int No = static_cast<unsigned int>(prm.No);
    const DoubleMatrix &ok = prm.opt.k[s];
    const DoubleMatrix &oz = prm.opt.z[s];

    double v = 0.0;
    for (unsigned int j=0; j<No; j++)
    {
        v += ok[i][j] * ( u_info[j].vl[ln] - oz[i][j] );
    }
    return v;
}

auto Problem2HDirichletDelta::boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double { return NAN; }

auto Problem2HDirichletDelta::mu(unsigned int, unsigned int) const -> double { return 1.0; }

auto Problem2HDirichletDelta::sign(double x) const -> double
{
    if (x < 0.0)       return -1.0;
    else if (x > 0.0)  return +1.0;
    else               return  0.0;
}

auto Problem2HDirichletDelta::printOptimalParameters() const -> void
{

}
