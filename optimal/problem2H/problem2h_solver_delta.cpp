#include "problem2h_solver_delta.h"

void Problem2HDirichletDelta::Main(int argc, char* argv[])
{
    QGuiApplication app(argc, argv);
    example2();
}

auto Problem2HDirichletDelta::example1() -> void
{
    // Equation parameters
    EquaParameter2H equaPrm;
    equaPrm.a = 1.0;
    equaPrm.alpha = +0.001;

    // Functional parameters
    FuncParameter2H funcPrm;
    funcPrm.Q1 << 0.120;// << 0.22 << 0.24;
    funcPrm.Q2 << 0.110;// << 0.23 << 0.25;

    equaPrm.Nt = 10;
    equaPrm.tm.resize(equaPrm.Nt);
    //    equaPrm.tm[0] = TimeNodePDE( 20, 0.2);
    //    equaPrm.tm[1] = TimeNodePDE( 40, 0.4);
    //    equaPrm.tm[2] = TimeNodePDE( 70, 0.7);
    //    equaPrm.tm[3] = TimeNodePDE(100, 1.0);
    //    equaPrm.tm[4] = TimeNodePDE(120, 1.2);
    //    equaPrm.tm[5] = TimeNodePDE(150, 1.5);
    //    equaPrm.tm[6] = TimeNodePDE(180, 1.8);
    //    equaPrm.tm[7] = TimeNodePDE(210, 2.1);
    //    equaPrm.tm[8] = TimeNodePDE(240, 2.4);
    //    equaPrm.tm[9] = TimeNodePDE(290, 2.9);

    for (unsigned int s=0; s<equaPrm.Nt; s++)
    {
        equaPrm.tm[s].i = (s+1)*30;
        equaPrm.tm[s].t = (s+1)*0.3;
    }

    equaPrm.initPulseParemeters(2);
    equaPrm.pulses[0] = InitialPulse2D(SpacePoint(0.2600, 0.2800), 0.120);
    equaPrm.pulses[1] = InitialPulse2D(SpacePoint(0.7300, 0.7200), 0.110);

    // Optimization parameters
    equaPrm.initParemeters(10, 2, 2);

    DoubleVector ox;
    ox << +0.005385 << -0.002199 << -0.006781 << +0.003281 << -0.002574 << -0.007911 << -0.002827 << +0.003458 << +0.002461 << -0.002868;
    ox << +0.006875 << +0.002599 << -0.009575 << -0.002249 << -0.001259 << +0.004423 << +0.008923 << +0.003974 << +0.004489 << +0.003762;
    ox << -0.009402 << -0.002957 << -0.006576 << +0.006415 << -0.008385 << +0.009921 << +0.001416 << +0.002321 << -0.009274 << -0.001151;
    ox << +0.001424 << +0.008344 << -0.005865 << -0.004854 <<  0.001312 << +0.004445 << -0.006242 << -0.006395 << -0.006137 << +0.004976;
    ox << -0.001543 << +0.004164 << -0.001115 << -0.005455 << -0.008401 << -0.004187 << +0.009259 << +0.005656 << -0.004587 << -0.002575;
    ox << -0.005981 << -0.009127 << -0.003624 << +0.002528 <<  0.002537 << -0.009252 << -0.003723 << -0.001212 << +0.004212 << -0.003419;
    ox <<  0.002306 << -0.002503 <<  0.006129 << -0.001589 <<  0.002329 << -0.002604 << -0.002893 << -0.003764 << +0.004645 << +0.002853;
    ox << -0.002551 << -0.004408 << -0.009114 << -0.009739 << -0.004921 << +0.001554 << +0.004218 << -0.005192 << -0.006846 << -0.008455;
    //ox <<  0.424800 <<  0.392500 <<  0.585300 <<  0.755400;
    //ox <<  0.352600 <<  0.738600 <<  0.859800 <<  0.155200;
    ox <<  0.421200 <<  0.392800 <<  0.584300 <<  0.751500;
    ox <<  0.352500 <<  0.733100 <<  0.852900 <<  0.150500;

    DoubleVector rx;
    rx << -0.125250 << -0.101923 << -0.122294 << -0.003009 <<  0.341065 <<  0.010695 << -0.142474 << -0.057123 <<  0.015842 <<  0.029494;
    rx <<  0.104959 <<  0.086728 <<  0.082632 <<  0.229110 <<  0.154363 <<  0.090279 << -0.147585 << -0.181149 <<  0.034003 <<  0.047790;
    rx <<  0.063607 <<  0.051835 << -0.150809 <<  0.021375 << -0.202155 <<  0.079491 <<  0.028452 <<  0.183324 << -0.067230 <<  0.010065;
    rx <<  0.072888 <<  0.172289 <<  0.009920 << -0.091460 <<  0.087764 <<  0.128026 << -0.175585 << -0.186896 << -0.073537 <<  0.076075;
    rx << -0.009327 << -0.004954 << -0.036839 << -0.042095 << -0.135228 <<  0.086103 <<  0.169792 <<  0.048349 << -0.007319 << -0.014170;
    rx << -0.098346 << -0.118459 << -0.152081 <<  0.003923 <<  0.055058 << -0.042668 << -0.030643 << -0.084040 <<  0.036591 << -0.037824;
    rx <<  0.039041 << -0.041061 << -0.008106 <<  0.007648 <<  0.130315 << -0.161222 << -0.037750 << -0.164920 <<  0.022565 <<  0.027091;
    rx <<  0.033127 <<  0.117900 << -0.141319 << -0.179106 << -0.006130 <<  0.095220 <<  0.010655 << -0.083137 << -0.044724 << -0.093025;
    rx <<  0.236424 << 0.408446  <<  0.695219 <<  0.940548;
    rx <<  0.400140 <<  0.601846 <<  0.729755 <<  0.268136;

    for (unsigned int i= 0; i<40; i++) ox[i] *= 10.0;
    for (unsigned int i=40; i<80; i++) ox[i] *= 10.0;

    //    DoubleVector ox;
    //    ox << -0.008415 << -0.009093 << -0.009581 << -0.008411 << -0.002504 << -0.001411 << -0.002457 << -0.003568 << -0.002481 << -0.001568;
    //    ox << -0.007568 << -0.009259 << -0.006758 << -0.002589 << -0.005989 << -0.002794 << -0.002389 << -0.003974 << -0.004279 << -0.006570;
    //    ox <<  0.001794 << -0.004570 << -0.002476 << -0.006325 << -0.004125 << -0.002354 << -0.001286 << -0.002321 << -0.004694 << -0.001281;
    //    ox << -0.004281 <<  0.003440 << -0.001575 <<  0.005440 <<  0.001403 <<  0.001845 <<  0.006942 <<  0.004512 <<  0.003267 <<  0.004456;
    //    ox <<  0.005403 << -0.004161 << -0.001461 << -0.009245 << -0.008121 << -0.001997 << -0.001299 << -0.005663 << -0.008647 << -0.002357;
    //    ox << -0.004593 <<  0.002837 << -0.006536 <<  0.002975 <<  0.008237 <<  0.009602 <<  0.002237 <<  0.009512 <<  0.008612 <<  0.007539;
    //    ox <<  0.003602 <<  0.005529 <<  0.006529 << -0.001455 <<  0.002539 << -0.002425 << -0.003428 << -0.009751 << -0.004145 << -0.002851;
    //    ox << -0.005251 << -0.008391 << -0.009111 << -0.002839 << -0.002831 <<  0.005444 << -0.001839 <<  0.002364 <<  0.006544 <<  0.008439;
    //    ox <<  0.344500 <<  0.382500 <<  0.826700 <<  0.916400;
    //    ox <<  0.228200 <<  0.636700 <<  0.718400 <<  0.373600;
    //    //ox <<  0.144500 <<  0.282500 <<  0.626700 <<  0.716400;
    //    //ox <<  0.328200 <<  0.686700 <<  0.618400 <<  0.473600;

    //    DoubleVector rx;
    //    rx << -0.107629 << -0.158044 << -0.108772 << -0.111354 << -0.087842 << -0.182702 << -0.178399 << -0.095018 <<  0.077698 <<  0.037216;
    //    rx << -0.029577 << -0.112046 << -0.001187 <<  0.164995 <<  0.045846 <<  0.077016 << -0.110522 << -0.062707 << -0.125581 << -0.139676;
    //    rx <<  0.015192 << -0.031486 << -0.055354 << -0.027126 << -0.169134 << -0.098646 <<  0.146419 <<  0.165281 << -0.146657 << -0.034489;
    //    rx << -0.066862 << -0.009731 <<  0.012910 <<  0.121064 <<  0.024004 <<  0.118247 <<  0.030304 <<  0.027209 <<  0.041410 <<  0.101089;
    //    rx <<  0.022128 << -0.077645 << -0.019998 << -0.089765 << -0.004332 <<  0.082125 <<  0.090855 <<  0.009114 << -0.094706 << -0.021583;
    //    rx <<  0.009802  << 0.115759 << -0.108259  << 0.032727 << -0.038252 << -0.036528 << -0.008582 <<  0.054865 << -0.001345 << -0.072882;
    //    rx <<  0.032744 <<  0.061813 <<  0.066693 <<  0.005660 <<  0.196938 <<  0.081715 << -0.114564 << -0.185153 << -0.147540 << -0.038948;
    //    rx << -0.073464 << -0.077934 << -0.116422 << -0.004037 << -0.045240 << -0.025326 <<  0.011253 <<  0.031801 <<  0.095917 <<  0.136037;
    //    rx <<  0.154108 <<  0.299017 <<  0.618269 <<  0.709296;
    //    rx <<  0.366528 <<  0.634731 <<  0.602577 <<  0.394226;
    
    //    for (unsigned int i= 0; i<40; i++) ox[i] *= 10.0;
    //    for (unsigned int i=40; i<80; i++) ox[i] *= 10.0;

    //    ox = rx;
    equaPrm.OptimalParameterFromVector(ox);
    equaPrm.RegularParameterFromVector(rx);

    // Penalty paramteres
    DoubleVector r; r << 0.0000;// << 0.0000 << 0.0000;
    // Regularization coefficients
    //    DoubleVector e; e << 0.0000 << 0.0000 << 0.0000;
    DoubleVector e; e << 0.0000 << 0.0100 << 0.0000;

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
        prob.optimizeZ = true;
        //prob.optimizeO = false;
        //prob.optimizeC = false;
        prob.funcPrm.vmin.resize(equaPrm.Nc, -0.099);
        prob.funcPrm.vmax.resize(equaPrm.Nc, +0.099);
        prob.LD = 30;
        prob.noise = 0.0;

        prob.funcPrm.regEpsilon = e[i];
        prob.funcPrm.r = r[i];

        if (i==0)
        {
            prob.equaPrm.OptimalParameterToVector(x);
            //            prob.checkGradient1(prob);
            //            prob.optimizeK = true;
            //            prob.optimizeZ = true;
            //            prob.optimizeO = false;
            //            prob.optimizeC = false;
            //            return;
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setProjection(new ProjectionEx1);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.0);
        g.setFunctionTolerance(0.0);
        g.setStepTolerance(0.0);
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
        DoubleVector y; prob.equaPrm.RegularParameterToVector(y);
        IPrinter::printSeperatorLine(nullptr, '*');
        printf("k : "); IPrinter::print(y.mid(00, 19), y.mid(00, 19).length(), 9, 6);
        printf("k : "); IPrinter::print(y.mid(20, 39), y.mid(20, 39).length(), 9, 6);
        printf("z : "); IPrinter::print(y.mid(40, 59), y.mid(40, 59).length(), 9, 6);
        printf("z : "); IPrinter::print(y.mid(60, 79), y.mid(60, 79).length(), 9, 6);
        printf("xy: "); IPrinter::print(y.mid(80, 87), y.mid(80, 87).length(), 9, 6);
        IPrinter::printSeperatorLine(nullptr, '*');

        g.calculate(x);

        IPrinter::printSeperatorLine(nullptr, '=');
        printf("k : "); IPrinter::print(x.mid(00, 19), x.mid(00, 19).length(), 9, 5);
        printf("k : "); IPrinter::print(x.mid(20, 39), x.mid(20, 39).length(), 9, 5);
        printf("z : "); IPrinter::print(x.mid(40, 59), x.mid(40, 59).length(), 9, 5);
        printf("z : "); IPrinter::print(x.mid(60, 79), x.mid(60, 79).length(), 9, 5);
        printf("xy: "); IPrinter::print(x.mid(80, 87), x.mid(80, 87).length(), 9, 5);
        IPrinter::printSeperatorLine(nullptr, '=');

        prob.equaPrm.OptimalParameterFromVector(x);
        prob.printLayers = true;
        if (prob.printLayers)
        {
            prob.setGridDimensions(Dimension(0.010, 0, 400), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
            std::vector<DoubleMatrix> u;
            spif_vectorH u_info, p_info;
            prob.solveForwardIBVP(u, u_info, true, x, 0.25);
        }

        puts("End");

    }
}

auto Problem2HDirichletDelta::example2() -> void
{
    // Equation parameters
    EquaParameter2H equaPrm;
    equaPrm.a = 1.0;
    equaPrm.alpha = +0.0;

    // Functional parameters
    FuncParameter2H funcPrm;
    funcPrm.Q1 << 0.120;// << 0.22 << 0.24;
    funcPrm.Q2 << 0.110;// << 0.23 << 0.25;

    equaPrm.Nt = 10;
    equaPrm.tm.resize(equaPrm.Nt);
    //equaPrm.tm[0] = TimeNodePDE( 20, 0.2);
    //equaPrm.tm[1] = TimeNodePDE( 40, 0.4);
    //equaPrm.tm[2] = TimeNodePDE( 70, 0.7);
    //equaPrm.tm[3] = TimeNodePDE(100, 1.0);
    //equaPrm.tm[4] = TimeNodePDE(120, 1.2);
    //equaPrm.tm[5] = TimeNodePDE(150, 1.5);
    //equaPrm.tm[6] = TimeNodePDE(180, 1.8);
    //equaPrm.tm[7] = TimeNodePDE(210, 2.1);
    //equaPrm.tm[8] = TimeNodePDE(240, 2.4);
    //equaPrm.tm[9] = TimeNodePDE(290, 2.9);

    for (unsigned int s=0; s<equaPrm.Nt; s++) equaPrm.tm[s] = TimeNodePDE((s+1)*30, (s+1)*0.3);

    equaPrm.initPulseParemeters(2);
    equaPrm.pulses[0] = InitialPulse2D(SpacePoint(0.2600, 0.2800), 0.120);
    equaPrm.pulses[1] = InitialPulse2D(SpacePoint(0.7300, 0.7200), 0.110);

    // Optimization parameters
    equaPrm.initParemeters(10, 2, 2);

    DoubleVector ox;
    ox << -0.053850 << -0.021990 << -0.067810 << -0.032810;
    ox << -0.015430 << +0.041640 << -0.011150 << -0.054550;
    ox <<  0.421200 <<  0.392800 <<  0.584300 <<  0.751500;
    ox <<  0.352500 <<  0.733100 <<  0.852900 <<  0.150500;

    //DoubleVector ox;
    //ox << -0.082745 << -0.032199 << -0.018485 << -0.094181;
    //ox << -0.008745 << +0.007857 << -0.007469 << -0.007415;
    //ox <<  0.252500 <<  0.333100 <<  0.493400 <<  0.350500;
    //ox <<  0.521200 <<  0.492800 <<  0.284300 <<  0.151500;

    DoubleVector rx;
    rx << -0.125250 << -0.101923 << -0.122294 << -0.003009;
    rx << -0.009327 << -0.004954 << -0.036839 << -0.042095;
    rx <<  0.236424 <<  0.408446 <<  0.695219 <<  0.940548;
    rx <<  0.400140 <<  0.601846 <<  0.729755 <<  0.268136;

    //for (unsigned int i=0;  i<4;  i++) ox[i] *= 1.05;
    //for (unsigned int i=4;  i<8;  i++) ox[i] *= 1.05;
    //for (unsigned int i=8;  i<12; i++) ox[i] *= 1.05;
    //for (unsigned int i=12; i<16; i++) ox[i] *= 1.05;

    //ox = rx;
    equaPrm.OptimalParameterFromVector(ox);
    equaPrm.RegularParameterFromVector(rx);

    // Penalty paramteres
    DoubleVector r; r << 0.0000;// << 0.0000 << 0.0000;
    // Regularization coefficients
    //DoubleVector e; e << 0.0000 << 0.0000 << 0.0000;
    DoubleVector e; e << 0.0000 << 0.0100 << 0.0000;

    DoubleVector e1; e1 << 0.1000 << 0.1000 << 0.1000;
    DoubleVector e2; e2 << 0.0100 << 0.0100 << 0.0100;

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem2HDirichletDelta prob;
        prob.setGridDimensions(Dimension(0.01, 0, 300), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
        prob.equaPrm = equaPrm;
        prob.funcPrm = funcPrm;
        prob.optimizeK = true;
        prob.optimizeZ = true;
        prob.optimizeO = true;
        prob.optimizeC = true;
        prob.funcPrm.vmin.resize(equaPrm.Nc, -0.099);
        prob.funcPrm.vmax.resize(equaPrm.Nc, +0.099);
        prob.LD = 30;
        prob.noise = 0.0;

        prob.funcPrm.regEpsilon = e[i];
        prob.funcPrm.r = r[i];

        if (i==0)
        {
            prob.equaPrm.OptimalParameterToVector(x);
            prob.checkGradient3(prob);
            //prob.optimizeK = true;
            //prob.optimizeZ = true;
            //prob.optimizeO = false;
            //prob.optimizeC = false;
            return;
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setProjection(new ProjectionEx1);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.0000001);
        g.setFunctionTolerance(0.0000001);
        g.setStepTolerance(0.0000001);
        g.setR1MinimizeEpsilon(e1[i], e2[i]);
        g.setMaxIterations(50);
        g.setNormalize(true);
        g.showExitMessage(true);
        //prob.gm = &g;

        IPrinter::printSeperatorLine(nullptr, '*');
        printf("k : "); IPrinter::print(x.mid( 0,  3), x.mid( 0,  3).length(), 9, 6);
        printf("z : "); IPrinter::print(x.mid( 4,  7), x.mid( 4,  7).length(), 9, 6);
        printf("o : "); IPrinter::print(x.mid( 8, 11), x.mid( 8, 11).length(), 9, 6);
        printf("c : "); IPrinter::print(x.mid(12, 15), x.mid(12, 15).length(), 9, 6);
        IPrinter::printSeperatorLine(nullptr, '*');
        DoubleVector y; prob.equaPrm.RegularParameterToVector(y);
        IPrinter::printSeperatorLine(nullptr, '*');
        printf("k : "); IPrinter::print(y.mid( 0,  3), y.mid( 0,  3).length(), 9, 6);
        printf("z : "); IPrinter::print(y.mid( 4,  7), y.mid( 4,  7).length(), 9, 6);
        printf("o : "); IPrinter::print(y.mid( 8, 11), y.mid( 8, 11).length(), 9, 6);
        printf("c : "); IPrinter::print(y.mid(12, 15), y.mid(12, 15).length(), 9, 6);
        IPrinter::printSeperatorLine(nullptr, '*');

        g.calculate(x);

        IPrinter::printSeperatorLine(nullptr, '=');
        printf("k : "); IPrinter::print(x.mid( 0,  3), x.mid( 0,  3).length(), 9, 6);
        printf("z : "); IPrinter::print(x.mid( 4,  7), x.mid( 4,  7).length(), 9, 6);
        printf("o : "); IPrinter::print(x.mid( 8, 11), x.mid( 8, 11).length(), 9, 6);
        printf("c : "); IPrinter::print(x.mid(12, 15), x.mid(12, 15).length(), 9, 6);
        IPrinter::printSeperatorLine(nullptr, '=');

//        prob.equaPrm.OptimalParameterFromVector(x);
//        prob.printLayers = true;
//        if (prob.printLayers)
//        {
//            prob.setGridDimensions(Dimension(0.010, 0, 800), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
//            std::vector<DoubleMatrix> u;
//            spif_vectorH u_info, p_info;
//            prob.solveForwardIBVP(u, u_info, true, x, 0.25);
//        }

//        puts("End");

//        prob.printLayers = false;
//        prob.setGridDimensions(Dimension(0.010, 0, 300), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
//        DoubleMatrix mx(10, 10, 0.0);
//        unsigned int i1 = 0;
//        for (unsigned int n=5; n<=95; n+=10)
//        {
//            x[14] = n*0.01;
//            unsigned int j1 = 0;
//            for (unsigned int m=5; m<=95; m+=10)
//            {
//                x[15] = m*0.01;

//                double f = prob.fx(x);
//                mx[j1][i1] = f;

//                printf("%d %d %14.10f\n", i1, j1, f);
//                j1++;
//            }
//            i1++;
//        }

//        FILE* file = fopen("data/matrix10.txt", "w");
//        IPrinter::print(mx, mx.rows(), mx.cols(), 14, 10, file);
//        fclose(file);

//        puts("End");
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
    const_cast<Problem2HDirichletDelta&>(prob).printLayers = true;
    prob.gradient(pv, ag);
    const_cast<Problem2HDirichletDelta&>(prob).printLayers = false;
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

auto Problem2HDirichletDelta::checkGradient3(const Problem2HDirichletDelta &prob) -> void
{
    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prob.equaPrm.OptimalParameterToVector(pv);
    printf("ok: "); IPrinter::print(pv.mid(0,  3), pv.mid(0,  3).length(), 9, 6);
    printf("oz: "); IPrinter::print(pv.mid(4,  7), pv.mid(4,  7).length(), 9, 6);
    printf("xy: "); IPrinter::print(pv.mid(8, 15), pv.mid(8, 15).length(), 9, 6);
    IPrinter::printSeperatorLine();

    DoubleVector rv;
    prob.equaPrm.RegularParameterToVector(rv);
    printf("rk: "); IPrinter::print(rv.mid(0,  3), rv.mid(0,  3).length(), 9, 6);
    printf("rz: "); IPrinter::print(rv.mid(4,  7), rv.mid(4,  7).length(), 9, 6);
    printf("xy: "); IPrinter::print(rv.mid(8, 15), rv.mid(8, 15).length(), 9, 6);
    IPrinter::printSeperatorLine();

    DoubleVector ag(pv.length());
    double functional = prob.fx(pv);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    const_cast<Problem2HDirichletDelta&>(prob).printLayers = false;
    prob.gradient(pv, ag);
    const_cast<Problem2HDirichletDelta&>(prob).printLayers = false;
    puts("Gradients are calculated.");
    //    return;

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    const unsigned int Nc = prob.equaPrm.Nc;
    const unsigned int No = prob.equaPrm.No;
    const unsigned int offset = Nc*No;

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
    {
        IPrinter::printSeperatorLine("k");
        DoubleVector pk0 = pv.mid(0, offset-1);
        DoubleVector ak0 = ag.mid(0, offset-1);
        DoubleVector nk1 = ng1.mid(0, offset-1);
        DoubleVector nk2 = ng2.mid(0, offset-1);

        IPrinter::print(pk0,pk0.length(),W,P);
        IPrinter::print(ak0,ak0.length(),W,P); ak0.L2Normalize();
        IPrinter::print(nk1,nk1.length(),W,P); nk1.L2Normalize();
        IPrinter::print(nk2,nk2.length(),W,P); nk2.L2Normalize();
        IPrinter::print(ak0,ak0.length(),W,P);
        IPrinter::print(nk1,nk1.length(),W,P);
        IPrinter::print(nk2,nk2.length(),W,P);
    }
    //z------------------------------------------------------//
    {
        IPrinter::printSeperatorLine("z");
        DoubleVector pz0 = pv.mid(offset, 2*offset-1);
        DoubleVector az0 = ag.mid(offset, 2*offset-1);
        DoubleVector nz1 = ng1.mid(offset, 2*offset-1);
        DoubleVector nz2 = ng2.mid(offset, 2*offset-1);

        IPrinter::print(pz0,pz0.length(),W,P);
        IPrinter::print(az0,az0.length(),W,P); az0.L2Normalize();
        IPrinter::print(nz1,nz1.length(),W,P); nz1.L2Normalize();
        IPrinter::print(nz2,nz2.length(),W,P); nz2.L2Normalize();
        IPrinter::print(az0,az0.length(),W,P);
        IPrinter::print(nz1,nz1.length(),W,P);
        IPrinter::print(nz2,nz2.length(),W,P);
    }
    //xi------------------------------------------------------//
    {
        IPrinter::printSeperatorLine("xi");
        DoubleVector pe0 = pv.mid(2*offset, 2*offset+2*prob.equaPrm.No-1);
        DoubleVector ae0 = ag.mid(2*offset, 2*offset+2*prob.equaPrm.No-1);
        DoubleVector ne1 = ng1.mid(2*offset, 2*offset+2*prob.equaPrm.No-1);
        DoubleVector ne2 = ng2.mid(2*offset, 2*offset+2*prob.equaPrm.No-1);

        IPrinter::print(pe0,pe0.length(),W,P);
        IPrinter::print(ae0,ae0.length(),W,P); ae0.L2Normalize();
        IPrinter::print(ne1,ne1.length(),W,P); ne1.L2Normalize();
        IPrinter::print(ne2,ne2.length(),W,P); ne2.L2Normalize();
        IPrinter::print(ae0,ae0.length(),W,P);
        IPrinter::print(ne1,ne1.length(),W,P);
        IPrinter::print(ne2,ne2.length(),W,P);
    }

    //eta------------------------------------------------------//
    {
        IPrinter::printSeperatorLine("eta");
        DoubleVector px0 = pv.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
        DoubleVector ax0 = ag.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
        DoubleVector nx1 = ng1.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
        DoubleVector nx2 = ng2.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);

        IPrinter::print(px0,px0.length(),W,P);
        IPrinter::print(ax0,ax0.length(),W,P); ax0.L2Normalize();
        IPrinter::print(nx1,nx1.length(),W,P); nx1.L2Normalize();
        IPrinter::print(nx2,nx2.length(),W,P); nx2.L2Normalize();
        IPrinter::print(ax0,ax0.length(),W,P);
        IPrinter::print(nx1,nx1.length(),W,P);
        IPrinter::print(nx2,nx2.length(),W,P);
        IPrinter::printSeperatorLine();
    }
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
            double pnt = 0.0;
#ifdef USE_PENALTY
            pnt = penalty(u_info, equaPrm);
#endif
            double nrm = 0.0;
#ifdef USE_NORM
            nrm = norm(equaPrm);
#endif
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

#if defined (DISCRETE_DELTA_TIME_1)
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
#endif

#if defined (DISCRETE_DELTA_TIME_2)
    const DoubleMatrix &ok = prm.opt.k;
    const DoubleMatrix &rk = prm.reg.k;
    const DoubleMatrix &oz = prm.opt.z;
    const DoubleMatrix &rz = prm.reg.z;

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            _norm += (ok[i][j] - rk[i][j])*(ok[i][j] - rk[i][j]);
            _norm += (oz[i][j] - rz[i][j])*(oz[i][j] - rz[i][j]);
        }
    }
#endif

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
#if defined (DISCRETE_DELTA_TIME_1)
    for (unsigned int j=0; j<No; j++)
    {
        const DoubleMatrix &ok = prm.opt.k[s];
        const DoubleMatrix &oz = prm.opt.z[s];

        const SpacePointInfoH &u_xij = u_info[j];
        vi += ok[i][j] * ( u_xij.vl[ln] - oz[i][j] );
    }
#endif

#if defined (DISCRETE_DELTA_TIME_2)
    for (unsigned int j=0; j<No; j++)
    {
        const SpacePointInfoH &u_xij = u_info[j];
        vi += prm.opt.k[i][j] * ( u_xij.vl[ln] - prm.opt.z[i][j] );
    }
#endif

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
#if defined (DISCRETE_DELTA_TIME_1)
                for (unsigned int s=0; s<Nt; s++)
                {
                    const unsigned int ln = 2*equaPrm.tm[s].i;
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        const SpacePointInfoH &pi = p_info[i];
                        for (unsigned int j=0; j<No; j++)
                        {
                            const SpacePointInfoH &uj = u_info[j];
                            double grad_Kij = 0.0;
                            double zij = equaPrm.opt.z[s][i][j];
#ifdef USE_PENALTY
                            grad_Kij += -(uj.vl[ln] - zij) * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
#endif

#ifdef USE_NORM
                            grad_Kij += +2.0*regEpsilon*(equaPrm.opt.k[s][i][j] - equaPrm.reg.k[s][i][j]);
#endif
                            g[gi++] += grad_Kij * (1.0/(double(Q1.length())*double(Q2.length())));
                        }
                    }
                }
#endif

#if defined (DISCRETE_DELTA_TIME_2)
                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfoH &pi = p_info[i];
                    for (unsigned int j=0; j<No; j++)
                    {
                        const SpacePointInfoH &uj = u_info[j];
                        double zij = equaPrm.opt.z[i][j];

                        double grad_Kij = 0.0;
                        for (unsigned int s=0; s<Nt; s++)
                        {
                            const unsigned int ln = 2*equaPrm.tm[s].i;
                            grad_Kij += -(uj.vl[ln] - zij) * pi.vl[ln];
#ifdef USE_PENALTY
                            grad_Kij += -(uj.vl[ln] - zij) * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
#endif
                        }


#ifdef USE_NORM
                        grad_Kij += +2.0*regEpsilon*(equaPrm.opt.k[s][i][j] - equaPrm.reg.k[s][i][j]);
#endif
                        g[gi++] += grad_Kij * (1.0/(double(Q1.length())*double(Q2.length())));
                    }
                }
#endif
            }
            else
            {
#if defined (DISCRETE_DELTA_TIME_1)
                for (unsigned int s=0; s<Nt; s++)
                {
#endif
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        for (unsigned int j=0; j<No; j++)
                        {
                            g[gi++] = 0.0;
                        }
                    }
#if defined (DISCRETE_DELTA_TIME_1)
                }
#endif
            }

            // z
            if (optimizeZ)
            {
                //puts("Calculating z gradients...");
#if defined (DISCRETE_DELTA_TIME_1)
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
#ifdef USE_PENALTY
                            grad_Zij += kij * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
#endif
#ifdef USE_NORM
                            grad_Zij += +2.0*regEpsilon*(equaPrm.opt.z[s][i][j] - equaPrm.opt.z[s][i][j]);
#endif
                            g[gi++] += grad_Zij * (1.0/(double(Q1.length())*double(Q2.length())));
                        }
                    }
                }
#endif

#if defined (DISCRETE_DELTA_TIME_2)
                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfoH &pi = p_info[i];
                    for (unsigned int j=0; j<No; j++)
                    {
                        double grad_Zij = 0.0;

                        double kij = equaPrm.opt.k[i][j];
                        for (unsigned int s=0; s<Nt; s++)
                        {
                            const unsigned int ln = 2*equaPrm.tm[s].i;
                            grad_Zij += kij * pi.vl[ln];
#ifdef USE_PENALTY
                            grad_Zij += kij * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
#endif
                        }
#ifdef USE_NORM
                        grad_Zij += +2.0*regEpsilon*(equaPrm.opt.z[s][i][j] - equaPrm.opt.z[s][i][j]);
#endif
                        g[gi++] = grad_Zij * (1.0/(double(Q1.length())*double(Q2.length())));
                    }
                }
#endif
            }
            else
            {
#if defined (DISCRETE_DELTA_TIME_1)
                for (unsigned int s=0; s<Nt; s++)
                {
#endif
                    for (unsigned int i=0; i<Nc; i++)
                    {
                        for (unsigned int j=0; j<No; j++)
                        {
                            g[gi++] = 0.0;
                        }
                    }
#if defined (DISCRETE_DELTA_TIME_1)
                }
#endif
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
#if defined (DISCRETE_DELTA_TIME_1)
                            double kij = equaPrm.opt.k[s][i][j];
#endif
#if defined (DISCRETE_DELTA_TIME_2)
                            double kij = equaPrm.opt.k[i][j];
#endif
                            gradXijX += -kij * uj.dx[ln] * p_info[i].vl[ln];
                            gradXijY += -kij * uj.dy[ln] * p_info[i].vl[ln];
#ifdef USE_PENALTY
                            gradXijX += -kij * uj.dx[ln] * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
                            gradXijY += -kij * uj.dy[ln] * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
#endif
                        }
                    }
#ifdef USE_NORM
                    gradXijX += 2.0*regEpsilon*(equaPrm.opt.ksi[j].x - equaPrm.reg.ksi[j].x);
                    gradXijY += 2.0*regEpsilon*(equaPrm.opt.ksi[j].y - equaPrm.reg.ksi[j].y);
#endif

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
#if defined (DISCRETE_DELTA_TIME_1)
                            double kij = equaPrm.opt.k[s][i][j];
                            double zij = equaPrm.opt.z[s][i][j];
#endif
#if defined (DISCRETE_DELTA_TIME_2)
                            double kij = equaPrm.opt.k[i][j];
                            double zij = equaPrm.opt.z[i][j];
#endif
                            gradEtaiX += -pi.dx[ln] * kij * (u_info[j].vl[ln] - zij);
                            gradEtaiY += -pi.dy[ln] * kij * (u_info[j].vl[ln] - zij);
                        }
                    }
#ifdef USE_NORM
                    gradEtaiX += 2.0*regEpsilon*(equaPrm.opt.eta[i].x - equaPrm.reg.eta[i].x);
                    gradEtaiY += 2.0*regEpsilon*(equaPrm.opt.eta[i].y - equaPrm.reg.eta[i].y);
#endif
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
    //std::vector<DeltaGrid2D> measuremntGirdList(No);
    //std::vector<DeltaGrid2D> cntrlDeltaGridList(Nc);
    const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList.resize(No);
    const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList.resize(Nc);
    for (unsigned int j=0; j<No; j++)
    {
        const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList[j].initGrid(N, hx, M, hy);
        const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList[j].distributeGauss(equaPrm.opt.ksi[j], MSRMT_SIGMA, MSRMT_SIGMA);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList[i].initGrid(N, hx, M, hy);
        const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList[i].distributeGauss(equaPrm.opt.eta[i], CNTRL_SIGMA, CNTRL_SIGMA);
    }
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) prepareInfo(No,  equaPrm.opt.ksi, u_info, 2*LLD+1);
    //----------------------------------------------------------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    initPulseWeightMatrix(equaPrm.pulses);
    //----------------------------------------------------------------------------------------------//
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
    //----------------------------------------------------------------------------------------------//
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*ht;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = 1.0*ht;
    SpaceNodePDE sn05, sn10;
    sn05.i = 0; sn05.x = 0.0; sn10.i = static_cast<int>(N); sn10.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        u05[m][0] = u10[m][0] = u05[m][N] = u10[m][N] = 0.0;
        //sn05.j = static_cast<int>(m); sn05.y = m*hy; u10[m][0] = f_boundary(sn05, tn10); u05[m][0] = f_boundary(sn00, tn05);
        //sn10.j = static_cast<int>(m); sn10.y = m*hy; u10[m][N] = f_boundary(sn10, tn10); u05[m][N] = f_boundary(sn10, tn05);
    }
    sn05.j = 0; sn05.y = 0.0; sn10.j = static_cast<int>(M); sn10.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {

        u05[0][n] = u10[0][n] = u05[M][n] = u10[M][n] = 0.0;
        //sn05.i = static_cast<int>(n); sn05.x = n*hx; u10[0][n] = f_boundary(sn05, tn10); u05[0][n] = f_boundary(sn05, tn05);
        //sn10.i = static_cast<int>(n); sn10.x = n*hx; u10[M][n] = f_boundary(sn10, tn10); u05[M][n] = f_boundary(sn10, tn05);
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
    //----------------------------------------------------------------------------------------------//
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
            u15[m][0] = u20[m][0] = u15[m][N] = u20[m][N] = 0.0;
            //sn0.j = static_cast<int>(m); sn0.y = m*hy; u15[m][0] = f_boundary(sn0, tn15); u20[m][0] = f_boundary(sn0, tn20);
            //sn1.j = static_cast<int>(m); sn1.y = m*hy; u15[m][N] = f_boundary(sn1, tn15); u20[m][N] = f_boundary(sn1, tn20);
        }
        sn0.j = static_cast<int>(0); sn0.y = 0*hy; sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            u15[0][n] = u20[0][n] = u15[M][n] = u20[M][n] = 0.0;
            //sn0.i = static_cast<int>(n); sn0.x = n*hx; u15[0][n] = f_boundary(sn0, tn15); u20[0][n] = f_boundary(sn0, tn20);
            //sn1.i = static_cast<int>(n); sn1.x = n*hx; u15[M][n] = f_boundary(sn1, tn15); u20[M][n] = f_boundary(sn1, tn20);
        }
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        //fxMatrixForward(u10, cntrlDeltaGridList, measuremntGirdList, 2*ln-2);
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
                dx[n-1] += ht_ht_025 * mfxMatrixForward[m][n];
            }
            dx[0]   -= u15[m][0]*m_aa_htht__hxhx_025_lambda;
            dx[N-2] -= u15[m][N]*m_aa_htht__hxhx_025_lambda;
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }
        if (use == true) add2Info(u15, u_info, 2*ln-1, hx, hy, measuremntGirdList); f_layerInfo(u15, 2*ln-1);
        /**************************************************** x direction apprx ***************************************************/
        /**************************************************** y direction apprx ***************************************************/
        //fxMatrixForward(u15, cntrlDeltaGridList, measuremntGirdList, 2*ln-1);
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
                dy[m-1] += ht_ht_025 * mfxMatrixForward[m][n];
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

    for (unsigned int j=0; j<No; j++) const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList[j].cleanGrid();
    const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList.clear();
    for (unsigned int i=0; i<Nc; i++) const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList[i].cleanGrid();
    const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList.clear();

    u00.clear();
    u05.clear();
    u10.clear();
    u15.clear();
    u20.clear();
    //puts("+++ solveForwardIBVP");
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
    //std::vector<DeltaGrid2D> measuremntGirdList(No);
    //std::vector<DeltaGrid2D> cntrlDeltaGridList(Nc);
    const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList.resize(No);
    const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList.resize(Nc);
    for (unsigned int j=0; j<No; j++)
    {
        const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList[j].initGrid(N, hx, M, hy);
        const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList[j].distributeGauss(equaPrm.opt.ksi[j], MSRMT_SIGMA, MSRMT_SIGMA);
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList[i].initGrid(N, hx, M, hy);
        const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList[i].distributeGauss(equaPrm.opt.eta[i], CNTRL_SIGMA, CNTRL_SIGMA);
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
        fxMatrixBackward(p10, cntrlDeltaGridList, measuremntGirdList, 2*ln+2, u_info);
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
                dx[n-1] += ht_ht_025 * mfxMatrixBackward[m][n];
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
        fxMatrixBackward(p15, cntrlDeltaGridList, measuremntGirdList, 2*ln+1, u_info);
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
                dy[m-1] += ht_ht_025 * mfxMatrixBackward[m][n];
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

    for (unsigned int j=0; j<No; j++) const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList[j].cleanGrid();
    const_cast<Problem2HDirichletDelta*>(this)->measuremntGirdList.clear();
    for (unsigned int i=0; i<Nc; i++) const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList[i].cleanGrid();
    const_cast<Problem2HDirichletDelta*>(this)->cntrlDeltaGridList.clear();

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

    DoubleMatrix &mForwardInitialMatrix = const_cast<Problem2HDirichletDelta*>(this)->mInitialMatrixForward;

    const unsigned int Nb = static_cast<unsigned int>(pulses.size());
    DeltaGrid2D *deltaGrids = new DeltaGrid2D[Nb];
    for (unsigned int s=0; s<Nb; s++)
    {
        deltaGrids[s].initGrid(N, hx, M, hy);
        deltaGrids[s].distributeGauss(pulses[s].theta, PULSE_SIGMA, PULSE_SIGMA);
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            mForwardInitialMatrix[m][n] = 0.0;
            for (unsigned int s=0; s<Nb; s++)
            {
                mForwardInitialMatrix[m][n] += pulses[s].q*deltaGrids[s].weight(n, m);
            }
        }
    }

    for (unsigned int s=0; s<Nb; s++) deltaGrids[s].cleanGrid();
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
    return mInitialMatrixForward[m][n];
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

    printf("I[%3d]: F:%.6f I:%.6f P:%.6f N:%.5f R:%.3f e:%.3f a:%10.6f\n", i, f, ing, pnt, nrm, funcPrm.r, funcPrm.regEpsilon, alpha);
    printf("min:%10.6f max:%10.6f U:%.8f\n", u.front().min(), u.front().max(), integralU(u.front()));
    printf("min:%10.6f max:%10.6f U:%.8f\n", u.back().min(),  u.back().max(),  integralU(u.back()));
    printf("k:%10.6f %10.6f %10.6f %8.4f z:%10.6f %10.6f %10.6f %10.6f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]);
    printf("o:%10.6f %10.6f %10.6f %8.4f c:%10.6f %10.6f %10.6f %10.6f\n", x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    IPrinter::printSeperatorLine("-");

    //unsigned int s = 2*(equaPrm.Nc*equaPrm.No);
    //printf("o: %6.4f %6.4f %6.4f %6.4f c: %6.4f %6.4f %6.4f %6.4f\n", x[s+0], x[s+1], x[s+2], x[s+3], x[s+4], x[s+5], x[s+6], x[s+7]);
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

    //    prob->optimizeK = (i%3 == 2);
    //    prob->optimizeZ = (i%3 == 0);
    //    prob->optimizeO = (i%3 == 1);
    //    prob->optimizeC = (i%3 == 1);

    //    prob->optimizeK = (i%2 == 0);
    //    prob->optimizeZ = (i%2 == 0);
    //    prob->optimizeO = (i%2 == 1);
    //    prob->optimizeC = (i%2 == 1);
}

auto Problem2HDirichletDelta::fxMatrixForward(const DoubleMatrix &u, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids, unsigned int ln) const -> void
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

    for (unsigned int m=0; m<=M; m++) for (unsigned int n=0; n<=N; n++) const_this->mfxMatrixForward[m][n] = 0.0;

    //bool first = true;
    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        const unsigned int sln = equaPrm.tm[s].i;
        if (ln == 2*sln) { wt = 2.0/ht; } else { continue; }

        //double ts = equaPrm.tm[s].t;
        //double sigma = ht;
        //double t = ln*ht*0.5;

        //if (first) { for (unsigned int m=0; m<=M; m++) for (unsigned int n=0; n<=N; n++) const_this->mfxMatrixForward[m][n] = 0.0; }
        //first = false;

        //if (ln == 2*sln || ln == 2*sln-1) { wt = 1.0/ht; } else { continue; }
        //printf("%d\n", ln);
        //if (2*sln <= ln && ln <= 2*sln+16) { wt = 2.0*(1.0/(sqrt(2.0*M_PI)*sigma)) * exp(((ts-t)*(ts-t))/(-2.0*sigma*sigma)); } else { continue; }
        //if (2*sln-8 <= ln && ln <= 2*sln+8) { wt = (1.0/(sqrt(2.0*M_PI)*sigma)) * exp(((ts-t)*(ts-t))/(-2.0*sigma*sigma)); } else { continue; }
        //if (wt != 0.0 && s == 0)
        //{
        //   printf("%d %d %.10f %.10f %.10f\n", ln, 2*sln, wt, ts, t);
        //}

        double* _u = new double[No];
        for (unsigned int j=0; j<No; j++)
        {
            _u[j] = measurementDeltaGrids[j].lumpPointGauss(u);
            _u[j] *= (1.0 + noise * ( rand() % 2==0 ? +1.0 : -1.0));
        }

        double *_v = new double[Nc];
        for (unsigned int i=0; i<Nc; i++)
        {
            _v[i] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
#if defined (DISCRETE_DELTA_TIME_1)
                _v[i] += equaPrm.opt.k[s][i][j] * (_u[j] - equaPrm.opt.z[s][i][j]);
#endif
#if defined (DISCRETE_DELTA_TIME_2)
                _v[i] += equaPrm.opt.k[i][j] * (_u[j] - equaPrm.opt.z[i][j]);
#endif
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
                    const_this->mfxMatrixForward[m][n] += _v[i] * dg.weight(n,m) * wt;
                }
            }
        }
        delete [] _v;
    }
}

auto Problem2HDirichletDelta::fxMatrixBackward(const DoubleMatrix &p, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids, unsigned int ln, const spif_vectorH &u_info) const -> void
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

    for (unsigned int m=0; m<=M; m++) for (unsigned int n=0; n<=N; n++) const_this->mfxMatrixBackward[m][n] = 0.0;

    //bool first = true;
    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        const unsigned int sln = equaPrm.tm[s].i;
        if (ln == 2*sln) { wt = 2.0/ht; } else { continue; }

        //double ts = equaPrm.tm[s].t;
        //double sigma = ht;
        //double t = ln*ht*0.5;

        //if ((ln == 2*sln-0) || (ln == 2*sln-1)) { wt = 1.0/ht; } else { continue; }
        //if ((ln == 2*sln+2) || (ln == 2*sln+1)) { wt = 1.0/ht; } else { continue; }
        //if ((ln == 2*sln-0) || (ln == 2*sln+1)) { wt = 1.0/ht; } else { continue; }
        //if (ln == 2*sln+1) { wt = 2.0/ht; } else { continue; }
        //printf("ln: %d %d\n", ln, sln);

        //if (first) for (unsigned int m=0; m<=M; m++) for (unsigned int n=0; n<=N; n++) const_this->mfxMatrixBackward[m][n] = 0.0;
        //first = false;
        //if (ln == 2*sln || ln == 2*sln+1) { wt = 1.0/ht; } else { continue; }
        //if (2*sln-16 <= ln && ln <= 2*sln) { wt = 2.0*(1.0/(sqrt(2.0*M_PI)*sigma)) * exp(((ts-t)*(ts-t))/(-2.0*sigma*sigma)); } else { continue; }
        //if (2*sln-8 <= ln && ln <= 2*sln+8) { wt = (1.0/(sqrt(2.0*M_PI)*sigma)) * exp(((ts-t)*(ts-t))/(-2.0*sigma*sigma)); } else { continue; }

        double* _p = new double[Nc];
        for (unsigned int i=0; i<Nc; i++) _p[i] = controlDeltaGrids[i].lumpPointGauss(p);

        double *_w = new double[No];
        for (unsigned int j=0; j<No; j++)
        {
            _w[j] = 0.0;
            for (unsigned int i=0; i<Nc; i++)
            {

#if defined (DISCRETE_DELTA_TIME_1)
                _w[j] += equaPrm.opt.k[s][i][j] * _p[i];
#endif

#if defined (DISCRETE_DELTA_TIME_2)
                _w[j] += equaPrm.opt.k[i][j] * _p[i];
#endif

#ifdef USE_PENALTY
                _w[j] += equaPrm.opt.k[s][i][j] * 2.0*r*gpi(i, s, u_info, equaPrm)*sgn(g0i(i, s, u_info, equaPrm));
#endif
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
                    const_this->mfxMatrixBackward[m][n] += _w[j] * dg.weight(n,m) * wt;
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

auto Problem2HDirichletDelta::add2Info(const DoubleMatrix &u, spif_vectorH &info, unsigned int ln, double /*hx*/, double /*hy*/, const std::vector<DeltaGrid2D> &deltaList) const -> void
{
    const unsigned int N = static_cast<unsigned int>(deltaList.size());

    for (unsigned int i=0; i<N; i++)
    {
        const DeltaGrid2D &deltagrid = deltaList[i];
        SpacePointInfoH &ui = info[i];
        ui.vl[ln] = deltagrid.lumpPointGauss(u, ui.dx[ln], ui.dy[ln]);
    }
}

auto Problem2HDirichletDelta::project(DoubleVector &, unsigned int) -> void {}

auto Problem2HDirichletDelta::project(DoubleVector &pv) const -> void
{
    unsigned int Nc = equaPrm.Nc;
    unsigned int No = equaPrm.No;
    unsigned int Nt = equaPrm.Nt;

    unsigned int start = 2*Nc*No;//*Nt;
    unsigned int end = 2*Nc*No + 2*No + 2*Nc - 1;

    for (unsigned int index = start; index <= end; index++)
    {
        if (pv[index] <= 0.05) pv[index] = 0.05;
        if (pv[index] >= 0.95) pv[index] = 0.95;
    }
    return;

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
    fxMatrixForward(u, cntrlDeltaGridList, measuremntGirdList, ln);

    //return;
//    if (ln%2==0 && printLayers)
//    {
//        QPixmap pxm;
//        visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
//        pxm.save(QString("data/images/all/%1_f.png").arg(ln/2, 4, 10, QChar('0')));

//        std::string filename = "data/txt/f/" + std::to_string(ln/2) + ".txt";
//        IPrinter::print(u, filename.data(), u.rows(), u.cols(), 14, 10);
//        return;
//    }

    //const Dimension time = timeDimension();
    //const unsigned int L = static_cast<const unsigned int> ( time.size() );
    //IPrinter::printVector(u, nullptr, 100);
    //printf("%4d %.8f %.8f\n", ln, u.min(), u.max());
    //return;

    if (printLayers)
    {
        FILE *file;
        if (ln == 0) file = fopen("data/data1.txt", "w"); else file = fopen("data/data1.txt", "a");

        Problem2HDirichletDelta* tmp = const_cast<Problem2HDirichletDelta*>(this);
        std::vector<DoubleMatrix> &rvu = tmp->vu;

        rvu.push_back(u);
        if (rvu.size() > 2*LD+1) rvu.erase(rvu.begin());

        if (ln%2==0)
        {
            //std::string filename = "data/txt/f/" + std::to_string(ln/2) + ".txt";
            //IPrinter::print(u, filename.data(), u.rows(), u.cols(), 14, 10);


            if (rvu.size() == 2*LD+1)
            {
                double fx = integral(rvu);
                //fprintf(file, "%4d %4d %14.10f\n", ln, ln-2*LD+1, fx);
                //printf("%d %d %.10f\n", ln, ln-LD, fx);
                fprintf(file, "%4d %6.4f %14.10f %14.10f %14.10f\n", ln, (ln-2*LD)*0.005, fx, u.min(), u.max());
            }
            else
            {
                //fprintf(file, "%.10f %.10f %.10f\n", 0.0, u.min(), u.max());
            }
        }

        fclose(file);

        //printf("%d,%.10f,%.10f\n", ln, u.min(), u.max());
        //visualString1(u, -1.00, +1.00, 100, 100, Qt::white, Qt::blue, QString("d:/img/string/%1.png").arg(ln,5));
    }
}

auto Problem2HDirichletDelta::b_layerInfo(const DoubleMatrix &p UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const -> void
{
    //fxMatrixBackward(p, cntrlDeltaGridList, measuremntGirdList, ln, u_info);

    return;

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

    mInitialMatrixForward.resize(M+1, N+1);
    mfxMatrixForward.resize(M+1, N+1);
    mfxMatrixBackward.resize(M+1, N+1);
}

auto Problem2HDirichletDelta::v(unsigned int i, unsigned int s, const EquaParameter2H &prm, const spif_vectorH &u_info) const -> double
{
    const unsigned int ln = 2*prm.tm[s].i;
    const unsigned int No = static_cast<unsigned int>(prm.No);
#if defined (DISCRETE_DELTA_TIME_1)
    const DoubleMatrix &ok = prm.opt.k[s];
    const DoubleMatrix &oz = prm.opt.z[s];
#endif
#if defined (DISCRETE_DELTA_TIME_2)
    const DoubleMatrix &ok = prm.opt.k;
    const DoubleMatrix &oz = prm.opt.z;
#endif

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
