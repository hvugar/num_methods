#include "hyperboliccontrolh.h"

void HyperbolicControlH::main()
{
    HyperbolicControlH hcx;

    //double a = 0.6;
    //double b = 1.1;
    double t = 0.92;
    //0.9544858420

    //    printf("%.8f %.16f\n", a, hcx.fx(a));
    //    printf("%.8f %.16f\n", b, hcx.fx(b));

    //    goldenSectionSearch(a, b, t, &hcx, 0.001);
    //    R1Minimize::HalphIntervalMethod(a, b, t, &hcx, 0.01);
    //    stranghLineSearch(0.8, 0.5, a, b, &hcx);
    //    printf("%.10f %.10f\n", a, b);

    //    for (double t=0.80; t<1.10; t+=0.01)
    //    {
    //        hcx.fx(t);
    //    }

    hcx.fx(t);
    printf("optimal t: %.10f\n", t);
}

HyperbolicControlH::HyperbolicControlH()
{
    U = 0.0;
    lamda = 0.25;
    a = 1.0;
    L = 1;
}

double HyperbolicControlH::fx(double t)
{
    t0 = 0.0;
    t1 = t;
    x0 = 0.0;
    x1 = 1.0;
    N = 100;
    hx = 0.01;

    t1 = t;
    ht = 0.01;
    M = (unsigned int) round((t1-t0)/ht);
    D = 10;
    xi = 0.2;
    Xi = 20;

    printf("%d %d %f %f\n", M, N, ht, hx);

    DoubleVector v((L+2)*(M+D-1));
    for (unsigned int j=0; j<=(M+D-2); j++)
    {
        v[0*(M+D-1)+j] = sin(j*ht);
        v[1*(M+D-1)+j] = cos(j*ht);
        v[2*(M+D-1)+j] = sin(2.0*j*ht);
    }

    double v_1[] = {0.9597836792,3.2300669991,0.8322871093,1.7038695673,2.4238508089,1.5468238795,1.3363084770,2.0028298435,1.9416819887,1.3855359528,1.5185076162,1.9184514233,1.6900852564,1.3395703375,1.5539820091,1.7940493048,1.5286390115,1.3030976966,1.5273323529,1.6643448935,1.4015846439,1.2609797910,1.4757775081,1.5419421621,1.2910589354,1.2084462389,1.4104626258,1.4318085685,1.1975947235,1.1503682876,1.3301306714,1.3188211557,1.1071450089,1.0884720840,1.2448872413,1.2072704875,1.0175656769,1.0213526071,1.1561256568,1.1003374132,0.9350283800,0.9572897610,1.0682247145,0.9942777199,0.8510454688,0.8896441091,0.9793446343,0.8906620145,0.7677168552,0.8186799577,0.8875030650,0.7869147674,0.6829376539,0.7417469540,0.7898415289,0.6844242692,0.6060793551,0.6723161855,0.6945608006,0.5810619892,0.5290817828,0.6054703165,0.6057104681,0.4886334863,0.4648494121,0.5481058040,0.5209463612,0.3957274773,0.3942782333,0.4790823478,0.4305212362,0.3156996210,0.3475233393,0.4240209920,0.3432756539,0.2475617744,0.3211870516,0.3731337213,0.2435433052,0.1804732672,0.3068330370,0.3058387210,0.1116937702,0.1245365379,0.3147615277,0.1813397567,-0.0603099505,0.2210762890,0.3667339721,-0.4167541316,0.6438095583,-0.2755322745,0.2410625112,0.0959413958,0.0768579965,0.0301410573,-0.0742400869,0.2963880924,-0.4548065146,0.8360259786,0.8414709848,1.0497962039,3.1710479537,0.9368006339,1.7340839237,2.4255098859,1.6165650983,1.3992843657,2.0193913838,1.9840701058,1.4638680628,1.5667959178,1.9418258752,1.7468009322,1.4149545430,1.5958279181,1.8217298077,1.5880693137,1.3713308281,1.5626681353,1.6904354936,1.4551648994,1.3174024137,1.5016994632,1.5614407433,1.3336386444,1.2503203805,1.4252989892,1.4430126342,1.2292717053,1.1806328148,1.3378634443,1.3251013459,1.1309879418,1.1104007048,1.2485908991,1.2110245575,1.0349058096,1.0346283899,1.1540941596,1.0999754088,0.9450443543,0.9604393248,1.0574825005,0.9874982653,0.8548049138,0.8863069899,0.9628440299,0.8786032136,0.7664678836,0.8113277602,0.8671278378,0.7690747489,0.6750879674,0.7304112201,0.7675279228,0.6620780434,0.5902848065,0.6543317838,0.6695743322,0.5571342116,0.5104056953,0.5851310928,0.5807025442,0.4661081170,0.4480384435,0.5294416388,0.4950963057,0.3691658793,0.3748746222,0.4626281569,0.4063055547,0.2840823989,0.3214678523,0.4081352184,0.3232814155,0.2132814404,0.2872819781,0.3556897788,0.2286991341,0.1467766273,0.2666915101,0.2848306773,0.0984067644,0.0909571157,0.2740874756,0.1642612197,-0.0707041278,0.1874682332,0.3403440088,-0.4090941901,0.6193831782,-0.3037904219,0.2522412601,0.0982132719,0.0735850640,0.0479587166,-0.0221863758,0.2357989368,-0.2459409271,0.5486898606,0.5403023059,-2.2136852291,-2.1621324249,-2.1105934031,-2.0591437954,-2.0079258470,-1.9571626539,-1.9070568348,-1.8577662152,-1.8093219706,-1.7615504905,-1.7141047904,-1.6664929504,-1.6180259742,-1.5678271436,-1.5149709461,-1.4585939833,-1.3979097398,-1.3322963852,-1.2615390713,-1.1860421843,-1.1068482271,-1.0255078347,-0.9438285515,-0.8634058758,-0.7850232449,-0.7083490780,-0.6322904616,-0.5558233496,-0.4786885322,-0.4014537101,-0.3249437446,-0.2495196928,-0.1748130515,-0.1001512132,-0.0252612265,0.0494761536,0.1234969670,0.1967092754,0.2695638322,0.3424673883,0.4152791007,0.4875258596,0.5590393648,0.6301601464,0.7012442106,0.7721709905,0.8425405760,0.9122564428,0.9816427613,1.0509483688,1.1199935615,1.1884679694,1.2563932278,1.3240540488,1.3915530092,1.4586929714,1.5253189061,1.5915348573,1.6574819827,1.7231026657,1.7882903863,1.8531289574,1.9177646646,1.9821042710,2.0458791357,2.1090409095,2.1718898832,2.2347280368,2.2975495512,2.3601653051,2.4224334298,2.4841527900,2.5448158339,2.6036764417,2.6600365958,2.7132953460,2.7627749308,2.8077795242,2.8479413570,2.8833809252,2.9145086233,2.9418827524,2.9663535421,2.9891212820,3.0114311389,3.0342354623,3.0581581512,3.0835347881,3.1102644662,3.1377265308,3.1650057223,3.0606231116,2.8529213573,2.5705746868,2.2427515164,1.8982924625,1.5656705118,1.2729583315,1.0476795818,0.9174379553,0.9092974268};
    for (unsigned int j=0; j<v.size(); j++)
    {
        v[j] = v_1[j];
    }

//    double min_step = 1.0;
//    double gold_eps = 0.001;

//    ConjugateGradient cg;
//    cg.setFunction(this);
//    cg.setGradient(this);
//    cg.setEpsilon1(0.001);
//    cg.setEpsilon2(0.001);
//    cg.setEpsilon3(0.001);
//    cg.setR1MinimizeEpsilon(min_step, gold_eps);
//    cg.setPrinter(this);
//    cg.setNormalize(true);
//    cg.showEndMessage(false);
//    cg.calculate(v);

    DoubleVector gr2(v.size());
    gradient(v, gr2);
    gr2.L2Normalize();

    double h = 0.01;
    DoubleVector gr1(v.size());
    IGradient::Gradient(this, h, v, gr1);
    gr1.L2Normalize();



    FILE* file = fopen("20160131_h.txt", "a");
    fprintf(file, "------------------------------------------------------------------------------------------------------------------------\n");
    fprintf(file, "T:%f hx:%f ht:%f M:%d N:%d x:%f X:%d J[v]:%.20f Numeric h: %f\n", t, hx, ht, M, N, xi, Xi, fx(v), h);
    unsigned int vc = (M+D-1);
    fprintf(file, "Controls. Count:%d\n", vc);
    //IPrinter::printVector(v, "v: ", v.size(), 0, v.size()-1, file);
    IPrinter::printVector(v, "v1: ", vc, 0*vc, 0*vc+(vc-1), file);
    IPrinter::printVector(v, "v2: ", vc, 1*vc, 1*vc+(vc-1), file);
    IPrinter::printVector(v, "v3: ", vc, 2*vc, 2*vc+(vc-1), file);
    fprintf(file, "Numerical gradients. Count:%d\n", vc);
    IPrinter::printVector(gr1, "gr1:", vc, 0*vc, 0*vc+(vc-1), file);
    IPrinter::printVector(gr1, "gr2:", vc, 1*vc, 1*vc+(vc-1), file);
    IPrinter::printVector(gr1, "gr3:", vc, 2*vc, 2*vc+(vc-1), file);
    fprintf(file, "Analytic gradient. Count:%d\n", vc);
    IPrinter::printVector(gr2, "gr1:", vc, 0*vc, 0*vc+(vc-1), file);
    IPrinter::printVector(gr2, "gr2:", vc, 1*vc, 1*vc+(vc-1), file);
    IPrinter::printVector(gr2, "gr3:", vc, 2*vc, 2*vc+(vc-1), file);
    fputs("Amplitudes:\n", file);
    DoubleMatrix u;
    pv = &v;
    IHyperbolicEquation::calculateU(u, hx, ht, M+D, N);
    for (unsigned int j=M; j<=M+D; j++)
    {
        char buffer[20];
        int n = sprintf(buffer, "u[%d]: ", j);
        buffer[n] = 0;
        IPrinter::printVector(u[j], buffer, u[j].size(), 0, 0, file);
    }
    fputs("------------------------------------------------------------------------------------------------------------------------\n", file);
    fclose(file);

    double rf = fx(v);
    printf("%.8f %.16f\n", t, rf);
    return rf;
}

double HyperbolicControlH::fx(const DoubleVector &v)
{
    pv = &v;
    DoubleMatrix u;
    IHyperbolicEquation::calculateU(u, hx, ht, M+D, N);

    double sum = 0.0;
    for (unsigned int j=M; j<=M+D; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==0 || i==N || j==M+D || j==0) alpha = 0.5;
            if (i==0 && j==M+D) alpha = 0.25;
            if (i==N && j==M+D) alpha = 0.25;
            sum += alpha*(u[j][i]-U)*(u[j][i]-U);
        }
    }
    sum = hx*ht*sum;

    return sum;
}

void HyperbolicControlH::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleMatrix u;
    IHyperbolicEquation::calculateU(u, hx, ht, M+D, N);

    pu = &u;
    DoubleMatrix p;
    IBackwardHyperbolicEquation::calculateU(p, hx, ht, M+D, N);

    for (unsigned j=2; j<=M+D; j++)
    {
        g[0*(M+D-1)+(j-2)] = -(p[j][1]-p[j][0])/hx;
        g[1*(M+D-1)+(j-2)] = +(p[j][N]-p[j][N-1])/hx;

        double sum = 0.0;
        for (unsigned int i=Xi; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==Xi || i==N) alpha = 0.5;
            sum += alpha*p[j][i];
        }
        g[2*(M+D-1)+(j-2)] = -hx*sum;
    }
    //    IGradient::Gradient(this, 0.01, v, g);
}

double HyperbolicControlH::f(unsigned int i, unsigned int j) const
{
    double sum  = 0.0;
    //double x = i*hx;
    const DoubleVector &v = *pv;
    double v3 = v[2*(M+D-1)+(j-2)];
    //version 1
    if (i>=Xi)
    {
        sum = v3;
    }
    // version 2
    //    if (fabs(x-xi) < (hx+0.000001))
    //    {
    //        sum = (1.0/hx) * v3 * ((hx-fabs(x-xi))/hx);
    //    }
    //version 3
    //    double sgm = 3.0*hx;
    //    double a = 1.0/(sgm*sqrt(2.0*M_PI));
    //    double b = 2.0*sgm*sgm;
    //    double g = a * exp(-((x-e[0])*(x-e[0]))/b);
    //    sum += v3 * g;
    return sum;
}

double HyperbolicControlH::bf(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u = *pu;
    if (M<=j)
    {
        return -(2.0*(u[j][i]-U));
    }
    return 0.0;
}

double HyperbolicControlH::fi1(unsigned int i) const
{
    return 2.0;
}

double HyperbolicControlH::fi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControlH::m1(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v1 = v[0*(M+D-1)+(j-2)];
    return v1;
}

double HyperbolicControlH::m2(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v2 = v[1*(M+D-1)+(j-2)];
    return v2;
}

double HyperbolicControlH::bfi1(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControlH::bfi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControlH::bm1(unsigned int j) const
{
    return 0.0;
}

double HyperbolicControlH::bm2(unsigned int j) const
{
    return 0.0;
}

void HyperbolicControlH::print(unsigned int iteration, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    printf("J[%d]: %.16f\n", iteration, fn->fx(v));
}
