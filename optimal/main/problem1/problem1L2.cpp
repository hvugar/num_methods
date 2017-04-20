#include "problem1L2.h"

void Problem1L2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem1L2 p;
    p.initialize();
    p.startOptimize();
}

Problem1L2::Problem1L2() {}

void Problem1L2::initialize()
{
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;
    withError = false;

    L = 2;

    N = 100;
    hx = 0.001;

    M = 1000;
    ht = 0.1;

    // initial temperatures
    vfi << +0.0;// << -5.0 << +5.0;
    // environment temperatures
    vtt << +0.0;// << -4.0 << +4.0;

    // initial temperature
    fi = 0.0;
    // environment temperature
    tt = 0.0;

    // золото
    {
        const double r0 = 19320.0; // kg/m^3   // плотность
        const double c0 = 130.0;   // C/(kg*S) // удельная теплоемкость
        const double k0 = 312.0;   //          // коэффициент теплопроводности

        const double h1 = 1000.0;      // коэффициент теплообмена ->
        const double h2 = 10.0;        // коэффициент теплообмена ->

        a = sqrt((k0/(c0*r0)));        // коэффициент температуропроворности
        lambda0 = h2/(c0*r0);          // коэффициент теплообмена ->
        lambda1 = h1/k0;               // коэффициент теплообмена ->
        lambda2 = h2/k0;               // коэффициент теплообмена ->
    }

    // мед
    {
        const double r0 = 8900.0;  // kg/m^3     // плотность
        const double c0 = 400.0;   // C/(kg*S)   // удельная теплоемкость
        const double k0 = 380.0;   // Vt/(m*S)   // коэффициент теплопроводности

        const double h1 = 1000.0;      // коэффициент теплообмена ->
        const double h2 = 10.0;        // коэффициент теплообмена ->

        a = sqrt((k0/(c0*r0)));        // коэффициент температуропроворности
        lambda0 = h2/(c0*r0);          // коэффициент теплообмена ->
        lambda1 = h1/k0;               // коэффициент теплообмена ->
        lambda2 = h2/k0;               // коэффициент теплообмена ->

        //printf("%.10f %.10f %.10f %.10f\n", lambda0, lambda1, lambda2, a);
    }

    /* коэффициенты регуляризации */
    alpha0 = 1.0;
    alpha1 = 0.0001;
    alpha2 = 0.0001;
    alpha3 = 0.0001;

    if (!optimizeK) alpha1 = 0.0;
    if (!optimizeZ) alpha2 = 0.0;
    if (!optimizeE) alpha3 = 0.0;

    k0 << 0.0 << 0.0;//-20.57571017653454 << -30.63314593795166;
    z0 << 0.0 << 0.0;//+10.33818417154749 << +10.47968970008047;
    e0 << 0.0 << 0.0;// +0.04500000000000 <<  +0.09500000000000;

    /* шаги числовых производных */
    hk = 0.001;
    hz = 0.001;
    he = 0.001;

    R = 0.0;

    /* пределы z параметров */
    vmin = -100.0;
    vmax = +100.0;
    d0 = (vmax+vmin)/2.0;
    d1 = (vmax-vmin)/2.0;

    /* температура стержня */
    V.resize(N+1);
    for (unsigned int n=0; n<=N; n++) V[n] = 10.0;
}

void Problem1L2::startOptimize()
{
    DoubleVector x0;
    if (optimizeK)
    {
        x0 << -8.5000 << -2.7000; //k
        //x0 << +1.0000 << +1.0000; //k
    }
    else
    {
        K.clear();
        K << -8.5000 << -2.7000; //k
    }

    if (optimizeZ)
    {
        x0 << +2.1000 << +4.9000; //z
        //x0 << +1.0000 << +1.0000; //z
    }
    else
    {
        Z.clear();
        Z << +2.1000 << +4.9000; //z
    }

    if (optimizeE)
    {
        x0 << +0.02000 << +0.08000; //e
    }
    else
    {
        E.clear();
        E << +0.02000 << +0.08000; //e
    }

    R = 1.0;
    optimize(x0);

        while (R < 10000000000.0)
        {
            IPrinter::printSeperatorLine();
            optimize(x0);
            R *= 10.0;
        }
}

void Problem1L2::optimize(DoubleVector &x0) const
{
    Problem1L2* p = const_cast<Problem1L2*>(this);

    ConjugateGradient g;
    g.setFunction(p);
    g.setGradient(p);
    g.setPrinter(p);
    g.setProjection(p);
    g.setEpsilon1(0.0001);//0.00000001
    g.setEpsilon2(0.0001);//0.00000001
    g.setEpsilon3(0.0001);//0.00000001
    g.setR1MinimizeEpsilon(1.0, 0.00001); //0.00000001
    g.setNormalize(true);
    g.showEndMessage(true);
    g.setResetIteration(false);
    g.calculate(x0);

    DoubleMatrix u;
    DoubleVector k,z,e;
    getParameters(k,z,e,x0);

    IPrinter::printSeperatorLine();
    printf("Optimal k: %20.14f %20.14f\n", k[0], k[1]);
    printf("Optimal z: %20.14f %20.14f\n", z[0], z[1]);
    printf("Optimal e: %20.14f %20.14f\n", e[0], e[1]);
    IPrinter::printSeperatorLine();

//    calculateU(u, k, z, e);
//    double aa1 = fx(x0);
//    printf("J: %.10f: ", aa1);
//    IPrinter::printVector(14, 10, u.row(M));
//    IPrinter::printSeperatorLine();

//    p->error = true;
//    for (unsigned int test=0; test<10; test++)
//    {
//        calculateU(u, k, z, e);
//        double aa1 = fx(x0);
//        printf("J: %16.10f: ", aa1);
//        IPrinter::printVector(14, 10, u.row(M));
//        //IPrinter::printSeperatorLine();
//    }
//    for (unsigned int i=0; i<100; i++)
//    {
//    p->M = 10000;
//    calculateU(u,k,z,e);
//    double _MAX_ = 0.0;
//    for (unsigned int m=2000; m<=M; m++)
//    {
//        double sum = 0.0;
//        sum += 0.5*mu(0)*(u[m][0]-V[0])*(u[m][0]-V[0]);
//        for (unsigned int n=1; n<=N-1; n++)
//        {
//            sum += mu(n)*(u[m][n]-V[n])*(u[m][n]-V[n]);
//        }
//        sum += 0.5*mu(N)*(u[m][N]-V[N])*(u[m][N]-V[N]);
//        sum *= hx;

//        if (sum > _MAX_) _MAX_ = sum;
//        //printf("M: %u %.20f %.20f\n",m,sum,_MAX_);
//    }
//    printf("%.20f\n",_MAX_);
//    }
}

void Problem1L2::print(unsigned int i, const DoubleVector &prm, const DoubleVector &g, double r, GradientMethod::MethodResult result) const
{
    C_UNUSED(result);
    //if (i>0) return;
    //if (i % 20 != 0 && result < 3) return;

    Problem1L2 *pm = const_cast<Problem1L2*>(this);
    pm->py = &prm;

    DoubleVector k,z,e;
    getParameters(k,z,e,prm);

    DoubleMatrix u;
    calculateU(u,k,z,e);

    //    if (result == GradientMethod::BREAK_DISTANCE_LESS || result == GradientMethod::BREAK_DISTANCE_LESS || result == GradientMethod::BREAK_GRADIENT_NORM_LESS)
    //    {
    //        DoubleVector v(M+1);
    //        for (unsigned int m=0; m<=M; m++) v[m] = vf(m,k,z,e,u);
    //        FILE *file1 = fopen("v_0.txt", "w");
    //        IPrinter::printVector(14,10,v,NULL,v.size(),0,0,file1);
    //        fclose(file);
    //        IPrinter::printVector(14,10,u.row(u.rows()-1),"u: ", 10, 0, 0, stdout);
    //    }

    //

    if (result == GradientMethod::FIRST_ITERATION)
    {
        FILE *file = fopen("control_v.txt", "w");
        fprintf(file, "%u", i);
        for (unsigned int m=0; m<=M; m++)
        {
            fprintf(file, ",%.10f",vf(m,k,z,e,u));
        }
        fprintf(file, "\n");
        fclose(file);
    }
    else
    {
        FILE *file = fopen("control_v.txt", "a");
        fprintf(file, "%u", i);
        for (unsigned int m=0; m<=M; m++)
        {
            fprintf(file, ",%.10f",vf(m,k,z,e,u));
        }
        fprintf(file, "\n");
        fclose(file);
    }

    double v = 0.0;//vf(M, k, z, e, u);
    //IPrinter::printSeperatorLine(NULL,'-',file);
    //IPrinter::printSeperatorLine(NULL,'-',stdout);

    //fprintf(file,"\n");
    //fprintf(file, "J[%d]: %.10f v: %.10f\n", i, r, v);

    unsigned int p=0;
    //fprintf(stdout, "---\n");
    //fprintf(file, "k: %14.10f %14.10f\n", k[0], k[1]);

    //if ( result == GradientMethod::FIRST_ITERATION) pm->R = RMAX;

    //printf("R: %.14f\n", pm->R);

    printf("%d,%.6f",i,r);

    //DoubleVector a1(prm.size());
    //pm->gradient(prm, a1);

    if (optimizeK)
    {
        DoubleVector a = g.mid(p,p+1);
        DoubleVector na = a;
        na.L2Normalize();

        DoubleVector n = a;

        DoubleVector cx = prm;
        double f1,f2;

        double x0 = prm[p];
        cx[p] = x0 - hk; f1 = fx(cx);
        cx[p] = x0 + hk; f2 = fx(cx);
        n[0] = (f2-f1)/(2.0*hk);
        cx[p] = x0;

        double x1 = prm[p+1];
        cx[p+1] = x1 - hk; f1 = fx(cx);
        cx[p+1] = x1 + hk; f2 = fx(cx);
        n[1] = (f2-f1)/(2.0*hk); cx[p+1] = x1;

        DoubleVector nn = n;
        nn.L2Normalize();

        //fprintf(file, "a: %14.10f %14.10f | %14.10f %14.10f\n", a[0], a[1], na[0], na[1]);
        //fprintf(file, "n: %14.10f %14.10f | %14.10f %14.10f\n", n[0], n[1], nn[0], nn[1]);

        p+=2;

        printf("|%.6f,%.6f,%.6f",k[0],a[0],n[0]);
        printf("|%.6f,%.6f,%.6f",k[1],a[1],n[1]);
    }

    //fprintf(file, "---\n");
    //fprintf(file, "z: %14.10f %14.10f\n", z[0], z[1]);
    if (optimizeZ)
    {
        DoubleVector a = g.mid(p,p+1);
        DoubleVector na = a;
        na.L2Normalize();

        DoubleVector n = a;

        DoubleVector cx = prm;
        double f1,f2;

        double x0 = prm[p];
        cx[p] = x0 - hz; f1 = fx(cx);
        cx[p] = x0 + hz; f2 = fx(cx);
        n[0] = (f2-f1)/(2.0*hz); cx[p] = x0;

        double x1 = prm[p+1];
        cx[p+1] = x1 - hz; f1 = fx(cx);
        cx[p+1] = x1 + hz; f2 = fx(cx);
        n[1] = (f2-f1)/(2.0*hz); cx[p+1] = x1;

        DoubleVector nn = n;
        nn.L2Normalize();

        //fprintf(file, "a: %14.10f %14.10f | %14.10f %14.10f\n", a[0], a[1], na[0], na[1]);
        //fprintf(file, "n: %14.10f %14.10f | %14.10f %14.10f\n", n[0], n[1], nn[0], nn[1]);
        \
        p+=2;
        printf("|%.6f,%.6f,%.6f",z[0],a[0],n[0]);
        printf("|%.6f,%.6f,%.6f",z[1],a[1],n[1]);
    }

    //fprintf(file, "---\n");
    //fprintf(file, "e: %14.10f %14.10f\n", e[0], e[1]);
    if (optimizeE)
    {
        DoubleVector a = g.mid(p,p+1);
        DoubleVector na = a;
        na.L2Normalize();

        DoubleVector n = a;

        DoubleVector cx = prm;
        double f1,f2;

        double x0 = prm[p];
        cx[p] = x0 - he; f1 = fx(cx);
        cx[p] = x0 + he; f2 = fx(cx);
        n[0] = (f2-f1)/(2.0*he); cx[p] = x0;

        double x1 = prm[p+1];
        cx[p+1] = x1 - he; f1 = fx(cx);
        cx[p+1] = x1 + he; f2 = fx(cx);
        n[1] = (f2-f1)/(2.0*he); cx[p+1] = x1;

        DoubleVector nn = n;
        nn.L2Normalize();

        //fprintf(file, "a: %14.10f %14.10f | %14.10f %14.10f\n", a[0], a[1], na[0], na[1]);
        //fprintf(file, "n: %14.10f %14.10f | %14.10f %14.10f\n", n[0], n[1], nn[0], nn[1]);
        p+=2;
        printf("|%.6f,%.6f,%.6f",e[0],a[0],n[0]);
        printf("|%.6f,%.6f,%.6f",e[1],a[1],n[1]);
    }

    printf("\n");

    //IPrinter::printVector(14,10,u.row(u.rows()-1),"u: ", 10, 0, 0, stdout);
}

void Problem1L2::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
{
    unsigned int p = 0;
    if (optimizeK) p+=2;

    if (optimizeZ) p+=2;

    /* z lower/upper limits */
    //    if (x.at(2) < zmin) x.at(2) = zmin;
    //    if (x.at(3) < zmin) x.at(3) = zmin;
    //    if (x.at(2) > zmax) x.at(2) = zmax;
    //    if (x.at(3) > zmax) x.at(3) = zmax;

    //    if (x.at(4) < 0.10) x.at(4) = 0.10;
    //    //if (x.at(4) > 0.90) x.at(4) = 0.90;
    //    if (x.at(4) > 0.40) x.at(4) = 0.40;

    //    //if (x.at(5) < 0.10) x.at(5) = 0.10;
    //    if (x.at(5) > 0.90) x.at(5) = 0.90;
    //    if (x.at(5) < 0.60) x.at(5) = 0.60;

    /* e lower/upper limits */
    if (optimizeE)
    {
        if (x.at(p) < 5*hx)  x.at(p) = 5*hx;
        if (x.at(p) > 0.045) x.at(p) = 0.045;

        if (x.at(p+1) < 0.055)     x.at(p+1) = 0.055;
        if (x.at(p+1) > (N-5)*hx)  x.at(p+1) = (N-5)*hx;
    }
}
