#include "art_problem1L2.h"

#define RMAX 0.0

void ArtProblem1L2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ArtProblem1L2 p;
    p.initialize();
    p.startOptimize();
}

ArtProblem1L2::ArtProblem1L2() {}

void ArtProblem1L2::initialize()
{
    optimizeK = true;
    optimizeZ = false;
    optimizeE = false;
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
    alpha1 = 0.0000;
    alpha2 = 0.0000;
    alpha3 = 0.0000;

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
    zmin = 0.0;
    zmax = 0.0;

    /* пределы z параметров */
    vmin = -25.0;
    vmax = +20.0;
    d0 = (vmax+vmin)/2.0;
    d1 = (vmax-vmin)/2.0;

    /* температура стержня */
    V.resize(N+1);
    for (unsigned int n=0; n<=N; n++) V[n] = 10.0;
}

void ArtProblem1L2::startOptimize()
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

    R = 100000000.0;
    optimize(x0);
//    while (R < 10000000000.0)
//    {
//        IPrinter::printSeperatorLine();
//        optimize(x0);
//        R *= 10.0;
//    }

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

void ArtProblem1L2::optimize(DoubleVector &x0) const
{
    ArtProblem1L2* p = const_cast<ArtProblem1L2*>(this);
    ConjugateGradient g;
    g.setFunction(p);
    g.setGradient(p);
    g.setPrinter(p);
    g.setProjection(p);
    g.setEpsilon1(0.0000000001);//0.00000001
    g.setEpsilon2(0.0000000001);//0.00000001
    g.setEpsilon3(0.0000000001);//0.00000001
    g.setR1MinimizeEpsilon(1.0, 0.00001); //0.00000001
    g.setNormalize(true);
    g.showEndMessage(true);
    g.setResetIteration(false);
    g.calculate(x0);
}

double ArtProblem1L2::fx(const DoubleVector &prm) const
{
    ArtProblem1L2* pm = const_cast<ArtProblem1L2*>(this);
    pm->px = &prm;

    DoubleVector k,z,e;
    getParameters(k,z,e,prm);

    unsigned int N1 = vfi.size();
    unsigned int N2 = vtt.size();
    double SUM = 0.0;
    for (unsigned int n1=0; n1<N1; n1++)
    {
        for (unsigned int n2=0; n2<N2; n2++)
        {
            pm->fi = vfi[n1];
            pm->tt = vtt[n2];

            DoubleMatrix u;
            calculateU(u, k, z, e);

            double sum = 0.0;
            sum += 0.5*mu(0)*(u[M][0]-V[0])*(u[M][0]-V[0]);
            for (unsigned int n=1; n<=N-1; n++)
            {
                sum += mu(n)*(u[M][n]-V[n])*(u[M][n]-V[n]);
            }
            sum += 0.5*mu(N)*(u[M][N]-V[N])*(u[M][N]-V[N]);
            sum *= hx;

            double pnlt = 0.0;
            double min = 0.0;

            min = fmin(0.0, gf(0, k, z, e, u));
            pnlt += 0.5*min*min;
            for (unsigned int m=1; m<=M-1; m++)
            {
                min = fmin(0.0, gf(m, k, z, e, u));
                pnlt += min*min;
            }
            min = fmin(0.0, gf(M, k, z, e, u));
            pnlt += 0.5*min*min;

            pnlt *= ht;

            SUM += alpha0*sum + R*pnlt;
        }
    }
    SUM *= ((1.0/N1)*(1.0/N2));

    double norm1 = 0.0;
    double norm2 = 0.0;
    double norm3 = 0.0;

    norm1 = (k[0] - k0[0])*(k[0] - k0[0]) + (k[1] - k0[1])*(k[1] - k0[1]);
    norm2 = (z[0] - z0[0])*(z[0] - z0[0]) + (z[1] - z0[1])*(z[1] - z0[1]);
    norm3 = (e[0] - e0[0])*(e[0] - e0[0]) + (e[1] - e0[1])*(e[1] - e0[1]);

    SUM += alpha1*norm1 + alpha2*norm2 + alpha3*norm3;

    return SUM;
}

void ArtProblem1L2::gradient(const DoubleVector &prm, DoubleVector &g)
{
    px = &prm;

    DoubleVector k,z,e;
    getParameters(k,z,e,prm);

    for (unsigned int i=0; i<g.size(); i++) g[i] = 0.0;

    unsigned int N1 = vfi.size();
    unsigned int N2 = vtt.size();

    for (unsigned int n1=0; n1<N1; n1++)
    {
        for (unsigned int n1=0; n1<N2; n1++)
        {
            fi = vfi[n1];
            tt = vtt[n1];

            DoubleMatrix u;
            calculateU(u, k, z, e);

            DoubleMatrix p;
            calculateP(p, u, k, z, e);

            unsigned int i = 0;

            // k gradient
            if (optimizeK)
            {
                for (unsigned int s=0; s<L; s++)
                {
                    unsigned int xi = (unsigned int)round(e[s] * N*10);

                    // Integral part of gradient
                    double sum = 0.0;
                    sum += 0.5*p.at(0, 0)*(u.at(0, xi) - z[s]);
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        sum += p[m][0]*(u[m][xi] - z[s]);
                    }
                    sum += 0.5*p[M][0]*(u[M][xi] - z[s]);
                    sum *= ht;

                    // Penalty part of gradient
                    double pnlt = 0.0;
                    pnlt += 0.5 * (u[0][xi]-z[s]) * sgn_min(0, k, z, e, u);
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        pnlt += (u[m][xi]-z[s]) * sgn_min(m, k, z, e, u);
                    }
                    pnlt += 0.5 * (u[M][xi]-z[s]) * sgn_min(M, k, z, e, u);
                    pnlt *= ht;

                    g[i] += -lambda1*a*a*sum + 2.0*R*pnlt;
                    i++;
                }
            }

            // z gradient
            if (optimizeZ)
            {
                for (unsigned int s=0; s<L; s++)
                {
                    // Integral part of gradient
                    double sum = 0.0;
                    sum += 0.5*p[0][0];
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        sum += p[m][0];
                    }
                    sum += 0.5*p[M][0];
                    sum *= ht;

                    // Penalty part of gradient
                    double pnlt = 0.0;
                    pnlt += 0.5 * sgn_min(0, k, z, e, u);
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        pnlt += sgn_min(m, k, z, e, u);
                    }
                    pnlt += 0.5 * sgn_min(M, k, z, e, u);
                    pnlt *= ht;

                    g[i] += lambda1*a*a*k[s]*sum - 2.0*R*k[s]*pnlt;
                    i++;
                }
            }

            // e gradient
            if (optimizeE)
            {
                for (unsigned int s=0; s<L; s++)
                {
                    unsigned int xi = (unsigned int)round(e.at(s) * N*10);

                    // Integral part of gradient
                    double sum = 0.0;
                    sum += 0.5 * p[0][0] * ((u.at(0, xi+1) - u.at(0, xi-1))/(2.0*hx));
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        sum += p[m][0] * ((u.at(m, xi+1) - u.at(m, xi-1))/(2.0*hx));
                    }
                    sum += 0.5 * p[M][0] * ((u.at(M, xi+1) - u.at(M, xi-1))/(2.0*hx));
                    sum *= ht;

                    // Penalty part of gradient
                    double pnlt = 0.0;
                    pnlt += 0.5 * k[s] * ((u.at(0, xi+1) - u.at(0, xi-1))/(2.0*hx)) * sgn_min(0, k, z, e, u);
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        pnlt += ((u.at(m, xi+1) - u.at(m, xi-1))/(2.0*hx)) * sgn_min(m, k, z, e, u);
                    }
                    pnlt += 0.5 * ((u.at(M, xi+1) - u.at(M, xi-1))/(2.0*hx)) * sgn_min(M, k, z, e, u);
                    pnlt *= ht;

                    g[i] += -lambda1*a*a*k[s]*sum + 2.0*R*k[s]*pnlt;
                    i++;
                }
            }
        }
    }

    for (unsigned int i=0; i<g.size(); i++) g[i] *= (1.0/N1)*(1.0/N2);

    unsigned int i = 0;
    if (optimizeK)
    {
        for (unsigned int s=0; s<L; s++)
        {
            g[i] += 2.0*alpha1*(k[s]-k0[s]);
            i++;
        }
    }
    if (optimizeZ)
    {
        for (unsigned int s=0; s<L; s++)
        {
            g[i] += 2.0*alpha2*(z[s]-z0[s]);
            i++;
        }
    }
    if (optimizeE)
    {
        for (unsigned int s=0; s<L; s++)
        {
            g[i] += 2.0*alpha3*(e[s]-e0[s]);
            i++;
        }
    }
}

void ArtProblem1L2::calculateU(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const
{
    u.clear();
    u.resize(M+1, N+1);

    unsigned int E0 = (unsigned int)round(e[0] * N*10);
    unsigned int E1 = (unsigned int)round(e[1] * N*10);

    double *da = (double*) malloc(sizeof(double)*(N+1));
    double *db = (double*) malloc(sizeof(double)*(N+1));
    double *dc = (double*) malloc(sizeof(double)*(N+1));
    double *dd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *de = (double*) malloc(sizeof(double)*(N+1));

    for (unsigned int n=0; n<=N; n++) u[0][n] = initial(n);

    for (unsigned int m=1; m<=M; m++)
    {
        // n = 0
        da[0] = 0.0;
        db[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambda1*a*a*ht)/hx + lambda0*ht;
        dc[0] = -(a*a*ht)/(hx*hx);
        dd[0] = u.at(m-1,0) + lambda0*ht*tt - ((lambda1*a*a*ht)/hx)*(k[0]*z[0] + k[1]*z[1]);

        // n = 1,...,N-1
        for (unsigned int n=1; n<=N-1; n++)
        {
            da[n] = -(a*a*ht)/(hx*hx);
            db[n] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
            dc[n] = -(a*a*ht)/(hx*hx);
            dd[n] = u.at(m-1,n) + lambda0*ht*tt;
        }

        // n = N
        da[N] = -(a*a*ht)/(hx*hx);
        db[N] = 1.0 + (a*a*ht)/(hx*hx) + (lambda2*a*a*ht)/hx + lambda0*ht;
        dc[N] = 0.0;
        dd[N] = u.at(m-1,N)  + lambda0*ht*tt + (lambda2*a*a*ht*tt)/hx;

        de[0]  =de[1] = 0.0;
        de[N-1]=de[N] = 0.0;
        for (unsigned int n=2; n<=N-2; n++)
        {
            de[n] = 0.0;

            double dif0 = fabs(n*hx - e[0]);
            if (dif0 <= hx)
            {
                de[n] = -((lambda1*(a*a*ht)/hx) * k[0]) * (1.0 - dif0/hx);
            }

            double dif1 = fabs(n*hx - e[1]);
            if (dif1 <= hx)
            {
                de[n] = -((lambda1*(a*a*ht)/hx) * k[1]) * (1.0 - dif1/hx);
            }
        }

        qovmaFirstRowM(da, db, dc, dd, rx, N+1, de);

        for (unsigned int i=0; i<=N; i++)
        {
            u[m][i] = rx[i];

            if (withError)
            {
                u[m][E0] += (((rand()%RAND_MAX) % 2 == 0) ? +0.1 : -0.1) * u[m][E0];
                u[m][E1] += (((rand()%RAND_MAX) % 2 == 0) ? +0.1 : -0.1) * u[m][E1];
            }
        }
    }

    free(de);
    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

void ArtProblem1L2::calculateP(DoubleMatrix &p, const DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e)
{
    p.clear();
    p.resize(M+1, N+1);

    double *da = (double*) malloc(sizeof(double)*(N+1));
    double *db = (double*) malloc(sizeof(double)*(N+1));
    double *dc = (double*) malloc(sizeof(double)*(N+1));
    double *dd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *de = (double*) malloc(sizeof(double)*(N+1));

    for (unsigned int n=0; n<=N; n++)
    {
        p[M][n] = -2.0*alpha0*mu(n)*(u[M][n] - V[n]);
    }

    for (unsigned int m=M-1; m != UINT32_MAX; m--)
    {
        // n = 0
        da[0] = 0.0;
        db[0] = -1.0 - (a*a*ht)/(hx*hx) - (lambda1*a*a*ht)/hx - lambda0*ht;
        dc[0] = (a*a*ht)/(hx*hx);
        dd[0] = -p[m+1][0];

        // n = 1,...,N-1
        for (unsigned int n=1; n<=N-1; n++)
        {
            da[n] = (a*a*ht)/(hx*hx);
            db[n] = -1.0-(2.0*a*a*ht)/(hx*hx) - lambda0*ht;
            dc[n] = (a*a*ht)/(hx*hx);
            dd[n] = -p[m+1][n] + R * ht * sgn_min(m, k, z, e, u);
            for (unsigned int i=0; i<L; i++)
            {
                dd[n] +=  R * ht * 2.0*sgn_min(m, k, z, e, u) * k[i] * delta(n,e,i);
            }
        }

        // n = N
        da[N] = (a*a*ht)/(hx*hx);
        db[N] = -1.0-(a*a*ht)/(hx*hx) - (lambda2*a*a*ht)/hx - lambda0*ht;
        dc[N] = 0.0;
        dd[N] = -p[m+1][N];

        de[0] = de[1] = 0.0;
        de[N] = de[N-1] = 0.0;
        for (unsigned int n=2; n<=N-2; n++)
        {
            de[n] = 0.0;
            for (unsigned int i=0; i<L; i++)
            {
                de[n] += k[i] * delta(n,e,i);
            }
            de[n] *= (lambda1*a*a*ht);
        }

        qovmaFirstColM(da, db, dc, dd, rx, N+1, de);

        for (unsigned int i=0; i<=N; i++) p[m][i] = rx[i];
    }

    free(de);
    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

double ArtProblem1L2::delta(unsigned int n, const DoubleVector &e, unsigned int i) const
{
    double sigma = 3.0*hx;
    double x = n*hx;
    return 1.0/(sqrt(2.0*M_PI)*sigma) * exp(-((x-e[i])*(x-e[i]))/(2.0*sigma*sigma));

    return 0.0;
}

double ArtProblem1L2::initial(unsigned int n UNUSED_PARAM) const
{
    return fi;
}

double ArtProblem1L2::vf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const
{
    unsigned int E0 = (unsigned int)round(e[0] * N*10);
    unsigned int E1 = (unsigned int)round(e[1] * N*10);

    return k[0]*(u[m][E0]-z[0]) + k[1]*(u[m][E1]-z[1]);
}

void ArtProblem1L2::print(unsigned int i, const DoubleVector &prm, const DoubleVector &g, double r, GradientMethod::MethodResult ) const
{
    ArtProblem1L2* pm = const_cast<ArtProblem1L2*>(this);

    DoubleVector k,z,e;
    getParameters(k,z,e,prm);

    IPrinter::printSeperatorLine();
    printf("J[%d]: %.10f\n", i, r);

    DoubleVector cx = prm;

    unsigned int p = 0;
    if (optimizeK)
    {
        DoubleVector w(3*L);
        pm->gradient(prm, w);

//        if (i==1)
//        {
//            FILE *file = fopen("sample1.txt", "w");
//            DoubleMatrix u;
//            calculateU(u, k, z, e);
//            for (unsigned int m=0; m<=M; m++) fprintf(file,"%.10f ", vf(m,k,z,e,u));
//            fprintf(file,"\n");

//            double temp = k[0];
//            k[0] = temp + hk;
//            calculateU(u, k, z, e);
//            for (unsigned int m=0; m<=M; m++) fprintf(file,"%.10f ", vf(m,k,z,e,u));
//            fprintf(file,"\n");

//            k[0] = temp - hk;
//            calculateU(u, k, z, e);
//            for (unsigned int m=0; m<=M; m++) fprintf(file,"%.10f ", vf(m,k,z,e,u));
//            fprintf(file,"\n");
//            k[0] = temp;
//            fclose(file);
//        }

        printf("%d ----------------------------------------------------------------------------------------------------------------------------------\n", i);

        DoubleVector a = w.mid(p,p+1);
        DoubleVector n = a;

        double f1,f2;
        cx = prm;
        double x0 = cx[p];

        printf("%d ----------------------------------------------------------------------------------------------------------------------------------\n", i);

        cx[p] = x0 - hk;
        f1 = fx(cx);

        cx[p] = x0 + hk;
        f2 = fx(cx);

        n[0] = (f2-f1)/(2.0*hk);


        printf("%d ----------------------------------------------------------------------------------------------------------------------------------\n", i);

        cx = prm;
        double x1 = cx[p+1];
        cx[p+1] = x1 - hk; f1 = fx(cx);
        cx[p+1] = x1 + hk; f2 = fx(cx);
        n[1] = (f2-f1)/(2.0*hk);

        printf("k0: %.10f a: %.10f n: %.10f\n", k[0], a[0], n[0]);
        printf("k1: %.10f a: %.10f n: %.10f\n", k[1], a[1], n[1]);
        a.L2Normalize();
        n.L2Normalize();
        printf("k0: %.10f a: %.10f n: %.10f\n", k[0], a[0], n[0]);
        printf("k1: %.10f a: %.10f n: %.10f\n", k[1], a[1], n[1]);
        p+=2;
    }
}
/*
void ArtProblem1L2::print(unsigned int i, const DoubleVector &prm, const DoubleVector &g, double r, GradientMethod::MethodResult result) const
{
    C_UNUSED(result);
    //if (i>0) return;
    //if (i % 20 != 0 && result < 3) return;

    ArtProblem1L2 *pm = const_cast<ArtProblem1L2*>(this);
    pm->px = &prm;

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
*/

void ArtProblem1L2::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
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

void ArtProblem1L2::getParameters(DoubleVector &k, DoubleVector &z, DoubleVector &e, const DoubleVector &prm) const
{
    unsigned int p = 0;

    if (optimizeK)
    {
        k = prm.mid(p,p+1);
        p+=L;
    }
    else
    {
        k = this->K;
    }

    if (optimizeZ)
    {
        z = prm.mid(p,p+1);
        p+=L;
    }
    else
    {
        z = this->Z;
    }

    if (optimizeE)
    {
        e = prm.mid(p,p+1);
        p+=L;
    }
    else
    {
        e = this->E;
    }
}

void ArtProblem1L2::qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);
    double *k = (double*)malloc(sizeof(double)*n);

    for (unsigned int i=n-1; i != UINT32_MAX; i--)
    {
        if (i == n-1)
        {
            p[i] = -a[i]/b[i];
            q[i] = +d[i]/b[i];
            k[i] = -e[i]/b[i];
        }
        else if (i == 1)
        {
            double m = b[i]+c[i]*p[i+1];
            p[i] = -(a[i]+c[i]*k[i+1])/m;
            q[i] = +(d[i]-c[i]*q[i+1])/m;
            k[i] = 0.0;
        }
        else if (i == 0)
        {
            double m = b[i]+c[i]*p[i+1];
            p[i] = 0.0;
            q[i] = +(d[i]-c[i]*q[i+1])/m;
            k[i] = 0.0;
        }
        else
        {
            double m = b[i]+c[i]*p[i+1];
            p[i] = -a[i]/m;
            q[i] = +(d[i]-c[i]*q[i+1])/m;
            k[i] = -(e[i]+c[i]*k[i+1])/m;
        }
    }

    for (unsigned int i=0; i<n; i++)
    {
        if (i==0)
        {
            x[i] = q[i];
        }
        else if (i==1)
        {
            x[i] = p[i]*x[i-1] + q[i];
        }
        else
        {
            x[i] = p[i]*x[i-1] + q[i] + k[i]*x[0];
        }
    }

    free(k);
    free(p);
    free(q);
}

void ArtProblem1L2::qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);

    unsigned int L = 0;
    for (unsigned int s=0; s<n; s++)
    {
        if (fabs(e[s]) != 0.0)
        {
            L+=1;
        }
    }
    unsigned int *E = (unsigned int *)malloc(sizeof(unsigned int)*L);

    unsigned int i = 0;
    for (unsigned int s=0; s<n; s++)
    {
        if (fabs(e[s]) != 0.0)
        {
            E[i++] = s;
        }
    }

    double **k = (double**) malloc(sizeof(double*)*L);
    for (unsigned int s=0; s<L; s++) k[s] = (double*)malloc(sizeof(double)*n);

    for (unsigned int i=0; i<n; i++)
    {
        if (i == 0)
        {
            p[0] = +d[0]/b[0];
            q[0] = -c[0]/b[0];

            for (unsigned int s=0; s<L; s++)
            {
                k[s][0] = -e[E[s]]/b[0];
            }
        }
        else if (i == (n-1))
        {
            p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
            q[i] = 0.0;

            for (unsigned int s=0; s<L; s++) k[s][i] = 0.0;
        }
        else
        {
            double m = b[i]+a[i]*q[i-1];
            p[i] = +(d[i]-a[i]*p[i-1])/m;
            q[i] = -c[i]/m;

            for (unsigned int s=0; s<L; s++)
            {
                if (i<(E[s]-1))
                    k[s][i] = -(a[i]*k[s][i-1])/m;
                else
                    k[s][i] = 0.0;

                if (i==E[s]-1) q[i] += -(a[i]*k[s][i-1])/m;
            }

            //            for (unsigned int s=0; s<L; s++)
            //            {
            //                if (i==E[s]-1) q[i] += -(a[i]*k[s][i-1])/m;
            //            }
        }
    }

    for (unsigned int i=n-1; i != UINT_MAX; i--)
    {
        if (i==(n-1))
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];

            for (unsigned int s=0; s<L; s++)
            {
                if (i<=E[s]-1)
                {
                    x[i] = x[i] + k[s][i]*x[E[s]];
                }
            }
        }
    }

    for (unsigned int s=0; s<L; s++) free(k[s]);
    free(k);
    free(E);
    free(q);
    free(p);
}

double ArtProblem1L2::gf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const
{
    return d1 - fabs(vd0(m, k, z, e, u));
}

double ArtProblem1L2::vd0(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const
{
    return d0 - vf(m, k, z, e, u);
}

double ArtProblem1L2::sgn_min(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const
{

    return sign(vd0(m, k, z, e, u)) * fmin(0.0, gf(m, k, z, e, u));
}

double ArtProblem1L2::sign(double x) const
{
    if (x < 0.0) return -1.0;
    if (x > 0.0) return +1.0;
    return 0.0;
}
