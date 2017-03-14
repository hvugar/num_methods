#include "problem1L2.h"

void Problem1L2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    //    for (unsigned int i=500; i<=3000; i+=100)
    //    {
    //        Problem1L2 p;
    //        double r = p.fx(i);
    //        printf("%d %18.10f\n", i, r);
    //    }


    Problem1L2 p;
    p.fx(500);
}

Problem1L2::Problem1L2()
{
    optimizeK = true;
    optimizeZ = true;
    optimizeE = false;

    L = 2;

    N = 100;
    hx = 0.001;

    ht = 0.1;

    //    N = 100;
    //    M = 50*100;
    //    hx = 0.01;
    //    ht = 0.01;
    //    h  = 0.001;

    // initial temperature
    fi = 0.0;
    // environment temperature
    tt = 0.0;

    //    const double r0 = 19320.0; // kg/m^3   // плотность
    //    const double c0 = 130.0;   // C/(kg*S) // удельная теплоемкость
    //    const double k0 = 312.0;   //          // коэффициент теплопроводности

    const double r0 = 8900.0; // kg/m^3      // плотность
    const double c0 = 400.0;   // C/(kg*S)   // удельная теплоемкость
    const double k0 = 380.0;   // Vt/(m*S)   // коэффициент теплопроводности

    const double h1 = 1000.0;      // коэффициент теплообмена ->
    const double h2 = 10.0;        // коэффициент теплообмена ->

    a = sqrt((k0/(c0*r0)));  // коэффициент температуропроворности
    lambda0 = h2/(c0*r0);          // коэффициент теплообмена ->
    lambda1 = h1/k0;               // коэффициент теплообмена ->
    lambda2 = h2/k0;               // коэффициент теплообмена ->

    //    a = 1.0;
    //    alpha = 0.01;
    //    lambda0 = 10.0;
    //    lambdal = 1.0;

    alpha0 = 1.0;
    alpha1 = 0.001;
    alpha2 = 0.001;
    alpha3 = 0.000;
    R = 1.0;

    hk = 0.0001;
    hz = 0.0001;
    he = 0.001;

    zmin = 9.5;
    zmax = 10.5;

    vmin = -10.0;
    vmax = 1000.0;
    d0 = (vmax-vmin)/2.0;
    d1 = (vmax-vmin)/2.0;

    V.resize(N+1);

    /* init V */
    for (unsigned int n=0; n<=N; n++)
    {
        //double h1 = 0.4/N;
        V[n] = 10.0;//4.2 - n*h1;
    }
    IPrinter::printVector(14, 10, V,"V: ");

    fprintf(stdout, "alpha0: %f alpha1: %f alpha2: %f alpha3: %f\n", alpha0, alpha1, alpha2, alpha3);
    fprintf(stdout, "ro: %.10f, c: %.10f k %.10f h1 %.10f h2 %.10f a %.10f lambda0 %.10f lambda1 %.10f lambda2 %.10f\n", r0, c0, k0, h1, h2, a, lambda0, lambda1, lambda2);
}

double Problem1L2::fx(double t) const
{
    Problem1L2* p = const_cast<Problem1L2*>(this);
    p->file = fopen("problem1_test2.txt", "w");
    p->M = (unsigned int)(t);

    fprintf(file, "hx: %f N: %d ht: %f M: %d L: %d\n", hx, N, ht, M, L);

    DoubleVector x0;
    if (optimizeK)
    {
        fprintf(file, "Optimizing k parameters\n");
        x0 << -8.5000 << -2.7000; //k
        //x0 << -4.8819170104 << 3.0421205175;
    }
    else
    {
        p->K.clear();
        p->K << -8.5000 << -2.7000; //k
    }

    if (optimizeZ)
    {
        fprintf(file, "Optimizing z parameters\n");
        x0 << +2.1000 << +4.9000; //z
        //x0 << 8.2918856107 << -0.4531024453;
    }
    else
    {
        p->z.clear();
        p->z << +2.1000 << +4.9000; //z
    }

    if (optimizeE)
    {
        fprintf(file, "Optimizing e parameters\n");
        p->e.clear();
        x0 << +0.02000 << +0.08000; //e
    }
    else
    {
        p->e.clear();
        p->e << +0.02000 << +0.08000; //e
    }

    //k << 3.50 << 3.70;
    //z << 4.20 << 4.20;
    //e << 0.40 << 0.90;

    ConjugateGradient g;
    g.setFunction(p);
    g.setGradient(p);
    g.setPrinter(p);
    g.setProjection(p);
    g.setEpsilon1(0.0001);//0.00000001
    g.setEpsilon2(0.0001);//0.00000001
    g.setEpsilon3(0.0001);//0.00000001
    g.setR1MinimizeEpsilon(0.1, 0.0001); //0.00000001
    g.setNormalize(true);
    g.showEndMessage(false);
    g.calculate(x0);

    //    FILE *file = fopen("data_z.txt", "w");
    //    for (unsigned int n=2; n<=N; n++)
    //    {
    //        DoubleVector x1 = x0;
    //        x1[1] = 1.0+0.2*n;
    //        double rx = fx(x1);
    //        fprintf(file, "%d,%18.14f\n", n, rx);
    //    }
    //    fclose(file);

    //    DoubleVector v(M+1);
    //    DoubleVector k,z,e;
    //    getComponents(k,z,e,x0);
    //    p->px = &x0;

    //    DoubleMatrix u;
    //    calculateU(u);
    //    unsigned int e0 = (unsigned int)round(e[0] * p->N*10);
    //    unsigned int e1 = (unsigned int)round(e[1] * p->N*10);
    //    FILE *file = fopen("data_v.txt", "w");

    //    for (unsigned int m=0; m<=p->M; m++)
    //    {
    //        v[m] = k[0]*(u[m][e0]-z[0])+k[1]*(u[m][e1]-z[1]);
    //        fprintf(file, "%18.14f\n", v[m]);
    //    }

    fclose(p->file);
    return fx(x0);
}

double Problem1L2::fx(const DoubleVector &prm) const
{
    const_cast<Problem1L2*>(this)->px = &prm;
    DoubleMatrix u;
    calculateU(u);

    return integral(prm, u) + norm(prm) + R*penalty(prm, u);
}

double Problem1L2::integral(const DoubleVector &prm, const DoubleMatrix &u) const
{
    double sum = 0.0;
    sum += 0.5*mu(0)*(u[M][0]-V[0])*(u[M][0]-V[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += mu(n)*(u[M][n]-V[n])*(u[M][n]-V[n]);
    }
    sum += 0.5*mu(N)*(u[M][N]-V[N])*(u[M][N]-V[N]);
    sum *= hx;
    return alpha0*sum;
}

double Problem1L2::norm(const DoubleVector &prm) const
{
    const_cast<Problem1L2*>(this)->px = &prm;
    DoubleVector k,z,e;
    getComponents(k,z,e,prm);

    double norm1 = 0.0;
    double norm2 = 0.0;
    double norm3 = 0.0;

    norm1 = k.at(0)*k.at(0) + k.at(1)*k.at(1);
    norm2 = z.at(0)*z.at(0) + z.at(1)*z.at(1);
    norm3 = e.at(0)*e.at(0) + e.at(1)*e.at(1);

    return alpha1*norm1 + alpha2*norm2 + alpha3*norm3;
}

double Problem1L2::penalty(const DoubleVector &prm, const DoubleMatrix &u) const
{
    double sum = 0.0;

    double max = fmax(0.0, g(prm, 0, u));
    sum += 0.5*max*max;
    for (unsigned int m=1; m<=M-1; m++)
    {
        max = fmax(0.0, g(prm, m, u));
        sum += max*max;
    }
    max = fmax(0.0, g(prm, M, u));
    sum += 0.5*max*max;
    sum *= ht;
    return sum;

    return 0.0;
}

double Problem1L2::g(const DoubleVector &prm, unsigned int m, const DoubleMatrix &u) const
{
    return vd0(prm, m, u) - d1;
}

double Problem1L2::vd0(const DoubleVector &prm, unsigned int m, const DoubleMatrix &u) const
{
    DoubleVector k,z,e;
    getComponents(k,z,e,prm);

    unsigned int xi0 = (unsigned int)round(e[0] * N*10);
    unsigned int xi1 = (unsigned int)round(e[1] * N*10);

    return fabs(d0 - ((k[0]*u[m][xi0]-z[0])+(k[1]*u[m][xi1]-z[1])));
}

double sign(double x)
{
    if (x < 0.0) return -1.0;
    if (x > 0.0) return +1.0;
    return 0.0;
}

void Problem1L2::gradient(const DoubleVector &prm, DoubleVector &g)
{
    px = &prm;
    DoubleVector k,z,e;
    getComponents(k,z,e,prm);

    DoubleMatrix u;
    calculateU(u);

    DoubleMatrix p;
    calculateP(p, u);

    unsigned int i = 0;
    if (optimizeK)
    {
        // k gradient
        for (unsigned int s=0; s<L; s++)
        {
            unsigned int xi = (unsigned int)round(e.at(s) * N*10);

            double sum = 0.0;
            sum += 0.5*p.at(0, 0)*(u.at(0, xi) - z[s]);
            for (unsigned int m=1; m<=M-1; m++)
            {
                sum += p.at(m, 0)*(u.at(m, xi) - z[s]);
            }
            sum += 0.5*p.at(M, 0)*(u.at(M, xi) - z[s]);
            sum *= ht;

            // Penalty
            double pnlt = 0.0;
            pnlt += 0.5*(u[0][xi]-z[s])*sign(vd0(prm,0,u))*fmax(0.0, g(prm, 0, u));
            for (unsigned int m=1; m<=M-1; m++)
            {
                pnlt += (u[m][xi]-z[s])*sign(vd0(prm,m,u))*fmax(0.0, g(prm, m, u));
            }
            pnlt += 0.5*(u[M][xi]-z[s])*sign(vd0(prm,M,u))*fmax(0.0, g(prm, M, u));
            pnlt *= ht;

            g.at(i) = -lambda1*a*a*sum + 2.0*alpha1*k.at(s) + 2.0*R*pnlt;
            i++;
        }
    }

    // z gradient
    if (optimizeZ)
    {
        for (unsigned int s=0; s<L; s++)
        {
            double sum = 0.0;
            sum += 0.5*p.at(0, 0);
            for (unsigned int m=1; m<=M-1; m++)
            {
                sum += p.at(m, 0);
            }
            sum += 0.5*p.at(M, 0);
            sum *= ht;

            // Penalty
            double pnlt = 0.0;
            pnlt += 0.5*k[s]*sign(vd0(prm,0,u))*fmax(0.0, g(rpm, 0, u));
            for (unsigned int m=1; m<=M-1; m++)
            {
                pnlt += k[s]*sign(vd0(prm,m,u))*fmax(0.0, g(rpm, m, u));
            }
            pnlt += 0.5*k[s]*sign(vd0(prm,M,u))*fmax(0.0, g(rpm, M, u));
            pnlt *= ht;

            g.at(i) = lambda1*a*a*k[s]*sum + 2.0*alpha2*z.at(s) - 2.0*R*pnlt;
            i++;
        }
    }

    if (optimizeE)
    {
        // e gradient
        for (unsigned int s=0; s<L; s++)
        {
            unsigned int xi = (unsigned int)round(e.at(s) * N*10);

            double sum = 0.0;
            sum += 0.5 * p.at(0, 0) * ((u.at(0, xi+1) - u.at(0, xi-1))/(2.0*hx));
            for (unsigned int m=1; m<=M-1; m++)
            {
                sum += p.at(m, 0) * ((u.at(m, xi+1) - u.at(m, xi-1))/(2.0*hx));
            }
            sum += 0.5 * p.at(M, 0) * ((u.at(M, xi+1) - u.at(M, xi-1))/(2.0*hx));
            sum *= ht;

            // Penalty
            double pnlt = 0.0;
            pnlt += 0.5*k[s]*((u.at(0, xi+1) - u.at(0, xi-1))/(2.0*hx))*sign(vd0(prm,0,u))*fmax(0.0, g(rpm, 0, u));
            for (unsigned int m=1; m<=M-1; m++)
            {
                pnlt += k[s]*((u.at(m, xi+1) - u.at(m, xi-1))/(2.0*hx))*sign(vd0(prm,m,u))*fmax(0.0, g(rpm, m, u));
            }
            pnlt += 0.5*k[s]*((u.at(M, xi+1) - u.at(M, xi-1))/(2.0*hx))*sign(vd0(prm,M,u))*fmax(0.0, g(rpm, M, u));
            pnlt *= ht;

            g.at(i) = -lambda1*a*a*k[s]*sum + 2.0*alpha3*e.at(s) + 2.0*R*pnlt;
            i++;
        }
    }
}

void Problem1L2::calculateU(DoubleMatrix &u) const
{
    const DoubleVector &x = *px;

    DoubleVector k,z,e;
    getComponents(k,z,e,x);

    u.clear();
    u.resize(M+1, N+1);

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
            //u[m][20] = u[m][20] + u[m][20]*0.05;
            //u[m][80] = u[m][80] - u[m][80]*0.05;
        }
    }

    free(de);
    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

void Problem1L2::calculateP(DoubleMatrix &p, const DoubleMatrix &u)
{
    const DoubleVector &x = *px;
    DoubleVector k,z,e;
    getComponents(k,z,e,x);

    p.clear();
    p.resize(M+1, N+1);

    double *da = (double*) malloc(sizeof(double)*(N+1));
    double *db = (double*) malloc(sizeof(double)*(N+1));
    double *dc = (double*) malloc(sizeof(double)*(N+1));
    double *dd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *de = (double*) malloc(sizeof(double)*(N+1));

    for (unsigned int n=0; n<=N; n++) p.at(M,n) = -2.0*alpha0*mu(n)*(u.at(M, n) - V.at(n));

    for (unsigned int m=M-1; m != UINT32_MAX; m--)
    {
        // n = 0
        da[0] = 0.0;
        db[0] = -1.0 - (a*a*ht)/(hx*hx) - (lambda1*a*a*ht)/hx - lambda0*ht;
        dc[0] = (a*a*ht)/(hx*hx);
        dd[0] = -p.at(m+1,0);

        // n = 1,...,N-1
        for (unsigned int n=1; n<=N-1; n++)
        {
            da[n] = (a*a*ht)/(hx*hx);
            db[n] = -1.0-(2.0*a*a*ht)/(hx*hx) - lambda0*ht;
            dc[n] = (a*a*ht)/(hx*hx);
            dd[n] = -p.at(m+1,n);
        }

        // n = N
        da[N] = (a*a*ht)/(hx*hx);
        db[N] = -1.0-(a*a*ht)/(hx*hx) - (lambda2*a*a*ht)/hx - lambda0*ht;
        dc[N] = 0.0;
        dd[N] = -p.at(m+1,N);

        de[0] = de[1] = 0.0;
        de[N] = de[N-1] = 0.0;
        for (unsigned int n=2; n<=N-2; n++)
        {
            de[n] = 0.0;
            double dif0 = fabs(n*hx - e[0]);
            if (dif0 <= hx)
            {
                de[n] = k[0]*(lambda1*(a*a*ht)/hx) * (1.0 - dif0/hx);
                de[n] += 2.0*R*sign(vd0(prm,m,u))*fmax(0.0, g(rpm, m, u))/hx * (1.0 - dif0/hx);
            }
            double dif1 = fabs(n*hx - e[1]);
            if (dif1 <= hx)
            {
                de[n] = k[1]*(lambda1*(a*a*ht)/hx) * (1.0 - dif1/hx);
                de[n] += 2.0*R*sign(vd0(prm,m,u))*fmax(0.0, g(rpm, m, u))/hx * (1.0 - dif0/hx);
            }
        }

        qovmaFirstColM(da, db, dc, dd, rx, N+1, de);

        for (unsigned int i=0; i<=N; i++) p.at(m, i) = rx[i];
    }

    free(de);
    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

double Problem1L2::initial(unsigned int n UNUSED_PARAM) const
{
    return fi;
}

void Problem1L2::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double r, GradientMethod::MethodResult result) const
{
    if (i % 10 != 0 && result < 3) return;

    Problem1L2 *pm = const_cast<Problem1L2*>(this);
    pm->px = &x;

    DoubleVector k,z,e;
    getComponents(k,z,e,x);

    DoubleMatrix u;
    calculateU(u);

    unsigned int e1 = (unsigned int)round(e[0] * N*10);
    unsigned int e2 = (unsigned int)round(e[1] * N*10);
    double v = k[0]*(u[M][e1]-z[0]) + k[1]*(u[M][e2]-z[1]);
    //IPrinter::printSeperatorLine(NULL,'-',file);

    //fprintf(file,"\n");
    //fprintf(file, "J[%d]: %.10f v: %.10f\n", i, r, v);

    unsigned int p=0;
    //fprintf(file, "---\n");
    //fprintf(file, "k: %14.10f %14.10f\n", k[0], k[1]);

    printf("%d,%10.6f,",i,r);

    if (optimizeK)
    {
        DoubleVector a = g.mid(p,p+1);
        DoubleVector na = a;
        na.L2Normalize();

        DoubleVector n = a;

        DoubleVector cx = x;
        double f1,f2;

        double x0 = x[p];
        cx[p] = x0 - hk; f1 = fx(cx);
        cx[p] = x0 + hk; f2 = fx(cx);
        n[0] = (f2-f1)/(2.0*hk); cx[p] = x0;

        cx = x;
        double x1 = x[p+1];
        cx[p+1] = x1 - hk; f1 = fx(cx);
        cx[p+1] = x1 + hk; f2 = fx(cx);
        n[1] = (f2-f1)/(2.0*hk); cx[p+1] = x1;

        DoubleVector nn = n;
        nn.L2Normalize();

        //fprintf(file, "a: %14.10f %14.10f | %14.10f %14.10f\n", a[0], a[1], na[0], na[1]);
        //fprintf(file, "n: %14.10f %14.10f | %14.10f %14.10f\n", n[0], n[1], nn[0], nn[1]);

        p+=2;

        printf("|%10.6f,%10.6f,%10.6f,",k[0],a[0],n[0]);
        printf("|%10.6f,%10.6f,%10.6f,",k[1],a[1],n[1]);
    }

    //fprintf(file, "---\n");
    //fprintf(file, "z: %14.10f %14.10f\n", z[0], z[1]);
    if (optimizeZ)
    {
        DoubleVector a = g.mid(p,p+1);
        DoubleVector na = a;
        na.L2Normalize();

        DoubleVector n = a;

        DoubleVector cx = x;
        double f1,f2;

        double x0 = x[p];
        cx[p] = x0 - hz; f1 = fx(cx);
        cx[p] = x0 + hz; f2 = fx(cx);
        n[0] = (f2-f1)/(2.0*hz); cx[p] = x0;

        double x1 = x[p+1];
        cx[p+1] = x1 - hz; f1 = fx(cx);
        cx[p+1] = x1 + hz; f2 = fx(cx);
        n[1] = (f2-f1)/(2.0*hz); cx[p+1] = x1;

        DoubleVector nn = n;
        nn.L2Normalize();

        //fprintf(file, "a: %14.10f %14.10f | %14.10f %14.10f\n", a[0], a[1], na[0], na[1]);
        //fprintf(file, "n: %14.10f %14.10f | %14.10f %14.10f\n", n[0], n[1], nn[0], nn[1]);
        \
        p+=2;
        printf("|%10.6f,%10.6f,%10.6f,",z[0],a[0],n[0]);
        printf("|%10.6f,%10.6f,%10.6f\n",z[1],a[1],n[1]);
    }

    //fprintf(file, "---\n");
    //fprintf(file, "e: %14.10f %14.10f\n", e[0], e[1]);
    if (optimizeE)
    {
        DoubleVector a = g.mid(p,p+1);
        DoubleVector na = a;
        na.L2Normalize();

        DoubleVector n = a;

        DoubleVector cx = x;
        double f1,f2;

        double x0 = x[p];
        cx[p] = x0 - he; f1 = fx(cx);
        cx[p] = x0 + he; f2 = fx(cx);
        n[0] = (f2-f1)/(2.0*he); cx[p] = x0;

        double x1 = x[p+1];
        cx[p+1] = x1 - he; f1 = fx(cx);
        cx[p+1] = x1 + he; f2 = fx(cx);
        n[1] = (f2-f1)/(2.0*he); cx[p+1] = x1;

        DoubleVector nn = n;
        nn.L2Normalize();

        //fprintf(file, "a: %14.10f %14.10f | %14.10f %14.10f\n", a[0], a[1], na[0], na[1]);
        //fprintf(file, "n: %14.10f %14.10f | %14.10f %14.10f\n", n[0], n[1], nn[0], nn[1]);
        p+=2;
        printf("%.6f,%.6f,%.6f,",e[0],a[0],n[0]);
        printf("%.6f,%.6f,%.6f\n",e[1],a[1],n[1]);
    }

    IPrinter::printVector(14,10,u.row(u.rows()-1),"u: ", 10, 0, 0, stdout);
    fflush(file);
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

void Problem1L2::getComponents(DoubleVector &k, DoubleVector &z, DoubleVector &e, const DoubleVector &x) const
{
    unsigned int p = 0;

    if (optimizeK)
    {
        k = x.mid(p,p+1);
        p+=L;
    }
    else
    {
        k = this->K;
    }

    if (optimizeZ)
    {
        z = x.mid(p,p+1);
        p+=L;
    }
    else
    {
        z = this->z;
    }

    if (optimizeE)
    {
        e = x.mid(p,p+1);
        p+=L;
    }
    else
    {
        e = this->e;
    }
}

void Problem1L2::qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
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

void Problem1L2::qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
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
