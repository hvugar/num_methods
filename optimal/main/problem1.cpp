#include "problem1.h"

void Problem1::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

#if !defined (_OPTIMIZE_K_) && !defined (_OPTIMIZE_Z_) && !defined (_OPTIMIZE_E_)
    puts("Nothing to optimize...");
    return;
#endif

    Problem1 e3;
    e3.initialize();
}

Problem1::Problem1()
{
}

void Problem1::initialize()
{
#ifndef _OPTIMIZE_K_
    k << 2.50 << 2.70;
#endif
#ifndef _OPTIMIZE_Z_
    z << 1.52 << 1.71;
#endif
#ifndef _OPTIMIZE_E_
    e << 0.25 << 0.75;
#endif

    /***************************************************************/
#ifdef _OPTIMIZE_K_
    xs << 2.50 << 2.70;
#endif
#ifdef _OPTIMIZE_Z_
    xs << 1.52 << 1.71;
#endif
#ifdef _OPTIMIZE_E_
    xs << 0.25 << 0.75;
#endif
    /***************************************************************/

    px = &xs;
    DoubleMatrix u;
    calculateU(u);
    V = u.row(M);

    DoubleVector x0;
#ifdef _OPTIMIZE_K_
    x0 << 3.50 << 3.70;
#endif
#ifdef _OPTIMIZE_Z_
    x0 << 2.52 << 2.71;
#endif
#ifdef _OPTIMIZE_E_
    x0 << 0.30 << 0.60;
#endif
    px = &x0;

    printNAGradinets(x0);

    ConjugateGradient g;
    g.setFunction(this);
    g.setGradient(this);
    g.setPrinter(this);
    g.setProjection(this);
    g.setEpsilon1(0.00000001);
    g.setEpsilon2(0.00000001);
    g.setEpsilon3(0.00000001);
    g.setR1MinimizeEpsilon(2.0, 0.00000001);
    g.setNormalize(true);
    g.calculate(x0);
}

void Problem1::printNAGradinets(const DoubleVector &x0)
{
    puts("------------------------------------------");

    // Calculating numerical gradients
    DoubleVector gn(x0.size());
    IGradient::Gradient(this, h, x0, gn);
    gn.L2Normalize();

    // Calculating analitical gradients
    DoubleVector ga(x0.size());
    gradient(x0, ga);
    ga.L2Normalize();

    unsigned int i;
    printf("Optimal:   ");
    i = 0;
#if defined(_OPTIMIZE_K_)
    printf("%12.8f %12.8f ", xs.at(i), xs.at(i+1));
    i+=2;
#endif
#if defined(_OPTIMIZE_Z_)
    printf("%12.8f %12.8f ", xs.at(i), xs.at(i+1));
    i+=2;
#endif
#if defined(_OPTIMIZE_E_)
    printf("%12.8f %12.8f ", xs.at(i), xs.at(i+1));
    i+=2;
#endif
    puts("");

    printf("Initial:   ");
    i = 0;
#if defined(_OPTIMIZE_K_)
    printf("%12.8f %12.8f ", x0.at(i), x0.at(i+1));
    i+=2;
#endif
#if defined(_OPTIMIZE_Z_)
    printf("%12.8f %12.8f ", x0.at(i), x0.at(i+1));
    i+=2;
#endif
#if defined(_OPTIMIZE_E_)
    printf("%12.8f %12.8f ", x0.at(i), x0.at(i+1));
    i+=2;
#endif
    puts("");

    printf("Numerical: ");
    i = 0;
#if defined(_OPTIMIZE_K_)
    printf("%12.8f %12.8f ", gn.at(i), gn.at(i+1));
    i+=2;
#endif
#if defined(_OPTIMIZE_Z_)
    printf("%12.8f %12.8f ", gn.at(i), gn.at(i+1));
    i+=2;
#endif
#if defined(_OPTIMIZE_E_)
    printf("%12.8f %12.8f ", gn.at(i), gn.at(i+1));
    i+=2;
#endif
    puts("");

    printf("Analytic:  ");
    i = 0;
#if defined(_OPTIMIZE_K_)
    printf("%12.8f %12.8f ", ga.at(i), ga.at(i+1));
    i+=2;
#endif
#if defined(_OPTIMIZE_Z_)
    printf("%12.8f %12.8f ", ga.at(i), ga.at(i+1));
    i+=2;
#endif
#if defined(_OPTIMIZE_E_)
    printf("%12.8f %12.8f ", ga.at(i), ga.at(i+1));
    i+=2;
#endif
    puts("");

    //#if defined (_OPTIMIZE_K_) && defined (_OPTIMIZE_Z_) && defined (_OPTIMIZE_E_)
    //    printf("Optimal:   %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", xs.at(0), xs.at(1), xs.at(2), xs.at(3), xs.at(4), xs.at(5));
    //    printf("Initial:   %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", x0.at(0), x0.at(1), x0.at(2), x0.at(3), x0.at(4), x0.at(5));
    //    printf("Numerical: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", g2.at(0), g2.at(1), g2.at(2), g2.at(3), g2.at(4), g2.at(5));
    //    printf("Analytic:  %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", g1.at(0), g1.at(1), g1.at(2), g1.at(3), g1.at(4), g1.at(5));
    //#else
    //    printf("Optimal:   %12.8f %12.8f\n", xs.at(0), xs.at(1));
    //    printf("Initial:   %12.8f %12.8f\n", x0.at(0), x0.at(1));
    //    printf("Numerical: %12.8f %12.8f\n", gn.at(0), gn.at(1));
    //    printf("Analytic:  %12.8f %12.8f\n", ga.at(0), ga.at(1));
    //#endif
    puts("------------------------------------------");
}

double Problem1::initial(unsigned int i UNUSED_PARAM) const
{
    return Ti;
}

double Problem1::fx(const DoubleVector &x)
{
    px = &x;

    DoubleVector k,z,e;
    getComponents(k,z,e,x);

    DoubleMatrix u;
    calculateU(u);

    double sum = 0.0;
    for (unsigned int n=0; n<N; n++)
    {
        unsigned int n1 = n+0;
        unsigned int n2 = n+1;
        double f1 = mu(n1)*(u.at(M, n1)-V.at(n1))*(u.at(M, n1)-V.at(n1));
        double f2 = mu(n2)*(u.at(M, n2)-V.at(n2))*(u.at(M, n2)-V.at(n2));
        sum = sum + (f1 + f2);
    }
    sum = 0.5*hx*sum;

#ifndef __V_NORM__
    double norm1 = 0.0;
    double norm2 = 0.0;
    double norm3 = 0.0;

    int i = 0;
#ifdef _OPTIMIZE_K_
    norm1 = sqrt((k.at(0)-xs.at(i))*(k.at(0)-xs.at(i))+(k.at(1)-xs.at(i+1))*(k.at(1)-xs.at(i+1)));
    i+=2;
#endif
#ifdef _OPTIMIZE_Z_
    norm2 = sqrt((z.at(0)-xs.at(i))*(z.at(0)-xs.at(i))+(z.at(1)-xs.at(i+1))*(z.at(1)-xs.at(i+1)));
    i+=2;
#endif
#ifdef _OPTIMIZE_E_
    norm3 = sqrt((e.at(0)-xs.at(i))*(e.at(0)-xs.at(i))+(e.at(1)-xs.at(i+1))*(e.at(1)-xs.at(i+1)));
    i+=2;
#endif
    return alpha0*sum + alpha1*norm1 + alpha2*norm2 + alpha3*norm3;
#else
    double norm = 0.0;
    for (unsigned int m=0; m<M; m++)
    {
        unsigned int m1 = m + 0;
        unsigned int m2 = m + 1;
        double f1 = v(m1, u);
        double f2 = v(m2, u);
        norm += (f1*f1 + f2*f2);
    }
    return alpha0*sum + alpha4*norm*0.5*ht;
#endif
}

void Problem1::gradient(const DoubleVector &x, DoubleVector &g)
{
    px = &x;

    DoubleVector k,z,e;
    getComponents(k,z,e,x);

    DoubleMatrix u;
    calculateU(u);

    DoubleMatrix psi;
    calculateP(psi, u);

    unsigned int i = 0;
    // k gradient
#ifdef _OPTIMIZE_K_
    for (unsigned int s=0; s<L; s++)
    {
        unsigned int xi = (unsigned int)round(e.at(s) * N);
        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = psi.at(m1, 0)*(u.at(m1, xi) - z[s]);
            double g2 = psi.at(m2, 0)*(u.at(m2, xi) - z[s]);
            sum = sum + (g1 + g2);
        }
        g.at(i) = 0.5*ht*(-lambda0*a*a)*sum;
#ifndef __V_NORM__
        g.at(i) += + 2.0*alpha1*(k.at(s)-xs.at(i));
#else
        double sumk = 0.0;
        for (unsigned int m=0; m<M; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double f1 = (u.at(m1,xi)-z.at(s))*v(m1, u);
            double f2 = (u.at(m2,xi)-z.at(s))*v(m2, u);
            sumk += (f1 + f2);
        }
        g.at(i) += 2.0*alpha4*sumk*0.5*ht;
#endif
        i++;
    }
#endif

    // z gradient
#ifdef _OPTIMIZE_Z_
    for (unsigned int s=0; s<L; s++)
    {
        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = psi.at(m1, 0);
            double g2 = psi.at(m2, 0);
            sum = sum + (g1 + g2);
        }
        g.at(i) = 0.5*ht*(lambda0*a*a)*k[s]*sum;
#ifndef __V_NORM__
        g.at(i) += 2.0*alpha2*(z.at(s)-xs.at(i));
#else
        double sumz = 0.0;
        for (unsigned int m=0; m<M; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double f1 = k.at(s)*v(m1, u);
            double f2 = k.at(s)*v(m2, u);
            sumz += (f1 + f2);
        }
        g.at(i) += 2.0*alpha4*sumz*0.5*ht;
#endif
        i++;
    }
#endif

    // e gradient
#ifdef _OPTIMIZE_E_
    for (unsigned int s=0; s<L; s++)
    {
        unsigned int xi = (unsigned int)round(e.at(s) * N);
        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = psi.at(m1, 0) * ((u.at(m1, xi+1) - u.at(m1, xi-1))/(2.0*hx));
            double g2 = psi.at(m2, 0) * ((u.at(m2, xi+1) - u.at(m2, xi-1))/(2.0*hx));
            sum = sum + (g1 + g2);
        }
        g.at(i) = 0.5*ht*(-lambda0*a*a)*k.at(s)*sum;
#ifndef __V_NORM__
        g.at(i) += 2.0*alpha3*(e.at(s)-xs.at(i));
#else
        double sume = 0.0;
        for (unsigned int m=0; m<M; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double f1 = k.at(s)*((u.at(m1,xi)-u.at(m1,xi-1))/hx)*v(m1, u);
            double f2 = k.at(s)*((u.at(m2,xi)-u.at(m2,xi-1))/hx)*v(m2, u);
            sume += (f1 + f2);
        }
        g.at(i) += 2.0*alpha4*sume*0.5*ht;
#endif
        i++;
    }
#endif
}

void Problem1::print(unsigned int i UNUSED_PARAM, const DoubleVector &x UNUSED_PARAM, const DoubleVector &g UNUSED_PARAM, double alpha UNUSED_PARAM, RnFunction *fn UNUSED_PARAM) const
{
    //    C_UNUSED(alpha);
    //    printf("J[%d]: %18.14f %16.12f ", i, fn->fx(x), alpha);

    //    int j = 0;
    //#ifdef _OPTIMIZE_K_
    //    printf("%14.10f %14.10f ", x.at(j), x.at(j+1));
    //    j+=2;
    //#endif
    //#ifdef _OPTIMIZE_Z_
    //    printf("%14.10f %14.10f ", x.at(j), x.at(j+1));
    //    j+=2;
    //#endif
    //#ifdef _OPTIMIZE_E_
    //    printf("%14.10f %14.10f ", x.at(j), x.at(j+1));
    //    j+=2;
    //#endif
    //    puts("");

    //#ifdef OPTIMIZE_REPLACEMENT
    //    printf("k:  %14.10f, %14.10f %14.10f, %14.10f %14.10f, %14.10f\n", x.at(0), x.at(1), x.at(2), x.at(3), x.at(4), x.at(5));
    //    printf("g:  %14.10f, %14.10f %14.10f, %14.10f %14.10f, %14.10f\n", g.at(0), g.at(1), g.at(2), g.at(3), g.at(4), g.at(5));
    //    DoubleVector w = g;
    //    w.L2Normalize();
    //    printf("g:  %14.10f, %14.10f %14.10f, %14.10f %14.10f, %14.10f\n", w.at(0), w.at(1), w.at(2), w.at(3), w.at(4), w.at(5));
    //#else
    //    printf("x:  %14.10f, %14.10f %14.10f, %14.10f\n", x.at(0), x.at(1), x.at(2), x.at(3));
    //    printf("g:  %14.10f, %14.10f %14.10f, %14.10f\n", g.at(0), g.at(1), g.at(2), g.at(3));
    //#endif
}

void Problem1::print(const DoubleVector &x, const DoubleVector &g UNUSED_PARAM, unsigned int i) const
{
    printf("J[%d]: %18.14f \n", i, const_cast<Problem1*>(this)->fx(x));

    int j = 0;
    printf("x: ");
#ifdef _OPTIMIZE_K_
    printf("%14.10f %14.10f ", x.at(j), x.at(j+1));
    j+=2;
#endif
#ifdef _OPTIMIZE_Z_
    printf("%14.10f %14.10f ", x.at(j), x.at(j+1));
    j+=2;
#endif
#ifdef _OPTIMIZE_E_
    printf("%14.10f %14.10f ", x.at(j), x.at(j+1));
    j+=2;
#endif
    puts("");

    j = 0;
    DoubleVector n = g;
    n.L2Normalize();
    printf("n: ");
#ifdef _OPTIMIZE_K_
    printf("%14.10f %14.10f ", n.at(j), n.at(j+1));
    j+=2;
#endif
#ifdef _OPTIMIZE_Z_
    printf("%14.10f %14.10f ", n.at(j), n.at(j+1));
    j+=2;
#endif
#ifdef _OPTIMIZE_E_
    printf("%14.10f %14.10f ", n.at(j), n.at(j+1));
    j+=2;
#endif
    puts("");

    j = 0;
    printf("g: ");
#ifdef _OPTIMIZE_K_
    printf("%14.10f %14.10f ", g.at(j), g.at(j+1));
    j+=2;
#endif
#ifdef _OPTIMIZE_Z_
    printf("%14.10f %14.10f ", g.at(j), g.at(j+1));
    j+=2;
#endif
#ifdef _OPTIMIZE_E_
    printf("%14.10f %14.10f ", g.at(j), g.at(j+1));
    j+=2;
#endif
    puts("");

    DoubleVector gn(x.size());
    IGradient::Gradient(const_cast<Problem1*>(this), h, x, gn);
    gn.L2Normalize();
    j = 0;
    gn.L2Normalize();
    printf("gn:");
#ifdef _OPTIMIZE_K_
    printf("%14.10f %14.10f ", gn.at(j), gn.at(j+1));
    j+=2;
#endif
#ifdef _OPTIMIZE_Z_
    printf("%14.10f %14.10f ", gn.at(j), gn.at(j+1));
    j+=2;
#endif
#ifdef _OPTIMIZE_E_
    printf("%14.10f %14.10f ", gn.at(j), gn.at(j+1));
    j+=2;
#endif

    puts("\n-----------------------------------------");
}

void Problem1::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
{
    int j = 0;
#ifdef _OPTIMIZE_K_
    j+=2;
#endif
#ifdef _OPTIMIZE_Z_
    j+=2;
#endif
#ifdef _OPTIMIZE_E_
    if (x.at(j) < 0.10) x.at(j) = 0.10;
    if (x.at(j) > 0.90) x.at(j) = 0.90;

    if (x.at(j+1) < 0.10) x.at(j+1) = 0.10;
    if (x.at(j+1) > 0.90) x.at(j+1) = 0.90;
#endif
}

void Problem1::calculateU(DoubleMatrix &u)
{
    DoubleVector x = *px;

    DoubleVector k,z,e;
    getComponents(k,z,e,x);

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);
    DoubleVector de(N+1);

    for (unsigned int m=0; m<=M; m++)
    {
        if (m==0) for (unsigned int n=0; n<=N; n++) u.at(0,n) = initial(n);
        else
        {
            // n = 0
            da[0] = 0.0;
            db[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambda0*a*a*ht)/hx + alpha*ht;
            dc[0] = -(a*a*ht)/(hx*hx);
            dd[0] = u.at(m-1,0) + alpha*ht*Te - ((lambda0*a*a*ht)/hx)*(k[0]*z[0] + k[1]*z[1]);

            // n = 1,...,N-1
            for (unsigned int n=1; n<=N-1; n++)
            {
                da[n] = -(a*a*ht)/(hx*hx);
                db[n] = 1.0 + (2.0*a*a*ht)/(hx*hx) + alpha*ht;
                dc[n] = -(a*a*ht)/(hx*hx);
                dd[n] = u.at(m-1,n) + alpha*ht*Te;
            }

            // n = N
            da[N] = -(a*a*ht)/(hx*hx);
            db[N] = 1.0 + (a*a*ht)/(hx*hx) + (lambdal*a*a*ht)/hx + alpha*ht;
            dc[N] = 0.0;
            dd[N] = u.at(m-1,N) + (lambdal*a*a*ht*Te)/hx + alpha*ht*Te;

            for (unsigned int n=0; n<=N; n++)
            {
                de[n] = 0.0;

                if (fabs(n*hx - e.at(0)) <= hx)
                {
                    de[n] = -k[0]*(lambda0*(a*a*ht)/hx) * (1.0 - fabs(n*hx - e.at(0))/hx);
                }
                if (fabs(n*hx - e.at(1)) <= hx)
                {
                    de[n] = -k[1]*(lambda0*(a*a*ht)/hx) * (1.0 - fabs(n*hx - e.at(1))/hx);
                }

                //if (fabs(n*hx - e.at(0)) <= DBL_EPSILON) { de[n] = -k[0]*(lambda0*(a*a*ht)/hx); }
                //if (fabs(n*hx - e.at(1)) <= DBL_EPSILON) { de[n] = -k[1]*(lambda0*(a*a*ht)/hx); }
            }

            qovmaFirstRow(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size(), de.data());

            for (unsigned int i=0; i<=N; i++) u.at(m, i) = rx[i];
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
    de.clear();
}

void Problem1::calculateP(DoubleMatrix &p, const DoubleMatrix &u)
{
    DoubleVector k,z,e;
    getComponents(k,z,e,*px);

    p.clear();
    p.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);
    DoubleVector de(N+1);

    for (unsigned int m=M; m != UINT32_MAX; m--)
    {
        if (m==M)
        {
            for (unsigned int n=0; n<=N; n++) p.at(M,n) = -2.0*alpha0*mu(n)*(u.at(M, n) - V.at(n));
        }
        else
        {
            // n = 0
            da[0] = 0.0;
            db[0] = -1.0 - (a*a*ht)/(hx*hx) - (lambda0*a*a*ht)/hx - alpha*ht;
            dc[0] = (a*a*ht)/(hx*hx);
            dd[0] = -p.at(m+1,0);

            // n = 1,...,N-1
            for (unsigned int n=1; n<=N-1; n++)
            {
                da[n] = (a*a*ht)/(hx*hx);
                db[n] = -1.0-(2.0*a*a*ht)/(hx*hx) - alpha*ht;
                dc[n] = (a*a*ht)/(hx*hx);
                dd[n] = -p.at(m+1,n);

#ifdef __V_NORM__
                if (fabs(n*hx - e.at(0)) <= hx)
                {
                    dd[n] += 2.0*alpha4*(1.0 - fabs(n*hx - e.at(0))/hx)*v(m,u);
                }
                if (fabs(n*hx - e.at(1)) <= hx)
                {
                    dd[n] += 2.0*alpha4*(1.0 - fabs(n*hx - e.at(1))/hx)*v(m,u);
                }
#endif
            }

            // n = N
            da[N] = (a*a*ht)/(hx*hx);
            db[N] = -1.0-(a*a*ht)/(hx*hx) - (lambdal*a*a*ht)/hx - alpha*ht;
            dc[N] = 0.0;
            dd[N] = -p.at(m+1,N);

            for (unsigned int n=0; n<=N; n++)
            {
                de[n] = 0.0;
                if (fabs(n*hx - e.at(0)) <= hx)
                {
                    de[n] = k.at(0)*(lambda0*(a*a*ht)/hx) * (1.0 - fabs(n*hx - e.at(0))/hx);
                }
                if (fabs(n*hx - e.at(1)) <= hx)
                {
                    de[n] = k.at(1)*(lambda0*(a*a*ht)/hx) * (1.0 - fabs(n*hx - e.at(1))/hx);
                }
            }

            qovmaFirstCol(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size(), de.data());

            for (unsigned int i=0; i<=N; i++) p.at(m, i) = rx[i];
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
    de.clear();
}

//void Problem1::calculateU1(DoubleMatrix &u)
//{
//    DoubleVector x = *px;

//    DoubleVector k,z,e;
//    getComponents(k,z,e,x);

//    u.clear();
//    u.resize(M+1, N+1);

//    for (unsigned int m=0; m<=M; m++)
//    {
//        if (m==0)
//        {
//            for (unsigned int n=0; n<=N; n++) u.at(0,n) = initial(n);
//        }
//        else
//        {
//            DoubleMatrix betta(2, N+1, 0.0);
//            DoubleVector qamma(2);

//            betta.at(0,0) = 1.0 + (a*a*ht)/(hx*hx) + (lambda0*a*a*ht)/hx + alpha*ht;
//            betta.at(0,1) = -(a*a*ht)/(hx*hx);
//            for (unsigned int n=2; n<=N; n++)
//            {
//                if (n == (unsigned int)(e.at(0)*N)) { betta.at(0,n) = -k.at(0)*(lambda0*(a*a*ht)/hx); }
//                if (n == (unsigned int)(e.at(1)*N)) { betta.at(0,n) = -k.at(1)*(lambda0*(a*a*ht)/hx); }
//            }
//            qamma.at(0) = u.at(m-1,0) + alpha*ht*Te - ((lambda0*a*a*ht)/hx)*(k.at(0)*z.at(0) + k.at(1)*z.at(1));

//            betta.at(1,N-1) = -(a*a*ht)/(hx*hx);
//            betta.at(1,N-0) = 1.0 + (a*a*ht)/(hx*hx) + (lambdal*a*a*ht)/hx + alpha*ht;
//            qamma.at(1) = u.at(m-1,N) + (lambdal*a*a*ht*Te)/hx + alpha*ht*Te;

//            DoubleMatrix alfa(N-1,3);
//            for (unsigned int n=1; n<=N-1; n++)
//            {
//                double A = -(a*a*ht)/(hx*hx);
//                double B = 1.0 + (2.0*a*a*ht)/(hx*hx) + alpha*ht;
//                double C = -(a*a*ht)/(hx*hx);
//                double D = u.at(m-1,n) + alpha*ht*Te;

//                alfa.at(n-1,0) = +D/A;
//                alfa.at(n-1,1) = -B/A;
//                alfa.at(n-1,2) = -C/A;
//            }

//            for (unsigned int i=0; i<=N-2; i++)
//            {
//                //double A = -(a*a*ht)/(hx*hx);
//                //double B = 1.0 + (2.0*a*a*ht)/(hx*hx) + alpha*ht;
//                //double C = -(a*a*ht)/(hx*hx);
//                //double D = u.at(m-1,i+1) + alpha*ht*Te;

//                //betta.at(0,i+1) = betta.at(0,i+1) + betta.at(0,i)*(-B/A);//alfa.at(i,1);
//                //betta.at(0,i+2) = betta.at(0,i+2) + betta.at(0,i)*(-C/A);//alfa.at(i,2);
//                //qamma.at(0)     = qamma.at(0)     - betta.at(0,i)*(+D/A);//alfa.at(i,0);

//                betta.at(0,i+1) = betta.at(0,i+1) + betta.at(0,i)*alfa.at(i,1);
//                betta.at(0,i+2) = betta.at(0,i+2) + betta.at(0,i)*alfa.at(i,2);
//                qamma.at(0)     = qamma.at(0)     - betta.at(0,i)*alfa.at(i,0);
//            }

//            DoubleMatrix m1(2,2);
//            m1.at(0,0) = betta.at(0,N-1);
//            m1.at(0,1) = betta.at(0,N-0);
//            m1.at(1,0) = betta.at(1,N-1);
//            m1.at(1,1) = betta.at(1,N-0);

//            DoubleVector b1(2);
//            b1.at(0) = qamma.at(0);
//            b1.at(1) = qamma.at(1);
//            //puts("---");
//            //IPrinter::printMatrix(m1,m1.rows(),m1.cols());
//            //puts("---");

//            DoubleVector u1(2);
//            GaussianElimination(m1, b1, u1);

//            m1.clear();
//            b1.clear();

//            u.at(m, N-1) = u1.at(0);
//            u.at(m, N-0) = u1.at(1);
//            puts("---");
//            printf("%4d %12.8f %12.8f %12.8f\n", m, u.at(m,N-2), u.at(m,N-1), u.at(m,N));

//            for (unsigned int i=N-2; i != UINT_MAX; i--)
//            {
//                u.at(m,i) = alfa.at(i,1)*u.at(m,i+1) + alfa.at(i,2)*u.at(m,i+2) + alfa.at(i,0);
//            }
//            printf("%4d %12.8f %12.8f %12.8f\n", m, u.at(m,N-2), u.at(m,N-1), u.at(m,N));

//            if (m==M)
//            {
//                //    printf("%4d %12.8f %12.8f %12.8f\n", m, u.at(m,N-2), u.at(m,N-1), u.at(m,N));
//            }

//            u1.clear();
//            betta.clear();
//            alfa.clear();
//        }
//    }
//}

void Problem1::getComponents(DoubleVector &k1, DoubleVector &z1, DoubleVector &e1, const DoubleVector &x) const
{
    unsigned int is = 0;
    unsigned int ie = 1;
#ifdef _OPTIMIZE_K_
    k1 = x.mid(is, ie);
    is += L;
    ie += L;
#else
    k1 = k;
#endif

#ifdef _OPTIMIZE_Z_
    z1 = x.mid(is, ie);
    is += L;
    ie += L;
#else
    z1 = z;
#endif

#ifdef _OPTIMIZE_E_
    e1 = x.mid(is, ie);
    is += L;
    ie += L;
#else
    e1 = e;
#endif
}

void qovmaFirstCol(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e)
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

void qovmaFirstRow(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);

    unsigned int L = 0;
    for (unsigned int s=2; s<n; s++) if (e[s] != 0.0) L+=1;
    unsigned int *E = (unsigned int *)malloc(sizeof(unsigned int)*L);

    unsigned int i = 0;
    for (unsigned int s=2; s<n; s++)
    {
        if (e[s] != 0.0)
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
            p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
            q[i] = -c[i]/(b[i]+a[i]*q[i-1]);

            for (unsigned int s=0; s<L; s++)
            {
                if (i<(E[s]-1))
                    k[s][i] = -(a[i]*k[s][i-1])/(b[i]+a[i]*q[i-1]);
                else
                    k[s][i] = 0.0;
            }

            for (unsigned int s=0; s<L; s++) if (i==E[s]-1) q[i] += -(a[i]*k[s][i-1])/(b[i]+a[i]*q[i-1]);
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

double Problem1::v(unsigned int j, const DoubleMatrix &u) const
{
    const DoubleVector &x = *px;
    DoubleVector k,z,e;
    getComponents(k,z,e,x);

    double v1 = k.at(0)*(u.at(j, (unsigned int)round(e.at(0)*N)) - z.at(0)) + k.at(1)*(u.at(j, (unsigned int)round(e.at(1)*N)) - z.at(1));

    unsigned int i = 0;
    DoubleVector ks = k;
    DoubleVector zs = z;
    DoubleVector es = e;
#ifdef _OPTIMIZE_K_
    ks.at(0) = xs.at(i);
    ks.at(1) = xs.at(i+1);
    i += 2;
#endif
#ifdef _OPTIMIZE_Z_
    zs.at(0) = xs.at(i);
    zs.at(1) = xs.at(i+1);
    i += 2;
#endif
#ifdef _OPTIMIZE_E_
    es.at(0) = xs.at(i);
    es.at(1) = xs.at(i+1);
    i += 2;
#endif

    double vs = ks.at(0)*(u.at(j, (unsigned int)round(es.at(0)*N)) - zs.at(0)) + ks.at(1)*(u.at(j, (unsigned int)round(es.at(1)*N)) - zs.at(1));

    return v1 - vs;

}

/*
void Example3::calculateU1(DoubleMatrix &u, const DoubleVector &x)
{
    DoubleVector k = x.mid(0, 1);
    DoubleVector z = x.mid(2, 3);
    DoubleVector e = x.mid(4, 5);

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector rb(N+1);
    DoubleVector rx(N+1);

    for (unsigned int m=0; m<=M; m++)
    {
        if (m==0)
        {
            for (unsigned int n=0; n<=N; n++) u.at(0,n) = initial(n);
        }
        else
        {
            DoubleMatrix ra(N+1, N+1, 0.0);

            ra(0,0) = 1.0 + (a*a*ht)/(hx*hx) + (lambda0*a*a*ht)/hx + alpha*ht;
            ra(0,1) = -(a*a*ht)/(hx*hx);
            rb[0] = u.at(m-1,0) + alpha*ht*Te - ((lambda0*a*a*ht)/(hx))*(k[0]*z[0]+k[1]*z[1]);

            ra(0, 40) = -k[0] * ((lambda0*a*a*ht)/(hx));
            ra(0, 70) = -k[1] * ((lambda0*a*a*ht)/(hx));
            //rb[0] = rb[0] - ((lambda0*a*a*ht)/(hx))*((k[0]*z[0]+k[1]*z[1]));

            for (unsigned int i=1; i<=N-1; i++)
            {
                ra(i,i-1) = -(a*a*ht)/(hx*hx);
                ra(i,i) = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + alpha*ht;
                ra(i,i+1) = -(a*a*ht)/(hx*hx);
                rb[i] = u.at(m-1, i) + alpha*ht*Te;
            }

            ra(N,N-1) = -(a*a*ht)/(hx*hx);
            ra(N,N)   = 1.0 + (a*a*ht)/(hx*hx) + (lambdal*a*a*ht)/hx + alpha*ht;
            rb[N] = u.at(m-1,N) + ((lambdal*a*a*ht)/(hx))*Te + alpha*ht*Te;

            GaussianElimination(ra, rb, rx);

            for (unsigned int i=0; i<=N; i++) u.at(m,i) = rx[i];
        }
    }
}

void Example3::calculateU2(DoubleMatrix &u, const DoubleVector &x)
{
    DoubleVector k = x.mid(0, 1);
    DoubleVector z = x.mid(2, 3);
    DoubleVector e = x.mid(4, 5);

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);

    DoubleVector de(L);

    unsigned int E[2];
    E[0] = (unsigned int)round(e[0]/hx);
    E[1] = (unsigned int)round(e[1]/hx);

    for (unsigned int m=0; m<=M; m++)
    {
        if (m==0) for (unsigned int n=0; n<=N; n++) u.at(0,n) = initial(n);
        else
        {
            // n = 0
            da[0] = 0.0;
            db[0] = 1.0+(a*a*ht)/(hx*hx) + (lambda0*a*a*ht)/hx + alpha*ht;
            dc[0] = -(a*a*ht)/(hx*hx);
            dd[0] = u.at(m-1,0) + alpha*ht*Te - ((lambda0*a*a*ht)/hx)*(k[0]*z[0] + k[1]*z[1]);

            de[0] = -k[0]*(lambda0*(a*a*ht)/hx);
            de[1] = -k[1]*(lambda0*(a*a*ht)/hx);

            // n = 1,...,N-1
            for (unsigned int n=1; n<=N-1; n++)
            {
                da[n] = -(a*a*ht)/(hx*hx);
                db[n] = 1.0+(2.0*a*a*ht)/(hx*hx) + alpha*ht;
                dc[n] = -(a*a*ht)/(hx*hx);
                dd[n] = u.at(m-1,n) + alpha*ht*Te;
            }

            // n = N
            da[N] = -(a*a*ht)/(hx*hx);
            db[N] = 1.0+(a*a*ht)/(hx*hx) + (lambdal*a*a*ht)/hx + alpha*ht;
            dc[N] = 0.0;
            dd[N] = u.at(m-1,N) + (lambdal*a*a*ht*Te)/hx + alpha*ht*Te;

            qovmaE(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size(), de.data(), E, 2);

            for (unsigned int i=0; i<=N; i++) u.at(m, i) = rx[i];
        }
        //printf("v: %10.6f |", k[0]*(u.at(m,40)-z[0])+k[1]*(u.at(m,70)-z[1]));
        //IPrinter::printVector(u.row(m));
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
    de.clear();
}

void Example3::calculateU3(DoubleMatrix &u, const DoubleVector &x)
{
    DoubleVector k = x.mid(0, 1);
    DoubleVector z = x.mid(2, 3);
    DoubleVector e = x.mid(4, 5);

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);

    for (unsigned int m=0; m<=1; m++)
    {
        if (m==0) for (unsigned int n=0; n<=N; n++) u.at(0,n) = initial(n);
        else
        {
            // n = 0
            da[0] = 0.0;
            db[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambda0*a*a*ht)/hx + alpha*ht;
            dc[0] = -(a*a*ht)/(hx*hx);
            dd[0] = u.at(m-1,0) + alpha*ht*Te - ((lambda0*a*a*ht)/hx)*(k[0]*z[0] + k[1]*z[1]);

            // n = 1,...,N-1
            for (unsigned int n=1; n<=N-1; n++)
            {
                da[n] = -(a*a*ht)/(hx*hx);
                db[n] = 1.0 + (2.0*a*a*ht)/(hx*hx) + alpha*ht;
                dc[n] = -(a*a*ht)/(hx*hx);
                dd[n] = u.at(m-1,n) + alpha*ht*Te;
            }

            // n = N
            da[N] = -(a*a*ht)/(hx*hx);
            db[N] = 1.0 + (a*a*ht)/(hx*hx) + (lambdal*a*a*ht)/hx + alpha*ht;
            dc[N] = 0.0;
            dd[N] = u.at(m-1,N) + (lambdal*a*a*ht*Te)/hx + alpha*ht*Te;

            unsigned int K = 2;
            DoubleMatrix ma(N+1-K, K+1);
            for (unsigned int i=0; i<ma.rows(); i++)
            {
                ma.at(i,1) = -db.at(i+1) / da.at(i+1);
                ma.at(i,2) = -dc.at(i+1) / da.at(i+1);
                ma.at(i,0) = +dd.at(i+1) / da.at(i+1);
            }
            DoubleVector qamma(K);
            qamma[0] = dd.at(0);
            qamma[1] = dd.at(N);

            DoubleMatrix beta(K, N+1, 0.0);

            calculate1(N, K, ma, beta, qamma, rx);

            for (unsigned int i=0; i<=N; i++) u.at(m, i) = rx[i];
        }
        IPrinter::printVector(u.row(m));
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();

}

void calculate1(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x)
{
    for (unsigned int eq=0; eq<K; eq++)
    {
        for (unsigned int i=0; i<=N-K; i++)
        {
            for (unsigned int j=1; j<=K; j++)
            {
                beta.at(eq,i+j) = beta.at(eq,i+j) + beta.at(eq,i)*a.at(i,j);
            }
            qamma[eq] = qamma[eq] - beta.at(eq,i)*a.at(i,0);
        }
    }

    DoubleMatrix m(K,K);
    for (unsigned int j=0; j<K; j++)
    {
        for (unsigned int i=K-1; i != UINT_MAX; i--) m.at(j,(K-1)-i) = beta.at(j,N-i);
    }
    DoubleVector b(K);
    for (unsigned int i=0; i<K; i++) b.at(i) = qamma.at(i);

    DoubleVector x1(K);
    GaussianElimination(m, b, x1);

    for (unsigned int i=0; i<K; i++) x[N-i] = x1.at((K-1)-i);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x[i] += a.at(i,j)*x[i+j];
    }
}


void Example3::calculate1(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x)
{
    for (unsigned int eq=0; eq<K; eq++)
    {
        for (unsigned int i=0; i<=N-K; i++)
        {
            for (unsigned int j=1; j<=K; j++)
            {
                beta.at(eq,i+j) = beta.at(eq,i+j) + beta.at(eq,i)*a.at(i,j);
            }
            qamma[eq] = qamma[eq] - beta.at(eq,i)*a.at(i,0);
        }
    }

    DoubleMatrix m(K,K);
    for (unsigned int j=0; j<K; j++)
    {
        for (unsigned int i=K-1; i != UINT_MAX; i--) m.at(j,(K-1)-i) = beta.at(j,N-i);
    }
    DoubleVector b(K);
    for (unsigned int i=0; i<K; i++) b.at(i) = qamma.at(i);

    DoubleVector x1(K);
    GaussianElimination(m, b, x1);

    //printf("%12.8f %12.8f\n", x1[0], x1[1]);

    for (unsigned int i=0; i<K; i++) x[N-i] = x1.at((K-1)-i);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x[i] += a.at(i,j)*x[i+j];
    }
}
*/
