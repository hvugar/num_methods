#include "problem1kzx.h"

void Problem1KZX::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    Problem1KZX p;

    DoubleVector k(p.L);
    k[0] = 2.0;
    k[1] = 2.0;

//    printf("Optimal:   %.10f %.10f\n", p.ks[0], p.ks[1]);
//    printf("Initial:   %.10f %.10f\n", k[0], k[1]);
//    double h = 0.001;
//    DoubleVector gn1(p.L,0.0);
//    IGradient::Gradient(&p, h, k, gn1);
//    DoubleVector gn2 = gn1;
//    gn2.L2Normalize();
//    printf("Numerical: %12.8f %12.8f %12.8f %12.8f %12.8f\n", gn1[0], gn1[1], gn1.L2Norm(), gn2[0], gn2[1]);

//    DoubleVector ga1(p.L,0.0);
//    p.gradient(k, ga1);
//    DoubleVector ga2 = ga1;
//    ga2.L2Normalize();
//    printf("Analytic:  %12.8f %12.8f %12.8f %12.8f %12.8f\n", ga1[0], ga1[1], ga1.L2Norm(), ga2[0], ga2[1]);
//    puts("------------------------------------------");

//    DoubleVector g1(p.L);
//    p.gradient(k, g1);
//    p.print(0, k, g1, 0.0, &p);

//    ConjugateGradient g;
//    g.setFunction(&p);
//    g.setGradient(&p);
//    g.setPrinter(&p);
//    g.setEpsilon1(0.00000001);
//    g.setEpsilon2(0.00000001);
//    g.setEpsilon3(0.00000001);
//    g.setR1MinimizeEpsilon(0.1, 0.00000001);
//    g.setNormalize(true);
//    g.calculate(k);

//    DoubleVector g2(p.L);
//    p.gradient(k, g2);
//    p.print(0, k, g2, 0.0, &p);
}

Problem1KZX::Problem1KZX()
{
    hx = 0.001;
    ht = 0.001;
    N = 1000;
    M = 1000;
    a = 1.0;
    alpha = 1.0;
    lambda0 = 1.0;
    lambdal = 1.0;
    L = 2;
    alpha1 = 10.0;
    alpha2 = 1.0;
    Te = 3.0;
    Ti = 2.0;

    // initialize --------------------------------------------------------
    ks << 1.50 << 1.70;
    zs << 2.89 << 2.81;
    es << 0.40 << 0.70;

    pk = &ks;
    pz = &zs;
    pe = &es;

    DoubleMatrix u;
    calculateU3(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
    V = u.row(M);
    u.clear();

    pk = NULL;
    pz = NULL;
    pe = NULL;
    // initialize --------------------------------------------------------

    pk = new DoubleVector;
    (*pk) << 2.5 << 2.6;

    pz = new DoubleVector;
    (*pz) << 2.89 << 2.81;

    pe = new DoubleVector;
    (*pe) << 0.40 << 0.70;
}

double Problem1KZX::fx(const DoubleVector &x)
{
    pk = const_cast<DoubleVector*>(&x);

    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
    pu = &u;

    double sum = 0.0;
    for (unsigned int n=0; n<=N-1; n++)
    {
        unsigned int n1 = n;
        unsigned int n2 = n + 1;
        double f1 = mu(n1) * (u.at(M, n1) - V[n1])*(u.at(M, n1) - V[n1]);
        double f2 = mu(n2) * (u.at(M, n2) - V[n2])*(u.at(M, n2) - V[n2]);
        sum = sum + (f1 + f2);
    }
    sum = 0.5*hx*sum;

    double norm = 0.0;
    for (unsigned int m=0; m<=M-1; m++)
    {
        unsigned int m1 = m + 0;
        unsigned int m2 = m + 1;

        double vm1 = vl(m1);
        double vm2 = vl(m2);

        double vs1 = vs(m1);
        double vs2 = vs(m2);

        norm = norm + ((vm1-vs1)*(vm1-vs1) + (vm2-vs2)*(vm2-vs2));
    }
    norm = 0.5*ht*norm;

    return alpha1*sum + alpha2*norm;
}

void Problem1KZX::gradient(const DoubleVector &x UNUSED_PARAM, DoubleVector &g UNUSED_PARAM)
{
    pk = const_cast<DoubleVector*>(&x);
    DoubleVector &z = *pz;
    DoubleVector &e = *pe;

    for (unsigned int s = 0; s<L; s++) g[s] = 0.0;

    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
    pu = &u;

    DoubleMatrix psi;
    calculateP3(psi, ht, hx, M, N, alpha, lambda0, lambdal, a);
    pp = &psi;

    DoubleVector E(L);
    E[0] = (unsigned int) round(e[0]/hx);
    E[1] = (unsigned int) round(e[1]/hx);
    for (unsigned int s = 0; s<L; s++)
    {
        g[s] = 0.0;
        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = (u.at(m1, E[s]) - z[s]) * ( -alpha*a*a*psi.at(m1, 0) + 2.0*alpha2*(vl(m1)-vs(m1)) );
            double g2 = (u.at(m2, E[s]) - z[s]) * ( -alpha*a*a*psi.at(m2, 0) + 2.0*alpha2*(vl(m2)-vs(m2)) );
            sum = sum + (g1 + g2);
        }
        sum = 0.5 * ht * sum;
        g[s] = sum;
    }
}

void Problem1KZX::print(unsigned int i, const DoubleVector &k UNUSED_PARAM, const DoubleVector &g UNUSED_PARAM, double alpha UNUSED_PARAM, RnFunction *fn UNUSED_PARAM) const
{
    C_UNUSED(alpha);
    printf("J[%d]: %.16f\n", i, fn->fx(k));

    DoubleVector gr1(L);
    const_cast<Problem1KZX*>(this)->gradient(k, gr1);
    gr1.L2Normalize();

    DoubleVector gr2(L);
    IGradient::Gradient(const_cast<Problem1KZX*>(this), 0.001, k, gr2);
    gr2.L2Normalize();


    printf("k:  %14.10f, %14.10f\n", k[0], k[1]);
    printf("g:  %14.10f, %14.10f\n", g[0], g[1]);
    printf("g:  %14.10f, %14.10f\n", gr1[0], gr1[1]);
    printf("g:  %14.10f, %14.10f\n", gr2[0], gr2[1]);
    puts("------------------------------------------");
}

void Problem1KZX::calculateU(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    calculateU3(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
}

void Problem1KZX::calculateU1(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    DoubleVector &k = *pk;
    DoubleVector &z = *pz;
    DoubleVector &e = *pe;

    std::vector<unsigned int>E(L);
    E[0] = (unsigned int)round(e[0]/hx);
    E[1] = (unsigned int)round(e[1]/hx);

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
            rb[0] = u.at(m-1,0) + alpha*ht*vm(m) - ((lambda0*a*a*ht)/(hx))*(k[0]*z[0]+k[1]*z[1]);

            ra(0, E[0]) = -k[0] * ((lambda0*a*a*ht)/(hx));
            ra(0, E[1]) = -k[1] * ((lambda0*a*a*ht)/(hx));
            //rb[0] = rb[0] - ((lambda0*a*a*ht)/(hx))*((k[0]*z[0]+k[1]*z[1]));

            for (unsigned int i=1; i<=N-1; i++)
            {
                ra(i,i-1) = -(a*a*ht)/(hx*hx);
                ra(i,i) = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + alpha*ht;
                ra(i,i+1) = -(a*a*ht)/(hx*hx);
                rb[i] = u.at(m-1, i) + alpha*ht*vm(m);
            }

            ra(N,N-1) = -(a*a*ht)/(hx*hx);
            ra(N,N)   = 1.0 + (a*a*ht)/(hx*hx) + (lambda0*a*a*ht)/hx + alpha*ht;
            rb[N] = u.at(m-1,N) + ((lambdal*a*a*ht)/(hx))*vr(m) + alpha*ht*vm(m);

            GaussianElimination(ra, rb, rx);

            for (unsigned int i=0; i<=N; i++) u.at(m,i) = rx[i];
        }
    }
}

void Problem1KZX::calculateU3(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    DoubleVector &k = *pk;
    DoubleVector &z = *pz;
    DoubleVector &e = *pe;

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);

    DoubleVector de(L);

//    std::vector<unsigned int> E(L);
//    E[0] = (unsigned int)round(e[0]/hx);
//    E[1] = (unsigned int)round(e[1]/hx);

    for (unsigned int m=0; m<=M; m++)
    {
        if (m==0)
        {
            for (unsigned int i=0;i<=N; i++) u.at(0,i) = initial(i);
        }
        else
        {
            da[0] = 0.0;
            db[0] = 1.0+(a*a*ht)/(hx*hx) + (lambda0*a*a*ht)/hx + alpha*ht;
            dc[0] = -(a*a*ht)/(hx*hx);
            dd[0] = u.at(m-1,0) + alpha*ht*Te - ((lambda0*a*a*ht)/hx)*(k[0]*z[0] + k[1]*z[1]);

            de[0] = (-lambda0*a*a*ht*k[0])/hx;
            de[1] = (-lambda0*a*a*ht*k[1])/hx;

            for (unsigned int n=1; n<=N-1; n++)
            {
                da[n] = -(a*a*ht)/(hx*hx);
                db[n] = 1.0+(2.0*a*a*ht)/(hx*hx) + alpha*ht;
                dc[n] = -(a*a*ht)/(hx*hx);
                dd[n] = u.at(m-1,n) + alpha*ht*Te;
            }

            da[N] = -(a*a*ht)/(hx*hx);
            db[N] = 1.0+(a*a*ht)/(hx*hx) + (lambdal*a*a*ht)/hx + alpha*ht;
            dc[N] = 0.0;
            dd[N] = u.at(m-1,N) + (lambdal*a*a*ht*Te)/hx + alpha*ht*Te;

            qovmaFirstRow(da.data(), db.data(), dc.data(), dd.data(), rx.data(), N+1, de.data());
            //qovmaE(da.data(), db.data(), dc.data(), dd.data(), rx.data(), N+1, de.data(), E.data(), L);

            for (unsigned int i=0; i<=N; i++) u.at(m, i) = rx[i];
        }
    }
    puts("1");
}

void Problem1KZX::calculateP3(DoubleMatrix &p, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    DoubleVector &k = *pk;
    DoubleVector &z = *pz;
    DoubleVector &e = *pe;
    DoubleMatrix &u = *pu;

    p.clear();
    p.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);
    DoubleVector de(N+1);

    std::vector<unsigned int> E(L);
    E[0] = (unsigned int)round(e[0]/hx);
    E[1] = (unsigned int)round(e[1]/hx);

    for (unsigned int m1=0; m1<=M; m1++)
    {
        unsigned int m = M-m1;

        if (m == M)
        {
           for (unsigned int n=0; n<=N; n++) p.at(M,n) = -2.0*alpha1*mu(n)*(u.at(M,n)-V.at(n));
        }
        else
        {
            da[0] = 0.0;
            db[0] = -1.0-(a*a*ht)/(hx*hx) - (lambda0*a*a*ht)/hx - alpha*ht;
            dc[0] = +(a*a*ht)/(hx*hx);
            dd[0] = -p.at(m+1,0);

            de[0] = 0.0;
            de[1] = 0.0;

            for (unsigned int n=1; n<=N-1; n++)
            {
                da[n] = +(a*a*ht)/(hx*hx);
                db[n] = -1.0-(2.0*a*a*ht)/(hx*hx) - alpha*ht;
                dc[n] = +(a*a*ht)/(hx*hx);
                dd[n] = -p.at(m+1,n) + alpha*ht*Te;

                double sum = k[0]*(u.at(m, E[0]) - z[0]) + k[1]*(u.at(m, E[1]) - z[1]);

                if (n==E[0]) dd[n] += 2*alpha2*ht*k[0]*sum*(1.0/hx);
                if (n==E[1]) dd[n] += 2*alpha2*ht*k[1]*sum*(1.0/hx);

                if (n==E[0]) de[n] = lambda0*a*a*ht*k[0]*(1.0/hx);
                if (n==E[1]) de[n] = lambda0*a*a*ht*k[1]*(1.0/hx);
            }

            da[N] = +(a*a*ht)/(hx*hx);
            db[N] = -1.0-(a*a*ht)/(hx*hx) - (lambdal*a*a*ht)/hx - alpha*ht;
            dc[N] = 0.0;
            dd[N] = -p.at(m+1,N);
            de[N] = 0.0;

            qovma2(da.data(), db.data(), dc.data(), dd.data(), rx.data(), N+1, de.data());

            for (unsigned int i=0; i<=N; i++) p.at(m, i) = rx[i];
        }
    }
}

double Problem1KZX::initial(unsigned int i UNUSED_PARAM) const
{
    return Ti;
}

double Problem1KZX::vm(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

double Problem1KZX::vr(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

double Problem1KZX::vl(unsigned int j UNUSED_PARAM) const
{
    const DoubleVector &k = *pk;
    const DoubleVector &z = *pz;
    const DoubleVector &e = *pe;
    const DoubleMatrix &u = *pu;

    DoubleVector E;
    E << (unsigned int)round(e[0]/hx);
    E << (unsigned int)round(e[1]/hx);

    return k[0]*(u.at(j, E[0])-z[0]) + k[1]*(u.at(j, E[1])-z[1]);
}

double Problem1KZX::vs(double j) const
{
    const DoubleMatrix &u = *pu;

    DoubleVector E;
    E << (unsigned int)round(es[0]/hx);
    E << (unsigned int)round(es[1]/hx);

    return ks[0]*(u(j, E[0])-zs[0]) + ks[1]*(u(j, E[1])-zs[1]);
}
