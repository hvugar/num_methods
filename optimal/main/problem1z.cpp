#include "problem1z.h"

void Problem1Z::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem1Z p;

    DoubleVector z(p.L);
    z[0] = 1.6;
    z[1] = 1.4;

    printf("Optimal:   %.10f %.10f\n", p.zs[0], p.zs[1]);
    printf("Initial:   %.10f %.10f\n", z[0], z[1]);
    double h = 0.001;
    DoubleVector gn1(p.L,0.0);
    IGradient::Gradient(&p, h, z, gn1);
    DoubleVector gn2 = gn1;
    gn2.L2Normalize();
    printf("Numerical: %12.8f %12.8f %12.8f %12.8f %12.8f\n", gn1[0], gn1[1], gn1.L2Norm(), gn2[0], gn2[1]);

    DoubleVector ga1(p.L,0.0);
    p.gradient(z, ga1);
    DoubleVector ga2 = ga1;
    ga2.L2Normalize();
    printf("Analytic:  %12.8f %12.8f %12.8f %12.8f %12.8f\n", ga1[0], ga1[1], ga1.L2Norm(), ga2[0], ga2[1]);
    puts("------------------------------------------");


    DoubleVector g1(p.L);
    p.gradient(z, g1);
    p.print(0, z, g1, 0.0, &p);

    ConjugateGradient g;
    g.setFunction(&p);
    g.setGradient(&p);
    g.setPrinter(&p);
    g.setEpsilon1(0.00000001);
    g.setEpsilon2(0.00000001);
    g.setEpsilon3(0.00000001);
    g.setR1MinimizeEpsilon(1.0, 0.00000001);
    g.setNormalize(true);
    g.calculate(z);
}

Problem1Z::Problem1Z()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    hx = 0.01;
    ht = 0.00005;
    N = 100;
    M = 20000;
    a = 1.0;

    alpha   = 1.0;
    lambda0 = 1.0;
    lambdal = 1.0;

    L = 2;

    xi.resize(L);
    xi[0] = 0.4;
    xi[1] = 0.7;

    Xi.resize(L);
    Xi[0] = 40;
    Xi[1] = 70;

    alpha1 = 1.0;
    alpha2 = 1.0;

    Te = 3.0;

    k.resize(L);
    k[0] = 1.1;
    k[1] = 1.2;

    zs.resize(L);
    zs[0] = 2.89;
    zs[1] = 2.81;
    pz = &zs;
    calculateV(zs);
}

Problem1Z::~Problem1Z()
{}

void Problem1Z::calculateU(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    calculateU2(m, ht, hx, M, N, alpha, lambda0, lambdal, a);
}

void Problem1Z::calculateP(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    calculateP2(m, ht, hx, M, N, alpha, lambda0, lambdal, a);
}

double Problem1Z::fx(const DoubleVector &z)
{
    pz = &z;
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

    //printf("%.10f %.10f\n", sum, norm);

    return alpha1*sum + alpha2*norm;
}

void Problem1Z::gradient(const DoubleVector &z, DoubleVector &g)
{
    pz = &z;
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
    pu = &u;

    DoubleMatrix psi;
    calculateP(psi, ht, hx, M, N, alpha, lambda0, lambdal, a);
    pp = &psi;

    for (unsigned int s = 0; s<L; s++)
    {
        g[s] = 0.0;
        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = k[s] * ( alpha*a*a*psi.at(m1, 0) - 2.0*alpha2*(vl(m1)-vs(m1)) );
            double g2 = k[s] * ( alpha*a*a*psi.at(m2, 0) - 2.0*alpha2*(vl(m2)-vs(m2)) );
            sum = sum + (g1 + g2);
        }
        sum = 0.5 * ht * sum;
        g[s] = sum;
    }
}

void Problem1Z::print(unsigned int i, const DoubleVector &z, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(alpha);
    C_UNUSED(alpha);
    printf("J[%d]: %.16f\n", i, fn->fx(z));
    printf("k:  %14.10f, %14.10f\n", z[0], z[1]);
    printf("g:  %14.10f, %14.10f\n", g[0], g[1]);
    puts("------------------------------------------");
}

double Problem1Z::initial(unsigned int i UNUSED_PARAM) const
{
    return Te;
}

double Problem1Z::vm(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

double Problem1Z::vl(unsigned int j UNUSED_PARAM) const
{
    const DoubleVector &z = *pz;
    const DoubleMatrix &u = *pu;
    return k[0]*(u.at(j, Xi[0])-z[0]) + k[1]*(u.at(j, Xi[1])-z[1]);
}

double Problem1Z::vr(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

struct CauchyProblemZU : public CauchyProblem
{
    CauchyProblemZU(double a, double hx, double alpha, double lambda0, double lambdal, double Te, unsigned int i, unsigned int N,
                   const DoubleVector &k, const DoubleVector &z, const std::vector<unsigned int> &Xi) :
        a(a), hx(hx), alpha(alpha), lambda0(lambda0), lambdal(lambdal), Te(Te), i(i), N(N), k(k), z(z), Xi(Xi) {}

    double f(double t UNUSED_PARAM, const DoubleVector &u) const
    {
        double a1 = -(a*a)/(hx*hx) - (lambda0*a*a)/(hx) - alpha;
        double b1 = +(a*a)/(hx*hx);
        double b2 = -(2.0*a*a)/(hx*hx) - alpha;
        double a2 = -(a*a)/(hx*hx) - (lambdal*a*a)/(hx) - alpha;
        double c1 = (lambda0*a*a)/hx;
        double c2 = (lambda0*a*a)/hx;

        double res = 0.0;
        if (i==0)
        {
            res = a1*u[0] + b1*u[1] + alpha*Te + c1*( k[0]*(u[Xi[0]]-z[0]) + k[1]*(u[Xi[1]]-z[1]) );
        }
        else if (i==N)
        {
            res = b1*u[N-1] + a2*u[N] + alpha*Te + c2*Te;
        }
        else
        {
            res = b1*u[i+1] + b2*u[i] + b1*u[i-1] + alpha*Te;
        }

        return res;
    }
    double a;
    double hx;
    double alpha;
    double lambda0;
    double lambdal;
    double Te;
    unsigned int i;
    unsigned int N;
    const DoubleVector &k;
    const DoubleVector &z;
    const std::vector<unsigned int> &Xi;
};

void Problem1Z::calculateU2(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    const DoubleVector &z = *pz;

    std::vector<CauchyProblem*> cps(N+1);
    for (unsigned int n=0; n<=N; n++)
    {
        CauchyProblemZU *cp = new CauchyProblemZU(a,hx, alpha, lambda0, lambdal, Te, n, N, k, z, Xi);
        cp->x0 = 0.0;
        cp->y0 = initial(n);
        cps[n] = cp;
    }
    CauchyProblem::rungeKutta(cps, 0.0, ht, M, u);
    u.transpose();
}

struct CauchyProblemZP : public CauchyProblem
{
    CauchyProblemZP(double a, double hx, double ht, double alpha, double alpha2, double lambda0, double lambdal, double Te, unsigned int i, unsigned int N,
                   const DoubleVector &k, const DoubleVector &z, const std::vector<unsigned int> &Xi, const DoubleMatrix &u) :
        a(a), hx(hx), ht(ht), alpha(alpha), alpha2(alpha2), lambda0(lambda0), lambdal(lambdal), Te(Te), i(i), N(N), k(k), z(z), Xi(Xi), u(u)
    { }

    double f(double t UNUSED_PARAM, const DoubleVector &p) const
    {
        double a1 = +(a*a)/(hx*hx) + (lambda0*a*a)/(hx) + alpha;
        double a2 = +(a*a)/(hx*hx) + (lambdal*a*a)/(hx) + alpha;
        double b1 = -(a*a)/(hx*hx);
        double b2 = +(2.0*a*a)/(hx*hx) + alpha;

        double res = 0.0;
        if (i==0)
        {
            res = a1*p[0] + b1*p[1];
        }
        else if (i==N)
        {
            res = b1*p[N-1] + a2*p[N];
        }
        else
        {
            res = b1*p[i+1] + b2*p[i] + b1*p[i-1];
        }

        if (i==Xi[0]) { res += -lambda0*a*a*k[0]*(1.0/hx)*p[0]; }
        if (i==Xi[1]) { res += -lambda0*a*a*k[1]*(1.0/hx)*p[0]; }


        unsigned int j = (unsigned int)(t/ht);
        double dd = p2->vl(j) - p2->vs(j);

        if (i==Xi[0]) { res += 2.0 * alpha2 * dd * k[0] * (1.0/hx); }
        if (i==Xi[1]) { res += 2.0 * alpha2 * dd * k[1] * (1.0/hx); }

        return res;
    }
    double a;
    double hx;
    double ht;
    double alpha;
    double alpha2;
    double lambda0;
    double lambdal;
    double Te;
    unsigned int i;
    unsigned int N;
    const DoubleVector &k;
    const DoubleVector &z;
    const std::vector<unsigned int> &Xi;
    const DoubleMatrix &u;
    Problem1Z *p2;
};

void Problem1Z::calculateP2(DoubleMatrix &psi, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    const DoubleVector &z = *pz;
    const DoubleMatrix &u = *pu;

    std::vector<CauchyProblem*> cps(N+1);

    for (unsigned int n=0; n<=N; n++)
    {
        CauchyProblemZP *cp = new CauchyProblemZP(a,hx, ht, alpha, alpha2, lambda0, lambdal, Te, n, N, k, z, Xi, u);
        cp->x0 = 1.0;
        cp->y0 = -2.0*alpha1*mu(n)*(u.at(M,n)-V[n]);
        cp->p2 = this;
        cps[n] = cp;
    }
    CauchyProblem::rungeKutta(cps, 1.0, -ht, M, psi);
    psi.transpose();
}

double Problem1Z::mu(unsigned int i UNUSED_PARAM) const
{
    return 1.0;
}

void Problem1Z::calculateV(const DoubleVector &z)
{
    //pz = &z;
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
    V = u.row(M);
}

double Problem1Z::vs(double j) const
{
    const DoubleMatrix &u = *pu;
    return k[0]*(u(j, Xi[0])-zs[0]) + k[1]*(u(j, Xi[1])-zs[1]);
}
