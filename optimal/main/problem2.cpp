#include "problem2.h"

void Problem2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem2 p;

    DoubleVector k(p.L);
    k[0] = 1.5;
    k[1] = 1.6;

    //    DoubleMatrix m1;
    //    p.calculateU(m1, p.ht, p.hx, p.M, p.N, p.alpha, p.lambda0, p.lambdal, p.a);
    //    IPrinter::printMatrix(m1);
    //    puts("---");
    //    DoubleMatrix m2;
    //    p.calculateU1(m2, p.ht, p.hx, p.M, p.N, p.alpha, p.lambda0, p.lambdal, p.a);
    //    m2.transpose();
    //    IPrinter::printMatrix(m2);
    //    return;

    //k[0] = -0.96916037;
    //k[1] = +1.07008630;

    //    puts("Numerical");
    //    double h = 0.001;
    //    DoubleVector gn(p.L,0.0);
    //    IGradient::Gradient(&p, h, k, gn);
    //    printf("%.10f %.10f %.10f %.10f\n", k[0], k[1], p.ks[0], p.ks[1]);
    //    printf("%.10f %.10f %.10f\n", gn[0], gn[1], gn.L2Norm());
    //    gn.L2Normalize();
    //    printf("%.10f %.10f\n", gn[0], gn[1]);

    //    puts("Analitic");
    //    DoubleVector ga(p.L,0.0);
    //    p.gradient(k, ga);
    //    printf("%.10f %.10f %.10f %.10f\n", k[0], k[1], p.ks[0], p.ks[1]);
    //    printf("%.10f %.10f %.10f\n", ga[0], ga[1], ga.L2Norm());
    //    ga.L2Normalize();
    //    printf("%.10f %.10f\n", ga[0], ga[1]);


//    printf("J[0]: %.10f\n", p.fx(k));
    DoubleVector g1(p.L);
    p.gradient(k, g1);
    p.print(0, k, g1, 0.0, &p);

    ConjugateGradient g;
    g.setFunction(&p);
    g.setGradient(&p);
    g.setPrinter(&p);
    g.setEpsilon1(0.00000001);
    g.setEpsilon2(0.00000001);
    g.setEpsilon3(0.00000001);
    g.setR1MinimizeEpsilon(0.001, 0.00000001);
    g.setNormalize(true);
    g.calculate(k);

}

Problem2::Problem2()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    hx = 0.005;
    ht = 0.0000125;
    N = 200;
    M = 80000;
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

    z.resize(L);
    z[0] = 2.89;
    z[1] = 2.81;

    alpha1 = 1.0;
    alpha2 = 0.0;

    Te = 3.0;

    //DoubleVector ks;
    ks.resize(L);
    ks[0] = 1.1;
    ks[1] = 1.2;
    calculateV(ks);
    //IPrinter::printVector(V);
}

Problem2::~Problem2()
{}

void Problem2::calculateU0(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    calculateU1(m, ht, hx, M, N, alpha, lambda0, lambdal, a);
}

void Problem2::calculateP0(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    calculateP1(m, ht, hx, M, N, alpha, lambda0, lambdal, a);
}

double Problem2::fx(const DoubleVector &k)
{
    pk = &k;
    DoubleMatrix u;
    calculateU0(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
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
        double vm1 = k[0]*(u.at(m1,Xi[0])-z[0]) + k[1]*(u.at(m1,Xi[1])-z[1]);
        double vm2 = k[0]*(u.at(m2,Xi[0])-z[0]) + k[1]*(u.at(m2,Xi[1])-z[1]);

        double vs1 = vs(m1);
        double vs2 = vs(m2);

        norm = norm + ((vm1-vs1)*(vm1-vs1) + (vm2-vs2)*(vm2-vs2));
    }
    norm = 0.5*ht*norm;

    return alpha1*sum + alpha2*norm;
}

void Problem2::gradient(const DoubleVector &k, DoubleVector &g)
{
    pk = &k;
    DoubleMatrix u;
    calculateU0(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
    pu = &u;

    DoubleMatrix psi;
    calculateP0(psi, ht, hx, M, N, alpha, lambda0, lambdal, a);
    pp = &psi;

    for (unsigned int s = 0; s<L; s++)
    {
        g[s] = 0.0;
        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = (u.at(m1, Xi[s]) - z[s]) * (-alpha*a*a*psi.at(m1, 0) + 2.0*alpha2*( k[0]*(u.at(m1,Xi[0])-z[0])+k[1]*(u.at(m1,Xi[1])-z[1]) - vs(m1) ));
            double g2 = (u.at(m2, Xi[s]) - z[s]) * (-alpha*a*a*psi.at(m2, 0) + 2.0*alpha2*( k[0]*(u.at(m2,Xi[0])-z[0])+k[1]*(u.at(m2,Xi[1])-z[1]) - vs(m2) ));
            sum = sum + (g1 + g2);
        }
        sum = 0.5 * ht * sum;
        g[s] = sum;
    }
}

void Problem2::print(unsigned int i, const DoubleVector &k, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(alpha);

    Problem2 *hc = dynamic_cast<Problem2*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(k));
    //printf("Norm: %.16f Alpha: %.16f %.16f\n", hc->norm(x), hc->alpha, alpha);
    //printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", O[0], O[1], O[2], O[3], O[4], O[5]);
    printf("k: %12.8f, %12.8f\n", k[0], k[1]);
    printf("g: %12.8f, %12.8f\n", g[0], g[1]);
}

double Problem2::initial(unsigned int i UNUSED_PARAM) const
{
    return 1.0;
}

double Problem2::vm(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

double Problem2::vl(unsigned int j UNUSED_PARAM) const
{
    const DoubleVector &k = *pk;
    const DoubleMatrix &u = *pu;
    return k[0]*(u.at(j, Xi[0])-z[0]) + k[1]*(u.at(j, Xi[1])-z[1]);
}

double Problem2::vr(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

void Problem2::calculateU(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    const DoubleVector &k = *pk;

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

            ra(0, Xi[0]) = -k[0] * ((lambda0*a*a*ht)/(hx));
            ra(0, Xi[1]) = -k[1] * ((lambda0*a*a*ht)/(hx));
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

struct CauchyProblemU : public CauchyProblem
{
    CauchyProblemU(double a, double hx, double alpha, double lambda0, double lambdal, double Te, unsigned int i, unsigned int N,
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

void Problem2::calculateU1(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    const DoubleVector &k = *pk;

    std::vector<CauchyProblem*> cps(N+1);
    for (unsigned int n=0; n<=N; n++)
    {
        CauchyProblem *cp = new CauchyProblemU(a,hx, alpha, lambda0, lambdal, Te, n, N, k, z, Xi);
        cp->x0 = 0.0;
        cp->y0 = initial(n);
        cps[n] = cp;
    }
    CauchyProblem::rungeKutta(cps, 0.0, ht, M, u);
    u.transpose();
}

void Problem2::calculateP(DoubleMatrix &psi, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR UNUSED_PARAM, double a)
{
    const DoubleVector &k = *pk;
    const DoubleMatrix &u = *pu;

    psi.clear();
    psi.resize(M+1, N+1, 0.0);

    DoubleVector d1(N+1);
    DoubleVector rx(N+1);

    for (unsigned int m1=0; m1<=M; m1++)
    {
        unsigned int m = M-m1;

        if (m==M)
        {
            for (unsigned int n=0; n<=N; n++)
                psi.at(M,n) = -2.0*alpha1*mu(n)*(u.at(M,n) - V[n]);
        }
        else
        {
            DoubleMatrix p(N+1, N+1, 0.0);

            p.at(0,0) = -1.0 - (a*a*ht)/(hx*hx) - (lambdaL*a*a*ht)/hx - lambdaM*ht;
            p.at(0,1) = (a*a*ht)/(hx*hx);
            d1[0]     = -psi.at(m+1,0);

            //p(0, Xi[0]) = -k[0] * ((lambdaL*a*a*ht)/(hx));
            //p(0, Xi[1]) = -k[1] * ((lambdaL*a*a*ht)/(hx));

            for (unsigned int n=1; n<=N-1; n++)
            {
                p.at(n,n-1) = (a*a*ht)/(hx*hx);
                p.at(n,n)   = -1.0 - 2.0*((a*a)*ht)/(hx*hx) - lambdaM*ht;
                p.at(n,n+1) = (a*a*ht)/(hx*hx);

                if (n>1) p.at(n,0) = 0.0;
                if (n==Xi[0]) p.at(n,0) = p.at(n,0) + (1/hx)*lambdaM*a*a*ht*k[0];
                if (n==Xi[1]) p.at(n,0) = p.at(n,0) + (1/hx)*lambdaM*a*a*ht*k[1];

                d1[n] = -psi.at(m+1, n);

                //                if (n==Xi[0])
                //                {
                //                    double A = 0.0;
                //                    for (unsigned int i=0; i<L; i++)
                //                    {
                //                        A += k[i] * (u.at(m, Xi[i]) - z[i]);
                //                    }
                //                    d1[n] += 2*alpha2*ht*(1.0/hx)*k[0]*A;
                //                }

                //                if (n==Xi[1])
                //                {
                //                    double A = 0.0;
                //                    for (unsigned int i=0; i<L; i++)
                //                    {
                //                        A += k[i] * (u.at(m, Xi[i]) - z[i]);
                //                    }
                //                    d1[n] += 2*alpha2*ht*(1.0/hx)*k[1]*A;
                //                }

                for (unsigned int j=0; j<L; j++)
                {
                    if (n==Xi[j])
                    {
                        double A = 0.0;
                        for (unsigned int i=0; i<L; i++)
                        {
                            A = A + k[i] * (u.at(m, Xi[i]) - z[i]);
                        }
                        d1[n] += 2*alpha2*ht*(1.0/hx)*k[j]*A;
                    }
                }

                //                for (unsigned int j=0; j<L; j++)
                //                {
                //                    if (n==Xi[j])
                //                    {
                //                        for (unsigned i=0; i<L; i++)
                //                        {
                //                            if (i!=j)
                //                            {
                //                                d1[n] += (1/hx)*2*alpha2*ht*k[j]*k[i]*(u.at(m,Xi[i]) - z[i]);
                //                            }
                //                        }
                //                    }
                //                }

                //                for (unsigned int i=0; i<L; i++)
                //                {
                //                    if (n==Xi[i])
                //                    {
                //                        d1[n] += (1/hx)*2*alpha2*ht*k[i]*k[i]*(u.at(m,n) - z[i]);
                //                    }
                //                }

            }

            p.at(N,N-1) = (a*a*ht)/(hx*hx);
            p.at(N,N)   = -1.0 - (a*a*ht)/(hx*hx) - (lambdaL*a*a*ht)/hx - lambdaM*ht;
            d1[N]       = -psi.at(m+1,N);

            GaussianElimination(p, d1, rx);

            for (unsigned int n=0; n<=N; n++) psi.at(m,n) = rx[n];
        }
    }
}

struct CauchyProblemP : public CauchyProblem
{
    CauchyProblemP(double a, double hx, double ht, double alpha, double alpha2, double lambda0, double lambdal, double Te, unsigned int i, unsigned int N,
                   const DoubleVector &k, const DoubleVector &z, const std::vector<unsigned int> &Xi, const DoubleMatrix &u) :
        a(a), hx(hx), ht(ht), alpha(alpha), alpha2(alpha2), lambda0(lambda0), lambdal(lambdal), Te(Te), i(i), N(N), k(k), z(z), Xi(Xi), u(u)  {}

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
        double dd = k[0]*(u.at(j, Xi[0])-z.at(0)) + k[1]*(u.at(j, Xi[1])-z.at(1)) - p2->vs(j);

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
    Problem2 *p2;
};

void Problem2::calculateP1(DoubleMatrix &psi, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    const DoubleVector &k = *pk;
    const DoubleMatrix &u = *pu;

    std::vector<CauchyProblem*> cps(N+1);

    for (unsigned int n=0; n<=N; n++)
    {
        CauchyProblemP *cp = new CauchyProblemP(a,hx, ht, alpha, alpha2, lambda0, lambdal, Te, n, N, k, z, Xi, u);
        cp->x0 = 1.0;
        cp->y0 = -2.0*alpha1*mu(n)*(u.at(M,n)-V[n]);
        cp->p2 = this;
        cps[n] = cp;
    }
    CauchyProblem::rungeKutta(cps, 1.0, -ht, M, psi);
    psi.transpose();
}

double Problem2::mu(unsigned int i UNUSED_PARAM) const
{
    return 1.0;
}

void Problem2::calculateV(const DoubleVector &k)
{
    pk = &k;
    DoubleMatrix u;
    calculateU0(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
    V = u.row(M);
}

double Problem2::vs(double j) const
{
    const DoubleMatrix &u = *pu;
    return ks[0]*(u(j, Xi[0])-z[0]) + ks[1]*(u(j, Xi[1])-z[1]);
}

