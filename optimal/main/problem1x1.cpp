#include "problem1x1.h"

void qovma(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e, unsigned int *E, unsigned int L)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);
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
                k[s][0] = -e[s]/b[0];
                //if (i%2==0) k[s][i] *= -1.0;
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

    const unsigned int j = (unsigned)0-1;
    for (unsigned int i=n-1; i != j; i--)
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

    free(q);
    free(q);
    for (unsigned int s=0; s<L; s++) free(k[s]);
    free(k);
}

void Problem1X1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem1X1 p;

    // 0.45 0.75
    // 0.45 0.65
    // 0.35 0.65
    DoubleVector xi(p.L);
    xi[0] = 0.45;
    xi[1] = 0.65;

//    DoubleMatrix m(p.N+1, p.N+1,0.0);
//    DoubleVector x(p.L);
//    for (unsigned int i1=0; i1<=p.N; i1++)
//    {
//        x[0] = i1*p.hx;
//        for (unsigned int i2=0; i2<=p.N; i2++)
//        {
//            x[1] = i2*p.hx;
//            m.at(i1,i2) = p.fx(x);
//            printf("%f %f %.10f\n", x[0], x[1], m.at(i1,i2));
//        }
//    }
//    FILE *f = fopen("d:\\data1.txt", "w");
//    IPrinter::printMatrix(m, 100, 100, NULL,f);
//    fclose(f);
//    return;

    printf("Optimal:   %.10f %.10f\n", p.xis[0], p.xis[1]);
    printf("Initial:   %.10f %.10f\n", xi[0], xi[1]);

    DoubleVector ga1(p.L,0.0);
    p.gradient(xi, ga1);
    DoubleVector ga2 = ga1;
    ga2.L2Normalize();
    printf("Analytic:  %12.8f %12.8f %12.8f %12.8f %12.8f\n", ga1[0], ga1[1], ga1.L2Norm(), ga2[0], ga2[1]);

    double h = 0.01;
    DoubleVector gn1(p.L,0.0);
    IGradient::Gradient(&p, h, xi, gn1);
    DoubleVector gn2 = gn1;
    gn2.L2Normalize();
    printf("Numerical: %12.8f %12.8f %12.8f %12.8f %12.8f\n", gn1[0], gn1[1], gn1.L2Norm(), gn2[0], gn2[1]);
    puts("------------------------------------------");

    DoubleVector g1(p.L);
    p.gradient(xi, g1);
    p.print(0, xi, g1, 0.0, &p);

    ConjugateGradient g;
    g.setFunction(&p);
    g.setGradient(&p);
    g.setPrinter(&p);
    g.setProjection(&p);
    g.setEpsilon1(0.00000001);
    g.setEpsilon2(0.00000001);
    g.setEpsilon3(0.00000001);
    g.setR1MinimizeEpsilon(1.0, 0.00000001);
    g.setNormalize(true);
    g.calculate(xi);

    DoubleVector g2(p.L);
    p.gradient(xi, g2);
    p.print(0, xi, g2, 0.0, &p);
}

Problem1X1::Problem1X1()
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

    k.resize(L);
    k[0] = 1.1;
    k[1] = 1.2;

    z.resize(L);
    z[0] = 2.89;
    z[1] = 2.81;

    alpha1 = 1.0;
    alpha2 = 0.0;

    Te = 3.0;
    Ti = 2.0;

    xis.resize(L);
    xis[0] = 0.4;
    xis[1] = 0.7;
    calculateV(xis);
}

Problem1X1::~Problem1X1()
{}

double Problem1X1::fx(const DoubleVector &xi)
{
    pxi = &xi;

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
//    for (unsigned int m=0; m<=M-1; m++)
//    {
//        unsigned int m1 = m + 0;
//        unsigned int m2 = m + 1;

//        double vm1 = vl(m1);
//        double vm2 = vl(m2);

//        double vs1 = vs(m1);
//        double vs2 = vs(m2);

//        norm = norm + ((vm1-vs1)*(vm1-vs1) + (vm2-vs2)*(vm2-vs2));
//    }
//    norm = 0.5*ht*norm;

    return alpha1*sum + alpha2*norm;
}

void Problem1X1::gradient(const DoubleVector &xi, DoubleVector &g)
{
    pxi = &xi;
//    std::vector<unsigned int> Xi(L);
//    Xi[0] = (unsigned int) round(xi[0] * N);
//    Xi[1] = (unsigned int) round(xi[1] * N);

    for (unsigned int s = 0; s<L; s++) g[s] = 0.0;

    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
    pu = &u;

    DoubleMatrix psi;
    calculateP(psi, ht, hx, M, N, alpha, lambda0, lambdal, a);
    pp = &psi;

    for (unsigned int s = 0; s<L; s++)
    {
        unsigned int r = (unsigned int)floor(xi[s]*N);

//        double h1 = fabs(xi[s] - (r+1)*hx)/hx;
//        double h2 = fabs(xi[s] - (r+0)*hx)/hx;

        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            //double g1 = k[s] * ((h2*u.at(m1, r+1) - h1*u.at(m1, r+0))/(hx)) * ( -alpha*a*a*psi.at(m1, 0) /*+ 2.0*alpha2*(vl(m1)-vs(m1))*/ );
            //double g2 = k[s] * ((h2*u.at(m2, r+1) - h1*u.at(m2, r+0))/(hx)) * ( -alpha*a*a*psi.at(m2, 0) /*+ 2.0*alpha2*(vl(m2)-vs(m2))*/ );
//            double g1 = k[s] * ((u.at(m1, Xi[s]+1) - u.at(m1, Xi[s]-1))/(2.0*hx)) * ( -alpha*a*a*psi.at(m1, 0) /*+ 2.0*alpha2*(vl(m1)-vs(m1))*/ );
//            double g2 = k[s] * ((u.at(m2, Xi[s]+1) - u.at(m2, Xi[s]-1))/(2.0*hx)) * ( -alpha*a*a*psi.at(m2, 0) /*+ 2.0*alpha2*(vl(m2)-vs(m2))*/ );
            double g1 = k[s] * ((u.at(m1, r+1) - u.at(m1, r+0))/(hx)) * ( -lambda0*a*a*psi.at(m1, 0) /*+ 2.0*alpha2*(vl(m1)-vs(m1))*/ );
            double g2 = k[s] * ((u.at(m2, r+1) - u.at(m2, r+0))/(hx)) * ( -lambda0*a*a*psi.at(m2, 0) /*+ 2.0*alpha2*(vl(m2)-vs(m2))*/ );
            sum = sum + (g1 + g2);
        }
        g[s] = 0.5 * ht * sum;
    }
}

void Problem1X1::print(unsigned int i, const DoubleVector &x, const DoubleVector &g UNUSED_PARAM, double alpha, RnFunction *fn) const
{
    C_UNUSED(alpha);
    printf("J[%d]: %.16f\n", i, fn->fx(x));
//    DoubleVector gr(x.size());
//    const_cast<Problem1X1*>(this)->gradient(x, gr);
    printf("k:  %14.10f, %14.10f\n", x[0], x[1]);
//    printf("g:  %14.10f, %14.10f\n", gr[0], gr[1]);
//    puts("------------------------------------------");

    DoubleVector ga1(L,0.0);
    const_cast<Problem1X1*>(this)->gradient(x, ga1);
    DoubleVector ga2 = ga1;
    ga2.L2Normalize();
    printf("Analytic:  %12.8f %12.8f %12.8f %12.8f %12.8f\n", ga1[0], ga1[1], ga1.L2Norm(), ga2[0], ga2[1]);

    double h = 0.01;
    DoubleVector gn1(L,0.0);
    IGradient::Gradient(const_cast<Problem1X1*>(this), h, x, gn1);
    DoubleVector gn2 = gn1;
    gn2.L2Normalize();
    printf("Numerical: %12.8f %12.8f %12.8f %12.8f %12.8f\n", gn1[0], gn1[1], gn1.L2Norm(), gn2[0], gn2[1]);
    puts("------------------------------------------");

}

void Problem1X1::project(DoubleVector &x, int index)
{
    if (x[index] <= 0.0) x[index] = 2.0*hx;
    if (x[index] >= 1.0) x[index] = 1.0 - 2.0*hx;
}

double Problem1X1::initial(unsigned int i UNUSED_PARAM) const
{
    return Ti;
}

double Problem1X1::vm(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

//double Problem1X1::vl(unsigned int j UNUSED_PARAM) const
//{
//    const DoubleMatrix &u  = *pu;

//    const DoubleVector &xi = *pxi;
//    std::vector<unsigned int> Xi(L);
//    Xi[0] = (unsigned int) round(xi[0] * N);
//    Xi[1] = (unsigned int) round(xi[1] * N);

//    return k[0]*(u.at(j, Xi[0])-z[0]) + k[1]*(u.at(j, Xi[1])-z[1]);
//}

double Problem1X1::vr(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

struct CauchyProblemX1U : public CauchyProblem
{
    CauchyProblemX1U(double a, double hx, double alpha, double lambda0, double lambdal, double Te, unsigned int i, unsigned int N,
                    const DoubleVector &k, const DoubleVector &z, const DoubleVector &x) :
        a(a), hx(hx), alpha(alpha), lambda0(lambda0), lambdal(lambdal), Te(Te), i(i), N(N), k(k), z(z), x(x) {}

    virtual ~CauchyProblemX1U() {}

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
            res = a1*u[0] + b1*u[1] + alpha*Te;
            double v = 0.0;
            for (unsigned int s=0; s<2; s++)
            {
                unsigned int r = (unsigned int)floor(x[s]*N);
                //printf("%d %f %f\n", r, x[s], t);

                double h1 = fabs(x[s] - (r+1)*hx);
                double h2 = fabs(x[s] - (r+0)*hx);
                //printf("%f %f\n", h1, h2);

                v += k[s] * ((h1*u[r+0] + h2*u[r+1])/hx - z[s]);
            }
            res += c1 * v;
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
    const DoubleVector &x;
};

void Problem1X1::calculateU(DoubleMatrix &u UNUSED_PARAM, double ht UNUSED_PARAM, double hx UNUSED_PARAM, unsigned int M UNUSED_PARAM, unsigned int N UNUSED_PARAM, double alpha UNUSED_PARAM, double lambda0 UNUSED_PARAM, double lambdal UNUSED_PARAM, double a UNUSED_PARAM)
{
//    const DoubleVector &xi = *pxi;

//    std::vector<CauchyProblem*> cps(N+1);
//    for (unsigned int n=0; n<=N; n++)
//    {
//        CauchyProblem *cp = new CauchyProblemX1U(a,hx, alpha, lambda0, lambdal, Te, n, N, k, z, *pxi);
//        cp->x0 = 0.0;
//        cp->y0 = initial(n);
//        cps[n] = cp;
//    }
//    CauchyProblem::rungeKutta(cps, 0.0, ht, M, u);
//    u.transpose();

//    for (unsigned int n=0; n<=N; n++)
//    {
//        CauchyProblemX1U *cp = dynamic_cast<CauchyProblemX1U*>(cps[n]);
//        delete cp;
//    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

struct CauchyProblemX1P1 : public CauchyProblem
{
    CauchyProblemX1P1(double a, double hx, double ht, double alpha, double alpha2, double lambda0, double lambdal, double Te, unsigned int i, unsigned int N,
                     const DoubleVector &k, const DoubleVector &z, const std::vector<unsigned int> &Xi, const DoubleMatrix &u) :
        a(a), hx(hx), ht(ht), alpha(alpha), alpha2(alpha2), lambda0(lambda0), lambdal(lambdal), Te(Te), i(i), N(N), k(k), z(z), Xi(Xi), u(u)  {}

    virtual ~CauchyProblemX1P1() {}

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

        //unsigned int j = (unsigned int)(t/ht);
        //double dd = p2->vl(j) - p2->vs(j);

        //if (i==Xi[0]) { res += 2.0 * alpha2 * dd * k[0] * (1.0/hx); }
        //if (i==Xi[1]) { res += 2.0 * alpha2 * dd * k[1] * (1.0/hx); }

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
    Problem1X1 *p2;
};

void Problem1X1::calculateP(DoubleMatrix &psi, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a)
{
    const DoubleVector &xi = *pxi;
    std::vector<unsigned int> Xi(L);
    Xi[0] = (unsigned int) round(xi[0] * N);
    Xi[1] = (unsigned int) round(xi[1] * N);
    const DoubleMatrix &u = *pu;

    std::vector<CauchyProblem*> cps(N+1);
    for (unsigned int n=0; n<=N; n++)
    {
        CauchyProblemX1P1 *cp = new CauchyProblemX1P1(a,hx, ht, alpha, alpha2, lambda0, lambdal, Te, n, N, k, z, Xi, u);
        cp->x0 = 1.0;
        cp->y0 = -2.0*alpha1*mu(n)*(u.at(M,n)-V[n]);
        cp->p2 = this;
        cps[n] = cp;
    }
    CauchyProblem::rungeKutta(cps, 1.0, -ht, M, psi);
    psi.transpose();
    for (unsigned int n=0; n<=N; n++)
    {
        CauchyProblemX1P1 *cp = dynamic_cast<CauchyProblemX1P1*>(cps[n]);
        delete cp;
    }
}

double Problem1X1::mu(unsigned int i UNUSED_PARAM) const
{
    return 1.0;
}

void Problem1X1::calculateV(const DoubleVector &xi)
{
    pxi = &xi;
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, alpha, lambda0, lambdal, a);
    V = u.row(M);
}

//double Problem1X1::vs(double j) const
//{
//    const DoubleMatrix &u = *pu;
//    std::vector<unsigned int> Xi(L);
//    Xi[0] = (unsigned int) round(xis[0] * N);
//    Xi[1] = (unsigned int) round(xis[1] * N);

//    return k[0]*(u(j, Xi[0])-z[0]) + k[1]*(u(j, Xi[1])-z[1]);
//}
