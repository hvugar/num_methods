#include "problem1L3.h"

void Problem1L3::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem1L3 p;
}

Problem1L3::Problem1L3()
{
    L = 3;
    N = 1000;
    M = 1000;
    hx = 0.001;
    ht = 0.001;
    h  = 0.001;

    Ti = 2.0;
    Te = 3.0;
    alpha = 0.01;
    lambda0 = 1.0;
    lambdal = 1.0;

    alpha0 = 1.0;
    alpha1 = 0.01;
    alpha2 = 0.01;
    alpha3 = 0.01;
    a = 1.0;

    V.resize(N+1);

    /* init V */
    for (unsigned int n=0; n<=N; n++)
    {
        //double h1 = 0.4/N;
        V[n] = 4.0;//4.2 - n*h1;
    }
    IPrinter::printVector(14, 10, V,"V: ");

    DoubleVector x0;
    x0 << -3.5000 << -3.7000 << -3.4000; //k
    x0 << +4.1000 << +3.9000 << +3.7000; //z
    x0 << +0.2500 << +0.5000 << +0.7500; //e

    //k << 3.50 << 3.70;
    //z << 4.20 << 4.20;
    //e << 0.40 << 0.90;

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

double Problem1L3::fx(const DoubleVector &x) const
{
    return integral(x) + norm(x);
}

double Problem1L3::integral(const DoubleVector &x) const
{
    const_cast<Problem1L3*>(this)->px = &x;

    DoubleMatrix u;
    calculateU(u);

    double sum = 0.0;
    sum += 0.5*mu(0)*(u[M][0]-V[0])*(u[M][0]-V[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += mu(n)*(u[M][n]-V[n])*(u[M][n]-V[n]);
    }
    sum += 0.5*mu(N)*(u[M][N]-V[N])*(u[M][N]-V[N]);
    sum = hx*sum;

    return alpha0*sum;
}

double Problem1L3::norm(const DoubleVector &x) const
{
    const_cast<Problem1L3*>(this)->px = &x;
    DoubleVector k = x.mid(0,2);
    DoubleVector z = x.mid(3,5);
    DoubleVector e = x.mid(6,8);

    double norm1 = 0.0;
    double norm2 = 0.0;
    double norm3 = 0.0;

    norm1 = k.at(0)*k.at(0) + k.at(1)*k.at(1) + k.at(2)*k.at(2);
    norm2 = z.at(0)*z.at(0) + z.at(1)*z.at(1) + z.at(2)*z.at(2);
    norm3 = e.at(0)*e.at(0) + e.at(1)*e.at(1) + e.at(2)*e.at(2);
    return alpha1*norm1 + alpha2*norm2 + alpha3*norm3;
}

void Problem1L3::gradient(const DoubleVector &x, DoubleVector &g)
{
    px = &x;
    DoubleVector k = x.mid(0,2);
    DoubleVector z = x.mid(3,5);
    DoubleVector e = x.mid(6,8);

    DoubleMatrix u;
    calculateU(u);

    DoubleMatrix p;
    calculateP(p, u);

    unsigned int i = 0;
    // k gradient
    for (unsigned int s=0; s<L; s++)
    {
        unsigned int xi = (unsigned int)round(e.at(s) * N);
        double sum = 0.0;
        sum += 0.5*p.at(0, 0)*(u.at(0, xi) - z[s]);
        for (unsigned int m=1; m<=M-1; m++)
        {
            sum += p.at(m, 0)*(u.at(m, xi) - z[s]);
        }
        sum += 0.5*p.at(M, 0)*(u.at(M, xi) - z[s]);

        g.at(i) = -lambda0*a*a*ht*sum + 2.0*alpha1*k.at(s);
        i++;
    }

    // z gradient
    for (unsigned int s=0; s<L; s++)
    {
        double sum = 0.0;
        sum += 0.5*p.at(0, 0);
        for (unsigned int m=1; m<=M-1; m++)
        {
            sum += p.at(m, 0);
        }
        sum += 0.5*p.at(M, 0);
        g.at(i) = lambda0*a*a*ht*k[s]*sum + 2.0*alpha2*z.at(s);
        i++;
    }

    // e gradient
    for (unsigned int s=0; s<L; s++)
    {
        unsigned int xi = (unsigned int)round(e.at(s) * N);
        double sum = 0.0;
        sum += 0.5 * p.at(0, 0) * ((u.at(0, xi+1) - u.at(0, xi-1))/(2.0*hx));
        for (unsigned int m=1; m<=M-1; m++)
        {
            sum += p.at(m, 0) * ((u.at(m, xi+1) - u.at(m, xi-1))/(2.0*hx));
        }
        sum += 0.5 * p.at(M, 0) * ((u.at(M, xi+1) - u.at(M, xi-1))/(2.0*hx));
        g.at(i) = -lambda0*a*a*ht*k[s]*sum + 2.0*alpha3*e.at(s);
        i++;
    }
}

void Problem1L3::calculateU(DoubleMatrix &u) const
{
    const DoubleVector &x = *px;
    DoubleVector k = x.mid(0,2);
    DoubleVector z = x.mid(3,5);
    DoubleVector e = x.mid(6,8);

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);
    DoubleVector de(N+1);

    for (unsigned int n=0; n<=N; n++) u[0][n] = initial(n);

    for (unsigned int m=1; m<=M; m++)
    {
        // n = 0
        da[0] = 0.0;
        db[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambda0*a*a*ht)/hx + alpha*ht;
        dc[0] = -(a*a*ht)/(hx*hx);
        dd[0] = u.at(m-1,0) + alpha*ht*Te - ((lambda0*a*a*ht)/hx)*(k[0]*z[0] + k[1]*z[1] + k[2]*z[2]);

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

        for (unsigned int n=2; n<=N-2; n++)
        {
            de[n] = 0.0;

            if (fabs(n*hx - e[0]) <= hx)
            {
                de[n] = -((lambda0*(a*a*ht)/hx) * k[0]) * (1.0 - fabs(n*hx - e[0])/hx);
            }
            if (fabs(n*hx - e[1]) <= hx)
            {
                de[n] = -((lambda0*(a*a*ht)/hx) * k[1]) * (1.0 - fabs(n*hx - e[1])/hx);
            }
            if (fabs(n*hx - e[2]) <= hx)
            {
                de[n] = -((lambda0*(a*a*ht)/hx) * k[2]) * (1.0 - fabs(n*hx - e[2])/hx);
            }
        }

        qovmaFirstRowM(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size(), de.data());

        for (unsigned int i=0; i<=N; i++)
        {
            u.at(m, i) = rx[i];
            //if (fabs(i*hx - e.at(0)) <= hx) { u[m][i] += u[m][i]*0.01; }
            //if (fabs(i*hx - e.at(1)) <= hx) { u[m][i] -= u[m][i]*0.01; }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
    de.clear();
}

void Problem1L3::calculateP(DoubleMatrix &p, const DoubleMatrix &u)
{
    const DoubleVector &x = *px;
    DoubleVector k = x.mid(0,2);
    DoubleVector z = x.mid(3,5);
    DoubleVector e = x.mid(6,8);

    p.clear();
    p.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);
    DoubleVector de(N+1);

    for (unsigned int n=0; n<=N; n++) p.at(M,n) = -2.0*alpha0*mu(n)*(u.at(M, n) - V.at(n));

    for (unsigned int m=M-1; m != UINT32_MAX; m--)
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
        }

        // n = N
        da[N] = (a*a*ht)/(hx*hx);
        db[N] = -1.0-(a*a*ht)/(hx*hx) - (lambdal*a*a*ht)/hx - alpha*ht;
        dc[N] = 0.0;
        dd[N] = -p.at(m+1,N);

        for (unsigned int n=2; n<=N-2; n++)
        {
            de[n] = 0.0;
            if (fabs(n*hx - e[0]) <= hx)
            {
                de[n] = k[0]*(lambda0*(a*a*ht)/hx) * (1.0 - fabs(n*hx - e[0])/hx);
            }
            if (fabs(n*hx - e[1]) <= hx)
            {
                de[n] = k[1]*(lambda0*(a*a*ht)/hx) * (1.0 - fabs(n*hx - e[1])/hx);
            }
            if (fabs(n*hx - e[2]) <= hx)
            {
                de[n] = k[2]*(lambda0*(a*a*ht)/hx) * (1.0 - fabs(n*hx - e[2])/hx);
            }
        }

        qovmaFirstColM(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size(), de.data());

        for (unsigned int i=0; i<=N; i++) p.at(m, i) = rx[i];
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
    de.clear();
}

double Problem1L3::initial(unsigned int n UNUSED_PARAM) const
{
    return Ti;
}

void Problem1L3::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double r, GradientMethod::MethodResult) const
{
    Problem1L3 *pm = const_cast<Problem1L3*>(this);
    DoubleVector k = x.mid(0,2);
    DoubleVector z = x.mid(3,5);
    DoubleVector e = x.mid(6,8);

    IPrinter::printSeperatorLine();

    DoubleVector ng(x.size());
    IGradient::Gradient(pm, h, x, ng);

    DoubleMatrix u;
    pm->px = &x;
    pm->calculateU(u);

    unsigned int e1 = (unsigned int)round(e[0]*N);
    unsigned int e2 = (unsigned int)round(e[1]*N);
    unsigned int e3 = (unsigned int)round(e[2]*N);
    double v = k[0]*(u[M][e1]-z[0]) + k[1]*(u[M][e2]-z[1]) + k[2]*(u[M][e3]-z[2]);

    DoubleVector ag = g;

    DoubleVector nag = ag;
    DoubleVector nng = ng;

    nag.L2Normalize();
    nng.L2Normalize();

    printf("J[%d]: %.10f v: %.10f\n", i, r, v);
    printf("k: %14.10f %14.10f %14.10f\n", k[0], k[1], k[2]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ag[0], ag[1], ag[2], nag[0], nag[1], nag[2]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ng[0], ng[1], ng[2], nng[0], nng[1], nng[2]);
    printf("z: %14.10f %14.10f %14.10f\n", z[0], z[1], z[2]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ag[3], ag[4], ag[5], nag[3], nag[4], nag[5]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ng[3], ng[4], ng[5], nng[3], nng[4], nng[5]);
    printf("e: %14.10f %14.10f %14.10f\n", e[0], e[1], e[2]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ag[6], ag[7], ag[8], nag[6], nag[7], nag[8]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ng[6], ng[7], ng[8], nng[6], nng[7], nng[8]);
    IPrinter::printVector(14,10,u.row(u.rows()-1),"u: ");
}

void Problem1L3::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
{
    /* z lower/upper limits */
    if (x.at(3) < 3.80) x.at(3) = 3.80;
    if (x.at(3) > 4.20) x.at(3) = 4.20;

    if (x.at(4) < 3.80) x.at(4) = 3.80;
    if (x.at(4) > 4.20) x.at(4) = 4.20;

    if (x.at(5) < 3.80) x.at(5) = 3.80;
    if (x.at(5) > 4.20) x.at(5) = 4.20;

    /* e lower/upper limits */
    if (x.at(6) < 5*hx) x.at(6) = 5*hx;
    if (x.at(6) > 0.30) x.at(6) = 0.30;

    if (x.at(7) < 0.35) x.at(7) = 0.35;
    if (x.at(7) > 0.65) x.at(7) = 0.65;

    if (x.at(8) < 0.70)  x.at(8) = 0.70;
    if (x.at(8) > (N-5)*hx) x.at(8) = (N-5)*hx;
}

void Problem1L3::getComponents(DoubleVector &, DoubleVector &, DoubleVector &, const DoubleVector &) const
{
}

void Problem1L3::qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
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

void Problem1L3::qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);

    //printf("e: %f %f %f %f %f %f %f %f\n", e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]);
    unsigned int L = 0;
    for (unsigned int s=0; s<n; s++)
    {
        if (fabs(e[s]) > 0.000001)
        {
            //printf("s: %4d %.14f\n", s, e[s]);
            L+=1;
        }
    }
    //printf("L %d\n", L);
    unsigned int *E = (unsigned int *)malloc(sizeof(unsigned int)*L);

    unsigned int i = 0;
    for (unsigned int s=0; s<n; s++)
    {
        if (fabs(e[s]) > 0.000001)
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
