#include "problem1m.h"

void Problem1M::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem1M p;
}

Problem1M::Problem1M()
{
    L = 2;
    N = 1000;
    M = 1000;
    hx = 0.001;
    ht = 0.001;
    h  = 0.001;

    Ti = 2.0;
    Te = 3.0;
    alpha = 1.0;
    lambda0 = 1.0;
    lambdal = 1.0;

    alpha0 = 1.0;
    alpha1 = 1.0;
    alpha2 = 1.0;
    alpha3 = 1.0;
    a = 1.0;

    V.resize(N+1);

    /* init V */
    for (unsigned int n=0; n<=N; n++) V[n] = 4.0;

    DoubleVector x0;
    x0 << 3.50 << 3.70; //k
    x0 << 4.10 << 3.71; //z
    x0 << 0.25 << 0.75; //e
    px = &x0;

    DoubleVector ag(x0.size());
    gradient(x0, ag);
    ag.L2Normalize();

    DoubleVector ng(x0.size());
    IGradient::Gradient(this, h, x0, ng);
    ng.L2Normalize();

    printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", ag[0], ag[1], ag[2], ag[3], ag[4], ag[5]);
    printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", ng[0], ng[1], ng[2], ng[3], ng[4], ng[5]);

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

double Problem1M::fx(const DoubleVector &x)
{
    px = &x;
    DoubleVector k = x.mid(0,1);
    DoubleVector z = x.mid(2,3);
    DoubleVector e = x.mid(4,5);

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

    double norm1 = 0.0;
    double norm2 = 0.0;
    double norm3 = 0.0;

    norm1 = sqrt(k.at(0)*k.at(0) + k.at(1)*k.at(1));
    norm2 = sqrt(z.at(0)*z.at(0) + z.at(1)*z.at(1));
    norm3 = sqrt(e.at(0)*e.at(0) + e.at(1)*e.at(1));
    return alpha0*sum + alpha1*norm1 + + alpha2*norm2 + alpha3*norm3;
}

void Problem1M::gradient(const DoubleVector &x, DoubleVector &g)
{
    px = &x;
    DoubleVector k = x.mid(0,1);
    DoubleVector z = x.mid(2,3);
    DoubleVector e = x.mid(4,5);

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
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = p.at(m1, 0)*(u.at(m1, xi) - z[s]);
            double g2 = p.at(m2, 0)*(u.at(m2, xi) - z[s]);
            sum = sum + (g1 + g2);
        }
        g.at(i) = 0.5*ht*(-lambda0*a*a)*sum;
        g.at(i) += 2.0*alpha1*k.at(s);
        i++;
    }

    // z gradient
    for (unsigned int s=0; s<L; s++)
    {
        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = p.at(m1, 0);
            double g2 = p.at(m2, 0);
            sum = sum + (g1 + g2);
        }
        g.at(i) = 0.5*ht*(lambda0*a*a)*k[s]*sum;
        g.at(i) += 2.0*alpha2*z.at(s);
        i++;
    }

    // e gradient
    for (unsigned int s=0; s<L; s++)
    {
        unsigned int xi = (unsigned int)round(e.at(s) * N);
        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = p.at(m1, 0) * ((u.at(m1, xi+1) - u.at(m1, xi-1))/(2.0*hx));
            double g2 = p.at(m2, 0) * ((u.at(m2, xi+1) - u.at(m2, xi-1))/(2.0*hx));
            sum = sum + (g1 + g2);
        }
        g.at(i) = 0.5*ht*(-lambda0*a*a)*k.at(s)*sum;
        g.at(i) += 2.0*alpha3*e.at(s);
        i++;
    }
}

void Problem1M::calculateU(DoubleMatrix &u)
{
    const DoubleVector &x = *px;
    DoubleVector k = x.mid(0,1);
    DoubleVector z = x.mid(2,3);
    DoubleVector e = x.mid(4,5);

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);
    DoubleVector de(N+1);

    for (unsigned int n=0; n<=N; n++) u.at(0,n) = initial(n);

    for (unsigned int m=1; m<=M; m++)
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
        }

        qovmaFirstRowM(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size(), de.data());

        for (unsigned int i=0; i<=N; i++)
        {
            u.at(m, i) = rx[i];
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
    de.clear();
}

void Problem1M::calculateP(DoubleMatrix &p, const DoubleMatrix &u)
{
    const DoubleVector &x = *px;
    DoubleVector k = x.mid(0,1);
    DoubleVector z = x.mid(2,3);
    DoubleVector e = x.mid(4,5);

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

            qovmaFirstColM(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size(), de.data());

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

double Problem1M::initial(unsigned int n UNUSED_PARAM) const
{
    return Ti;
}

void Problem1M::print(unsigned int i UNUSED_PARAM, const DoubleVector &x UNUSED_PARAM, const DoubleVector &g UNUSED_PARAM, double alpha UNUSED_PARAM, RnFunction *fn UNUSED_PARAM) const
{
    IPrinter::printSeperatorLine();
    double f = fn->fx(x);
    printf("%14.10f\n", f);
    printf("x: %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", x[0], x[1], x[2], x[3], x[4], x[5]);

    DoubleVector ag(x.size());
    const_cast<Problem1M*>(this)->gradient(x, ag);
    ag.L2Normalize();

    DoubleVector ng(x.size());
    IGradient::Gradient(const_cast<Problem1M*>(this), h, x, ng);
    ng.L2Normalize();

    printf("a: %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", ag[0], ag[1], ag[2], ag[3], ag[4], ag[5], ag.L2Norm());
    printf("n: %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", ng[0], ng[1], ng[2], ng[3], ng[4], ng[5], ng.L2Norm());

    //const_cast<DoubleVector*>(px) = &x;
    DoubleMatrix u;
    const_cast<Problem1M*>(this)->calculateU(u);
    IPrinter::printVector(14,10,u.row(u.rows()-1));
}

void Problem1M::print(const DoubleVector &x, const DoubleVector &g UNUSED_PARAM, unsigned int i) const
{
}

void Problem1M::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
{
    if (x.at(4) < 0.10) x.at(4) = 0.10;
    if (x.at(4) > 0.90) x.at(4) = 0.90;

    if (x.at(5) < 0.10) x.at(5) = 0.10;
    if (x.at(5) > 0.90) x.at(5) = 0.90;
}

void qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e)
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

void qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e)
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
