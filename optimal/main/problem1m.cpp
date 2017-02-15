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
    x0 << +3.5000 << +3.7000; //k
    //x0 << +4.1000 << +3.9000; //z
    x0 << +4.0000 << +4.0000;
    x0 << +0.2000 << +0.8000; //e
    //k << 3.50 << 3.70;
    //z << 4.10 << 3.90;
    //e << 0.2 << 0.8;

    //x0 << 3.4004970612 << 3.6671254827 << 4.0801398905 << 3.8790071382 << 0.9299999944 << 0.9500000000;

//    px = &x0;

//    DoubleVector ag(x0.size());
//    gradient(x0, ag);
//    ag.L2Normalize();

//    DoubleVector ng(x0.size());
//    IGradient::Gradient(this, h, x0, ng);
//    ng.L2Normalize();

//    DoubleMatrix u;
//    calculateU(u);

//    printf("x: %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], integral(x0), norm(x0), fx(x0));
//    printf("a: %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", ag[0], ag[1], ag[2], ag[3], ag[4], ag[5]);
//    printf("n: %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", ng[0], ng[1], ng[2], ng[3], ng[4], ng[5]);
//    IPrinter::printVector(18, 10, u.row(u.rows()-1),"u: ");
    //printf("%18.10f %18.10f %18.10f %18.10f\n", ag[0], ag[1], ag[2], ag[3]);
    //printf("%18.10f %18.10f %18.10f %18.10f\n", ng[0], ng[1], ng[2], ng[3]);
    //printf("%18.10f %18.10f %18.10f\n",x0[0], x0[1], fx(x0));
    //printf("%18.10f %18.10f\n", ag[0], ag[1]);
    //printf("%18.10f %18.10f\n", ng[0], ng[1]);

    ConjugateGradient g;
    g.setFunction(this);
    g.setGradient(this);
    g.setPrinter(this);
    g.setProjection(this);
    g.setEpsilon1(0.00000001);
    g.setEpsilon2(0.00000001);
    g.setEpsilon3(0.00000001);
    g.setR1MinimizeEpsilon(0.01, 0.00000001);
    g.setNormalize(true);
    g.calculate(x0);
}

double Problem1M::fx(const DoubleVector &x) const
{
    return integral(x) + norm(x);
//    px = &x;
//    DoubleVector k = x.mid(0,1);
//    DoubleVector z = x.mid(2,3);
//    DoubleVector e = x.mid(4,5);

//    DoubleMatrix u;
//    calculateU(u);

//    double sum = 0.0;
//    sum += 0.5*mu(0)*(u.at(M, 0)-V.at(0))*(u.at(M, 0)-V.at(0));
//    for (unsigned int n=1; n<=N-1; n++)
//    {
//        sum += mu(n)*(u.at(M, n)-V.at(n))*(u.at(M, n)-V.at(n));
//    }
//    sum += 0.5*mu(N)*(u.at(M, N)-V.at(N))*(u.at(M, N)-V.at(N));
//    sum = hx*sum;

//    double norm1 = 0.0;
//    double norm2 = 0.0;
//    double norm3 = 0.0;

//    norm1 = sqrt(k.at(0)*k.at(0) + k.at(1)*k.at(1));
//    norm2 = sqrt(z.at(0)*z.at(0) + z.at(1)*z.at(1));
//    norm3 = sqrt(e.at(0)*e.at(0) + e.at(1)*e.at(1));
//    return alpha0*sum + alpha1*norm1 + alpha2*norm2 + alpha3*norm3;
}

double Problem1M::integral(const DoubleVector &x) const
{
    const_cast<Problem1M*>(this)->px = &x;
    DoubleVector k = x.mid(0,1);
    DoubleVector z = x.mid(2,3);
    DoubleVector e = x.mid(4,5);

    DoubleMatrix u;
    calculateU(u);

    double sum = 0.0;
    sum += 0.5*mu(0)*(u.at(M, 0)-V.at(0))*(u.at(M, 0)-V.at(0));
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += mu(n)*(u.at(M, n)-V.at(n))*(u.at(M, n)-V.at(n));
    }
    sum += 0.5*mu(N)*(u.at(M, N)-V.at(N))*(u.at(M, N)-V.at(N));
    sum = hx*sum;

    return alpha0*sum;
}

double Problem1M::norm(const DoubleVector &x) const
{
    const_cast<Problem1M*>(this)->px = &x;
    DoubleVector k = x.mid(0,1);
    DoubleVector z = x.mid(2,3);
    DoubleVector e = x.mid(4,5);

    double norm1 = 0.0;
    double norm2 = 0.0;
    double norm3 = 0.0;

    norm1 = sqrt(k.at(0)*k.at(0) + k.at(1)*k.at(1));
    norm2 = sqrt(z.at(0)*z.at(0) + z.at(1)*z.at(1));
    norm3 = sqrt(e.at(0)*e.at(0) + e.at(1)*e.at(1));
    return alpha1*norm1 + alpha2*norm2 + alpha3*norm3;
}

void Problem1M::gradient(const DoubleVector &x, DoubleVector &g)
{
    px = &x;
    DoubleVector k = x.mid(0,1);
    DoubleVector z = x.mid(2,3);
    DoubleVector e = x.mid(4,5);

//    IGradient::Gradient(this,h,x,g);
//    return;

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

void Problem1M::calculateU(DoubleMatrix &u) const
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

            //if (fabs(i*hx - e.at(0)) <= hx) { u[m][i] *= 1.05; }
            //if (fabs(i*hx - e.at(1)) <= hx) { u[m][i] *= 1.05; }

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

void Problem1M::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double fx) const
{
    IPrinter::printSeperatorLine();
    Problem1M *pm = const_cast<Problem1M*>(this);
    printf("Iteration %d Integral: %18.10f Norm: %18.10f Funksional: %18.10f Funksional: %18.10f\n", i, pm->integral(x), pm->norm(x), pm->fx(x), fx);

    DoubleVector ag(x.size());
    pm->gradient(x, ag);
    ag.L2Normalize();

    DoubleVector ng(x.size());
    IGradient::Gradient(pm, h, x, ng);
    ng.L2Normalize();

    DoubleMatrix u;
    pm->px = &x;
    pm->calculateU(u);

    DoubleVector gg = g;
    gg.L2Normalize();

    printf("x: %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", x[0], x[1], x[2], x[3], x[4], x[5]);
    printf("g: %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", gg[0], gg[1], gg[2], gg[3], gg[4], gg[5], gg.L2Norm());
    printf("a: %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", ag[0], ag[1], ag[2], ag[3], ag[4], ag[5], ag.L2Norm());
    printf("n: %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", ng[0], ng[1], ng[2], ng[3], ng[4], ng[5], ng.L2Norm());
    IPrinter::printVector(18,10,u.row(u.rows()-1),"u: ");
}

void Problem1M::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
{
    if (x.at(4) < 0.05) x.at(4) = 0.05;
    if (x.at(4) > 0.95) x.at(4) = 0.95;

    if (x.at(5) < 0.05) x.at(5) = 0.05;
    if (x.at(5) > 0.95) x.at(5) = 0.95;
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
