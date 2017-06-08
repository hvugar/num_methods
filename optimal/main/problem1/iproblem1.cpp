#include "iproblem1.h"

double IProblem1::fx(const DoubleVector &y) const
{
    IProblem1* pm = const_cast<IProblem1*>(this);
    pm->py = &y;

    DoubleVector k,z,e;
    getParameters(k,z,e,y);

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
            SUM += alpha0*integral(u);
            SUM += R*penalty(k, z, e, u);
        }
    }
    SUM *= ((1.0/N1)*(1.0/N2));
    SUM += norm(k,z,e);
    return SUM;
}

double IProblem1::integral(const DoubleMatrix &u) const
{
    double sum = 0.0;
    sum += 0.5*mu(0)*(u[M][0]-V[0])*(u[M][0]-V[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += mu(n)*(u[M][n]-V[n])*(u[M][n]-V[n]);
    }
    sum += 0.5*mu(N)*(u[M][N]-V[N])*(u[M][N]-V[N]);
    sum *= hx;
    return sum;
}

double IProblem1::norm(const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const
{
    double norm1 = 0.0;
    double norm2 = 0.0;
    double norm3 = 0.0;

    for (unsigned int s=0; s<L; s++)
    {
        norm1 += (k[s] - k0[s])*(k[s] - k0[s]);
        norm2 += (z[s] - z0[s])*(z[s] - z0[s]);
        norm3 += (e[s] - e0[s])*(e[s] - e0[s]);
    }
    return alpha1*norm1 + alpha2*norm2 + alpha3*norm3;
}

double IProblem1::penalty(const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const
{
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

    return pnlt;
}

void IProblem1::gradient(const DoubleVector &y, DoubleVector &g)
{
    py = &y;

    DoubleVector k,z,e;
    getParameters(k,z,e,y);

    for (unsigned int i=0; i<g.size(); i++) g[i] = 0.0;

    unsigned int N1 = vfi.size();
    unsigned int N2 = vtt.size();

    for (unsigned int n1=0; n1<N1; n1++)
    {
        for (unsigned int n2=0; n2<N2; n2++)
        {
            fi = vfi[n1];
            tt = vtt[n2];

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
                    // Integral part of gradient
                    double sum = 0.0;
                    sum += 0.5*p.at(0, 0)*(u_xi(0,e[s],u) - z[s]);
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        sum += p[m][0]*(u_xi(m,e[s],u) - z[s]);
                    }
                    sum += 0.5*p[M][0]*(u_xi(M,e[s],u) - z[s]);
                    sum *= ht;

                    // Penalty part of gradient
                    double pnlt = 0.0;
                    pnlt += 0.5 * (u_xi(0,e[s],u)-z[s]) * sgn_min(0, k, z, e, u);
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        pnlt += (u_xi(m,e[s],u)-z[s]) * sgn_min(m, k, z, e, u);
                    }
                    pnlt += 0.5 * (u_xi(M,e[s],u)-z[s]) * sgn_min(M, k, z, e, u);
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
                    // Integral part of gradient
                    double sum = 0.0;
                    sum += 0.5 * p[0][0] * u_xi_d(0,e[s],u);
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        sum += p[m][0] * u_xi_d(m,e[s],u);
                    }
                    sum += 0.5 * p[M][0] * u_xi_d(M,e[s],u);
                    sum *= ht;

                    // Penalty part of gradient
                    double pnlt = 0.0;
                    pnlt += 0.5 * u_xi_d(0,e[s],u) * sgn_min(0, k, z, e, u);
                    for (unsigned int m=1; m<=M-1; m++)
                    {
                        pnlt += u_xi_d(m,e[s],u) * sgn_min(m, k, z, e, u);
                    }
                    pnlt += 0.5 * u_xi_d(M,e[s],u) * sgn_min(M, k, z, e, u);
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

void IProblem1::calculateU(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const
{
    calculateU2(u, k, z, e);
}

void IProblem1::calculateU1(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const
{
    u.clear();
    u.resize(M+1, N+1);

    double *da = (double*) malloc(sizeof(double)*(N+1));
    double *db = (double*) malloc(sizeof(double)*(N+1));
    double *dc = (double*) malloc(sizeof(double)*(N+1));
    double *dd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *de = (double*) malloc(sizeof(double)*(N+1));

    double a_a_ht_hx_hx = (a*a*ht)/(hx*hx);
    double lambda1_a_a_ht_hx = (lambda1*a*a*ht)/hx;
    double lambda2_a_a_ht_hx = (lambda2*a*a*ht)/hx;
    double lambda2_a_a_ht_tt_hx = (lambda2*a*a*ht*tt)/hx;
    double lambda0_ht = lambda0*ht;
    double lambda0_ht_tt = lambda0*ht*tt;

    for (unsigned int n=0; n<=N; n++) u[0][n] = initial(n);

    for (unsigned int m=1; m<=M; m++)
    {
        double kzs = 0.0;
        for (unsigned int s=0; s<L; s++) kzs += k[s]*z[s];

        // n = 0
        da[0] = 0.0;
        db[0] = 1.0 + a_a_ht_hx_hx + lambda1_a_a_ht_hx + lambda0_ht;
        dc[0] = -a_a_ht_hx_hx;
        dd[0] = u.at(m-1,0) + lambda0_ht_tt - lambda1_a_a_ht_hx*kzs;

        //        da[0] = 0.0;
        //        db[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambda1*a*a*ht)/hx + lambda0*ht;
        //        dc[0] = -(a*a*ht)/(hx*hx);
        //        dd[0] = u.at(m-1,0) + lambda0_ht_tt - lambda1_a_a_ht_hx*kzs;

        // n = 1,...,N-1
        for (unsigned int n=1; n<=N-1; n++)
        {
            da[n] = -a_a_ht_hx_hx;
            db[n] = 1.0 + (2.0*a_a_ht_hx_hx) + lambda0_ht;
            dc[n] = -a_a_ht_hx_hx;
            dd[n] = u.at(m-1,n) + lambda0_ht_tt;

            //            da[n] = -(a*a*ht)/(hx*hx);
            //            db[n] = 1.0 + 2.0*(a*a*ht)/(hx*hx) + lambda0*ht;
            //            dc[n] = -(a*a*ht)/(hx*hx);
            //            dd[n] = u.at(m-1,n) + lambda0*ht*tt;
        }

        // n = N
        da[N] = -a_a_ht_hx_hx;
        db[N] = 1.0 + a_a_ht_hx_hx + lambda2_a_a_ht_hx + lambda0_ht;
        dc[N] = 0.0;
        dd[N] = u.at(m-1,N)  + lambda0_ht_tt + lambda2_a_a_ht_tt_hx;

        //        da[N] = -(a*a*ht)/(hx*hx);
        //        db[N] = 1.0 + (a*a*ht)/(hx*hx) + (lambda2*a*a*ht)/hx + lambda0*ht;
        //        dc[N] = 0.0;
        //        dd[N] = u.at(m-1,N)  + lambda0*ht*tt + (lambda2*a*a*ht*tt)/hx;

        de[0]  =de[1] = 0.0;
        de[N-1]=de[N] = 0.0;
        for (unsigned int n=2; n<=N-2; n++)
        {
            de[n] = 0.0;
            for (unsigned int s=0; s<L; s++)
            {
                double diff = fabs(n*hx - e[s]);
                if (diff <= hx)
                {
                    de[n] += k[s] * (1.0 - diff/hx);
                }
            }

            de[n] *= -lambda1_a_a_ht_hx;

            //if (de[n]!=0.0) printf("%.10f %d\n", de[n], n);
        }

        for (unsigned int i=0; i<=N; i++) rx[i] = 0.0;

        qovmaFirstRowM(da, db, dc, dd, rx, N+1, de);

        for (unsigned int i=0; i<=N; i++)
        {
            u[m][i] = rx[i];
        }
        //        IPrinter::printVector(u.row(m));

        //        if (withError)
        //        {
        //            for (unsigned int s=0; s<L; s++)
        //            {
        //                unsigned int E = (unsigned int) round(e[s]);
        //                double w = (rand() % 2000)*0.001 - 1.0;
        //                u[m][E] += sign(w)*persent * u[m][E];
        //            }
        //        }
    }

    free(de);
    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

void IProblem1::calculateU2(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const
{
    u.clear();
    u.resize(M+1, N+1);

    double *da = (double*) malloc(sizeof(double)*(N+1));
    double *db = (double*) malloc(sizeof(double)*(N+1));
    double *dc = (double*) malloc(sizeof(double)*(N+1));
    double *dd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *de = (double*) malloc(sizeof(double)*(N+1));

    double a_a_ht_hx_hx = (a*a*ht)/(hx*hx);
    double lambda1_a_a_ht_hx = (lambda1*a*a*ht)/hx;
    double lambda2_a_a_ht_hx = (lambda2*a*a*ht)/hx;
    double lambda2_a_a_ht_tt_hx = (lambda2*a*a*ht*tt)/hx;
    double lambda0_ht = lambda0*ht;
    double lambda0_ht_tt = lambda0*ht*tt;

    for (unsigned int n=0; n<=N; n++) u[0][n] = initial(n);

    for (unsigned int m=1; m<=M; m++)
    {
        da[0] = 0.0;
        db[0] = 0.0;
        dc[0] = 0.0;

        double kzs = 0.0;
        for (unsigned int s=0; s<L; s++) kzs += k[s]*z[s];

        // n = 0
        de[0] = 1.0 + a_a_ht_hx_hx + lambda1_a_a_ht_hx + lambda0_ht;
        de[1] = -a_a_ht_hx_hx;
        for (unsigned int n=2; n<=N; n++)
        {
            de[n] = 0.0;
            for (unsigned int s=0; s<L; s++)
            {
                double diff = fabs(n*hx - e[s]);
                if (diff <= hx)
                {
                    de[n] += k[s] * (1.0 - diff/hx);
                }
            }
            de[n] *= -lambda1_a_a_ht_hx;
        }

        dd[0] = u.at(m-1,0) + lambda0_ht_tt - lambda1_a_a_ht_hx*kzs;

        // n = 1,...,N-1
        for (unsigned int n=1; n<=N-1; n++)
        {
            da[n] = -a_a_ht_hx_hx;
            db[n] = 1.0 + (2.0*a_a_ht_hx_hx) + lambda0_ht;
            dc[n] = -a_a_ht_hx_hx;
            dd[n] = u.at(m-1,n) + lambda0_ht_tt;
        }

        // n = N
        da[N] = -a_a_ht_hx_hx;
        db[N] = 1.0 + a_a_ht_hx_hx + lambda2_a_a_ht_hx + lambda0_ht;
        dc[N] = 0.0;
        dd[N] = u.at(m-1,N)  + lambda0_ht_tt + lambda2_a_a_ht_tt_hx;

        qovmaFirstRowM2(de, NULL, 0, da, db, dc, dd, rx, N+1);

        for (unsigned int i=0; i<=N; i++)
        {
            u[m][i] = rx[i];
        }

        if (withError)
        {
            for (unsigned int s=0; s<L; s++)
            {
                unsigned int E = (unsigned int) ceil(e[s] * N * DD);
                double w = (rand() % 2000)*0.001 - 1.0;
                u[m][E] += sign(w)*persent * u[m][E];
                //printf("%d %d %f %f %d\n", m, s, w, e[s], E);
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

void IProblem1::qovmaFirstRowM2(double *e, unsigned int *I UNUSED_PARAM, unsigned int L1 UNUSED_PARAM, double *a, double *b, double *c, double *d, double *x, unsigned int N) const
{
    //if (fabs(e[0]) <= DBL_EPSILON) throw std::exception("first item is zero.");

    unsigned int CNT=0;
    for (unsigned int i=2; i<N; i++)
    {
        if (fabs(e[i]) > DBL_EPSILON) CNT++;
    }

    double *q = (double*) malloc(sizeof(double) * N);
    double *p = (double*) malloc(sizeof(double) * N);
    double **k = (double**) malloc(sizeof(double*) * CNT);
    for (unsigned int i=0; i<CNT; i++)
    {
        k[i] = (double*) malloc(sizeof(double) * N);
    }

    unsigned int *E = (unsigned int *)malloc(sizeof(unsigned int)*CNT);
    unsigned int s = 0;
    for (unsigned int i=2; i<N; i++)
    {
        if (fabs(e[i]) > DBL_EPSILON)
        {
            E[s] = i;
            s++;
        }
    }

    //printf("%d %d %f\n", CNT, E[0], e[E[0]]);

    for (unsigned int i=0; i<N; i++)
    {
        if (i == 0)
        {
            q[0] = -e[1]/e[0];
            for (unsigned int s=0; s<CNT; s++)
            {
                k[s][0] = -e[E[s]]/e[0];
            }
            p[0] = d[0]/e[0];
            e[0] = 1.0;
        }
        else if (i == N-1)
        {
            p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
            q[i] = 0.0;
            for (unsigned int s=0; s<CNT; s++)
            {
                k[s][i] = 0.0;
            }
        }
        else
        {
            double m = b[i]+a[i]*q[i-1];
            q[i] = -c[i]/m;
            p[i] = (d[i] - a[i]*p[i-1])/m;
            for (unsigned int s=0; s<CNT; s++)
            {
                if (i<(E[s]-1))
                {
                    k[s][i] = -(a[i]*k[s][i-1])/m;
                }
                else
                {
                    k[s][i] = 0.0;
                }
                if (i==E[s]-1)
                {
                    q[i] += -(a[i]*k[s][i-1])/m;
                }
            }
        }
    }

    for (unsigned int i=N-1; i != UINT_MAX; i--)
    {
        if (i==(N-1))
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];
            for (unsigned int s=0; s<CNT; s++)
            {
                if (i<=E[s]-1)
                {
                    x[i] = x[i] + k[s][i]*x[E[s]];
                }
            }
        }
    }

    for (unsigned int s=0; s<CNT; s++) free(k[s]);
    free(k);
    free(E);
    free(q);
    free(p);
}

void IProblem1::qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
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
        }
    }

    //    printf("%d %d %.14f\n", L, E[0], e[E[0]]);

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

    //    IPrinter::printVector(x, n);

    for (unsigned int s=0; s<L; s++) free(k[s]);
    free(k);
    free(E);
    free(q);
    free(p);
}

void IProblem1::calculateP(DoubleMatrix &p, const DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const
{
    p.clear();
    p.resize(M+1, N+1);

    double *da = (double*) malloc(sizeof(double)*(N+1));
    double *db = (double*) malloc(sizeof(double)*(N+1));
    double *dc = (double*) malloc(sizeof(double)*(N+1));
    double *dd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *de = (double*) malloc(sizeof(double)*(N+1));

    double a_a_ht__hx_hx = (a*a*ht)/(hx*hx);
    double lambda1_a_a_ht__hx = (lambda1*a*a*ht)/hx;
    double lambda2_a_a_ht__hx = (lambda2*a*a*ht)/hx;
    double lambda0_ht = lambda0*ht;
    double lambda1_a_a_ht = lambda1*a*a*ht;

    for (unsigned int n=0; n<=N; n++)
    {
        p[M][n] = -2.0*alpha0*mu(n)*(u[M][n] - V[n]);
    }

    for (unsigned int m=M-1; m != UINT32_MAX; m--)
    {
        // n = 0
        //        da[0] = 0.0;
        //        db[0] = -1.0 - (a*a*ht)/(hx*hx) - (lambda1*a*a*ht)/hx - lambda0*ht;
        //        dc[0] = (a*a*ht)/(hx*hx);
        //        dd[0] = -p[m+1][0];

        da[0] = 0.0;
        db[0] = -1.0 - a_a_ht__hx_hx - lambda1_a_a_ht__hx - lambda0_ht;
        dc[0] = a_a_ht__hx_hx;
        dd[0] = -p[m+1][0];

        // n = 1,...,N-1
        for (unsigned int n=1; n<=N-1; n++)
        {
            //            da[n] = (a*a*ht)/(hx*hx);
            //            db[n] = -1.0-(2.0*a*a*ht)/(hx*hx) - lambda0*ht;
            //            dc[n] = (a*a*ht)/(hx*hx);
            //            dd[n] = -p[m+1][n] + R * ht * sgn_min(m, k, z, e, u);

            da[n] = a_a_ht__hx_hx;
            db[n] = -1.0-2.0*a_a_ht__hx_hx - lambda0_ht;
            dc[n] = a_a_ht__hx_hx;
            dd[n] = -p[m+1][n] + R * ht * sgn_min(m, k, z, e, u);

            for (unsigned int s=0; s<L; s++)
            {
                double _delta_ = delta(n,e,s);
                dd[n] +=  R * ht * 2.0*sgn_min(m, k, z, e, u) * k[s] * _delta_;
            }
        }

        // n = N
        //        da[N] = (a*a*ht)/(hx*hx);
        //        db[N] = -1.0-(a*a*ht)/(hx*hx) - (lambda2*a*a*ht)/hx - lambda0*ht;
        //        dc[N] = 0.0;
        //        dd[N] = -p[m+1][N];

        da[N] = a_a_ht__hx_hx;
        db[N] = -1.0-a_a_ht__hx_hx - lambda2_a_a_ht__hx - lambda0_ht;
        dc[N] = 0.0;
        dd[N] = -p[m+1][N];

        de[0] = de[1] = 0.0;
        de[N] = de[N-1] = 0.0;
        for (unsigned int n=2; n<=N-2; n++)
        {
            de[n] = 0.0;
            for (unsigned int s=0; s<L; s++)
            {
                double _delta_ = delta(n,e,s);
                de[n] += k[s] * _delta_;
            }
            de[n] *= lambda1_a_a_ht;
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

void IProblem1::print(unsigned int i, const DoubleVector &y, const DoubleVector &g, double r, GradientMethod::MethodResult) const
{
    IProblem1 *pm = const_cast<IProblem1*>(this);
    DoubleVector k,z,e;
    getParameters(k,z,e,y);

    DoubleMatrix u;
    pm->fi = vfi.at(0);
    pm->tt = vtt.at(0);
    pm->calculateU(u, k, z, e);

    IPrinter::printSeperatorLine();

    DoubleVector n(y.size());
    IGradient::Gradient(pm, 0.001, y, n);

    DoubleVector ak = g.mid(0*L,1*L-1);
    DoubleVector az = g.mid(1*L,2*L-1);
    DoubleVector ae = g.mid(2*L,3*L-1);

    DoubleVector nk = n.mid(0*L,1*L-1);
    DoubleVector nz = n.mid(1*L,2*L-1);
    DoubleVector ne = n.mid(2*L,3*L-1);

    printf("J[%d]: %.10f R: %.2f  I: %.10f P: %.10f\n", i, r, R, pm->integral(u), pm->penalty(k,z,e,u));
    printf("k: "); for (unsigned int s=0; s<L; s++) printf("%18.10f", k[s]); printf("\n");
    //printf("a: "); for (unsigned int s=0; s<L; s++) printf("%18.10f", ak[s]); printf(" | "); ak.L2Normalize(); for (unsigned int s=0; s<L; s++) printf("%18.10f", ak[s]); printf("\n");
    //printf("n: "); for (unsigned int s=0; s<L; s++) printf("%18.10f", nk[s]); printf(" | "); nk.L2Normalize(); for (unsigned int s=0; s<L; s++) printf("%18.10f", nk[s]); printf("\n");

    printf("z: "); for (unsigned int s=0; s<L; s++) printf("%18.10f", z[s]); printf("\n");
    //printf("a: "); for (unsigned int s=0; s<L; s++) printf("%18.10f", az[s]); printf(" | "); az.L2Normalize(); for (unsigned int s=0; s<L; s++) printf("%18.10f", az[s]); printf("\n");
    //printf("n: "); for (unsigned int s=0; s<L; s++) printf("%18.10f", nz[s]); printf(" | "); nz.L2Normalize(); for (unsigned int s=0; s<L; s++) printf("%18.10f", nz[s]); printf("\n");

    printf("e: "); for (unsigned int s=0; s<L; s++) printf("%18.10f", e[s]); printf("\n");
    //printf("a: "); for (unsigned int s=0; s<L; s++) printf("%18.10f", ae[s]); printf(" | "); ae.L2Normalize(); for (unsigned int s=0; s<L; s++) printf("%18.10f", ae[s]); printf("\n");
    //printf("n: "); for (unsigned int s=0; s<L; s++) printf("%18.10f", ne[s]); printf(" | "); ne.L2Normalize(); for (unsigned int s=0; s<L; s++) printf("%18.10f", ne[s]); printf("\n");

    DoubleVector v(M+1); for (unsigned int m=0; m<=M; m++) v[m] = vf(m,k,z,e,u);

    IPrinter::printVector(14,10,V,"V: ");
    IPrinter::printVector(14,10,u.row(u.rows()-1),"U: ");
    IPrinter::printVector(14,10,v,"v: ");

    v.clear();
    u.clear();

    n.clear();
}

double IProblem1::initial(unsigned int n UNUSED_PARAM) const
{
    return fi;
}

double IProblem1::mu(unsigned int n UNUSED_PARAM) const
{
    return 1.0;
}

double IProblem1::delta(unsigned int n, const DoubleVector &e, unsigned int s) const
{
    double sigma = 3.0*hx;
    double x = n*hx;
    return 1.0/(sqrt(2.0*M_PI)*sigma) * exp(-((x-e[s])*(x-e[s]))/(2.0*sigma*sigma));

    //    double diff = fabs(n*hx - e[s]);
    //    if (diff <= hx)
    //    {
    //        return (1.0 - diff/hx)/hx;
    //    }
    //    return 0.0;
}

double IProblem1::vf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const
{
    double v = 0.0;
    for (unsigned int s=0; s<L; s++)
    {
        v += k[s]*(u_xi(m,e[s],u)-z[s]);
    }
    return v;
}

double IProblem1::sgn_min(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const
{
    double aa = sign(vd0(m, k, z, e, u));
    double bb = fmin(0.0, gf(m, k, z, e, u));;
    return sign(vd0(m, k, z, e, u)) * fmin(0.0, gf(m, k, z, e, u));
}

double IProblem1::gf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const
{
    return d1 - fabs(vd0(m, k, z, e, u));
}

double IProblem1::vd0(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const
{
    return d0 - vf(m, k, z, e, u);
}

double IProblem1::sign(double x) const
{
    if (x < 0.0) return -1.0;
    if (x > 0.0) return +1.0;
    return 0.0;
}

double IProblem1::u_xi(unsigned int m, double e, const DoubleMatrix &u) const
{
//    unsigned int xi = (unsigned int)ceil(e * N * DD);
//    return (u[m][xi]*((xi+1)*hx - e) + u[m][xi+1]*(e - (xi)*hx)) / hx;

    //printf("*** 1 %f \n", e);
    unsigned int xi = (unsigned int)ceil(e * N * DD);
    double aaa = (u[m][xi]*((xi+1)*hx - e) + u[m][xi+1]*(e - (xi)*hx)) / hx;
    //puts("*** 2");
    return aaa;
}

double IProblem1::u_xi_d(unsigned int m, double e, const DoubleMatrix &u) const
{
    unsigned int xi = (unsigned int)round(e * N*DD);
    return (u.at(m, xi+1) - u.at(m, xi-1))/(2.0*hx);
}

void IProblem1::qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
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

void IProblem1::getParameters(DoubleVector &k, DoubleVector &z, DoubleVector &e, const DoubleVector &y) const
{
    unsigned int p = 0;

    if (optimizeK)
    {
        k = y.mid(p,p+(L-1));
        p+=L;
    }
    else
    {
        k = this->K;
    }

    if (optimizeZ)
    {
        z = y.mid(p,p+(L-1));
        p+=L;
    }
    else
    {
        z = this->Z;
    }

    if (optimizeE)
    {
        e = y.mid(p,p+(L-1));
        p+=L;
    }
    else
    {
        e = this->E;
    }
}
