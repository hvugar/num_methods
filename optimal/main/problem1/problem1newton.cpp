#include "problem1newton.h"
#include "problem1L2.h"

Problem1NewtonF::Problem1NewtonF() {}

Problem1NewtonF::~Problem1NewtonF() {}

void Problem1NewtonF::calculate(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e)
{
    unsigned int M = timeDimension().sizeN();
    unsigned int N = spaceDimension(Dimension::Dim1).sizeN();
    double ht = timeDimension().step();
    double hx = spaceDimension(Dimension::Dim1).step();

    u.clear();
    u.resize(M+1, N+1);

    double *da = (double*) malloc(sizeof(double)*(N+1));
    double *db = (double*) malloc(sizeof(double)*(N+1));
    double *dc = (double*) malloc(sizeof(double)*(N+1));
    double *dd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *de = (double*) malloc(sizeof(double)*(N+1));

    for (unsigned int n=0; n<=N; n++) u[0][n] = 0.0;//initial(n);

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
        }
    }

    free(de);
    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

double Problem1NewtonF::initial(const SpaceNode &) const
{
    return fi;
}

double Problem1NewtonF::boundary(const SpaceNode &, const TimeNode &, BoundaryType) const
{
    return 0.0;
}

double Problem1NewtonF::f(const SpaceNode &, const TimeNode &) const
{
    return 0.0;
}

void Problem1NewtonF::qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
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

////////////////////////////////////////////////////////////////////////////

Problem1NewtonB::Problem1NewtonB() {}

Problem1NewtonB::~Problem1NewtonB() {}

void Problem1NewtonB::calculate(DoubleMatrix &p, const DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e)
{
    C_UNUSED(p);
    C_UNUSED(u);
    C_UNUSED(k);
    C_UNUSED(z);
    C_UNUSED(e);
//    unsigned int M = timeDimension().sizeN();
//    unsigned int N = spaceDimension(Dimension::Dim1).sizeN();
//    double ht = timeDimension().step();
//    double hx = spaceDimension(Dimension::Dim1).step();

//    p.clear();
//    p.resize(M+1, N+1);

//    double *da = (double*) malloc(sizeof(double)*(N+1));
//    double *db = (double*) malloc(sizeof(double)*(N+1));
//    double *dc = (double*) malloc(sizeof(double)*(N+1));
//    double *dd = (double*) malloc(sizeof(double)*(N+1));
//    double *rx = (double*) malloc(sizeof(double)*(N+1));
//    double *de = (double*) malloc(sizeof(double)*(N+1));

//    for (unsigned int n=0; n<=N; n++)
//    {
//        p[M][n] = -2.0*alpha0*mu(n)*(u[M][n] - V[n]);
//    }

//    for (unsigned int m=M-1; m != UINT32_MAX; m--)
//    {
//        // n = 0
//        da[0] = 0.0;
//        db[0] = -1.0 - (a*a*ht)/(hx*hx) - (lambda1*a*a*ht)/hx - lambda0*ht;
//        dc[0] = (a*a*ht)/(hx*hx);
//        dd[0] = -p[m+1][0];

//        // n = 1,...,N-1
//        for (unsigned int n=1; n<=N-1; n++)
//        {
//            da[n] = (a*a*ht)/(hx*hx);
//            db[n] = -1.0-(2.0*a*a*ht)/(hx*hx) - lambda0*ht;
//            dc[n] = (a*a*ht)/(hx*hx);
//            dd[n] = -p[m+1][n];

//            double dif0 = fabs(n*hx - e[0]);
//            if (dif0 <= hx)
//            {
//                dd[n] +=  ht * R * sgn_max(m, k, z, e, u) * k[0] * (1.0/hx);// * (1.0 - dif0/hx);
//            }
//            double dif1 = fabs(n*hx - e[1]);
//            if (dif1 <= hx)
//            {
//                dd[n] += ht * R * sgn_max(m, k, z, e, u) * k[1] * (1.0/hx);// * (1.0 - dif1/hx);
//            }
//        }

//        // n = N
//        da[N] = (a*a*ht)/(hx*hx);
//        db[N] = -1.0-(a*a*ht)/(hx*hx) - (lambda2*a*a*ht)/hx - lambda0*ht;
//        dc[N] = 0.0;
//        dd[N] = -p[m+1][N];

//        de[0] = de[1] = 0.0;
//        de[N] = de[N-1] = 0.0;
//        for (unsigned int n=2; n<=N-2; n++)
//        {
//            de[n] = 0.0;
//            double dif0 = fabs(n*hx - e[0]);
//            if (dif0 <= hx)
//            {
//                de[n] = k[0]*(lambda1*(a*a*ht)/hx) * (1.0 - dif0/hx);
//            }

//            double dif1 = fabs(n*hx - e[1]);
//            if (dif1 <= hx)
//            {
//                de[n] = k[1]*(lambda1*(a*a*ht)/hx) * (1.0 - dif1/hx);
//            }
//        }

//        qovmaFirstColM(da, db, dc, dd, rx, N+1, de);

//        for (unsigned int i=0; i<=N; i++) p[m][i] = rx[i];
//    }

//    free(de);
//    free(rx);
//    free(dd);
//    free(dc);
//    free(db);
//    free(da);
}

double Problem1NewtonB::initial(const SpaceNode &) const
{
    return 0.0;
}

double Problem1NewtonB::boundary(const SpaceNode &, const TimeNode &, BoundaryType) const
{
    return 0.0;
}

double Problem1NewtonB::f(const SpaceNode &, const TimeNode &) const
{
    return 0.0;
}

void Problem1NewtonB::qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const
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
