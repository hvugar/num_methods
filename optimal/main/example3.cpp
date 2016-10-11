#include "example3.h"
#include <cmethods.h>
#include <float.h>
#include <vector>

void Example3::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    Example3 e3;
}

Example3::Example3()
{
    z << 1.52 << 1.71;
    e << 0.40 << 0.70;

    DoubleVector x0;
    x0.resize(2);
    x0[0] = 2.50;
    x0[1] = 2.70; //k

//    x0 << 1.52 << 1.71; //z
//    x0 << 0.40 << 0.70; //e

    px = &x0;
    DoubleMatrix u;
    calculateU(u);
    IPrinter::printVector(u.row(M));
    V = u.row(M);

    DoubleVector x1;
    x1.resize(2);
    x1[0] = 3.50;
    x1[1] = 3.70; //k
//    x1 << 1.52 << 1.71; //z
//    x1 << 0.40 << 0.70; //x
    px = &x1;

    printf("Optimal:   %.10f %.10f %.10f %.10f %.10f %.10f\n", x0.at(0), x0.at(1), x0.at(2), x0.at(3), x0.at(4), x0.at(5));
    printf("Initial:   %.10f %.10f %.10f %.10f %.10f %.10f\n", x1.at(0), x1.at(1), x1.at(2), x1.at(3), x1.at(4), x1.at(5));

    double h = 0.001;
    DoubleVector g2(2, 0.0);
    IGradient::Gradient(this, h, x1, g2);
    //DoubleVector gn2 = g2;
    g2.L2Normalize();
    printf("Numerical: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", g2[0], g2[1], g2[2], g2[3], g2[4], g2[5]);

    DoubleVector g1(x0.size(),0.0);
    gradient(x1, g1);
    //DoubleVector gn1 = g1;
    g1.L2Normalize();
    printf("Analytic:  %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", g1[0], g1[1], g1[2], g1[3], g1[4], g1[5]);

    puts("------------------------------------------");

}

double Example3::fx(const DoubleVector &x)
{
    px = &x;

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

    DoubleVector k = x.mid(0, 1);
//    DoubleVector z = x.mid(2, 3);
//    DoubleVector e = x.mid(4, 5);

    return alpha0*sum + alpha1*x.EuclideanNorm() + alpha2*z.EuclideanNorm() + alpha3*e.EuclideanNorm();
}

void Example3::gradient(const DoubleVector &x, DoubleVector &g)
{
    px = &x;
//    DoubleVector z = x.mid(2, 3);
//    DoubleVector e = x.mid(4, 5);

    DoubleMatrix u;
    calculateU(u);

    DoubleMatrix psi;
    calculateP(psi, u);

    double sum = 0.0;

    // k gradient
    g[0] = 0.0;
    for (unsigned int m=0; m<=M-1; m++)
    {
        unsigned int m1 = m + 0;
        unsigned int m2 = m + 1;
        double g1 = -lambda0*a*a*psi.at(m1, 0)*(u.at(m1, 400) - z[0]);
        double g2 = -lambda0*a*a*psi.at(m2, 0)*(u.at(m2, 400) - z[0]);
        sum = sum + (g1 + g2);
    }
    sum = 0.5*ht*sum + 2.0*alpha1*x.at(0);
    g[0] += sum;

    sum = 0.0;
    g[1] = 0.0;
    for (unsigned int m=0; m<=M-1; m++)
    {
        unsigned int m1 = m + 0;
        unsigned int m2 = m + 1;
        double g1 = -lambda0*a*a*psi.at(m1, 0)*(u.at(m1, 700) - z[1]);
        double g2 = -lambda0*a*a*psi.at(m2, 0)*(u.at(m2, 700) - z[1]);
        sum = sum + (g1 + g2);
    }
    sum = 0.5*ht*sum + 2.0*alpha1*x.at(1);
    g[1] += sum;

//    for (unsigned int s = 0; s<L; s++)
//    {
//        g[s] = 0.0;
//        double sum = 0.0;
//        for (unsigned int m=0; m<=M-1; m++)
//        {
//            unsigned int m1 = m + 0;
//            unsigned int m2 = m + 1;
//            double g1 = -lambda0*a*a*psi.at(m1, 0)*(u.at(m1, Xi[s]) - z[s]);
//            double g2 = -lambda0*a*a*psi.at(m2, 0)*(u.at(m2, Xi[s]) - z[s]);
//            sum = sum + (g1 + g2);
//        }
//        sum = 0.5*ht*sum + 2.0*alpha1*k[s];
//        g[s] += sum;
//    }

    g[2] = g[3] = g[4] = g[5] = 0.0;
}

void Example3::print(unsigned int i UNUSED_PARAM, const DoubleVector &x UNUSED_PARAM, const DoubleVector &g UNUSED_PARAM, double alpha UNUSED_PARAM, RnFunction *fn UNUSED_PARAM) const
{}

void Example3::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
{}

void Example3::calculateU(DoubleMatrix &u)
{
    DoubleVector k = *px;
//    DoubleVector k = x.mid(0, 1);
//    DoubleVector z = x.mid(2, 3);
//    DoubleVector e = x.mid(4, 5);

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
            db[0] = 1.0+(a*a*ht)/(hx*hx) + (lambda0*a*a*ht)/hx + alpha*ht;
            dc[0] = -(a*a*ht)/(hx*hx);
            dd[0] = u.at(m-1,0) + alpha*ht*Te - ((lambda0*a*a*ht)/hx)*(k[0]*z[0] + k[1]*z[1]);

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

            for (unsigned int n=0; n<=N; n++)
            {
                de[n] = 0.0;
                if (fabs(n*hx - e.at(0)) <= DBL_EPSILON) de[n] = -k[0]*(lambda0*(a*a*ht)/hx);
                if (fabs(n*hx - e.at(1)) <= DBL_EPSILON) de[n] = -k[1]*(lambda0*(a*a*ht)/hx);
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

void Example3::calculateP(DoubleMatrix &p, const DoubleMatrix &u)
{
    DoubleVector x = *px;
    DoubleVector k = x.mid(0, 1);
//    DoubleVector z = x.mid(2, 3);
//    DoubleVector e = x.mid(4, 5);

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
            db[0] = -1.0-(a*a*ht)/(hx*hx) - (lambda0*a*a*ht)/hx - alpha*ht;
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
                if (fabs(n*hx - e.at(0)) <= DBL_EPSILON) de[n] = -k[0]*(lambda0*(a*a*ht)/hx);
                if (fabs(n*hx - e.at(1)) <= DBL_EPSILON) de[n] = -k[1]*(lambda0*(a*a*ht)/hx);
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

void Example3::qovmaFirstCol(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e)
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
            p[i] = +(d[i]-c[i]*q[i+1])/m;
            q[i] = -(a[i]+c[i]*k[i+1])/m;
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
            p[i] = +(d[i]-c[i]*q[i+1])/m;
            q[i] = -a[i]/m;
            k[i] = -(e[i]+c[i]*k[i+1])/m;
        }
    }

    for (unsigned int i=0; i<n; i++)
    {
        if (i==0)
        {
            x[i] = p[i];
        }
        else if (i==1)
        {
            x[i] = p[i] + q[i]*x[i-1];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i-1] + k[i]*x[0];
        }
    }

    free(k);
    free(p);
    free(q);
}

double Example3::initial(unsigned int i UNUSED_PARAM) const
{
    return Ti;
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
