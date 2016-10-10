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
//    k << 2.50 << 2.70;
//    z << 1.52 << 1.71;
//    e << 0.40 << 0.70;

    DoubleVector x;
    x << 2.50 << 2.70; //k
    x << 1.52 << 1.71; //z
    x << 0.40 << 0.70; //x

    //    sx = &x;

    DoubleMatrix u;
    calculateU(u, x);
    IPrinter::printVector(u.row(M));

//    DoubleMatrix u1;
//    calculateU1(u1, x);
//    IPrinter::printVector(u1.row(M));

//    DoubleMatrix u2;
//    calculateU2(u2, x);
//    IPrinter::printVector(u2.row(M));

//    DoubleMatrix u3;
//    calculateU3(u3, x);
//    IPrinter::printVector(u3.row(M));


    //V = u.row(M);

    //fx(x);
}

double Example3::fx(const DoubleVector &x)
{
    px = &x;

    DoubleMatrix u;
    calculateU(u, x);

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
    DoubleVector z = x.mid(2, 3);
    DoubleVector e = x.mid(4, 5);

    return alpha0*sum + alpha1*k.EuclideanNorm() + alpha2*z.EuclideanNorm() + alpha3*e.EuclideanNorm();
}

void Example3::gradient(const DoubleVector &x, DoubleVector &g)
{}

void Example3::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const
{}

void Example3::project(DoubleVector &x, int index)
{}

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

void Example3::calculateU(DoubleMatrix &u, const DoubleVector &x)
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

double Example3::initial(unsigned int i UNUSED_PARAM) const
{
    return Ti;
}

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

