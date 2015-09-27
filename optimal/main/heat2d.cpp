#include "heat2d.h"
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

void TomasAlgorithm(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x)
{
    if (x.size() != a.size() || x.size() != b.size() || x.size() != c.size() || x.size() != d.size())
        return;

    int n = x.size();
    DoubleVector p(n);
    DoubleVector q(n);

    for (int i=0; i<n; i++)
    {
        if (i==0)
        {
            p[0] = d[0]/b[0];
            q[0] = -c[0]/b[0];
        } else
            if(i==n-1)
            {
                p[n-1] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[n-1] = 0.0;
            }
            else
            {
                p[i] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[i] = -c[i]/(b[i]+a[i]*q[i-1]);
            }
    }

    for (int i=n-1; i>=0; i--)
    {
        if (i==n-1)
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];
        }
    }
}

void printLayer(const DoubleMatrix& x)
{
    int m = x.size()/10;
    for (int j=0; j<x.size(); j++)
    {
        if (j%m==0)
        {
            for (int i=0; i<x[j].size(); i++)
            {
                if (i%m==0) printf("%14.10f ", x[j][i]);
            }
            puts("");
        }
    }
}

Heat2DControl::Heat2DControl()
{
    N1 = 100;
    N2 = 100;
    M = 1000;
    C = (M+1)*(N2+1)*(N1+1);

    t0 = 0.0;
    t1 = 1.0;

    x11 = 0.0;
    x12 = 1.0;

    x21 = 0.0;
    x22 = 1.0;

    h1 = (x12 - x11)/N1;
    h2 = (x22 - x21)/N2;
    dt = (t1 - t0)/M;

    a1 = 1.0;
    a2 = 1.0;
}

Heat2DControl::~Heat2DControl()
{
}

double Heat2DControl::fx(const DoubleVector &f0)
{
    calculateU(f0);

    double sum = 0.0;

    //printLayer(mu);

    for (int j=0; j<N2; j++)
    {
        for (int i=0; i<N1; i++)
        {
            int j1 = j;
            int j2 = j+1;
            int i1 = i;
            int i2 = i+1;

            double f1 = mu[j1][i1] - U(i1*h1, j1*h2);
            double f2 = mu[j1][i2] - U(i2*h1, j1*h2);
            double f3 = mu[j2][i1] - U(i1*h1, j2*h2);
            double f4 = mu[j2][i2] - U(i2*h1, j2*h2);

            //printf("%d %d %.10f %.10f %.10f %.10f\n", j, i, f1, f2, f3, f4);

            sum = sum + (0.25*(h1*h2))*(f1*f1 + f2*f2 + f3*f3 + f4*f4);
        }
    }

    double norm = 0.0;
    for (int k=0; k<=M; k++)
    {
        for (int j=0; j<=N2; j++)
        {
            for (int i=0; i<=N1; i++)
            {
                double a = f0[k*(N2+1)*(N1+1) + j*(N1+1) + i];
                norm += (a - f(i*h1, j*h2, k*dt))*(a - f(i*h1, j*h2, k*dt));
            }
        }
    }
    sum = sum + sqrt(norm);

    //printf("sum : %.10f\n", sum);
    return sum;
}

void Heat2DControl::gradient(double step, const DoubleVector &f, DoubleVector &g)
{
    calculateU(f);
    calculateP(f, g);
}

double Heat2DControl::u(double x1, double x2, double t) { return x1*x1 + x2*x2 + t*t; }
double Heat2DControl::U(double x1, double x2) { return x1*x1 + x2*x2 + 1.0; }

double Heat2DControl::f(double x1, double x2, double t) { return 2.0*t - 4.0; }
double Heat2DControl::fi(double x1, double x2) { return x1*x1 + x2*x2; }
double Heat2DControl::m1(double x2, double t) { return x2*x2 + t*t; }
double Heat2DControl::m2(double x2, double t) { return x2*x2 + t*t + 1.0; }
double Heat2DControl::m3(double x1, double t) { return x1*x1 + t*t; }
double Heat2DControl::m4(double x1, double t) { return x1*x1 + t*t + 1.0; }

double Heat2DControl::psi_fi(int i, int j)
{
    return 0.0;//2*mu(0, 0, t1);// - U(x1, x2);
}
double Heat2DControl::psi_m1(double x2, double t) { return 0.0; }
double Heat2DControl::psi_m2(double x2, double t) { return 0.0; }
double Heat2DControl::psi_m3(double x1, double t) { return 0.0; }
double Heat2DControl::psi_m4(double x1, double t) { return 0.0; }

void Heat2DControl::calculateU(const DoubleVector& f)
{
    DoubleMatrix u0;
    DoubleMatrix u1;

    u0.resize(N2+1);
    u1.resize(N2+1);

    for (int j=0; j<N2+1; j++)
    {
        u0[j].resize(N1+1);
        u1[j].resize(N1+1);
        for (int i=0; i<N1+1; i++)
        {
            u0[j][i] = fi((h1*i), (h2*j));
        }
    }

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    double alpha1 = 0.0;
    double beta1  = 0.0;
    double alpha2 = 0.0;
    double beta2  = 0.0;

    for (int k=0; k<M; k++)
    {
        // Approximation on x direction
        alpha1 = -(a1*dt)/(2.0*h1*h1);
        beta1  = 1.0 + (a1*dt)/(h1*h1);
        alpha2 = (a2*dt)/(2.0*h2*h2);
        beta2  = 1.0 - (a2*dt)/(h2*h2);

        a.resize(N1-1);
        b.resize(N1-1);
        c.resize(N1-1);
        d.resize(N1-1);
        x.resize(N1-1);

        for (int j=1; j<N2; j++)
        {
            for (int i=1; i<N1; i++)
            {
                a[i-1] = alpha1;
                b[i-1] = beta1;
                c[i-1] = alpha1;
                d[i-1] = alpha2*u0[j-1][i] + beta2*u0[j][i] + alpha2*u0[j+1][i] + (dt/2.0)*f[(k*(N1+1)*(N2+1))+(j)*(N1+1)+i];
            }

            a[0]     = 0.0;
            c[N1-2]  = 0.0;
            d[0]    -= alpha1 * m1(h2*j, dt*(k+0.5));
            d[N1-2] -= alpha1 * m2(h2*j, dt*(k+0.5));
            TomasAlgorithm(a, b, c, d, x);
            for (int i=1; i<N1; i++)
            {
                u1[j][i] = x[i-1];
            }
            u1[j][0]  = m1(h2*j, dt*(k+0.5));
            u1[j][N1] = m2(h2*j, dt*(k+0.5));
        }

        for (int i=0; i<=N1; i++)
        {
            u1[0][i]  = m3(h1*i, dt*(k+0.5));
            u1[N2][i] = m4(h1*i, dt*(k+0.5));
        }

        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        //        printf("\nk = %f\n", k+0.5);
        //        printLayer(u1);

        ///////////////////////////////////////////////////////////////////////
        // Approximation on y direktion
        alpha1 = -(a2*dt)/(2.0*h2*h2);
        beta1  = 1.0 + (a2*dt)/(h2*h2);
        alpha2 = (a1*dt)/(2.0*h1*h1);
        beta2  = 1.0 - (a1*dt)/(h1*h1);

        a.resize(N2-1);
        b.resize(N2-1);
        c.resize(N2-1);
        d.resize(N2-1);
        x.resize(N2-1);

        for (int i=1; i<N1; i++)
        {
            for (int j=1; j<N2; j++)
            {
                a[j-1] = alpha1;
                b[j-1] = beta1;
                c[j-1] = alpha1;
                d[j-1] = alpha2*u1[j][i-1] + beta2*u1[j][i] + alpha2*u1[j][i+1] + (dt/2.0)*f[(k*(N1+1)*(N2+1))+(j)*(N1+1)+i];
            }
            a[0]     = 0.0;
            c[N2-2]  = 0.0;
            d[0]    -= alpha1 * m3(h1*i, dt*(k+1.0));
            d[N2-2] -= alpha1 * m4(h1*i, dt*(k+1.0));
            TomasAlgorithm(a, b, c, d, x);
            for (int j=1; j<N2; j++)
            {
                u0[j][i] = x[j-1];
            }
            u0[0][i]  = m3(h1*i, dt*(k+1));
            u0[N2][i] = m4(h1*i, dt*(k+1));
        }

        for (int j=0; j<=N2; j++)
        {
            u0[j][0]  = m1(h2*j, dt*(k+1));
            u0[j][N1] = m2(h2*j, dt*(k+1));
        }

        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        //        printf("k = %d\n", k);
        //        printLayer(u0);
    }

    mu = u0;
//    printLayer(mu);

    for (int j=0; j<N2+1; j++)
    {
        u1[j].clear();
        u0[j].clear();
    }
    u1.clear();
    u0.clear();
}

void Heat2DControl::calculateP(const DoubleVector& f0, DoubleVector& g)
{
    DoubleMatrix p0;
    DoubleMatrix p1;

    p0.resize(N2+1);
    p1.resize(N2+1);

    dt = -dt;

    for (int j=0; j<N2+1; j++)
    {
        p0[j].resize(N1+1);
        p1[j].resize(N1+1);
        for (int i=0; i<N1+1; i++)
        {
            p0[j][i] = -2.0*(mu[j][i] - U(i*h1, j*h2));
            g[M*(N2+1)*(N1+1) + j*(N1+1) + i] = -p0[j][i] + 2.0*(f0[M*(N2+1)*(N1+1) + j*(N1+1) + i] - f(i*h1, j*h2, t1));
        }
    }

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    double alpha1 = 0.0;
    double beta1  = 0.0;
    double alpha2 = 0.0;
    double beta2  = 0.0;

    for (int k=0; k<M; k++)
    {
        // Approximation on x direction
        alpha1 = (a1*dt)/(2.0*h1*h1);
        beta1  = 1.0 - (a1*dt)/(h1*h1);
        alpha2 = -(a2*dt)/(2.0*h2*h2);
        beta2  = 1.0 + (a2*dt)/(h2*h2);

        a.resize(N1-1);
        b.resize(N1-1);
        c.resize(N1-1);
        d.resize(N1-1);
        x.resize(N1-1);

        for (int j=1; j<N2; j++)
        {
            for (int i=1; i<N1; i++)
            {
                a[i-1] = alpha1;
                b[i-1] = beta1;
                c[i-1] = alpha1;
                d[i-1] = alpha2*p0[j-1][i] + beta2*p0[j][i] + alpha2*p0[j+1][i];
            }
            a[0]     = 0.0;
            c[N1-2]  = 0.0;
            d[0]    -= alpha1 * psi_m1(h2*j, dt*(k+0.5));
            d[N1-2] -= alpha1 * psi_m2(h2*j, dt*(k+0.5));
            TomasAlgorithm(a, b, c, d, x);
            for (int i=1; i<N1; i++)
            {
                p1[j][i] = x[i-1];
            }
            p1[j][0]  = psi_m1(h2*j, dt*(k-0.5));
            p1[j][N1] = psi_m2(h2*j, dt*(k-0.5));
        }

        for (int i=0; i<=N1; i++)
        {
            p1[0][i]  = psi_m3(h1*i, dt*(k-0.5));
            p1[N2][i] = psi_m4(h1*i, dt*(k-0.5));
        }
        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        ///////////////////////////////////////////////////////////////////////
        // Approximation on y direktion
        alpha1 = (a2*dt)/(2.0*h2*h2);
        beta1  = 1.0 - (a2*dt)/(h2*h2);
        alpha2 = -(a1*dt)/(2.0*h1*h1);
        beta2  = 1.0 + (a1*dt)/(h1*h1);

        a.resize(N2-1);
        b.resize(N2-1);
        c.resize(N2-1);
        d.resize(N2-1);
        x.resize(N2-1);
        for (int i=1; i<N1; i++)
        {
            for (int j=1; j<N2; j++)
            {
                a[j-1] = alpha1;
                b[j-1] = beta1;
                c[j-1] = alpha1;
                d[j-1] = alpha2*p1[j][i-1] + beta2*p1[j][i] + alpha2*p1[j][i+1];
            }
            a[0]     = 0.0;
            c[N2-2]  = 0.0;
            d[0]    -= alpha1 * psi_m3(h1*i, dt*(k+1.0));
            d[N2-2] -= alpha1 * psi_m4(h1*i, dt*(k+1.0));
            TomasAlgorithm(a, b, c, d, x);
            for (int j=1; j<N2; j++)
            {
                p0[j][i] = x[j-1];
            }
            p0[0][i]  = psi_m3(h1*i, dt*(k-1));
            p0[N2][i] = psi_m4(h1*i, dt*(k-1));
        }
        for (int j=0; j<=N2; j++)
        {
            p0[j][0]  = psi_m1(h2*j, dt*(k-1));
            p0[j][N1] = psi_m2(h2*j, dt*(k-1));
        }
        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        int main_idx = (M-k-1)*(N2+1)*(N1+1);
        for (int j=0; j<N2+1; j++)
        {
            p0[j].resize(N1+1);
            p1[j].resize(N1+1);
            for (int i=0; i<N1+1; i++)
            {
                p0[j][i] = -2.0*(mu[j][i] - U(i*h1, j*h2));
                g[main_idx + j*(N1+1) + i] = -p0[j][i] + 2.0*(f0[main_idx + j*(N1+1) + i]- f(i*h1, j*h2, (M-k-1)*dt));
            }
        }

    }
    dt = -dt;
    //puts("---");
    //printLayer(p0);

//    for (int k=0; k<=M; k++)
//    {
//        if (k%100==0)
//        {
//            printf("k=%d\n", k);
//            for (int j=0; j<=N2; j++)
//            {
//                if (j%10==0)
//                {
//                    for (int i=0; i<=N1; i++)
//                    {
//                        if (i%10==0) printf("%.10f ", mg[k*(N2+1)*(N1+1) + j*(N1+1) + i]);
//                    }
//                    puts("");
//                }
//            }
//        }
//    }
}

void Heat2DControl::main()
{
    /* Function */
    Heat2DControl hc;

    DoubleVector f0;
    f0.resize(hc.C);

    for (int i=0; i<f0.size(); i++)
        f0[i] = 2.01*(hc.dt*(i/((hc.N1+1)*(hc.N2+1))))-2.0*hc.a1-2.0*hc.a2;

    /* Minimization */
//    SteepestDescentGradient g1;
//    g1.setFunction(&hc);
//    g1.setEpsilon1(0.0000001);
//    g1.setEpsilon2(0.0000001);
//    g1.setGradientStep(0.000001);
//    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
//    g1.setPrinter(new Heat2DControlPrinter);
//    g1.calculate(hc.mf);

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.1, 0.0000001);
    g2.setPrinter(new Heat2DControlPrinter);
    g2.calculate(f0);

    for (int k=0; k<=hc.M; k++)
    {
        if (k%100==0)
        {
            printf("k=%d\n", k);
            for (int j=0; j<=hc.N2; j++)
            {
                if (j%10==0)
                {
                    for (int i=0; i<=hc.N1; i++)
                    {
                        if (i%10==0) printf("%.10f ", f0[k*(hc.N2+1)*(hc.N1+1) + j*(hc.N1+1) + i]);
                    }
                    puts("");
                }
            }
        }
    }
}

void Heat2DControlPrinter::print(unsigned int iterationCount, const DoubleVector& x, const DoubleVector &s, double alpha, RnFunction* f) const
{
    //Heat2DControl *c = dynamic_cast<CFunction3*>(f);

    printf("J[%2d]: %.16f  \n", iterationCount, f->fx(x));
}
