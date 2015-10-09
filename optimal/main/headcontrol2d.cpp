#include "headcontrol2d.h"
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

void TomasAlgorithm1(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x);
void printLayer1(const DoubleMatrix& x);

HeadControl2D::HeadControl2D()
{
    t0 = 0.0;
    t1 = 1.0;
    x10 = x20 = 0.0;
    x11 = x21 = 1.0;

    N1 = 100;
    N2 = 100;
    M  = 1000;

    L = 3;
    C = (M+1)*(N2+1)*(N1+1);

    h1 = (x11-x10) / N1;
    h2 = (x21-x20) / N2;
    ht  = (t1 - t0) / M;

    a1 = a2 = 1.0;

    U.resize(N2+1);
    for (unsigned j=0; j<=N2; j++) U[j].resize(N1+1);
}

void HeadControl2D::calculateU(const DoubleVector& E, DoubleMatrix& u)
{
    DoubleMatrix u0;
    DoubleMatrix u1;

    u0.resize(N2+1);
    u1.resize(N2+1);

    for (unsigned int j=0; j<=N2; j++) // x2-->
    {
        u0[j].resize(N1+1);
        u1[j].resize(N1+1);
        for (unsigned int i=0; i<=N1; i++) // x1-->
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

    for (unsigned int k=0; k<M; k++)
    {
        // Approximation on x direction
        alpha1 = -(a1*ht)/(2.0*h1*h1);
        beta1  = 1.0 + (a1*ht)/(h1*h1);
        alpha2 = (a2*ht)/(2.0*h2*h2);
        beta2  = 1.0 - (a2*ht)/(h2*h2);

        a.resize(N1-1);
        b.resize(N1-1);
        c.resize(N1-1);
        d.resize(N1-1);
        x.resize(N1-1);

        for (unsigned int j=1; j<N2; j++)
        {
            for (unsigned int i=1; i<N1; i++)
            {
                a[i-1] = alpha1;
                b[i-1] = beta1;
                c[i-1] = alpha1;
                d[i-1] = alpha2*u0[j-1][i] + beta2*u0[j][i] + alpha2*u0[j+1][i] + (ht/2.0) * f(E, i*h1, j*h2, k);
            }

            a[0]     = 0.0;
            c[N1-2]  = 0.0;
            d[0]    -= alpha1 * m1(h2*j, ht*(k+0.5));
            d[N1-2] -= alpha1 * m2(h2*j, ht*(k+0.5));
            TomasAlgorithm1(a, b, c, d, x);
            for (unsigned int i=1; i<N1; i++)
            {
                u1[j][i] = x[i-1];
            }
            u1[j][0]  = m1(h2*j, ht*(k+0.5));
            u1[j][N1] = m2(h2*j, ht*(k+0.5));
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            u1[0][i]  = m3(h1*i, ht*(k+0.5));
            u1[N2][i] = m4(h1*i, ht*(k+0.5));
        }

        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        ///////////////////////////////////////////////////////////////////////
        // Approximation on y direktion
        alpha1 = -(a2*ht)/(2.0*h2*h2);
        beta1  = 1.0 + (a2*ht)/(h2*h2);
        alpha2 = (a1*ht)/(2.0*h1*h1);
        beta2  = 1.0 - (a1*ht)/(h1*h1);

        a.resize(N2-1);
        b.resize(N2-1);
        c.resize(N2-1);
        d.resize(N2-1);
        x.resize(N2-1);

        for (unsigned int i=1; i<N1; i++)
        {
            for (unsigned int j=1; j<N2; j++)
            {
                a[j-1] = alpha1;
                b[j-1] = beta1;
                c[j-1] = alpha1;
                d[j-1] = alpha2*u1[j][i-1] + beta2*u1[j][i] + alpha2*u1[j][i+1] + (ht/2.0) * f(E, i*h1, j*h2, k);
            }
            a[0]     = 0.0;
            c[N2-2]  = 0.0;
            d[0]    -= alpha1 * m3(h1*i, ht*(k+1.0));
            d[N2-2] -= alpha1 * m4(h1*i, ht*(k+1.0));
            TomasAlgorithm1(a, b, c, d, x);
            for (unsigned int j=1; j<N2; j++)
            {
                u0[j][i] = x[j-1];
            }
            u0[0][i]  = m3(h1*i, ht*(k+1));
            u0[N2][i] = m4(h1*i, ht*(k+1));
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            u0[j][0]  = m1(h2*j, ht*(k+1));
            u0[j][N1] = m2(h2*j, ht*(k+1));
        }

        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();
    }

    u = u0;
    printLayer1(u);

    for (unsigned int j=0; j<=N2; j++)
    {
        u1[j].clear();
        u0[j].clear();
    }
    u1.clear();
    u0.clear();
}

void HeadControl2D::calculateP(const DoubleVector& E, DoubleVector& g)
{
    DoubleMatrix p0;
    DoubleMatrix p1;

    p0.resize(N2+1);
    p1.resize(N2+1);

    ht = -ht;

    for (unsigned int j=0; j<N2+1; j++)
    {
        p0[j].resize(N1+1);
        p1[j].resize(N1+1);
        for (unsigned int i=0; i<N1+1; i++)
        {
            p0[j][i] = -2.0*(mu[j][i] - U[j][i]);
            //g[M*(N2+1)*(N1+1) + j*(N1+1) + i] = -p0[j][i] + 2.0*(f0[M*(N2+1)*(N1+1) + j*(N1+1) + i] - f(i*h1, j*h2, t1));
        }
    }

    for (unsigned int i=0; i<=g.size(); i++) g[i] = 0.0;

    for (unsigned int j=0; j<=N2; j++)
    {
        g[2*M*L + 0] = -g1(t1) * ((p0[j][21] - p0[j][19]) / (2*h1));
        g[2*M*L + 2] = -g2(t1) * ((p0[j][51] - p0[j][49]) / (2*h1));
        g[2*M*L + 4] = -g3(t1) * ((p0[j][81] - p0[j][79]) / (2*h1));
    }

    for (unsigned int i=0; i<=N1; i++)
    {
        g[2*M*L + 1] = -g1(t1) * ((p0[41][i] - p0[39][i]) / (2*h2));
        g[2*M*L + 3] = -g2(t1) * ((p0[21][i] - p0[19][i]) / (2*h2));
        g[2*M*L + 5] = -g3(t1) * ((p0[61][i] - p0[59][i]) / (2*h2));
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

    for (unsigned int k=0; k<M; k++)
    {
        // Approximation on x direction
        alpha1 = (a1*ht)/(2.0*h1*h1);
        beta1  = 1.0 - (a1*ht)/(h1*h1);
        alpha2 = -(a2*ht)/(2.0*h2*h2);
        beta2  = 1.0 + (a2*ht)/(h2*h2);

        a.resize(N1-1);
        b.resize(N1-1);
        c.resize(N1-1);
        d.resize(N1-1);
        x.resize(N1-1);

        for (unsigned int j=1; j<N2; j++)
        {
            for (unsigned int i=1; i<N1; i++)
            {
                a[i-1] = alpha1;
                b[i-1] = beta1;
                c[i-1] = alpha1;
                d[i-1] = alpha2*p0[j-1][i] + beta2*p0[j][i] + alpha2*p0[j+1][i];
            }
            a[0]     = 0.0;
            c[N1-2]  = 0.0;
            d[0]    -= alpha1 * psi_m1(h2*j, ht*(k+0.5));
            d[N1-2] -= alpha1 * psi_m2(h2*j, ht*(k+0.5));
            TomasAlgorithm1(a, b, c, d, x);
            for (unsigned int i=1; i<N1; i++)
            {
                p1[j][i] = x[i-1];
            }
            p1[j][0]  = psi_m1(h2*j, ht*(k-0.5));
            p1[j][N1] = psi_m2(h2*j, ht*(k-0.5));
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            p1[0][i]  = psi_m3(h1*i, ht*(k-0.5));
            p1[N2][i] = psi_m4(h1*i, ht*(k-0.5));
        }
        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        ///////////////////////////////////////////////////////////////////////
        // Approximation on y direktion
        alpha1 = (a2*ht)/(2.0*h2*h2);
        beta1  = 1.0 - (a2*ht)/(h2*h2);
        alpha2 = -(a1*ht)/(2.0*h1*h1);
        beta2  = 1.0 + (a1*ht)/(h1*h1);

        a.resize(N2-1);
        b.resize(N2-1);
        c.resize(N2-1);
        d.resize(N2-1);
        x.resize(N2-1);
        for (unsigned int i=1; i<N1; i++)
        {
            for (unsigned int j=1; j<N2; j++)
            {
                a[j-1] = alpha1;
                b[j-1] = beta1;
                c[j-1] = alpha1;
                d[j-1] = alpha2*p1[j][i-1] + beta2*p1[j][i] + alpha2*p1[j][i+1];
            }
            a[0]     = 0.0;
            c[N2-2]  = 0.0;
            d[0]    -= alpha1 * psi_m3(h1*i, ht*(k+1.0));
            d[N2-2] -= alpha1 * psi_m4(h1*i, ht*(k+1.0));
            TomasAlgorithm1(a, b, c, d, x);
            for (unsigned int j=1; j<N2; j++)
            {
                p0[j][i] = x[j-1];
            }
            p0[0][i]  = psi_m3(h1*i, ht*(k-1));
            p0[N2][i] = psi_m4(h1*i, ht*(k-1));
        }
        for (unsigned int j=0; j<=N2; j++)
        {
            p0[j][0]  = psi_m1(h2*j, ht*(k-1));
            p0[j][N1] = psi_m2(h2*j, ht*(k-1));
        }
        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        //int main_idx = (M-k-1)*(N2+1)*(N1+1);

        for (unsigned int j=0; j<N2+1; j++)
        {
            p0[j].resize(N1+1);
            p1[j].resize(N1+1);
            for (unsigned int i=0; i<N1+1; i++)
            {
                p0[j][i] = -2.0*(mu[j][i] - U[j][i]);
                //g[main_idx + j*(N1+1) + i] = -p0[j][i] + 2.0*(f0[main_idx + j*(N1+1) + i]- f(i*h1, j*h2, (M-k-1)*dt));
            }
        }

        double t = (M-k-1)*ht;
        for (unsigned int j=0; j<=N2; j++)
        {
            g[2*(M-k-1)*L + 0] = -g1(t) * ((p0[j][21] - p0[j][19]) / (2*h1));
            g[2*(M-k-1)*L + 2] = -g2(t) * ((p0[j][51] - p0[j][49]) / (2*h1));
            g[2*(M-k-1)*L + 4] = -g3(t) * ((p0[j][81] - p0[j][79]) / (2*h1));
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            g[2*(M-k-1)*L + 1] = -g1(t) * ((p0[41][i] - p0[39][i]) / (2*h2));
            g[2*(M-k-1)*L + 3] = -g2(t) * ((p0[21][i] - p0[19][i]) / (2*h2));
            g[2*(M-k-1)*L + 5] = -g3(t) * ((p0[61][i] - p0[59][i]) / (2*h2));
        }

    }
    ht = -ht;
}

double HeadControl2D::fx(const DoubleVector& E)
{
    calculateU(E, mu);

    printLayer1(mu);

    double sum = 0.0;

    for (unsigned int j=0; j<N2; j++)
    {
        for (unsigned int i=0; i<N1; i++)
        {
            int j1 = j;
            int j2 = j+1;
            int i1 = i;
            int i2 = i+1;

            double f1 = mu[j1][i1] - U[j1][i1];
            double f2 = mu[j1][i2] - U[j1][i2];
            double f3 = mu[j2][i1] - U[j2][i1];
            double f4 = mu[j2][i2] - U[j2][i2];

            sum = sum + (0.25*(h1*h2))*(f1*f1 + f2*f2 + f3*f3 + f4*f4);
        }
    }

    //    double norm = 0.0;
    //    for (int k=0; k<=M; k++)
    //    {
    //        for (int j=0; j<=N2; j++)
    //        {
    //            for (int i=0; i<=N1; i++)
    //            {
    //                double a = f0[k*(N2+1)*(N1+1) + j*(N1+1) + i];
    //                norm += (a - f(i*h1, j*h2, k*dt))*(a - f(i*h1, j*h2, k*dt));
    //            }
    //        }
    //    }
    //    sum = sum + sqrt(norm);

    return sum;
}

void HeadControl2D::gradient(double step, const DoubleVector& E, DoubleVector& g)
{
    calculateU(E, mu);
    printLayer1(mu);
    calculateP(E, g);
}

double HeadControl2D::u(double x1, double x2, double t)
{
    return x1*x1 + x1*x2 + x2*x2*x2 + x1*x2*t + t*t;
}

double HeadControl2D::fi(double x1, double x2)
{
    return u(x1, x2, t0);
}

double HeadControl2D::m1(double x1, double t)
{
    return u(x1, x20, t);
}

double HeadControl2D::m2(double x1, double t)
{
    return u(x1, x21, t);
}

double HeadControl2D::m3(double x2, double t)
{
    return u(x10, x2, t);
}

double HeadControl2D::m4(double x2, double t)
{
    return u(x11, x2, t);
}

double HeadControl2D::psi_fi(int i, int j)
{
    return -2.0*(mu[j][i] - U[j][i]);
}
double HeadControl2D::psi_m1(double x2, double t) { return 0.0; }
double HeadControl2D::psi_m2(double x2, double t) { return 0.0; }
double HeadControl2D::psi_m3(double x1, double t) { return 0.0; }
double HeadControl2D::psi_m4(double x1, double t) { return 0.0; }

double HeadControl2D::f(const DoubleVector& E, double x1, double x2, unsigned int k)
{
    double sum = 0.0;

    if (fabs(x1 - E[2*L*k + 0]) < 0.000001 && fabs(x2 - E[2*L*k + 1]) < 0.000001) { sum += g1(k*ht); }
    if (fabs(x1 - E[2*L*k + 2]) < 0.000001 && fabs(x2 - E[2*L*k + 3]) < 0.000001) { sum += g2(k*ht); }
    if (fabs(x1 - E[2*L*k + 4]) < 0.000001 && fabs(x2 - E[2*L*k + 5]) < 0.000001) { sum += g3(k*ht); }

    return sum;// + x1*x2 - 2.0 - 6.0*x2 + 2.0*t;
}

void HeadControl2D::initialize()
{
    DoubleVector E;
    E.resize(2*(M+1)*L);

    for (unsigned int i=0; i<=M; i++)
    {
        E[2*L*i + 0] = 0.2;
        E[2*L*i + 1] = 0.4;
        E[2*L*i + 2] = 0.5;
        E[2*L*i + 3] = 0.2;
        E[2*L*i + 4] = 0.8;
        E[2*L*i + 5] = 0.6;
    }

    calculateU(E, U);
}

void HeadControl2D::main()
{
    /* Function */
    HeadControl2D hc;
    hc.initialize();

    DoubleVector E0;
    E0.resize(hc.C);

    for (unsigned int i=0; i<E0.size(); i++)
        E0[i] = 500000.0;

    printf("OK\n");

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
    g2.setPrinter(new HeadControl2DPrinter);
    g2.calculate(E0);

//    for (unsigned int k=0; k<=hc.M; k++)
//    {
//        if (k%100==0)
//        {
//            printf("k=%d\n", k);
//            for (unsigned int j=0; j<=hc.N2; j++)
//            {
//                if (j%10==0)
//                {
//                    for (unsigned int i=0; i<=hc.N1; i++)
//                    {
//                        if (i%10==0) printf("%.10f ", E0[k*(hc.N2+1)*(hc.N1+1) + j*(hc.N1+1) + i]);
//                    }
//                    puts("");
//                }
//            }
//        }
//    }
}

void TomasAlgorithm1(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x)
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

void printLayer1(const DoubleMatrix& x)
{
    int m = x.size()/10;
    for (unsigned int j=0; j<x.size(); j++)
    {
        if (j%m==0)
        {
            for (unsigned int i=0; i<x[j].size(); i++)
            {
                if (i%m==0) printf("%14.10f ", x[j][i]);
            }
            puts("");
        }
    }
}

void HeadControl2DPrinter::print(unsigned int iterationCount, const DoubleVector &m_x, const DoubleVector &s, double m_alpha, RnFunction *f) const
{

}

