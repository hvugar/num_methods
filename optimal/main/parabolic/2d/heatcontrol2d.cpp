#include "heatcontrol2d.h"

void HeatControl2D::main()
{
    /* Function */
    HeatControl2D hc(100, 10, 10);

    DoubleVector f((hc.M+1)*(hc.N2+1)*(hc.N1+1));
    for (unsigned int i=0; i<f.size(); i++)
    {
        f[i] = 2.0;//2.0*t - 4.0;
    }

    /* Minimization */
    ConjugateGradient g2;
    g2.setGradient(&hc);
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setEpsilon3(0.0000001);
    g2.setR1MinimizeEpsilon(0.1, 0.0000001);
    g2.setPrinter(&hc);
    g2.setNormalize(true);
    g2.calculate(f);

    IPrinter::printCube(f, hc.M, hc.N2, hc.N1);
}

HeatControl2D::HeatControl2D(unsigned int m, unsigned int n2, unsigned int n1)
{
    N1 = n1;
    N2 = n2;
    M  = m;

    t0 = 0.0;
    t1 = 1.0;

    x10 = 0.0;
    x11 = 1.0;

    x20 = 0.0;
    x21 = 1.0;

    a1 = a2 = 1.0;

    ht = (t1-t0)/M;
    h1 = (x11-x10)/N1;
    h2 = (x21-x20)/N2;

    U.resize(N2+1);
    for (unsigned int j=0; j<=N2; j++)
    {
        U[j].resize(N1+1);
        for (unsigned int i=0; i<=N1; i++) U[j][i] = u(i*h1, j*h2, M*ht);
    }
}

HeatControl2D::~HeatControl2D() {}

double HeatControl2D::fx(const DoubleVector &f)
{
    pf = &f;
    DoubleMatrix u;
    IParabolicEquation2D::caluclateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

    double sum = 0.0;
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            double alpha = 1.0;
            if (j==0 || j==N2 || i==0 || i==N1) alpha = 0.5;
            if (j==0 && i==0)   alpha = 0.25;
            if (j==0 && i==N1)  alpha = 0.25;
            if (j==N2 && i==0)  alpha = 0.25;
            if (j==N2 && i==N1) alpha = 0.25;
            sum += alpha*(u[j][i] - U[j][i])*(u[j][i] - U[j][i]);
        }
    }
    sum = h1*h2*sum;

    double norm = 0.0;
    for (unsigned int k=0; k<M; k++)
    {
        for (unsigned int j=0; j<N2; j++)
        {
            for (unsigned int i=0; i<N1; i++)
            {
                unsigned int i1 = i;
                unsigned int i2 = i+1;
                unsigned int j1 = j;
                unsigned int j2 = j+1;
                unsigned int k1 = k;
                unsigned int k2 = k+1;

                double f1 = f[k1*(N2+1)*(N1+1) + j1*(N1+1) + i1] - fxt(i1*h1, j1*h2, k1*ht);
                double f2 = f[k1*(N2+1)*(N1+1) + j1*(N1+1) + i2] - fxt(i2*h1, j1*h2, k1*ht);
                double f3 = f[k1*(N2+1)*(N1+1) + j2*(N1+1) + i1] - fxt(i1*h1, j2*h2, k1*ht);
                double f4 = f[k1*(N2+1)*(N1+1) + j2*(N1+1) + i2] - fxt(i2*h1, j2*h2, k1*ht);
                double f5 = f[k2*(N2+1)*(N1+1) + j1*(N1+1) + i1] - fxt(i1*h1, j1*h2, k2*ht);
                double f6 = f[k2*(N2+1)*(N1+1) + j1*(N1+1) + i2] - fxt(i2*h1, j1*h2, k2*ht);
                double f7 = f[k2*(N2+1)*(N1+1) + j2*(N1+1) + i1] - fxt(i1*h1, j2*h2, k2*ht);
                double f8 = f[k2*(N2+1)*(N1+1) + j2*(N1+1) + i2] - fxt(i2*h1, j2*h2, k2*ht);

                norm += f1*f1 + f2*f2 + f3*f3 + f4*f4 + f5*f5 + f6*f6 + f7*f7 + f8*f8;
            }
        }
    }
    norm = norm * (h1*h2*ht)*0.125;
    return sum + norm;
}

void HeatControl2D::gradient(const DoubleVector &f, DoubleVector &g)
{
    pf = &f;
    DoubleMatrix u;
    IParabolicEquation2D::caluclateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);
    pu = &u;
    DoubleCube psi;
    IBackwardParabolicEquation2D::caluclateMVD(psi, h1, h2, ht, N1, N2, M, a1, a2);
    //calculateP(f0, u, g);
    for (unsigned int k=0; k<=M; k++)
    {
        // calculating gradient
        for (unsigned int j=0; j<=N2; j++)
        {
            for (unsigned i=0; i<=N1; i++)
            {
                int index = k*(N1+1)*(N2+1)+j*(N1+1)+i;
                g[index] = -psi[k][j][i] + 2*(f[index] - fxt(i*h1, j*h2, k*ht));
            }
        }
    }
    psi.clear();
}

double HeatControl2D::u(double x1, double x2, double t) const
{
    return x1*x1 + x2*x2 + t*t;
}

double HeatControl2D::initial(unsigned int i, unsigned int j) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    return u(x1, x2, t0);
}

double HeatControl2D::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
   double x1 = i*h1;
   double x2 = j*h2;
   double t  = 0.5*k*ht;
   return u(x1, x2, t);
}

double HeatControl2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    unsigned int n = (k%2==0 ? k/2 : (k+1)/2)*(N1+1)*(N2+1)+j*(N1+1)+i;
    return (*pf)[n];
}

double HeatControl2D::binitial(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u = *pu;
    return -2.0*(u[j][i] - U[j][i]);
}

double HeatControl2D::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HeatControl2D::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HeatControl2D::fxt(double x1, double x2, double t)
{
    C_UNUSED(x1);
    C_UNUSED(x2);
    return 2.0*t - 2.0*a1 - 2.0*a2;
}

void HeatControl2D::print(unsigned int i, const DoubleVector &f0, const DoubleVector &s, double alpha, RnFunction *f) const
{
    C_UNUSED(i);
    C_UNUSED(s);
    C_UNUSED(alpha);

    HeatControl2D *hc = dynamic_cast<HeatControl2D*>(f);
    printf("J: %.16f\n", hc->fx(f0));
}

