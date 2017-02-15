#include "hyperboliccontrol2d21.h"

void HyperbolicControl2D21::Main(int argc, char* argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    HyperbolicControl2D21 hc;
    hc.fx(10.0);

    //hc.file = fopen("hyperboliccontrol2d21.txt", "w");
    //hc.file = stdout;
    //    for (double t=0.1; t<=10.1; t+=0.1)
    //    {
    //        printf("%f %.8f\n", t, hc.fx(t));
    //        //fprintf(hc.file, "%f %.10f\n", t, hc.fx(t));
    //        //fflush(hc.file);
    //    }
    //fclose(hc.file);
}

HyperbolicControl2D21::HyperbolicControl2D21()
{
    h1 = 0.005;
    h2 = 0.005;
    ht = 0.0025;
    N1 = (unsigned)ceil((x11 - x10)/h1);
    N2 = (unsigned)ceil((x21 - x20)/h2);
    M  = (unsigned)ceil((t1 - t0)/ht);
    L = 2;

    e.resize(2);
    e[0] = 0.2;
    e[1] = 0.2;

    alpha0 = 1.0;
    alpha1 = 10.0;
    alpha2 = 2.0;
    alpha3 = 1.0;
    qamma = 0.2;

    a1 = 1.0;
    a2 = 1.0;
}

double HyperbolicControl2D21::fx(double T) const
{
    const_cast<HyperbolicControl2D21*>(this)->t1 = T;
    const_cast<HyperbolicControl2D21*>(this)->M  = (unsigned)ceil((t1 - t0)/ht);

    DoubleMatrix m;
    IHyperbolicEquation2D::calculateU1(m, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    DoubleVector x;
    double rf = fx(x);
    return rf;
}

double HyperbolicControl2D21::fx(const DoubleVector &x) const
{
    C_UNUSED(x);

    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    const DoubleMatrix &u0 = c.matrix(M);
    const DoubleMatrix &u1 = c.matrix(M-2);
    double sum = 0.0;

    double sum1 = 0.0;
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            double k = 1.0;
            if (i==0 || i==N1) k *= 0.5;
            if (j==0 || j==N2) k *= 0.5;
            sum1 = sum1 + k * (u0[j][i]-U0) * (u0[j][i]-U0);

        }
    }
    sum1 = h1*h2*sum1;

    double sum2 = 0.0;
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            double k = 1.0;
            if (i==0 || i==N1) k *= 0.5;
            if (j==0 || j==N2) k *= 0.5;
            sum2 = sum2 + k * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1) * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
        }
    }
    sum2 = h1*h2*sum2;

    sum = sum1 + alpha0*sum2;
    return sum;
}

double HyperbolicControl2D21::initial1(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2D21::initial2(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2D21::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D21::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = fxt(i, j, k);
    return sum;
}

double HyperbolicControl2D21::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    sum += alpha1*exp(-alpha2*((x1-e[0])*(x1-e[0])+(x2-e[1])*(x2-e[1]))-alpha3*t);
    return sum;
}
