#include "cfunction3.h"

CFunction3::CFunction3(double t0, double t1, double h) : RnFunction()
{
    this->t0 = t0;
    this->t1 = t1;
    this->h = h;
    n = (int)ceil((t1-t0)/h) + 1;

    t = DoubleVector(n);
    x1 = DoubleVector(n);
    x2 = DoubleVector(n);
    psi1 = DoubleVector(n);
    psi2 = DoubleVector(n);

    for (int i=0; i<n; i++)
    {
        t[i] = i*h;
        x1[i] = 0.0;
        x2[i] = 0.0;
        psi1[i] = 0.0;
        psi2[i] = 0.0;
    }
}

CFunction3::~CFunction3()
{}

double CFunction3::fx(const DoubleVector& u)
{
    calculate_x(u);

    double sum = 0.0;
    int i=0;
    DoubleVector x(2);
    DoubleVector U(2);
    for (i=0; i<(n-1); i++)
    {
        int j=i+1;

        x[0] = x1[j];
        x[1] = x2[j];
        U[0] = u[j];
        U[1] = u[j+n];
        double fj = fx0(t[j], x, U);

        x[0] = x1[i];
        x[1] = x2[i];
        U[0] = u[i];
        U[1] = u[i+n];
        double fi = fx0(t[i], x, U);

        sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
    }

    x[0] = x1[n-1];
    x[1] = x2[n-1];
    sum = sum + T(t[n-1], x);
    return sum;
}

void CFunction3::gradient(const DoubleVector& u, DoubleVector &g, double gradient_step)
{
    DoubleVector x(2);
    DoubleVector psi(2);

    calculate_x(u);
    calculate_psi(u);

    DoubleVector g1(n);
    DoubleVector g2(n);

    for (int i=0; i<n; i++)
    {
        int j = i+n;

        x[0] = x1[i];
        x[1] = x2[i];

        psi[0] = psi1[i];
        psi[1] = psi2[i];

        g[i] = -2.0*(u[i]-2.0*t[i]) + psi1[i];
        g[j] = -2.0*(u[j]-3.0*t[i]) + psi2[i];

        g1[i] = g[i];
        g2[i] = g[j];

        //        DoubleVector up(2);
        //        DoubleVector ud(2);

        //        up[0] = u[i] + gradient_step;
        //        up[1] = u[j];

        //        ud[0] = u[i] - gradient_step;
        //        ud[1] = u[j];

        //        g[i] = (H(t[i], x, up, psi) - H(t[i], x, ud, psi)) / (2 * gradient_step);

        //        up[0] = u[i];
        //        up[1] = u[j] + gradient_step;

        //        ud[0] = u[i];
        //        ud[1] = u[j] - gradient_step;
        //        g[i+n] = (H(t[i], x, up, psi) - H(t[i], x, ud, psi)) / (2 * gradient_step);
    }
}

double CFunction3::fx0(double t, const DoubleVector& x, const DoubleVector& u) const
{
    double x1 = x[0];
    double x2 = x[1];
    double u1 = u[0];
    double u2 = u[1];
    return (x1-t*t)*(x1-t*t) + (x2-t*t*t)*(x2-t*t*t) + (u1-2.0*t)*(u1-2.0*t) + (u2-3.0*t)*(u2-3.0*t);
}

double CFunction3::T(double t, const DoubleVector& x) const
{
    double x1 = x[0];
    double x2 = x[1];
    return (x1 - 1.0)*(x1 - 1.0) + (x2 - 1.0)*(x2 - 1.0);
}

double CFunction3::fx1(double t, const DoubleVector& x, const DoubleVector& u) const
{
    //    double x1 = x[0];
    double x2 = x[1];
    double u1 = u[0];
    //    double u2 = u[1];
    return 2.0*t + (x2-t*t*t) + (u1-2.0*t);
}

double CFunction3::fx2(double t, const DoubleVector& x, const DoubleVector& u) const
{
    double x1 = x[0];
    //    double x2 = x[1];
    //    double u1 = u[0];
    double u2 = u[1];
    return 3.0*t*t + (x1-t*t) + (u2-3.0*t);
}

double CFunction3::fp1(double t, const DoubleVector& x, const DoubleVector& psi, const DoubleVector& u)
{
    double x1 = x[0];
    double psi2 = psi[1];
    return 2.0 * (x1-t*t) - psi2;
}

double CFunction3::fp2(double t, const DoubleVector& x, const DoubleVector& psi, const DoubleVector& u)
{
    double x2 = x[1];
    double psi1 = psi[0];
    return 2.0 * (x2-t*t*t) - psi1;
}

double CFunction3::H(double t, const DoubleVector& x, const DoubleVector& u, const DoubleVector& psi)
{
    return -1.0 * fx0(t, x, u) + psi[0] * fx1(t, x, u) + psi[1] * fx2(t, x, u);
}

void CFunction3::calculate_x(const DoubleVector& u)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    x1[0] = 0.0;
    x2[0] = 0.0;
    DoubleVector x(2);
    x[0] = x1[0];
    x[1] = x2[0];

    DoubleVector U(2);

    h = +fabs(h);
    for (int i=0; i<n-1; i++)
    {
        U[0] = u[i];
        U[1] = u[i+n];

        x[0] = x1[i];
        x[1] = x2[i];
        k1[0] = fx1(t[i], x, U);
        k1[1] = fx2(t[i], x, U);

        x[0] = x1[i] + (h/2.0) * k1[0];
        x[1] = x2[i] + (h/2.0) * k1[1];
        k2[0] = fx1(t[i]+h/2.0, x, U);
        k2[1] = fx2(t[i]+h/2.0, x, U);

        x[0] = x1[i] + (h/2.0) * k2[0];
        x[1] = x2[i] + (h/2.0) * k2[1];
        k3[0] = fx1(t[i]+h/2.0, x, U);
        k3[1] = fx2(t[i]+h/2.0, x, U);

        x[0] = x1[i] + h * k3[0];
        x[1] = x2[i] + h * k3[1];
        k4[0] = fx1(t[i]+h, x, U);
        k4[1] = fx2(t[i]+h, x, U);

        x1[i+1] = x1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        x2[i+1] = x2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
}

void CFunction3::calculate_psi(const DoubleVector& u)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    DoubleVector x(2);

    psi1[n-1] = -2.0*(x1[n-1] - 1.0);
    psi2[n-1] = -2.0*(x2[n-1] - 1.0);

    //    DoubleVector _x1(2);
    //    DoubleVector _x2(2);
    //    double dx = 0.000001;
    //    _x1[0] = x1[n-1] + dx;
    //    _x1[1] = x2[n-1];
    //    _x2[0] = x1[n-1] - dx;
    //    _x2[1] = x2[n-1];
    //    psi1[n-1] = -1.0*(T(0.0, _x1)-T(0.0, _x2))/(2.0*dx);
    //    _x1[0] = x1[n-1];
    //    _x1[1] = x2[n-1] + dx;
    //    _x2[0] = x1[n-1];
    //    _x2[1] = x2[n-1] - dx;
    //    psi2[n-1] = -1.0*(T(0.0, _x1)-T(0.0, _x2))/(2.0*dx);

    DoubleVector psi(2);
    psi[0] = psi1[n-1];
    psi[1] = psi2[n-1];

    DoubleVector U(2);

    h = -fabs(h);
    for (int i=n-1; i>0; i--)
    {
        x[0] = x1[i];
        x[1] = x2[i];
        psi[0] = psi1[i];
        psi[1] = psi2[i];

        U[0] = u[i];
        U[1] = u[n+i];

        k1[0] = fp1(t[i], x, psi, U);
        k1[1] = fp2(t[i], x, psi, U);
        psi[0] = psi1[i] + (h/2.0) * k1[0];
        psi[1] = psi2[i] + (h/2.0) * k1[1];
        k2[0] = fp1(t[i]+h/2.0, x, psi, U);
        k2[1] = fp2(t[i]+h/2.0, x, psi, U);
        psi[0] = psi1[i] + (h/2.0) * k2[0];
        psi[1] = psi2[i] + (h/2.0) * k2[1];
        k3[0] = fp1(t[i]+h/2.0, x, psi, U);
        k3[1] = fp2(t[i]+h/2.0, x, psi, U);
        psi[0] = psi1[i] + h * k3[0];
        psi[1] = psi2[i] + h * k3[1];
        k4[0] = fp1(t[i]+h, x, psi, U);
        k4[1] = fp2(t[i]+h, x, psi, U);
        psi1[i-1] = psi1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        psi2[i-1] = psi2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
}

void CFunction3::main()
{
    /* Function */
    CFunction3 c(0.0, 1.0, 0.0001);

    /* initial point */
    DoubleVector u0(2*c.n);
    puts("-----------------------------------------------------------------");
    for (int i=0; i<c.n; i++) u0[i]     = 2.0 * c.t[i];
    for (int i=0; i<c.n; i++) u0[i+c.n] = 3.0 * c.t[i];

    for (int i=0; i<2*c.n; i++) u0[i] = 0.1;

    /* Minimization */
    SteepestDescentGradient g1;
    g1.setFunction(&c);
    g1.setEpsilon1(0.0000001);
    g1.setEpsilon2(0.0000001);
    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
    g1.setPrinter(&c);
    g1.calculate(u0);


    /* Minimization */
//    ConjugateGradient g2;
//    g2.setFunction(&c);
//    g2.setEpsilon1(0.0000001);
//    g2.setEpsilon2(0.0000001);
//    g2.setGradientStep(0.0000001);
//    g2.setR1MinimizeEpsilon(0.01, 0.0000001);
//    g2.setPrinter(new CFunction3Printer);
//    g2.calculate(u0);
}

void CFunction3::print(unsigned int iterationCount, const DoubleVector& u, const DoubleVector &s, double alpha, RnFunction* f) const
{
    CFunction3 *c = dynamic_cast<CFunction3*>(f);

    printf("J[%2d]: %.10f  \n", iterationCount, f->fx(u));

    DoubleVector u1(c->n);
    DoubleVector u2(c->n);

    for (int i=0; i<c->n; i++)
    {
        u1[i] = u[i];
        u2[i] = u[i+c->n];
    }
    Printer::printVector(c->x1, "x1", 10);
    Printer::printVector(c->x2, "x2", 10);
    Printer::printVector(c->psi1, "p1", 10);
    Printer::printVector(c->psi2, "p2", 10);
    Printer::printVector(u1, "u1", 10);
    Printer::printVector(u2, "u2", 10);
}

void CFunction3::print(const char* s, const std::vector<double>& x) const
{
    unsigned int i;

    unsigned int n1 = x.size()/2.0;
    unsigned int n2 = x.size();

    printf("----------- %d %d\n", n1, n2);

    printf("%s1: ", s);
    for (i=0; i<n1; i++)
    {
        if ( i%((n1-1)/10) == 0 )
        {
            if (x[i] < 0)
            {
                printf("%10.8f", x[i]);
            }
            else
            {
                printf("%+10.8f", x[i]);
            }
        }
        if ( i%((n1-1)/10) == 0 && i != n1-1 )
        {
            printf(", ");
        }
    }
    printf("\n");

    printf("%s2: ", s);
    for (i=n1; i<n2; i++)
    {
        if ( i%((n1-1)/10) == 0 )
        {
            if (x[i] < 0)
            {
                printf("%10.8f", x[i]);
            }
            else
            {
                printf("%+10.8f", x[i]);
            }
        }
        if ( i%((n1-2)/10) == 0 && i != n1-1 )
        {
            printf(", ");
        }
    }
    printf("\n");

}

