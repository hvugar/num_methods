#include "cfunction2.h"

CFunction2::CFunction2(double t0, double t1, double h)
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

CFunction2::~CFunction2()
{}

double CFunction2::fx(const DoubleVector& u)
{
    calculate_x(u);

    double sum = 0.0;
    int i=0;
    DoubleVector x(2);
    for (i=0; i<(n-1); i++)
    {
        int j=i+1;
        x[0] = x1[j]; x[1] = x2[j];
        double fj = fx0(t[j], x, u[j]);
        x[0] = x1[i]; x[1] = x2[i];
        double fi = fx0(t[i], x, u[i]);
        sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
    }
    return sum;
}

void CFunction2::gradient(const DoubleVector& u, DoubleVector &g, double gradient_step)
{
    DoubleVector x(2);
    DoubleVector psi(2);
    calculate_x(u);
    calculate_psi(u);
    for (int i=0; i<n; i++)
    {
        x[0] = x1[i];
        x[1] = x2[i];

        psi[0] = psi1[i];
        psi[1] = psi2[i];

        //double u1 = u[i] + gradient_step;
        //double u2 = u[i] - gradient_step;
        //g[i] = (H(t[i], x, u1, psi) - H(t[i], x, u2, psi)) / (2 * gradient_step);
        g[i] = psi2[i];
    }
}

double CFunction2::fx0(double t, const DoubleVector& x, double u) const
{
    return (x[0]-(t*t)/2.0)*(x[0]-(t*t)/2.0) + (x[1] - t)*(x[1] - t);
}

double CFunction2::F(double t, const DoubleVector &x, double u) const
{
    return 0.0;
}

double CFunction2::fx1(double t, const DoubleVector& x, double u) const
{
    return x[1];
}

double CFunction2::fx2(double t, const DoubleVector& x, double u) const
{
    return -6.0*x[0] + x[1] + u - t + 1.0;
}

double CFunction2::fp1(double t, const DoubleVector& x, const DoubleVector& psi, double u) const
{
//    double h = 0.000001;
//    DoubleVector x1(2);
//    x1[0] = x[0] + h; x1[1] = x[1];
//    DoubleVector x2(2);
//    x2[0] = x[0] - h; x2[1] = x[1];
//    return -1.0 * (H(t, x1, u, psi) - H(t, x2, u, psi)) / (2 * h);

    return 2.0*x[0] - t*t + 6.0*psi[1];
}

double CFunction2::fp2(double t, const DoubleVector& x, const DoubleVector& psi, double u) const
{
//    double h = 0.000001;
//    DoubleVector x1(2);
//    x1[0] = x[0]; x1[1] = x[1] + h;
//    DoubleVector x2(2);
//    x2[0] = x[0]; x2[1] = x[1] - h;
//    return -1.0 * (H(t, x1, u, psi) - H(t, x2, u, psi)) / (2 * h);

    return 2.0*(x[1]-t) - psi[0] - psi[1];
}

double CFunction2::H(double t, const DoubleVector& x, double u, const DoubleVector& psi) const
{
    return -fx0(t, x, u) + psi[0] * fx1(t, x, u) + psi[1] * fx2(t, x, u);
}

void CFunction2::calculate_x(const DoubleVector& u)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    x1[0] = 0.0;
    x2[0] = 0.0;
    DoubleVector _x(2);
    _x[0] = x1[0];
    _x[1] = x2[0];

    h = +fabs(h);
    for (int i=0; i<n-1; i++)
    {
        _x[0] = x1[i];
        _x[1] = x2[i];
        k1[0] = fx1(t[i], _x, u[i]);
        k1[1] = fx2(t[i], _x, u[i]);

        _x[0] = x1[i] + (h/2.0) * k1[0];
        _x[1] = x2[i] + (h/2.0) * k1[1];
        k2[0] = fx1(t[i]+h/2.0, _x, u[i]);
        k2[1] = fx2(t[i]+h/2.0, _x, u[i]);

        _x[0] = x1[i] + (h/2.0) * k2[0];
        _x[1] = x2[i] + (h/2.0) * k2[1];
        k3[0] = fx1(t[i]+h/2.0, _x, u[i]);
        k3[1] = fx2(t[i]+h/2.0, _x, u[i]);

        _x[0] = x1[i] + h * k3[0];
        _x[1] = x2[i] + h * k3[1];
        k4[0] = fx1(t[i]+h, _x, u[i]);
        k4[1] = fx2(t[i]+h, _x, u[i]);

        x1[i+1] = x1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        x2[i+1] = x2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
}

void CFunction2::calculate_psi(const DoubleVector& u)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    DoubleVector _x(2);

    psi1[n-1] = 0.0;
    psi2[n-1] = 0.0;

//    DoubleVector _x1(2);
//    DoubleVector _x2(2);
//    double dx = 0.000001;
//    _x1[0] = x1[n-1] + dx;
//    _x1[1] = x2[n-1];
//    _x2[0] = x1[n-1] - dx;
//    _x2[1] = x2[n-1];
//    psi1[n-1] = -1.0*(F(0.0, _x1, 0.0)-F(0.0, _x2, 0.0))/(2.0*dx);
//    _x1[0] = x1[n-1];
//    _x1[1] = x2[n-1] + dx;
//    _x2[0] = x1[n-1];
//    _x2[1] = x2[n-1] - dx;
//    psi2[n-1] = -1.0*(F(0.0, _x1, 0.0)-F(0.0, _x2, 0.0))/(2.0*dx);

    DoubleVector _psi(2);
    _psi[0] = psi1[n-1];
    _psi[1] = psi2[n-1];

    h = -fabs(h);
    for (int i=n-1; i>0; i--)
    {
        _x[0] = x1[i];
        _x[1] = x2[i];
        _psi[0] = psi1[i];
        _psi[1] = psi2[i];

        k1[0] = fp1(t[i], _x, _psi, u[i]);
        k1[1] = fp2(t[i], _x, _psi, u[i]);
        _psi[0] = psi1[i] + (h/2.0) * k1[0];
        _psi[1] = psi2[i] + (h/2.0) * k1[1];
        k2[0] = fp1(t[i]+h/2.0, _x, _psi, u[i]);
        k2[1] = fp2(t[i]+h/2.0, _x, _psi, u[i]);
        _psi[0] = psi1[i] + (h/2.0) * k2[0];
        _psi[1] = psi2[i] + (h/2.0) * k2[1];
        k3[0] = fp1(t[i]+h/2.0, _x, _psi, u[i]);
        k3[1] = fp2(t[i]+h/2.0, _x, _psi, u[i]);
        _psi[0] = psi1[i] + h * k3[0];
        _psi[1] = psi2[i] + h * k3[1];
        k4[0] = fp1(t[i]+h, _x, _psi, u[i]);
        k4[1] = fp2(t[i]+h, _x, _psi, u[i]);
        psi1[i-1] = psi1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        psi2[i-1] = psi2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
}

void CFunction2::main()
{
    /* Function */
    CFunction2 c(0.0, 1.0, 0.0001);

    /* initial point */
    DoubleVector u0(c.n);

    for (int i=0; i<c.n; i++) u0[i] = 0.00001;
    /* Minimization */
    SteepestDescentGradient g1;
    g1.setFunction(&c);
    g1.setEpsilon1(0.0000001);
    g1.setEpsilon2(0.0000001);
    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
    g1.setPrinter(&c);
    g1.calculate(u0);

    puts("-----------------------------------------------------------------");
    for (int i=0; i<c.n; i++) u0[i] = 0.00001;
    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&c);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setR1MinimizeEpsilon(0.01, 0.0000001);
    g2.setPrinter(&c);
//    g2.calculate(u0);
}

void CFunction2::print(unsigned int iterationCount, const DoubleVector& u, const DoubleVector &s, double alpha, RnFunction* f) const
{
    CFunction2* f2 = dynamic_cast<CFunction2*>(f);
    printf("J[%2d]: %.10f\n", iterationCount, f->fx(u));
    print("psi1", f2->x1);
    print("psi2", f2->x2);
    print("psi1", f2->psi1);
    print("psi2", f2->psi2);
    print("u", u);
}

void CFunction2::print(const char* s, const std::vector<double>& x) const
{
    unsigned int i;
    unsigned int n = x.size();
    printf("%s:\t", s);
    for (i=0; i<n; i++)
    {
        if ( i%((n-1)/10) == 0 )
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
        if ( i%((n-1)/10) == 0 && i != n-1 )
        {
            printf(", ");
        }
    }
    printf("\n");
}
