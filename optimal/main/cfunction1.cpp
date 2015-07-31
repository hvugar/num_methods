#include "cfunction1.h"
#include "utils.h"

CFunction1::CFunction1(double t0, double t1, double h) : RnFunction()
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

CFunction1::~CFunction1()
{}

double CFunction1::fx(const DoubleVector& u)
{
    calculate_x(u);

    double sum = 0.0;
    int i=0;
    DoubleVector x(2);
    for (i=0; i<(n-1); i++)
    {
        int j=i+1;
        //double fj = (x1[j]-t[j]*t[j]*t[j])*(x1[j]-t[j]*t[j]*t[j]) + (x2[j] - t[j])*(x2[j] - t[j]) + (2*u[j] - t[j])*(2*u[j] - t[j]);
        //double fi = (x1[i]-t[i]*t[i]*t[i])*(x1[i]-t[i]*t[i]*t[i]) + (x2[i] - t[i])*(x2[i] - t[i]) + (2*u[i] - t[i])*(2*u[i] - t[i]);
        x[0] = x1[j]; x[1] = x2[j];
        double fj = fx0(t[j], x, u[j]);
        x[0] = x1[i]; x[1] = x2[i];
        double fi = fx0(t[i], x, u[i]);
        sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
    }
    sum = sum + (x2[n-1] - 1.0) * (x2[n-1] - 1.0);
    return sum;
}

void CFunction1::gradient(double gradient_step, const DoubleVector& u, DoubleVector &g)
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

        double u1 = u[i] + gradient_step;
        double u2 = u[i] - gradient_step;
        g[i] = (H(t[i], x, u1, psi) - H(t[i], x, u2, psi)) / (2 * gradient_step);
        //g[i] = -4.0*(2.0*u[i] - t[i]) - 2.0*psi2[i];
    }
}

double CFunction1::fx0(double t, const DoubleVector& x, double u) const
{
    double x1 = x[0];
    double x2 = x[1];
    double a1 = x1-t*t*t;
    double a2 = x2-t;
    double a3 = 2*u-t;
    return a1*a1 + a2*a2 + a3*a3;
}

double CFunction1::F(double t, const DoubleVector& x, double u) const
{
    return (x[1] - 1.0) * (x[1] - 1.0);
}

double CFunction1::fx1(double t, const DoubleVector& x, double u) const
{
    double x2 = x[1];
    return 3.0*x2*x2;
}

double CFunction1::fx2(double t, const DoubleVector& x, double u) const
{
    double x1 = x[0];
    double x2 = x[1];
    return x1 + x2 - 2.0*u - t*t*t + 1.0;
}

double CFunction1::fp1(double t, const DoubleVector& x, const DoubleVector& psi, double u)
{
    double h = 0.000001;
    DoubleVector x1(2);
    x1[0] = x[0] + h; x1[1] = x[1];
    DoubleVector x2(2);
    x2[0] = x[0] - h; x2[1] = x[1];
    return -1.0 * (H(t, x1, u, psi) - H(t, x2, u, psi)) / (2 * h);

    //double x1 = x[0];
    //double p2 = psi[1];
    //return 2.0 * (x1 - t*t*t) - p2;
}

double CFunction1::fp2(double t, const DoubleVector& x, const DoubleVector& psi, double u)
{
    double h = 0.000001;
    DoubleVector x1(2);
    x1[0] = x[0]; x1[1] = x[1] + h;
    DoubleVector x2(2);
    x2[0] = x[0]; x2[1] = x[1] - h;
    return -1.0 * (H(t, x1, u, psi) - H(t, x2, u, psi)) / (2 * h);

    //double x2 = x[1];
    //double p1 = psi[0];
    //double p2 = psi[1];
    //return 2.0 * (x2 - t) - 6.0 * x2 * p1 - p2;
}

double CFunction1::H(double t, const DoubleVector& x, double u, const DoubleVector& psi)
{
    return -1.0 * fx0(t, x, u) + psi[0] * fx1(t, x, u) + psi[1] * fx2(t, x, u);
}

void CFunction1::calculate_x(const DoubleVector& u)
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

void CFunction1::calculate_psi(const DoubleVector& u)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    DoubleVector _x(2);

    psi1[n-1] = 0.0;
    psi2[n-1] = -2.0 * (x2[n-1] - 1.0);

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

void CFunction1::main()
{
    /* Function */
    CFunction1 c(0.0, 1.0, 0.001);

    /* initial point */
    DoubleVector u0(c.n);

    for (int i=0; i<c.n; i++) u0[i] = 0.00001;
    /* Minimization */
    SteepestDescentGradient g1;
    g1.setFunction(&c);
    g1.setEpsilon(0.0000001);
    g1.setGradientStep(0.0000001);
    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
    g1.setX(u0);
    g1.setPrinter(new CFunction1Printer);
    g1.calculate();

    puts("-----------------------------------------------------------------");
    for (int i=0; i<c.n; i++) u0[i] = 0.00001;
    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&c);
    g2.setEpsilon(0.0000001);
    g2.setGradientStep(0.0000001);
    g2.setR1MinimizeEpsilon(0.01, 0.0000001);
    g2.setX(u0);
    g2.setPrinter(new CFunction1Printer);
    g2.calculate();
}

void CFunction1Printer::print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const
{
    printf("J[%2d]: %.10f  ", iterationCount, f->fx(m_x));
    print("u", m_x);
}

void CFunction1Printer::print(const char* s, const std::vector<double>& x) const
{
    unsigned int i;
    unsigned int n = x.size();
    printf("%s: ", s);
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

