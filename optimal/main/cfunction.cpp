#include "cfunction.h"

ControlFunction::ControlFunction(double t0, double t1, double h) : RnFunction()
{
    this->t0 = t0;
    this->t1 = t1;
    this->h = h;
    n = (int)ceil((t1-t0)/h) + 1;

    t = DoubleVector(n);
    //    x1 = DoubleVector(n);
    //    x2 = DoubleVector(n);
    //    psi1 = DoubleVector(n);
    //    psi2 = DoubleVector(n);

    for (int i=0; i<n; i++)
    {
        t[i] = i*h;
        //        x1[i] = 0.0;
        //        x2[i] = 0.0;
        //        psi1[i] = 0.0;
        //        psi2[i] = 0.0;
    }
}

ControlFunction::~ControlFunction()
{
}

double ControlFunction::fx(const DoubleVector &u)
{
    return Integral(u);
}

void ControlFunction::gradient(double gradient_step, const DoubleVector& u, DoubleVector &g)
{
    DoubleVector psi1(n);
    DoubleVector psi2(n);

    DoubleVector x1(n);
    DoubleVector x2(n);

    DoubleVector x(2);
    DoubleVector psi(2);
    calculate_x(u, x1, x2);
    calculate_psi(u, psi1, psi2, x1, x2);
    for (int i=0; i<n; i++)
    {
        x[0] = x1[i];
        x[1] = x2[i];

        psi[0] = psi1[i];
        psi[1] = psi2[i];

        double u1 = u[i] + gradient_step;
        double u2 = u[i] - gradient_step;
        g[i] = (H(t[i], x, u1, psi) - H(t[i], x, u2, psi)) / (2 * gradient_step);
    }
}

double ControlFunction::H(double t, const DoubleVector &x, double u, const DoubleVector &psi)
{
    return -1.0 * fx0->fx(t, x, u) + psi[0] * fx1->fx(t, x, u) + psi[1] * fx2->fx(t, x, u);
}

void ControlFunction::calculate_x(const DoubleVector& u, DoubleVector& x1, DoubleVector& x2)
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
        k1[0] = fx1->fx(t[i], _x, u[i]);
        k1[1] = fx2->fx(t[i], _x, u[i]);

        _x[0] = x1[i] + (h/2.0) * k1[0];
        _x[1] = x2[i] + (h/2.0) * k1[1];
        k2[0] = fx1->fx(t[i]+h/2.0, _x, u[i]);
        k2[1] = fx2->fx(t[i]+h/2.0, _x, u[i]);

        _x[0] = x1[i] + (h/2.0) * k2[0];
        _x[1] = x2[i] + (h/2.0) * k2[1];
        k3[0] = fx1->fx(t[i]+h/2.0, _x, u[i]);
        k3[1] = fx2->fx(t[i]+h/2.0, _x, u[i]);

        _x[0] = x1[i] + h * k3[0];
        _x[1] = x2[i] + h * k3[1];
        k4[0] = fx1->fx(t[i]+h, _x, u[i]);
        k4[1] = fx2->fx(t[i]+h, _x, u[i]);

        x1[i+1] = x1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        x2[i+1] = x2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
}

void ControlFunction::calculate_psi(const DoubleVector& u, DoubleVector& psi1, DoubleVector& psi2, DoubleVector& x1, DoubleVector& x2)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    DoubleVector _x(2);

    DoubleVector _x1(2);
    DoubleVector _x2(2);
    double dx = 0.000001;
    _x1[0] = x1[n-1] + dx;
    _x1[1] = x2[n-1];
    _x2[0] = x1[n-1] - dx;
    _x2[1] = x2[n-1];
    psi1[n-1] = -1.0*(T->fx(0.0, _x1) - T->fx(0.0, _x2))/(2.0*dx);
    _x1[0] = x1[n-1];
    _x1[1] = x2[n-1] + dx;
    _x2[0] = x1[n-1];
    _x2[1] = x2[n-1] - dx;
    psi2[n-1] = -1.0*(T->fx(0.0, _x1) - T->fx(0.0, _x2))/(2.0*dx);

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

        double h_grad = 0.000001;

        _x1[0] = x1[i] + h_grad;
        _x1[1] = x2[i];
        _x2[0] = x1[i] - h_grad;
        _x2[1] = x2[i];
        k1[0] = -1.0 * (H(t[i], _x1, u[i], _psi) - H(t[i], _x2, u[i], _psi)) / (2 * h);
        _x1[0] = x1[i];
        _x1[1] = x2[i] + h_grad;
        _x2[0] = x1[i];
        _x2[1] = x2[i] - h_grad;
        k1[1] = -1.0 * (H(t[i], _x1, u[i], _psi) - H(t[i], _x2, u[i], _psi)) / (2 * h);
        //k1[0] = fp1->fx(t[i], _x, _psi, u[i]);
        //k1[1] = fp2->fx(t[i], _x, _psi, u[i]);
        _psi[0] = psi1[i] + (h/2.0) * k1[0];
        _psi[1] = psi2[i] + (h/2.0) * k1[1];

        _x1[0] = x1[i] + h_grad;
        _x1[1] = x2[i];
        _x2[0] = x1[i] - h_grad;
        _x2[1] = x2[i];
        k2[0] = -1.0 * (H(t[i]+h/2.0, _x1, u[i], _psi) - H(t[i], _x2, u[i], _psi)) / (2 * h);
        _x1[0] = x1[i];
        _x1[1] = x2[i] + h_grad;
        _x2[0] = x1[i];
        _x2[1] = x2[i] - h_grad;
        k2[1] = -1.0 * (H(t[i]+h/2.0, _x1, u[i], _psi) - H(t[i], _x2, u[i], _psi)) / (2 * h);
        //k2[0] = fp1->fx(t[i]+h/2.0, _x, _psi, u[i]);
        //k2[1] = fp2->fx(t[i]+h/2.0, _x, _psi, u[i]);
        _psi[0] = psi1[i] + (h/2.0) * k2[0];
        _psi[1] = psi2[i] + (h/2.0) * k2[1];

        _x1[0] = x1[i] + h_grad;
        _x1[1] = x2[i];
        _x2[0] = x1[i] - h_grad;
        _x2[1] = x2[i];
        k3[0] = -1.0 * (H(t[i]+h/2.0, _x1, u[i], _psi) - H(t[i], _x2, u[i], _psi)) / (2 * h);
        _x1[0] = x1[i];
        _x1[1] = x2[i] + h_grad;
        _x2[0] = x1[i];
        _x2[1] = x2[i] - h_grad;
        k3[1] = -1.0 * (H(t[i]+h/2.0, _x1, u[i], _psi) - H(t[i], _x2, u[i], _psi)) / (2 * h);
        //k3[0] = fp1->fx(t[i]+h/2.0, _x, _psi, u[i]);
        //k3[1] = fp2->fx(t[i]+h/2.0, _x, _psi, u[i]);
        _psi[0] = psi1[i] + h * k3[0];
        _psi[1] = psi2[i] + h * k3[1];

        _x1[0] = x1[i] + h_grad;
        _x1[1] = x2[i];
        _x2[0] = x1[i] - h_grad;
        _x2[1] = x2[i];
        k4[0] = -1.0 * (H(t[i]+h, _x1, u[i], _psi) - H(t[i], _x2, u[i], _psi)) / (2 * h);
        _x1[0] = x1[i];
        _x1[1] = x2[i] + h_grad;
        _x2[0] = x1[i];
        _x2[1] = x2[i] - h_grad;
        k4[1] = -1.0 * (H(t[i]+h, _x1, u[i], _psi) - H(t[i], _x2, u[i], _psi)) / (2 * h);
        //k4[0] = fp1->fx(t[i]+h, _x, _psi, u[i]);
        //k4[1] = fp2->fx(t[i]+h, _x, _psi, u[i]);
        psi1[i-1] = psi1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        psi2[i-1] = psi2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
}

double ControlFunction::Integral(const DoubleVector &u)
{
    DoubleVector x1(n);
    DoubleVector x2(n);

    calculate_x(u, x1, x2);

    double sum = 0.0;
    int i=0;
    DoubleVector x(2);
    for (i=0; i<(n-1); i++)
    {
        int j=i+1;
        x[0] = x1[j];
        x[1] = x2[j];
        double fj = fx0->fx(t[j], x, u[j]);

        x[0] = x1[i];
        x[1] = x2[i];
        double fi = fx0->fx(t[i], x, u[i]);

        sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
    }

    x[0] = x1[n-1];
    x[1] = x2[n-1];
    sum = sum + T->fx(t[n-1], x);
    return sum;
}

void ControlFunction::main()
{
    struct t_fx0 : public CFunction { virtual double fx(double t, const DoubleVector &x, double u) { return (x[0] - t*t*t)*(x[0] - t*t*t)+(x[1]-t)*(x[1]-t)+(2*u-t)*(2*u-t); } };
    struct t_trm : public CFunction { virtual double fx(double t, const DoubleVector &x) { return (x[1] - 1.0) * (x[1] - 1.0); } };

    struct t_fx1 : public CFunction { virtual double fx(double t, const DoubleVector &x, double u) { return 3.0 * x[1] * x[1]; } };
    struct t_fx2 : public CFunction { virtual double fx(double t, const DoubleVector &x, double u) { return x[0] + x[1] - 2.0*u - t*t*t + 1.0; } };

    //struct t_px1 : public CFunction { virtual double fx(double t, const DoubleVector &x, const DoubleVector &psi, double u) { return 2.0 * (x[0] - t*t*t) - psi[1]; } };
    //struct t_px2 : public CFunction { virtual double fx(double t, const DoubleVector &x, const DoubleVector &psi, double u) { return 2.0 * (x[1] - t) - 6.0 * x[1] * psi[0] - psi[1]; } };

    ControlFunctionPrinter cfp;

    /* Function */
    ControlFunction c(0.0, 1.0, 0.001);
    c.fx0 = new t_fx0;
    c.T   = new t_trm;
    c.fx1 = new t_fx1;
    c.fx2 = new t_fx2;
    //c.fp1 = new t_px1;
    //c.fp2 = new t_px2;

    c.dx = new CFunction*[2];
    c.dx[0] = new t_fx1;
    c.dx[1] = new t_fx2;

    /* initial point */
    DoubleVector u0(c.n);

    for (int i=0; i<c.n; i++) u0[i] = 0.00001;
    /* Minimization */
    SteepestDescentGradient g1;
    g1.setFunction(&c);
    g1.setEpsilon(0.0000001);
    g1.setGradientStep(0.0000001);
    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
    g1.setPrinter(&cfp);
    g1.calculate(u0);

    puts("-----------------------------------------------------------------");
    for (int i=0; i<c.n; i++) u0[i] = 0.00001;
    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&c);
    g2.setEpsilon(0.0000001);
    g2.setGradientStep(0.0000001);
    g2.setR1MinimizeEpsilon(0.01, 0.0000001);
    g2.setPrinter(&cfp);
    g2.calculate(u0);
}

void ControlFunctionPrinter::print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const
{
    printf("J[%2d]: %.10f  ", iterationCount, f->fx(m_x));
    print("u", m_x);
}

void ControlFunctionPrinter::print(const char* s, const DoubleVector& x) const
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

