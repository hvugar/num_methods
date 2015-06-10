#include "cfunction2.h"

CFunction2::CFunction2(double t0, double t1, double h)
{
    this->t0 = t0;
    this->t1 = t1;
    this->h = h;
    n = (int)ceil((t1-t0)/h) + 1;

    t = std::vector<double>(n);
    x1 = std::vector<double>(n);
    x2 = std::vector<double>(n);
    psi1 = std::vector<double>(n);
    psi2 = std::vector<double>(n);

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

double CFunction2::fx(const std::vector<double>& u)
{
    double sum = 0.0;
    int i=0;
    for (i=0; i<(n-1); i++)
    {
        int j=i+1;
        double fj = fx0(t[j], x1[j], x2[j], u[j]);//(x1[j]-t[j]*t[j]*t[j])*(x1[j]-t[j]*t[j]*t[j]) + (x2[j] - t[j])*(x2[j] - t[j]) + (2*u[j] - t[j])*(2*u[j] - t[j]);
        double fi = fx0(t[i], x1[i], x2[i], u[i]);//(x1[i]-t[i]*t[i]*t[i])*(x1[i]-t[i]*t[i]*t[i]) + (x2[i] - t[i])*(x2[i] - t[i]) + (2*u[i] - t[i])*(2*u[i] - t[i]);
        sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
    }
    return sum;
}

void CFunction2::gradientJ(double grad_step, std::vector<double>& g, const std::vector<double>& u)
{
    std::vector<double> x(2);
    std::vector<double> psi(2);

    test(u);

    printX("t", t);
    printX("u", u);
    printX("x1", x1);
    printX("x2", x2);
    printX("p1", psi1);
    printX("p2", psi2);

    for (int i=0; i<n; i++)
    {
        x[0] = x1[i];
        x[1] = x2[i];

        psi[0] = psi1[i];
        psi[1] = psi2[i];

//        double u1 = u[i] + grad_step;
//        double u2 = u[i] - grad_step;
//        g[i] = (H(t[i], x1[i], x2[i], u1, psi1[i], psi2[i]) - H(t[i], x1[i], x2[i], u2, psi1[i], psi2[i])) / (2 * grad_step);
        g[i] = -psi2[i];
    }
    printX("gr", g);
}

double CFunction2::fx0(double t, double x1, double x2, double u)
{
    return (x1-(t*t)/2.0)*(x1-(t*t)/2.0) + (x2 - t)*(x2 - t);
}

double CFunction2::fx1(double t, double x1, double x2, double u)
{
    return x2;
}

double CFunction2::fx2(double t, double x1, double x2, double u)
{
    return -6.0*x1 + x2 + u - t + 1.0;
}

double CFunction2::fp1(double t, double x1, double x2, double p1, double p2, double u)
{
    return 2.0*x1 - t*t + 6.0*p2;
}

double CFunction2::fp2(double t, double x1, double x2, double p1, double p2, double u)
{
    return 2.0*(x2-t) - p1 - p2;
}

double CFunction2::H(double t, double x1, double x2, double u, double p1, double p2)
{
    return -fx0(t, x1, x2, u) + p1 * fx1(t, x1, x2, u) + p2 * fx2(t, x1, x2, u);
}

void CFunction2::test(const std::vector<double>& u)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    x1[0] = 0.0;
    x2[0] = 0.0;
    //std::vector<double> _x(2);
    //_x[0] = x1[0];
    //_x[1] = x2[0];
    double _x1, _x2;
    h = +fabs(h);
    for (int i=0; i<n-1; i++)
    {
        _x1 = x1[i];
        _x2 = x2[i];
        k1[0] = fx1(t[i], _x1, _x2, u[i]);
        k1[1] = fx2(t[i], _x1, _x2, u[i]);

        _x1 = x1[i] + (h/2.0) * k1[0];
        _x2 = x2[i] + (h/2.0) * k1[1];
        k2[0] = fx1(t[i]+h/2.0, _x1, _x2, u[i]);
        k2[1] = fx2(t[i]+h/2.0, _x1, _x2, u[i]);

        _x1 = x1[i] + (h/2.0) * k2[0];
        _x2 = x2[i] + (h/2.0) * k2[1];
        k3[0] = fx1(t[i]+h/2.0, _x1, _x2, u[i]);
        k3[1] = fx2(t[i]+h/2.0, _x1, _x2, u[i]);

        _x1 = x1[i] + h * k3[0];
        _x2 = x2[i] + h * k3[1];
        k4[0] = fx1(t[i]+h, _x1, _x2, u[i]);
        k4[1] = fx2(t[i]+h, _x1, _x2, u[i]);

        x1[i+1] = x1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        x2[i+1] = x2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }

    psi1[n-1] = 0.0;
    psi2[n-1] = 0.0;
    //std::vector<double> _psi(2);
    //_psi[0] = psi1[n-1];
    //_psi[1] = psi2[n-1];
    double _psi1, _psi2;
    h = -fabs(h);
    for (int i=n-1; i>0; i--)
    {
        //_x[0] = x1[i];
        //_x[1] = x2[i];
        //_psi[0] = psi1[i];
        //_psi[1] = psi2[i];
        _x1 = x1[i];
        _x2 = x2[i];
        _psi1 = psi1[i];
        _psi2 = psi2[i];

        k1[0] = fp1(t[i], x1[i], x2[i], psi1[i], psi2[i], u[i]);
        k1[1] = fp2(t[i], x1[i], x2[i], psi1[i], psi2[i], u[i]);

        _psi1 = psi1[i] + (h/2.0) * k1[0];
        _psi2 = psi2[i] + (h/2.0) * k1[1];
        k2[0] = fp1(t[i]+h/2.0, x1[i], x2[i], _psi1, _psi2, u[i]);
        k2[1] = fp2(t[i]+h/2.0, x1[i], x2[i], _psi1, _psi2, u[i]);

        _psi1 = psi1[i] + (h/2.0) * k2[0];
        _psi2 = psi2[i] + (h/2.0) * k2[1];
        k3[0] = fp1(t[i]+h/2.0, x1[i], x2[i], _psi1, _psi2, u[i]);
        k3[1] = fp2(t[i]+h/2.0, x1[i], x2[i], _psi1, _psi2, u[i]);

        _psi1 = psi1[i] + h * k3[0];
        _psi2 = psi2[i] + h * k3[1];
        k4[0] = fp1(t[i]+h, x1[i], x2[i], _psi1, _psi2, u[i]);
        k4[1] = fp2(t[i]+h, x1[i], x2[i], _psi1, _psi2, u[i]);

        psi1[i-1] = psi1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        psi2[i-1] = psi2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
}
