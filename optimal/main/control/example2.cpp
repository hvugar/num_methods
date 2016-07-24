#include "example2.h"
#include <math.h>

void Example2::main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    Example2 e2;
    C_UNUSED(e2);
}

Example2::Example2()
{
    t0 = 0.0;
    t1 = 1.0;
    h = 0.001;
    M = 1000;

    DoubleVector k0(8);

    x0 = 0.0;
    xT = 18.42060880;

    k0[0] = +40.0;
    k0[1] = +30.0;
    k0[2] = +50.0;
    k0[3] = +60.0;
    k0[4] = +80.0;
    k0[5] = +10.0;
    k0[6] = +5.0;
    k0[7] = +1.0;
    pK = &k0;

    //    ConjugateGradient cg;
    //    cg.setFunction(this);
    //    cg.setGradient(this);
    //    cg.setPrinter(this);
    //    cg.setProjection(NULL);
    //    cg.calculate(k0);

    double xT = 19.2116600200;

    DoubleVector x(M+1);
    calculateX(x);
}

double Example2::fx(double t, double x) const
{
    double K = getK(x0);
    return 3.0*t*x + K*x0 + (2.0*t - 3.0*t*t*t);
}

double Example2::getK(double x) const
{
    //    if (DBL_MIN <= x && x < -3.0000) (*pK)[0] = +40.0;
    //    if (-3.0000 <= x && x < -0.1000) (*pK)[1] = +30.0;
    //    if (-0.1000 <= x && x < +0.0050) (*pK)[2] = +50.0;
    //    if (+0.0050 <= x && x < +0.0300) (*pK)[3] = +60.0;
    //    if (+0.0300 <= x && x < +1.0000) (*pK)[4] = +80.0;
    //    if (+1.0000 <= x && x < +2.0000) (*pK)[5] = +10.0;
    //    if (+2.0000 <= x && x < +3.0000) (*pK)[6] = +5.0;
    //    if (+3.0000 <= x && x < DBL_MAX) (*pK)[7] = +1.0;

    double K = 0.0;
    if (DBL_MIN <= x && x < -3.0000) K = (*pK)[0];
    if (-3.0000 <= x && x < -0.1000) K = (*pK)[1];
    if (-0.1000 <= x && x < +0.0050) K = (*pK)[2];
    if (+0.0050 <= x && x < +0.0300) K = (*pK)[3];
    if (+0.0300 <= x && x < +1.0000) K = (*pK)[4];
    if (+1.0000 <= x && x < +2.0000) K = (*pK)[5];
    if (+2.0000 <= x && x < +3.0000) K = (*pK)[6];
    if (+3.0000 <= x && x < DBL_MAX) K = (*pK)[7];
    return K;
}

void Example2::calculateX(DoubleVector &x)
{
    double _x0 = x0;
    for (unsigned int i=0; i<10; i++)
    {
        t0 = (i*100)*h;
        double K = getK(_x0);
        printf("%.1f %.1f x0:%5.2f K:%5.2f ", t0, t0+0.1, _x0, K);

        DoubleVector _x;
        runge_kutta(this, t0, _x0, 100, _x, h);

        IPrinter::printVector(_x, NULL);
        for (unsigned int j=0; j<=100; j++) x[i*100+j] = _x[j];
        _x0 = _x[100];
    }

    IPrinter::printVector(x, "x:");
    IPrinter::printVector(x.data(), x.size(), "", x.size(), 0, 0, "data.txt");
}

double Example2::fx(const DoubleVector &x)
{
    return x[x.size()-1];
}

void Example2::gradient(const DoubleVector &K, DoubleVector &g)
{
    pK = &K;
    DoubleVector x(M+1);
    calculateX(x);
    px = &x;

    struct ODEPSi : public OrdDifEquation {

        Example2 *e;
        DoubleVector *psi;

        virtual double fx(double t, double psi) const {
            double sum = 0.0;
            unsigned int j = (unsigned int)(round(t/e->h));
            for (unsigned int i=0; i<1000; i+=100)
            {
                if (i == j)
                {
                    double int_sum = 0.0;
                    for (unsigned int p=i+1; p<=100; p++)
                    {
                        int_sum += (*this->psi)[p] * e->getK((*e->px)[p]);
                    }
                    int_sum *= e->h;
                    sum += int_sum;
                }
            }

            return -(3.0*t)*psi-sum;
        }
    };
    DoubleVector psi(M+1);

    ODEPSi ode;
    ode.e = this;
    ode.psi = &psi;

    psi[M] = -2.0*x[M];
    double dx = -h;
    double t1 = 1.0;
    for (int i=M; i>=1; i--)
    {
        double _t = t1 - (M-i)*h;
        double _x = x[i];

        double k1 = ode.fx(_t, _x);
        double k2 = ode.fx(_t+dx/2.0, _x+(dx/2.0)*k1);
        double k3 = ode.fx(_t+dx/2.0, _x+(dx/2.0)*k2);
        double k4 = ode.fx(_t+dx, _x+dx*k3);
        psi[i-1] = psi[i] + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
    }

    IPrinter::printVector(psi, "psi ");

    // s=3;
    for (unsigned int i=0; i<=M; i++)
    {
    }
}

void Example2::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double a, RnFunction *fn) const
{

}


void Example2::runge_kutta(OrdDifEquation *ode, double x0, double y0, unsigned int N, DoubleVector &y, double h)
{
    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    y.clear();
    y.resize(N+1,0.0);

    y[0] = y0;
    for (unsigned int i=0; i<N; i++)
    {
        double _x = x0 + i*h;
        double _y = y[i];

        k1 = ode->fx(_x, _y);
        k2 = ode->fx(_x+h/2.0, _y+(h/2.0)*k1);
        k3 = ode->fx(_x+h/2.0, _y+(h/2.0)*k2);
        k4 = ode->fx(_x+h, _y+h*k3);
        y[i+1] = y[i] + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);;
    }
}
