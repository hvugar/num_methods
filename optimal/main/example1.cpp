#include "example1.h"

Example1::Example1()
{
    N = 4;
}

struct CauchyProblemA : public CauchyProblem
{
    virtual double f(double t, const DoubleVector &x) const
    {
//        double a01 = x[0];
//        double a02 = x[1];
//        double a03 = x[2];
//        double a04 = x[3];
//        double b   = x[4];
//        double M   = x[5];

        DoubleVector a = x.mid(0, 3);
        double s0 = e->S0(a, x[4], t);

        if (tp == 0) return s0 * x[i] - (e->A(0,i,t)*x[0] + e->A(1,i,t)*x[1] + e->A(2,i,t)*x[2] + e->A(3,i,t)*x[3]);
        if (tp == 1) return s0 * x[4] + (e->B(0,t)  *x[0] + e->B(1,t)  *x[1] + e->B(2,t)  *x[2] + e->B(3,t)  *x[3]);
        if (tp == 2) return s0 * x[5];
        return 0.0;
    }
    int tp;
    int i;
    Example1 *e;
};

void Example1::calculate()
{
    DoubleVector a(4);
    a.at(0) = +1.0;
    a.at(1) = +0.0;
    a.at(2) = -2.0;
    a.at(3) = +0.0;

    std::vector<CauchyProblem*> cps(N+2);
    for (unsigned int i=0; i<N; i++)
    {
        CauchyProblemA* cp = new CauchyProblemA;
        cp->i  = i;
        cp->tp = 0;
        cp->x0 = 0.0;
        cp->y0 = a[i];
        cp->e = this;
        cps[i] = cp;
    }

    CauchyProblemA* cp1 = new CauchyProblemA;
    cp1->i  = 0;
    cp1->tp = 1;
    cp1->x0 = 0.0;
    cp1->y0 = 10.7551;
    cp1->e = this;
    cps[N] = cp1;

    CauchyProblemA* cp2 = new CauchyProblemA;
    cp2->i  = 0;
    cp2->tp = 2;
    cp2->x0 = 0.0;
    cp2->y0 = 1.0;
    cp2->e = this;
    cps[N+1] = cp2;

    double h = 0.025;
    DoubleMatrix m;
    CauchyProblem::rungeKutta(cps, 0.0, h, 80, m);

    IPrinter::printVector(m.row(0), NULL, 8);
    IPrinter::printVector(m.row(1), NULL, 8);
    IPrinter::printVector(m.row(2), NULL, 8);
    IPrinter::printVector(m.row(3), NULL, 8);
    IPrinter::printVector(m.row(4), NULL, 8);
    IPrinter::printVector(m.row(5), NULL, 8);
}

double Example1::A(unsigned int i, unsigned int j, double t UNUSED_PARAM) const
{
    double data[N][N] =
    {
        {+2.0, -1.0, +1.0, +0.0},
        {+0.0, +0.0, +1.0, +2.0},
        {+1.0, -3.0, +0.0, +1.0},
        {+0.0, +2.0, -3.0, +0.0}
    };
    return data[i][j];
}

double Example1::B(unsigned int i, double t) const
{
    if (i==0) return -6.0*t*t + 16.0*t + 3.0*cos(t) - sin(t) - 12.0;
    if (i==1) return -3.0*t*t + cos(t) + 2.0*sin(t) - 2.0;
    if (i==2) return t*t + 12.0*t + 5.0*sin(t) - cos(t) - 11.0;
    if (i==3) return t*t - 8.0*t - 3.0*sin(t) - 3.0*cos(t) + 7.0;
    return 0.0;
}

double Example1::S0(const DoubleVector &a, double b, double t) const
{
    double c[4];
    c[0] = a.at(0)*A(0,0,t) + a.at(1)*A(1,0,t) + a.at(2)*A(2,0,t) + a.at(3)*A(3,0,t);
    c[1] = a.at(0)*A(0,1,t) + a.at(1)*A(1,1,t) + a.at(2)*A(2,1,t) + a.at(3)*A(3,1,t);
    c[2] = a.at(0)*A(0,2,t) + a.at(1)*A(1,2,t) + a.at(2)*A(2,2,t) + a.at(3)*A(3,2,t);
    c[3] = a.at(0)*A(0,3,t) + a.at(1)*A(1,3,t) + a.at(2)*A(2,3,t) + a.at(3)*A(3,3,t);
    double sum1 = c[0]*a.at(0) + c[1]*a.at(1) + c[2]*a.at(2) + c[3]*a.at(3);

    double sum2 = a.at(0)*B(0,t) + a.at(1)*B(1,t) + a.at(2)*B(2,t) + a.at(3)*B(3,t);

    double sum3 = a.at(0)*a.at(0) + a.at(1)*a.at(1) + a.at(2)*a.at(2) + a.at(3)*a.at(3) + b*b;

    return (sum1-sum2*b)/sum3;
}


