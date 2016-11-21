#include "bordertest1.h"

void BorderTest1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderTest1 bt;
    bt.alpha = bt.X(0);
    bt.betta = bt.X(bt.N);
    DoubleVector x1;
    bt.calculate1(x1);
    IPrinter::printVector(14,10,x1);
    DoubleVector x2;
    bt.calculate2(x2);
    IPrinter::printVector(14,10,x2);
}

BorderTest1::BorderTest1()
{
}

double BorderTest1::a(unsigned int i) const
{
    double t = i*h;
#ifdef SAMPLE_1
    return 2.0*t;
#endif
#ifdef SAMPLE_2
    return -3.0;
#endif
}

double BorderTest1::b(unsigned int i) const
{
    double t = i*h;
#ifdef SAMPLE_1
    return 6.0*t - 2.0*t*t*t*t;
#endif
#ifdef SAMPLE_2
    return 4.0*exp(t)*cos(2.0*t);
#endif
}

double BorderTest1::X(unsigned int i) const
{
    double t = i*h;
#ifdef SAMPLE_1
    return t*t*t;
#endif
#ifdef SAMPLE_2
    return exp(t)*sin(2.0*t);
#endif
}

void BorderTest1::calculate1(DoubleVector &u)
{
    u.clear();
    u.resize(N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector dd(N-1);
    DoubleVector rx(N-1);

    for (unsigned int i=1; i<=N-1; i++)
    {
        da[i-1] = 1.0;
        db[i-1] = (-2.0-h*h*a(i));
        dc[i-1] = 1.0;
        dd[i-1] = h*h*b(i);
    }

    da[0]   = 0.0;
    dc[N-2] = 0.0;

    u.at(0) = X(0);
    u.at(N) = X(N);

    dd.at(0)   -= u.at(0);
    dd.at(N-2) -= u.at(N);

    tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

    for (unsigned int i=1; i<=N-1; i++)
    {
        u.at(i) = rx.at(i-1);
    }


    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();

}

void BorderTest1::calculate2(DoubleVector &u)
{
    u.clear();
    u.resize(N+1);

    u.at(0) = X(0);
    u.at(N) = X(N);

    DoubleVector betta(N-1, 0.0);
    betta.at(0) = (-2.0-h*h*a(1));
    betta.at(1) = 1.0;
    double qamma = h*h*b(1) - u.at(0);

    //    DoubleVector alpha1(N-1);
    //    DoubleVector alpha2(N-1);
    //    DoubleVector alpha3(N-1);

    //    for (unsigned int i=0; i<=N-2; i++)
    //    {
    //        alpha1.at(i) = 2.0 + a(i)*h*h;
    //        alpha2.at(i) = -1.0;
    //        alpha3.at(i) = b(i)*h*h;
    //    }

    for (unsigned int i=0; i<N-3; i++)
    {
        betta.at(i+1) = betta.at(i+1) + betta.at(i)*(2.0+a(i+2)*h*h);
        betta.at(i+2) = betta.at(i+2) + betta.at(i)*(-1.0);
        qamma = qamma - betta.at(i)*(b(i+2)*h*h);
    }

    DoubleMatrix m(2,2);
    m.at(0,0) = betta.at(N-3);
    m.at(0,1) = betta.at(N-2);
    m.at(1,0) = 1.0;
    m.at(1,1) = (-2.0-a(N-1)*h*h);
    DoubleVector b1(2);
    b1.at(0) = qamma;
    b1.at(1) = b(N-1)*h*h - u.at(N);

    DoubleVector x1(2);
    GaussianElimination(m,b1,x1);

    //printf("%18.14f %18.14f %18.14f\n", betta.at(N-4), betta.at(N-3), qamma);
    //printf("%18.14f %18.14f\n", x1.at(0), x1.at(1));

    u.at(N-2) = x1.at(0);
    u.at(N-1) = x1.at(1);
    for (unsigned int i=N-2; i!=1; i--)
    {
        u.at(i-1) = (2.0 + a(i)*h*h)*u.at(i)
                + (-1.0)*u.at(i+1)
                + (b(i)*h*h);
    }
}


