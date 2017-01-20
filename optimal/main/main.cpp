#include "problem1.h"
#include "loadedsystems.h"
#include "example1.h"
#include "example2.h"
#include "example3.h"
#include "example4.h"
#include "example5.h"
#include "problem1.h"
#include <cmethods.h>

#include <../border/borderparabolicd.h>
#include <../border/borderparabolicn.h>
#include <../border/borderparabolic2d.h>
#include <../border/borderhyperbolic2d.h>

#include "bordertest1.h"
#include "sampleboundaryproblem1.h"

#include <float.h>
#include <time.h>

#include "bordertest.h"

#include <grid/parabolicequationgird1d.h>

class A : public ParabolicEquationGird1D
{
public:
    inline double U(unsigned int n, unsigned int m) const;

protected:
    virtual double initial(unsigned int n) const;
    virtual double boundary(unsigned int m, Boundary boundary) const;
    virtual double f(unsigned int n, unsigned int m) const;
    virtual double a(unsigned int n, unsigned int m) const;
};

inline double A::U(unsigned int n, unsigned int m) const
{
    double x = n*setting.hx;
    double t = m*setting.hx;
    return x*x + t;
}

double A::initial(unsigned int n) const
{
    return U(n, 0);
}

double A::boundary(unsigned int m, Boundary boundary) const
{
    unsigned int N = setting.N;
    if (boundary == Left) return U(0,m);
    if (boundary == Right) return U(N,m);
    return 0.0;
}

double A::f(unsigned int n, unsigned int m) const
{
    return 1.0 - 2.0*a(n,m);
}

double A::a(unsigned int n UNUSED_PARAM, unsigned int m UNUSED_PARAM) const
{
    return 1.0;
}

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    srand(time(NULL));

//    A a;
//    DoubleMatrix m;
//    a.gridMethod(m);
//    IPrinter::printMatrix(14, 10, m);
//    IPrinter::printSeperatorLine();

//    DoubleMatrix u1;
//    a.calculateN4L2RD(u1);
//    IPrinter::printMatrix(14, 10, u1);

    //BorderTest::Main(argc, argv);
    //BorderTest1::Main(argc, argv);
    //Example2::Main(argc, argv);
    //Problem1::Main(argc, argv);
    //Example4::Main(argc, argv);
    BorderParabolicD::Main(argc, argv);
    //BorderParabolicN::Main(argc, argv);
    //BorderParabolic2D::Main(argc, argv);
    //BorderHyperbolic2D::Main(argc, argv);
    //Example5 e5;
    //BoundaryValueProblem1::Main(argc, argv);
    return 0;
}
