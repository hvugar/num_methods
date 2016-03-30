#include "borderparabolic2d.h"

#define SAMPLE2

void BorderParabolic2D::main()
{
    FILE* file = fopen("file1.txt", "w");
    BorderParabolic2D bp;
    DoubleMatrix m;
    bp.caluclateMVD1(m, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M);
    IPrinter::printMatrix(m, bp.N2, bp.N1, NULL, file);
    IPrinter::printMatrix(m, 10, 10, NULL, stdout);
    fclose(file);
}

BorderParabolic2D::BorderParabolic2D()
{
    a1 = a2 = 1.0;
    x10 = x20 = t0 = 0.0;
    x11 = x21 = t1 = 1.0;
    h1 = 0.01;
    h2 = 0.01;
    ht = 0.005;
    N1 = (unsigned int)(ceil(x11-x10)/h1);
    N2 = (unsigned int)(ceil(x21-x20)/h2);
    M  = (unsigned int)(ceil(t1-t0)/ht);
}

BorderParabolic2D::~BorderParabolic2D()
{}

double BorderParabolic2D::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;

    return x1*x1 + x2*x2 + t*t;
}

double BorderParabolic2D::initial(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);

#ifdef SAMPLE1
    return u(i, j, 0);
#endif
#ifdef SAMPLE2
    return 4.0;
#endif
}

double BorderParabolic2D::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    double t = k*ht;
#ifdef SAMPLE1
    return u(i, j, k);
#endif
#ifdef SAMPLE2
    return initial(i,j)+3.0*t;
#endif
}

double BorderParabolic2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;

    static double sgm1 = 3.0*h1;
    static double sgm2 = 3.0*h2;
    static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    static double gause_b = 2.0*sgm1*sgm2;

#ifdef SAMPLE1
    return 2.0*t - 2.0*a1 - 2.0*a2;
#endif
#ifdef SAMPLE2
    double sum = 0.0;
    sum += (10*t) * gause_a * exp(-((x1-0.5)*(x1-0.5) + (x2-0.5)*(x2-0.5))/gause_b);
    //if (i==50 && j==50) return (100*t)/(h1*h2);
#endif
    return sum;
}
