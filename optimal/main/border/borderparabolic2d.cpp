#include "borderparabolic2d.h"

#define SAMPLE3

void BorderParabolic2D::main(int argc, char **argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    BorderParabolic2D bp;
    DoubleMatrix m;
    bp.caluclateMVD(m, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, bp.a1, bp.a2);
    //IPrinter::printMatrix(m, 10, 10, NULL, stdout);

    //    FILE* file = fopen("20160516.txt", "w");
    //    IPrinter::printMatrix(m, bp.N2, bp.N1, NULL, file);
    //    fclose(file);
}

BorderParabolic2D::BorderParabolic2D()
{
    a1 = 0.38;
    a2 = 0.38;
    x10 = x20 = t0 = 0.0;
    x11 = x21 = 1.0;
    t1 = 1.0;

    h1 = 0.001;
    h2 = 0.001;
    ht = 0.0000005;

    N1 = 1000;//(unsigned int)(ceil((x11-x10)/h1));
    N2 = 1000;//(unsigned int)(ceil((x21-x20)/h2));
    M  = 2000000;//(unsigned int)(ceil((t1-t0)/ht));
}

double BorderParabolic2D::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = 0.5*k*ht;

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
#ifdef SAMPLE3
    //    if (i<=2 || i>=N1-2 || j<=2 || j>=N2-2) return 5.0;
    return 0.0;
#endif
}

double BorderParabolic2D::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    //    double t = 0.5*k*ht;

#ifdef SAMPLE1
    return u(i, j, k);
#endif
    //#ifdef SAMPLE2
    //    return initial(i,j)+3.0*t;
    //#endif
    //#ifdef SAMPLE3
    //    return 1.0;
    //#endif

    return 1.0;
}

double BorderParabolic2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    //    double x1 = i*h1;
    //    double x2 = j*h2;
    double t = 0.5*k*ht;

    //    static double sgm1 = 3.0*h1;
    //    static double sgm2 = 3.0*h2;
    //    static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    //    static double gause_b = 2.0*sgm1*sgm2;

#ifdef SAMPLE1
    return 2.0*t - 2.0*a1 - 2.0*a2;
#endif

#ifdef SAMPLE2
    double sum = 0.0;
    sum += (10*t) * gause_a * exp(-((x1-0.5)*(x1-0.5) + (x2-0.5)*(x2-0.5))/gause_b);
    //    if (i==N1/2 && j==N2/2) return (100*t)/(h1*h2);
#endif

#ifdef SAMPLE3
    double sum = 0.0;
    //    double power = 2000.0;
    //    double delta1 = (1.0/h1);//(1.0/(sqrt(2.0*M_PI)*sgm1)) * exp(-((x1-k*ht)*(x1-k*ht))/(2.0*sgm1*sgm1));
    //    double delta2 = (1.0/h2);//(1.0/(sqrt(2.0*M_PI)*sgm2)) * exp(-((x1-j*h2)*(x1-j*h2))/(2.0*sgm2*sgm2));
    //    if ( (i==N1/2) && (j==(N2/2)) && (k <= 20) ) sum = power * delta1 * delta2;
    //    if ( k <= 20) sum = 10000.0;
#endif
    return 0.0;
}
