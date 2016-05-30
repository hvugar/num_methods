#include "headers.h"

////double X(double t) { return t*t + 1.0; }
////double fx(double t) { return 2.0*t; }
////double fa(double t, double a) { return 4.0*a + (2.0*t-4.0*t*t-7.75); }
////double fb(double t, double b) { return 4.0*b + 3.0; }

////double fX(double t) { return t*t + sin(t); }
//////double fx(double t) { return 2.0*t; }
////double fa(double t, double a) { return 2.0*a + (2.0*t-2.0*t*t-2.0*sin(t)+cos(t)-0.3295); }
////double fb(double t, double b) { return 2.0*b + 3.0; }

//double fx(double t) { return 3.0*(1.35*t-0.3)*(1.35*t-0.3) - 0.8*(t-0.2)*(t-0.2)*(t-0.2) - 2.5*(t-0.05)*(t-0.05)*(t-0.05)*(t-0.05) +0.05; }
//double dx(double t) { return 2.0*1.35*3.0*(1.35*t-0.3) - 3.0*0.8*(t-0.2)*(t-0.2) - 4.0*2.5*(t-0.05)*(t-0.05)*(t-0.05); }
//double sum(double ht) {
//    double s = 0.0;
//    for (int i=1; i<10; i++) s += fx(i*ht);
//    return s;
//}

double f(double t, double x)
{
    return exp(t) + t*t*t*x + x*x*x - t*t*t*exp(t) - exp(3.0*t);
}

int main(int argc, char ** argv)
{
    //BorderParabolic::main(argc, argv);
    //BorderParabolic2D::main(argc, argv);
    //BorderHyperbolic::main();

    //printf("%f %f\n", f(0.0), f(1.0));

    double dx = 0.001;
    unsigned int N = 1001;
//    unsigned int P = 10000;

    DoubleVector a(N);
    RungeKutta::calculate(f, 0.0, 1.0, a, dx);
//    DoubleVector b(N);
//    RungeKutta::calculate(fb, 0.0, 0.0, b, dx);

    IPrinter::printVector(a, "a:");
//    IPrinter::printVector(b, "b:");

//    DoubleVector x(N);
    DoubleVector X(N);

//    x[P] = a[P]/(1.0-b[P]);

    for (unsigned int i=0; i<N; i++)
    {
        //x[i] = a[i] + b[i] * x[P];
        X[i] = exp(i*dx);
    }
//    IPrinter::printVector(x, "x:");
    IPrinter::printVector(X, "X:");

    return 0;
}
