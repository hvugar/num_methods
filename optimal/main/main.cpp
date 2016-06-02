#include "headers.h"

////double X(double t) { return t*t + 1.0; }
////double fx(double t) { return 2.0*t; }
////double fa(double t, double a) { return 4.0*a + (2.0*t-4.0*t*t-7.75); }
////double fb(double t, double b) { return 4.0*b + 3.0; }

<<<<<<< .mine
//double fx(double t) { return 3.0*(1.35*t-0.3)*(1.35*t-0.3) - 0.8*(t-0.2)*(t-0.2)*(t-0.2) - 2.5*(t-0.05)*(t-0.05)*(t-0.05)*(t-0.05) +0.05; }
//double dx(double t) { return 2.0*1.35*3.0*(1.35*t-0.3) - 3.0*0.8*(t-0.2)*(t-0.2) - 4.0*2.5*(t-0.05)*(t-0.05)*(t-0.05); }
//double sum(double ht) {
//    double s = 0.0;
//    for (int i=1; i<10; i++) s += fx(i*ht);
//    return s;
//}
=======
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
>>>>>>> .r671

double fx(double t) { return t*t + 1.0; }
double fa(double t, double a) { return 3.0 * a + (2.0*t - 3.0*t*t - 14.61); }
double fb1(double t, double b1) { return 3.0*b1 + 4.0; }
double fb2(double t, double b2) { return 3.0*b2 + 5.0; }

double ode1(double x, double y) { return x; }

int main(int argc, char ** argv)
{
    //BorderParabolic::main(argc, argv);
    //BorderParabolic2D::main(argc, argv);
    //BorderHyperbolic::main();

<<<<<<< .mine
    SampleBorderHyperBolic::main(argc, argv);
    return 0;

    for (int i=0; i<1; i++)
    {
        printf("OK\n");
    }

    unsigned int N1 = 100;
    double *px = (double*)malloc(sizeof(double)*(N1+1));
    double *py = (double*)malloc(sizeof(double)*(N1+1));
    runge_kutta_rk4(0.0, 0.0, 1.0, 0.0, N1, px, py, ode1);
    printf("%9.6f\n", px[0]);
    for (unsigned int i=0; i<=N1; i++) if (i%10==0) printf("%10.6f", px[i]);
    printf("\n");
    for (unsigned int i=0; i<=N1; i++) if (i%10==0) printf("%10.6f", py[i]);
    printf("\n");
    free(py);
    free(px);

    return 0;

    printf("%f %f\n", fx(0.0), fx(1.0));
=======
    //printf("%f %f\n", f(0.0), f(1.0));
>>>>>>> .r671

<<<<<<< .mine
    double dx = 0.01;
    unsigned int N = 101;
    unsigned int P1 = 20;
    unsigned int P2 = 70;
=======
    double dx = 0.001;
    unsigned int N = 1001;
//    unsigned int P = 10000;
>>>>>>> .r671

<<<<<<< .mine
    DoubleVector a(N);
    RungeKutta::calculate(fa, 0.0, 1.0, a, dx);
    DoubleVector b1(N);
    RungeKutta::calculate(fb1, 0.0, 0.0, b1, dx);
    DoubleVector b2(N);
    RungeKutta::calculate(fb2, 0.0, 0.0, b2, dx);
=======
    DoubleVector a(N);
    RungeKutta::calculate(f, 0.0, 1.0, a, dx);
//    DoubleVector b(N);
//    RungeKutta::calculate(fb, 0.0, 0.0, b, dx);
>>>>>>> .r671

<<<<<<< .mine
    IPrinter::printVector(a, "a: ");
    IPrinter::printVector(b1, "b1:");
    IPrinter::printVector(b2, "b2:");
=======
    IPrinter::printVector(a, "a:");
//    IPrinter::printVector(b, "b:");
>>>>>>> .r671

<<<<<<< .mine
    double _a1 = 1.0 - b1[P1];
    double _b1 = 0.0 - b2[P1];
    double _c1 = a[P1];
=======
//    DoubleVector x(N);
    DoubleVector X(N);
>>>>>>> .r671

    double _a2 = 0.0 - b1[P2];
    double _b2 = 1.0 - b2[P2];
    double _c2 = a[P2];

<<<<<<< .mine
    DoubleVector x(N);
    DoubleVector X(N);
=======
    for (unsigned int i=0; i<N; i++)
    {
        //x[i] = a[i] + b[i] * x[P];
        X[i] = exp(i*dx);
    }
//    IPrinter::printVector(x, "x:");
    IPrinter::printVector(X, "X:");
>>>>>>> .r671

    x[P2] = (_c1*_a2 - _c2*_a1)/(_b1*_a2 - _b2*_a1);
    x[P1] = (_c1*_b2 - _c2*_b1)/(_b2*_a1 - _b1*_a2);

    for (unsigned int i=0; i<N; i++)
    {
        x[i] = a[i] + b1[i] * x[P1] + b2[i] * x[P2];
        X[i] = fx(i*dx);
    }
    IPrinter::printVector(x, "x: ");
    IPrinter::printVector(X, "X: ");

    return 0;
}
