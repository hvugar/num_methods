#include "sampleexample1.h"

#define arguments1 \
    double a111 = x[0]; \
    double a112 = x[1]; \
    double a113 = x[2]; \
    double b111 = x[3]; \
    double b112 = x[4]; \
    double b113 = x[5]; \
    double b211 = x[6]; \
    double b212 = x[7]; \
    double b213 = x[8]; \
    double q1   = x[9]; \
    double M    = x[10];

#define arguments2 \
    double a121 = x[0]; \
    double a122 = x[1]; \
    double a123 = x[2];

#define arguments3 \
    double a131 = x[0]; \
    double a132 = x[1]; \
    double a133 = x[2];

double _A11(double t) { return 0.0; }
double _A12(double t) { return t; }
double _A13(double t) { return 2.0*t; }
double _A21(double t) { return 3.0*t; }
double _A22(double t) { return 0.0; }
double _A23(double t) { return -1.0; }
double _A31(double t) { return 1.0; }
double _A32(double t) { return 2.0; }
double _A33(double t) { return 0.0; }

double _B111(double t) { return 0.0; }
double _B112(double t) { return 1.0; }
double _B113(double t) { return 0.0; }
double _B121(double t) { return 0.0; }
double _B122(double t) { return 0.0; }
double _B123(double t) { return 0.0; }
double _B131(double t) { return 0.0; }
double _B132(double t) { return -1.0; }
double _B133(double t) { return 0.0; }

double _B211(double t) { return 0.0; }
double _B212(double t) { return 0.0; }
double _B213(double t) { return 0.0; }
double _B221(double t) { return 0.0; }
double _B222(double t) { return 0.0; }
double _B223(double t) { return 1.0; }
double _B231(double t) { return 0.0; }
double _B232(double t) { return 0.0; }
double _B233(double t) { return 0.0; }

double _C1(double t) { return -t*t*t - 6.0*t*t + 2.0*t*(cos(t)+sin(t)-1.0) + cos(t) + 3.87532; }
double _C2(double t) { return -6.0*t*t + t*(5.0-sin(t)) + sin(t) - 1.56836; }
double _C3(double t) { return -2.0*t*t - 2.0*t + 3.0*cos(t) - sin(t) + 1.12468; }

double R01(double t, double *x, unsigned int n) { arguments1 return (a111*a111 + a112*a112 + a113*a113) + (b111*b111 + b112*b112 + b113*b113) + (b211*b211 + b212*b212 + b213*b213); }

double S01(double t, double *x, unsigned int n)
{
    arguments1
    return ((a111*(a111*_A11(t)+a112*_A21(t)+a113*_A31(t)) + a112*(a111*_A12(t)+a112*_A22(t)+a113*_A32(t)) + a113*(a111*_A13(t)+a112*_A23(t)+a113*_A33(t)))
            + (a111*(b111*_B111(t)+b112*_B112(t)+a113*_B113(t)) + a112*(b111*_B121(t)+b112*_B122(t)+b113*_B123(t)) + a113*(b111*_B131(t)+b112*_B132(t)+b113*_B133(t)))
            + (a111*(b111*_B211(t)+b112*_B212(t)+a113*_B213(t)) + a112*(b111*_B221(t)+b112*_B222(t)+b113*_B223(t)) + a113*(b111*_B231(t)+b112*_B232(t)+b113*_B233(t)))
            - (a111 * _C1(t) + a112 * _C2(t) + a113 * _C3(t))) / R01(t, x, n);
}

double R02(double t, double *x, unsigned int n)
{
    return 0.0;
}

double S02(double t, double *x, unsigned int n)
{
    return 0.0;
}

double R03(double t, double *x, unsigned int n)
{
    return 0.0;
}

double S03(double t, double *x, unsigned int n)
{
    return 0.0;
}

double alpha111(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * a111 - (a111*_A11(t) + a112*_A21(t) + a113*_A31(t)); }
double alpha112(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * a112 - (a111*_A12(t) + a112*_A22(t) + a113*_A32(t)); }
double alpha113(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * a113 - (a111*_A13(t) + a112*_A23(t) + a113*_A33(t)); }

double betta111(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b111 - (a111*_B111(t) + a112*_B121(t) + a113*_B131(t)); }
double betta112(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b112 - (a111*_B112(t) + a112*_B122(t) + a113*_B132(t)); }
double betta113(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b113 - (a111*_B113(t) + a112*_B123(t) + a113*_B133(t)); }

double betta211(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b211 - (a111*_B211(t) + a112*_B121(t) + a113*_B231(t)); }
double betta212(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b212 - (a111*_B212(t) + a112*_B122(t) + a113*_B232(t)); }
double betta213(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b213 - (a111*_B213(t) + a112*_B123(t) + a113*_B233(t)); }

double qamma1(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * q1 - (a111*_C1(t) + a112*_C2(t) + a113*_C3(t)); }

double Mt(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * M; }

//double alpha121(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * a121 - (a121*_A11(t) + a122*_A21(t) + a123*_A31(t)); }
//double alpha122(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * a122 - (a121*_A12(t) + a122*_A22(t) + a123*_A32(t)); }
//double alpha123(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * a123 - (a121*_A13(t) + a122*_A23(t) + a123*_A33(t)); }

//double alpha131(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * a131 - (a131*_A11(t) + a132*_A21(t) + a133*_A31(t)); }
//double alpha132(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * a132 - (a131*_A12(t) + a132*_A22(t) + a133*_A32(t)); }
//double alpha133(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * a133 - (a131*_A13(t) + a132*_A23(t) + a133*_A33(t)); }

void SampleMain()
{
    double t0 = 0.0;
    double t1 = 0.5;
    unsigned int N = 60;
    double h = (t1 - t0) / N;
    unsigned int n = 11;

    double *x0 = (double*) malloc(sizeof(double) * n);
    // alpha
    x0[0] = 2.0;
    x0[1] = 0.0;
    x0[2] = 0.0;
    // beta1
    x0[3] = 0.0;
    x0[4] = 0.0;
    x0[5] = 0.0;
    // beta2
    x0[6] = 0.0;
    x0[7] = 0.0;
    x0[8] = 0.0;
    // qamma1
    x0[9] = -1.3569;
    // M
    x0[10] = 1.0;

    double **x1 = (double**)malloc(sizeof(double*) * n);
    for (unsigned int i=0; i<n; i++) x1[i] = (double*)malloc(sizeof(double) * N+1);

    ODE1stOrderEquationN eqs[11];
    eqs[0] = &alpha111;
    eqs[1] = &alpha112;
    eqs[2] = &alpha113;

    eqs[3] = &betta111;
    eqs[4] = &betta112;
    eqs[5] = &betta113;

    eqs[6] = &betta211;
    eqs[7] = &betta212;
    eqs[8] = &betta213;

    eqs[9] = &qamma1;
    eqs[10] = &Mt;

    runge_kutta_rk4_system(t0, t1, x0, x1, n, N, h, eqs);

    IPrinter::printVector(x1[0], N+1, "alpha1_11", 4);
    IPrinter::printVector(x1[1], N+1, "alpha1_12", 4);
    IPrinter::printVector(x1[2], N+1, "alpha1_13", 4);

    IPrinter::printVector(x1[3], N+1, "betta1_11", 4);
    IPrinter::printVector(x1[4], N+1, "betta1_12", 4);
    IPrinter::printVector(x1[5], N+1, "betta1_13", 4);

    IPrinter::printVector(x1[6], N+1, "betta2_11", 4);
    IPrinter::printVector(x1[7], N+1, "betta2_12", 4);
    IPrinter::printVector(x1[8], N+1, "betta2_13", 4);

    IPrinter::printVector(x1[9], N+1, "qamma_1  ", 4);
    IPrinter::printVector(x1[10], N+1, "M        ", 4);
}
