#include "headers.h"

#include "widget/qsimplewavewidget.h"
#include <matrix.h>
#include <time.h>

double fx1(double t, double *x, unsigned int n)
{
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    return x1 + 2.0*x2 + x3 - 2.0*sin(t) - cos(t) - t*t;
}

double fx2(double t, double *x, unsigned int n)
{
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    return 2.0*x1 + x2 + x3 - 2.0*t*t - t - sin(t) + 1.0;
}

double fx3(double t, double *x, unsigned int n)
{
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    return x1 - x3 + cos(t) - sin(t) - t*t;
}

int main(int argc, char ** argv)
{
    srand(time(NULL));

//    double y0[] = {0.0, 0.0, 1.0};
//    ODE1stOrderEquationN equations[3];
//    equations[0] = fx1;
//    equations[1] = fx2;
//    equations[2] = fx3;

//    unsigned int n = 3;
//    unsigned int N = 1000;
//    double h = 0.001;
//    double **y1 = (double **)malloc(sizeof(double*)*n);

//    for (unsigned int j=0; j<3; j++) y1[j] = (double*)malloc(sizeof(double)*(N+1));
//    runge_kutta_rk4_system(0.0, 1.0, y0, y1, 3, N, 0.001, equations);

//    DoubleVector v(N+1);
//    for (unsigned int i=0; i<=N; i++) v[i] = y1[2][i];
//    IPrinter::printVector(v, "x1");
    SampleLoaderBorder::main();

//    std::vector<RnFunction*> fs(3);
//    RnFunctionX1 f1;
//    RnFunctionX2 f2;
//    RnFunctionX3 f3;
//    fs[0] = &f1;
//    fs[1] = &f2;
//    fs[2] = &f3;

//    DoubleVector x0(3);
//    x0[0] = 0.0;
//    x0[1] = 0.0;
//    x0[2] = 1.0;

//    double t0 = 0.0;
//    double t1 = 1.0;
//    unsigned int N = 10;
//    double h = 0.1;

//    DoubleMatrix m;
//    CauchyProblemSystem(fs, t0, x0, m, t1, h, N);

//    IPrinter::printVector(m[0], "x1");
//    IPrinter::printVector(m[1], "x2");
//    IPrinter::printVector(m[2], "x3");

//    DoubleVector a1(N+1);
//    DoubleVector a2(N+1);
//    DoubleVector a3(N+1);
//    for (unsigned int i=0; i<=N; i++)
//    {
//        double t = i*h;
//        a1[i] = t*t;
//        a2[i] = t + sin(t);
//        a3[i] = cos(t);
//    }
//    IPrinter::printVector(a1, "x1");
//    IPrinter::printVector(a2, "x1");
//    IPrinter::printVector(a3, "x1");

//    SampleBorderHyperBolic::main(argc, argv);
//    QSimpleWaveWidget::main(argc, argv);

    return 0;
}
