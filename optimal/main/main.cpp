#include "headers.h"

<<<<<<< .mine
#include "border/sampleexample1.h"

||||||| .r690
#include "widget/qsimplewavewidget.h"
#include <matrix.h>
#include <time.h>
#include "border/sampleexample1.h"

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

=======
>>>>>>> .r691
int main(int argc, char ** argv)
{
<<<<<<< .mine
    SampleMain();
||||||| .r690
    srand(time(NULL));

    SampleMain();

=======
//    SampleMain();
//    SampleLoaderBorder::main();
    HyperbolicControl2D21::main(argc, argv);
>>>>>>> .r691
    return 0;
}
