#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gridmethod.h>
#include <rungekutta.h>

#include "cfunction1.h"
#include "cfunction2.h"
#include "cfunction.h"
#include "rosenbrock.h"
#include "bealesfunction.h"
#include "boothfunction.h"
#include "heatcontrol.h"
#include "pointcontrol.h"
#include "pointcontrol1.h"
#include "pointcontrol2.h"
#include "utils.h"

void write(const double *x, unsigned int n, const char* filename)
{
    FILE *f = fopen(filename, "w");
    for (unsigned int i=0; i<n; i++)
    {
        if (i%10==0)
            fprintf(f, "%.10f\n", x[i]);
    }
    fclose(f);
}

double f(double t, double x)
{
    return x + 2*t - t*t;
//    return 2*t - x + t*t;
}

void a()
{
    const unsigned int n =  10001;
    double *x = (double*) malloc(sizeof(double) * n);

    double dt = 0.0001;
    x[0] = 0.0;
    double t = 0.0;
    double _x0 = x[0];

    for (unsigned int i=1; i<n; i++)
    {
        if (fabs(t-0.2) < dt/10.0) _x0 = _x0 + 10.0;
        //if (fabs(t-0.5) < dt/10.0) _x0 = _x0 + 0.0;//10.6541501861;
        //if (fabs(t-0.8) < dt/10.0) _x0 = _x0 + 0.0;// 14.3815984660;

//        if (fabs(t-0.2) < dt/10.0) _x0 = _x0 + 0.3;
//        if (fabs(t-0.5) < dt/10.0) _x0 = _x0 + 0.4;
//        if (fabs(t-0.8) < dt/10.0) _x0 = _x0 + 0.5;

        double k1 = f(t, _x0);
        double k2 = f(t+dt/2.0, _x0+(dt/2.0)*k1);
        double k3 = f(t+dt/2.0, _x0+(dt/2.0)*k2);
        double k4 = f(t+dt,     _x0+dt*k3);

        _x0 = _x0 + (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        t = t + dt;
        x[i] = _x0;

        if (i >= 1998 && i <= 2002) printf("%d %.6f %.10f\n", i, t, _x0);
    }

    printf("%.6f %.10f\n", t, x[n-1]);

    write(x, n, "test.txt");

    free(x);
}

int main()
{
//    a();
    PointControl2::main();
    return 0;
}

