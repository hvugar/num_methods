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
#include "cfunction3.h"
#include "cfunction.h"
#include "rosenbrock.h"
#include "bealesfunction.h"
#include "boothfunction.h"
#include "heatcontrol.h"
#include "pointcontrol.h"
#include "pointcontrol1.h"
#include "pointcontrol2.h"
#include "utils.h"

double u(double x, double t) { return x*x*x + t*t; }
double fi(double x) { return u(x, 0.0); }
double psi1(double t) { return u(0.0, t); }
double psi2(double t) { return u(1.0, t); }
double f(double x, double t) { return 2.0*t - 6.0*x; }

int main()
{
    const unsigned int N = 1000;
    const unsigned int M = 1000;
    const double dx = 1.0/N;
    const double dt = 1.0/M;

    DoubleMatrix U;
    U.resize(M+1);

    U[0].resize(N+1);
    for (unsigned int i=0; i<=N; i++) U[0][i] = fi(dx*i);
    for (unsigned int j=0; j<=M; j++)
    {
        U[j].resize(N+1);
        U[j][0] = psi1(dt*j);
        U[j][N] = psi2(dt*j);
    }

    for (unsigned int j=0; j<M; j++)
    {
        for (unsigned int i=1; i<=N-1; i++)
        {
            U[j+1][i] = U[j][i] + dt*((U[j][i+1] - 2.0*U[j][i] + U[j][i-1])/(dx*dx) + f(i*dx, (j+1)*dt));
        }
    }

    unsigned int p=100;
    for (unsigned int j=0; j<=M; j++)
    {
        if (j%p == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                if (i%p == 0) printf("%14.8f ", U[j][i]);
            }
            puts("");
        }
    }

//    GridMethod::VariableDirectionsMethod(0, 0, 0, 0, 0, 0);
    return 0;
}



