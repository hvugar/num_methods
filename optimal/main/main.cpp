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
#include <doublevector.h>
#include <parabolicequation.h>
#include <hyperbolicequation.h>
#include <r1minimize.h>

#include <iostream>
#include <stdexcept>

#include "cfunction1.h"
#include "cfunction2.h"
#include "cfunction3.h"
#include "cfunction.h"
#include "rosenbrock.h"
#include "bealesfunction.h"
#include "boothfunction.h"

#include "parabolic/1d/heatcontrol.h"
#include "parabolic/1d/heatcontrol1.h"
#include "parabolic/1d/heatcontroldeltaf.h"
#include "parabolic/1d/heatcontroldeltax.h"
#include "parabolic/2d/heatcontrol2d.h"
#include "parabolic/2d/heatcontrol2delta.h"
#include "parabolic/2d/heatcontrol2deltaf.h"
#include "parabolic/2d/heatcontrol2deltax.h"

#include "hyperbolic/hyperbolic1dx.h"
#include "hyperbolic/hyperboliccontrol1d2.h"
#include "hyperbolic/hyperboliccontrol1d3.h"
#include "hyperbolic/hyperboliccontrol1d4.h"
#include "hyperbolic/hyperboliccontrol1dt.h"

#include "point/pointcontrol11.h"
#include "point/pointcontrol.h"
#include "point/pointcontrol1.h"
#include "point/pointcontrol2.h"

#include "discrete/discreteheat.h"

#include <cmethods.h>

double f(double x, double y)
{
    return x*x + y*y;
}

int main()
{
//    double Nx = 10;
//    double Ny = 10;
//    double hx = 0.1;
//    double hy = 0.1;

//    double sum = 0.25*hx*hy*(f(0*hx, 0*hy) + f(0*hx,Ny*hy) + f(Nx*hx, 0*hy) + f(Nx*hx, Ny*hy));
//    for (unsigned int i=1; i<=Nx-1; i++)
//    {
//        sum += 0.50*hx*hy*f(i*hx, 0*hy);
//        sum += 0.50*hx*hy*f(i*hx, Ny*hy);
//    }
//    for (unsigned int j=1; j<=Ny-1; j++)
//    {
//        sum += 0.50*hx*hy*f(0*hx, j*hy);
//        sum += 0.50*hx*hy*f(Nx*hx, j*hy);
//    }
//        for (unsigned int i=1; i<=Nx-1; i++)
//        {
//            for (unsigned int j=1; j<=Ny-1; j++)
//            {
//                sum += hx*hy*f(i*hx, j*hy);
//            }
//        }
//    printf("%.12f\n", sum);


//    sum = 0.0;
//    for (unsigned int i=0; i<=Nx-1; i++)
//    {
//        for (unsigned int j=0; j<=Ny-1; j++)
//        {
//            sum += f(i*hx, j*hy) + f(i*hx, (j+1)*hy) + f((i+1)*hx, j*hy) + f((i+1)*hx, (j+1)*hy);
//        }
//    }
//    sum = sum * 0.25 * hx *hy;
//    printf("%.12f\n", sum);


//    HyperbolicControl1DT::main();
    puts("------------------------------------------------------------------------------------------------");
//    HyperbolicControl1D4::main();
    DiscreteHeat::main();
//    HeatControl1::main();
//    HeatControl2Delta::main();
//    PointControl11::main();
//    Rosenbrock::main();
    return 0;
}



