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

#include "cfunction1.h"
#include "cfunction2.h"
#include "cfunction3.h"
#include "cfunction.h"
#include "rosenbrock.h"
#include "bealesfunction.h"
#include "boothfunction.h"
#include "pointcontrol.h"
#include "pointcontrol1.h"
#include "pointcontrol2.h"
#include "utils.h"

#include "heat/heatcontrol.h"
#include "heat/heatcontrol2d.h"
#include "heat/heatcontroldelta.h"
#include "heat/headcontrol2d.h"

double fxt(double x, double y, double z) { return x*y*z; }

int main()
{
    unsigned int M = 1000;
    unsigned int N1 = 1000;
    unsigned int N2 = 1000;

    double ht = 1.0/M;
    double h1 = 1.0/N1;
    double h2 = 1.0/N2;

    double norm = 0.0;

    for (unsigned int k=0; k<M; k++)
    {
        for (unsigned int j=0; j<N2; j++)
        {
            for (unsigned int i=0; i<N1; i++)
            {
                unsigned int i1 = i;
                unsigned int i2 = i+1;
                unsigned int j1 = j;
                unsigned int j2 = j+1;
                unsigned int k1 = k;
                unsigned int k2 = k+1;

                double f1 = fxt(i1*h1, j1*h2, k1*ht);
                double f2 = fxt(i2*h1, j1*h2, k1*ht);
                double f3 = fxt(i1*h1, j2*h2, k1*ht);
                double f4 = fxt(i2*h1, j2*h2, k1*ht);
                double f5 = fxt(i1*h1, j1*h2, k2*ht);
                double f6 = fxt(i2*h1, j1*h2, k2*ht);
                double f7 = fxt(i1*h1, j2*h2, k2*ht);
                double f8 = fxt(i2*h1, j2*h2, k2*ht);

                norm += f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
            }
        }
    }
    norm = norm * (h1*h2*ht)*0.125;

    printf("%.8f\n", norm);

//    HeatControl2D::main();
    return 0;
}



