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

#include "heat2d.h"

int main()
{
    DoubleMatrix u;
    DoubleMatrix f;
    DoubleMatrix x1;
    DoubleMatrix x2;

    Heat2DControl hc;
    hc.calculateX(u, f, x1, x2);
//    GridMethod::VariableDirectionsMethod(0, 0, 0, 0, 0, 0);
    return 0;
}



