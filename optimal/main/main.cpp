#include <stdio.h>
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

int main()
{
//    PointControl::main();
//    PointControl1::main();
    PointControl2::main();
    //        Rosenbrock::main();
    //    BealesFunction::main();
    //    BoothFunction::main();
    //    puts("*****************************************************************");
    //    CFunction1::main();
    //    CFunction2::main();
    //    ControlFunction::main();
    //    HeatControl::main();

    return 0;
}

