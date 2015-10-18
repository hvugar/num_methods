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
#include "heat/heatcontroldeltaf.h"
#include "heat/heatcontroldeltax.h"
#include "heat/heatcontrol2d.h"
#include "heat/heatcontrol2delta.h"
#include "heat/heatcontrol2deltaf.h"

double fxt(double x, double y, double z) { return x*y*z; }

int main()
{
    HeatControlDeltaX::main();
    return 0;
}



