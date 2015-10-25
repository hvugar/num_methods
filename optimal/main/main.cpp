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
#include <heatequation.h>

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

#include "heat/1d/heatcontrol.h"
#include "heat/1d/heatcontroldeltaf.h"
#include "heat/1d/heatcontroldeltax.h"
#include "heat/2d/heatcontrol2d.h"
#include "heat/2d/heatcontrol2delta.h"
#include "heat/2d/heatcontrol2deltaf.h"
#include "heat/2d/heatcontrol2deltax.h"

#include "hiperbolic/hiperbolic1dx.h"

int main()
{
//    HeatControl2DeltaX::main();
    Hiperbolic1DX::main();

    return 0;
}



