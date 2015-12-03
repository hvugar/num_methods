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
#include "pointcontrol.h"
#include "pointcontrol1.h"
#include "pointcontrol2.h"
#include "utils.h"

#include "heat/1d/heatcontrol.h"
#include "heat/1d/heatcontrol1.h"
#include "heat/1d/heatcontroldeltaf.h"
#include "heat/1d/heatcontroldeltax.h"
#include "heat/2d/heatcontrol2d.h"
#include "heat/2d/heatcontrol2delta.h"
#include "heat/2d/heatcontrol2deltaf.h"
#include "heat/2d/heatcontrol2deltax.h"

#include "hyperbolic/hyperbolic1dx.h"
#include "hyperbolic/hyperboliccontrol1d2.h"
#include "hyperbolic/hyperboliccontrol1d3.h"
#include "hyperbolic/hyperboliccontrol1d4.h"
#include "hyperbolic/hyperboliccontrol1dt.h"

int main()
{
    //HyperbolicControl1DT::main();
    //puts("------------------------------------------------------------------------------------------------");
    HyperbolicControl1D4::main();
    //HeatControl2DeltaX::main();
    return 0;
}



