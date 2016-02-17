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

#include "control/cfunction.h"
#include "control/cfunction1.h"
#include "control/cfunction2.h"
#include "control/cfunction3.h"

#include "rnfunction/rosenbrock.h"
#include "rnfunction/bealesfunction.h"
#include "rnfunction/boothfunction.h"

#include "parabolic/1d/heatcontrol.h"
#include "parabolic/1d/heatcontroldeltaf.h"
#include "parabolic/1d/heatcontroldeltax.h"

#include "parabolic/2d/heatcontrol2d.h"
#include "parabolic/2d/heatcontrol2delta.h"
#include "parabolic/2d/heatcontrol2deltaf.h"
#include "parabolic/2d/heatcontrol2deltax.h"

#include "hyperbolic/1d/hyperbolic1dx.h"
#include "hyperbolic/1d/hyperboliccontrol1d.h"
#include "hyperbolic/1d/hyperboliccontrolx.h"
#include "hyperbolic/1d/hyperboliccontrolh.h"
#include "hyperbolic/2d/hyperboliccontrol2d.h"
#include "hyperbolic/2d/hyperboliccontrol2dm.h"
#include "hyperbolic/2d/hyperboliccontrol2dmv.h"
#include "hyperbolic/2d/hyperboliccontrol2dmx.h"

#include "point/pointcontrol11.h"
#include "point/pointcontrol.h"
#include "point/pointcontrol1.h"
#include "point/pointcontrol2.h"

#include "discrete/discreteheat.h"
#include "discrete/discretehyperbolic.h"
#include "discrete/discretehyperbolic1.h"

#include "border/borderparabolic2d.h"
#include "border/borderhyperbolic2d.h"

int main()
{
//    DoubleMatrix u;
//    BorderHyperbolic2D s;
//    s.calculateU(u, s.h1, s.h2, s.ht, s.N1, s.N2, s.M);
//    //IPrinter::printMatrix(u3[u3.size()-1]);
//    IPrinter::printMatrix(u);
//    printf("---\n");
//    FILE *file = fopen("u.txt", "w");
//    IPrinter::printMatrix(u, s.N2, s.N1, NULL, file);
//    fclose(file);

//    HeatControl2DeltaX::main();
//    HeatControlDeltaX::main();
//    DiscreteHyperbolic1::main();
//    HyperbolicControlH::main();

    HyperbolicControl2DMV::main();

    return 0;
}
