#include <float.h>
#include <time.h>
#include <cmethods.h>
#include <grid/pibvp.h>
#include <grid/hibvp.h>

#include "problem1.h"
#include "loadedsystems.h"
#include "example1.h"
#include "example2.h"
#include "example3.h"
#include "example4.h"
#include "example5.h"
#include "problem1.h"

#include <../border/borderparabolicd.h>
#include <../border/borderparabolicn.h>
#include <../border/borderparabolic2d.h>
#include <../border/borderhyperbolic2d.h>

#include <../hyperbolic/2d/hyperboliccontrol2d.h>
#include <../hyperbolic/2d/hyperboliccontrol2d1.h>
#include <../hyperbolic/2d/hyperboliccontrol2d21.h>
#include <../hyperbolic/2d/hyperboliccontrol2d22.h>
#include <../hyperbolic/2d/hyperboliccontrol2d23.h>
#include <../hyperbolic/2d/hyperboliccontrol2d24.h>
#include <../hyperbolic/2d/hyperboliccontrol2d24.h>
#include <../hyperbolic/2d/hyperboliccontrol2dm.h>
#include <../hyperbolic/2d/hyperboliccontrol2dmv.h>
#include <../hyperbolic/2d/hyperboliccontrol2dmx.h>

#include "bordertest.h"
#include "bordertest1.h"
#include "sampleboundaryproblem1.h"

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    srand(time(NULL));

    //HyperbolicControl2DM::Main(argc, argv);

    //BorderTest::Main(argc, argv);
    //BorderTest1::Main(argc, argv);
    //Example2::Main(argc, argv);
    //Problem1::Main(argc, argv);
    //Example4::Main(argc, argv);
    //BorderParabolicD::Main(argc, argv);
    //BorderParabolicN::Main(argc, argv);
    //BorderParabolic2D::Main(argc, argv);
    //BorderHyperbolic2D::Main(argc, argv);
    //Example5 e5;
    BoundaryValueProblem1::Main(argc, argv);
    return 0;
}
