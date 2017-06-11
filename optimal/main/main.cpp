#include <float.h>
#include <time.h>
#include <cmethods.h>
#include <grid/pibvp.h>
#include <grid/hibvp.h>

#include "problem1/problem1L1.h"
#include "problem1/problem1L2.h"
#include "problem1/problem1L3.h"

#include "problem1/art_problem1.h"

//#include "loadedsystems.h"
#include "example1.h"
#include "example2.h"
#include "example3.h"
#include "example4.h"
//#include "example5.h"
#include "example6.h"
#include "example7.h"
#include "example.h"

//#include <../border/borderparabolicd.h>
//#include <../border/borderparabolicn.h>
//#include <../border/borderparabolic2d.h>
//#include <../border/borderhyperbolic2d.h>

#include <../border/grid/parabolicibvp1.h>
//#include <../border/grid/parabolicibvp2.h>
//#include <../border/grid/hyperbolicibvp1.h>

//#include <../hyperbolic/2d/hyperboliccontrol2d.h>
//#include <../hyperbolic/2d/hyperboliccontrol2d1.h>
//#include <../hyperbolic/2d/hyperboliccontrol2d21.h>
//#include <../hyperbolic/2d/hyperboliccontrol2d22.h>
//#include <../hyperbolic/2d/hyperboliccontrol2d23.h>
//#include <../hyperbolic/2d/hyperboliccontrol2d24.h>
//#include <../hyperbolic/2d/hyperboliccontrol2d24.h>
//#include <../hyperbolic/2d/hyperboliccontrol2dm.h>
//#include <../hyperbolic/2d/hyperboliccontrol2dmv.h>
//#include <../hyperbolic/2d/hyperboliccontrol2dmx.h>

//#include <../parabolic/1d/heatcontrol.h>
//#include <../parabolic/1d/heatcontrol1.h>

#include <../rnfunction/rosenbrock.h>
#include <../rnfunction/quadraticfunction.h>

#include "bordertest.h"
#include "bordertest1.h"
#include "sampleboundaryproblem1.h"

//#include <xlsxwriter.h>

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    //srand(time(NULL));

    //HyperbolicControl2DM::Main(argc, argv);
    //HyperbolicControl2DMX::Main(argc, argv);
    //HyperbolicControl2DMV::Main(argc, argv);

    //BorderTest::Main(argc, argv);
    //BorderTest1::Main(argc, argv);
    //Example2::Main(argc, argv);
    //Example4::Main(argc, argv);
    //Example6::Main(argc, argv);
    //Example7::Main(argc, argv);
    Example::Main(argc, argv);

//    ArtProblem1::Main(argc, argv);
//    ArticleProblem1L3::Main(argc, argv);

    //Example4::Main(argc, argv);
    //BorderParabolicD::Main(argc, argv);
    //BorderParabolicN::Main(argc, argv);
    //BorderParabolic2D::Main(argc, argv);
    //BorderHyperbolic2D::Main(argc, argv);
    //Example5 e5;
    //BoundaryValueProblem1::Main(argc, argv);

    //ParabolicIBVP1::Main(argc, argv);
    //ParabolicIBVP2::Main(argc, argv);
    //HyperbolicIBVP1::Main(argc, argv);

    //HeatControl::Main(argc, argv);
    //HeatControl1::Main(argc, argv);

    //Rosenbrock::Main(argc, argv);
    //QuadraticFunction::Main(argc,argv);

//    lxw_workbook  *workbook  = workbook_new("hello_world.xlsx");
//    lxw_worksheet *worksheet = workbook_add_worksheet(workbook, NULL);
//    worksheet_write_string(worksheet, 0, 0, "Hello", NULL);
//    worksheet_write_number(worksheet, 1, 0, 123, NULL);
//    workbook_close(workbook);

    return 0;
}
