#include <float.h>
#include <time.h>
#include <cmethods.h>
#include <grid/pibvp.h>
#include <grid/hibvp.h>
//#include <ode/cauchyp.h>
#include <ode/nlode1o.h>

#include "problem1/problem1L1.h"
#include "problem1/problem1L2.h"
#include "problem1/problem1L3.h"

#include "problem1/art_problem1.h"
#include "problem1/loadedheatequation.h"

#include "problem4/problem4ex1.h"

#include "high_order/singledifequ.h"
#include "high_order/systemdifequ.h"

//#include <../border/borderparabolicd.h>
//#include <../border/borderparabolicn.h>
//#include <../border/borderparabolic2d.h>
//#include <../border/borderhyperbolic2d.h>

#include <../border/grid/parabolicibvp1.h>
//#include <../border/grid/parabolicibvp2.h>
//#include <../border/grid/hyperbolicibvp1.h>
#include <../border/grid/newtonheatequationex1.h>

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

#include "ivp/nlode1oex1.h"

#include "load_sys/slodenlcsv.h"
#include "load_sys/slodenlcsm.h"
#include "load_sys/slodenlcsv2.h"

#include "numintegralexp1.h"

#include "utils/matrix.h"

#include "nonlinearfunctionsystem.h"
#include "load_sys/lode1oex1.h"

#include "matrixtest.h"

class A
{
public:
    A()
    {
        printf("A\n");
    }
};

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    //MatrixTest::Main(argc, argv);

    //srand(time(NULL));

    //NonLinearFunctionSystem::Main(argc, argv);

    //Example1::Main(argc, argv);
    //SystemLinearODENonLocalContionsV::Main(argc, argv);
    //SystemLinearODENonLocalContionsM::Main(argc, argv);
    Problem4Ex1::Main(argc, argv);
    //SystemLinearODENonLocalContionsV2::Main(argc, argv);
    //LinearODE1stOrderEx1::Main(argc, argv);

    //NumIntegralExp1::Main(argc, argv);

    //NonLinearODE1stOrderEx1::Main(argc, argv);
    //ParabolicIBVP1::Main(argc, argv);

    //HyperbolicControl2DM::Main(argc, argv);
    //HyperbolicControl2DMX::Main(argc, argv);
    //HyperbolicControl2DMV::Main(argc, argv);

    //BorderTest::Main(argc, argv);
    //BorderTest1::Main(argc, argv);
    //Example2::Main(argc, argv);
    //Example4::Main(argc, argv);
    //Example6::Main(argc, argv);
    //Example7::Main(argc, argv);

    //SingleDifEquation::Main(argc, argv);
    //SystemDifEquation::Main(argc, argv);

    //ArtProblem1::Main(argc, argv);
    //ArticleProblem1L3::Main(argc, argv);

    //BorderParabolicD::Main(argc, argv);
    //BorderParabolicN::Main(argc, argv);
    //BorderParabolic2D::Main(argc, argv);
    //BorderHyperbolic2D::Main(argc, argv);
    //Example5 e5;
    //BoundaryValueProblem1::Main(argc, argv);

    //ParabolicIBVP1::Main(argc, argv);
    //ParabolicIBVP2::Main(argc, argv);
    //HyperbolicIBVP1::Main(argc, argv);
    //NewtonHeatEquationEx1::Main(argc, argv);
    //Problem1L3::Main(argc, argv);
    //LoadedHeatEquation::Main(argc, argv);

    //HeatControl::Main(argc, argv);
    //HeatControl1::Main(argc, argv);

    //Rosenbrock::Main(argc, argv);
    //QuadraticFunction::Main(argc,argv);

    return 0;
}
