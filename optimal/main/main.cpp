#include <float.h>
#include <time.h>
#include <cmethods.h>
#include <grid/pibvp.h>
#include <grid/hibvp.h>

#include <ode/nlode1o.h>
#include <ode/lode1o.h>
#include <ode/firstorderlinearodeex1.h>
#include <ode/firstordernonlinearodeex1.h>
#include <ode/secondorderlinearodeex1.h>

#include "../problem2P/problem2p_solver.h"
#include "../problem2H/problem2h_solver_delta.h"

#include "../problem1H/problem1h_example.h"

#include <utils/matrix.h>
#include <utils/random.h>

#include <grid/hyperbolicibvp1.h>
#include <grid/parabolicibvp1.h>
#include <grid/newtonheatequationex1.h>
#include <deltagrid.h>

#include <QtGui>
#include <r1minimize.h>

#include <ode/firstorderlinearodeex1.h>
#include <ode/firstordernonlinearodeex1.h>
#include <ode/secondorderlinearodeex1.h>

#include "test/deltagrid2dext1.h"
#include <testwaveequation.h>

#include <problem0h_solver.h>

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    FirstOrderLinearODEEx1::Main(argc, argv);
    //FirstOrderNonLinearODErEx1::Main(argc, argv);
    //SecondOrderLinearODEEx1::Main(argc, argv);

    //LinearODE1stOrderEx1::Main(argc, argv);
    //NonLinearODE1stOrderEx2::Main(argc, argv);
    //SecondOrderLinearODEEx1::Main(argc, argv);
    //CdIHyperbolicIBVP1::Main(argc, argv);

    //Problem0HFunctional::Main(argc, argv);
    //TestWaveEquation::Main(argc, argv);
    //TestWaveEquationEx1::Main(argc, argv);

    //DeltaGrid2DExt1::Main(argc, argv);
    //DeltaGrid1DExt1::Main(argc, argv);

    //return 0;

    //srand(static_cast<unsigned int>(time(nullptr)));

    //NonLocalSystem::Main(argc, argv);


    //CcIHyperbolicIBVP1::Main(argc, argv);
    //CdIHyperbolicIBVP1::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //ConjugateCC1IHyperbolicIBVP1::Main(argc, argv);
    //Problem1HDirichlet1::Main(argc, argv);
    //Problem2HDirichlet::Main(argc, argv);
    //Problem2HDirichletDelta::Main(argc, argv);

    //CCParabolicIBVP1::Main(argc, argv);
    //CCIHyperbolicIBVP1::Main(argc, argv);
    //HyperbolicIBVP2::Main(argc, argv);
    //return 0;
    //QGuiApplication app(argc, argv);

    //Problem2Article::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //Problem2PNeumann::Main(argc, argv);
    //Problem2HNDirichletForward1 p;
    //Problem2HNDirichlet::Main(argc, argv);

    //BorderHyperbolic2D::Main(argc, argv);
    //BorderHyperbolic2DN::Main(argc, argv);
    //srand(time(NULL));

    return 0;
}
