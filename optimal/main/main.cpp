 #include <float.h>
#include <time.h>
#include <cmethods.h>
#include <grid/pibvp.h>
#include <grid/hibvp.h>
#include <utils/matrix.h>
#include <utils/random.h>
#include <deltagrid.h>
#include <r1minimize.h>

#include <ode/nlode1o.h>
#include <ode/lode1o.h>
//#include <ode/firstorderlinearodeex1.h>
//#include <ode/firstordernonlinearodeex1.h>
//#include <ode/secondorderlinearodeex1.h>
#include <heat_equation_ibvp.h>
#include <wave_equation_ibvp.h>

#include "../problem0H/problem0h_solver.h"
#include "../problem1H/problem1h_example.h"
#include "../problem2P/problem2p_solver.h"
#include "../problem2H/problem2h_solver.h"

#include "../problem1P/problem1p_solver.h"

#include "test/delta_grid_2d_ext1.h"
#include "test/nonlinear_equation_ex1.h"
#include "conjugate_gradinet_test.h"

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    std::cout << "__cplusplus: " << __cplusplus << std::endl;

    //ConjugateGradinetTest::Main(argc, argv);

    //FirstOrderLinearODEEx1::Main(argc, argv);
    //FirstOrderNonLinearODErEx1::Main(argc, argv);
    //SecondOrderLinearODEEx1::Main(argc, argv);

    //LinearODE1stOrderEx1::Main(argc, argv);
    //NonLinearODE1stOrderEx2::Main(argc, argv);
    //SecondOrderLinearODEEx1::Main(argc, argv);

    //WaveEquationIBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //WaveEquationFBVP::Main(argc, argv);

    //HeatEquationIBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    HeatEquationFBVP::Main(argc, argv);

    //p1p::ProblemSolver::Main(argc, argv);

    return 0;


    //DeltaGrid2DExt1::Main(argc, argv);
    //puts("-------------------------");
    //DeltaGrid1DExt1::Main(argc, argv);
    //puts("-------------------------");

    //return 0;

    //srand(static_cast<unsigned int>(time(nullptr)));

    //NonLocalSystem::Main(argc, argv);


    //IPrinter::printSeperatorLine();
    //ConjugateCC1IHyperbolicIBVP1::Main(argc, argv);
    //Problem0HFunctional::Main(argc, argv);
    //Problem1HDirichlet1::Main(argc, argv);

    Problem2HSolver::Main(argc, argv);
    //Problem2HDirichlet::Main(argc, argv);
    //Problem2HDirichletDelta::Main(argc, argv);

    //CCParabolicIBVP1::Main(argc, argv);
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
