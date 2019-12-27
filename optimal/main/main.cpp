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
#include <ode/first_order_linear_ode.h>
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
    srand(time(0));

//    DoubleMatrix m(5,5);
//    Random::fillMatrix(m, 1, 5, 2);
//    IPrinter::print(m, 5, 5, 6, 2);
//    IPrinter::printSeperatorLine();
//    DoubleMatrix i = m;
//    i.inverse();
//    IPrinter::print(i, 5, 5, 6, 2);
//    IPrinter::printSeperatorLine();
//    DoubleMatrix e1 = i*m;
//    IPrinter::print(e1, 5, 5, 6, 2);
//    IPrinter::printSeperatorLine();
//    DoubleMatrix e2 = m*i;
//    IPrinter::print(e2, 5, 5, 6, 2);
//    IPrinter::printSeperatorLine();
//    return 0;

    //ConjugateGradinetTest::Main(argc, argv);

    //FirstOrderLinearODEEx1::Main(argc, argv);
    //FirstOrderNonLinearODErEx1::Main(argc, argv);
    //SecondOrderLinearODEEx1::Main(argc, argv);

    FirstOrderLinearODEEx1::Main(argc, argv);
    //NonLinearODE1stOrderEx2::Main(argc, argv);
    //SecondOrderLinearODEEx1::Main(argc, argv);

    //WaveEquationIBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //WaveEquationFBVP::Main(argc, argv);

    //HeatEquationIBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //HeatEquationFBVP::Main(argc, argv);

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
    h0p::ProblemSolver::Main(argc, argv);
    //Problem1HDirichlet1::Main(argc, argv);

    //Problem2HSolver::Main(argc, argv);
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
