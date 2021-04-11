#include <float.h>
#include <time.h>
#include <cmath>
#include <iostream>

//#include <cmethods.h>
//#include <grid/pibvp.h>
//#include <grid/hibvp.h>
//#include <utils/matrix.h>
//#include <utils/random.h>
//#include <deltagrid.h>
//#include <r1minimize.h>

#include <ode/first_order_linear_ode.h>
#include <ode/second_order_linear_ode.h>
//#include <ode/firstordernonlinearodeex1.h>
//#include <ode/secondorderlinearodeex1.h>
#include <heat_equation_ibvp.h>
//#include <wave_equation_ibvp.h>

//#include "../problem0H/problem0h_solver.h"
//#include "../problem1H/problem1h_example.h"
//#include "../problem2P/problem2p_solver.h"
//#include "../problem2H/problem2h_solver.h"

//#include "../problem1P/problem1p_solver.h"
//#include "../problem1H/problem1h_solver.h"
#include "../problem3P/solver.h"
#include "../problem3P/solver1.h"
#include "../problem3P/solver2.h"
#include "../problem3P/solver3.h"
//#include "../problem3P/heat_equation_ibvp.h"

#include "test/delta_grid_2d_ext1.h"
//#include "test/nonlinear_equation_ex1.h"
//#include "conjugate_gradinet_test.h"
#include <QGuiApplication>

#include <chrono>
#include <thread>
#include "../border/test_function.h"

int main(int argc, char *argv[])
{
    //srand(static_cast<unsigned int>(time(nullptr)));
    QGuiApplication app(argc, argv);

    //test1();
    //return 0;

    //ConjugateGradinetTest::Main(argc, argv);

    FirstOrderLinearODEIVP::Main(argc, argv);
    //FirstOrderLinearODEFVP::Main(argc, argv);
    //SecondOrderLinearODEIBVP::Main(argc, argv);
    //SecondOrderLinearODEFBVP::Main(argc, argv);
    //FirstOrderNonLinearODErEx1::Main(argc, argv);
    //SecondOrderLinearODEEx1::Main(argc, argv);

    //HeatEquationIBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //HeatEquationFBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //LoadedHeatEquationIBVP::Main(argc, argv);
    //LoadedHeatEquationFBVP::Main(argc, argv);

    //WaveEquationIBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //WaveEquationFBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();

    //WaveEquationIBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //WaveEquationFBVP::Main(argc, argv);
    //IPrinter::printSeperatorLine();

    //h0p::ProblemSolver::Main(argc, argv);
    //p1p::ProblemSolver::Main(argc, argv);
    //h1p::ProblemSolver::Main(argc, argv);
    //p3p::Solver::Main(argc, argv);
    //p3p1::Solver1::Main(argc, argv);
    //p3p0::Functional::Main(argc, argv);
    ///p3p3::Functional::Main(argc, argv);
    //p3p2::Functional::Main(argc, argv);
    //p3p::HeatEquationIBVP1::Main(argc, argv);

    //return 0;

    //DeltaGrid1DExt1::Main(argc, argv);
    //puts("-------------------------");
    //DeltaGrid2DExt1::Main(argc, argv);
    //puts("-------------------------");

    //return 0;

    //srand(static_cast<unsigned int>(time(nullptr)));

    //NonLocalSystem::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //ConjugateCC1IHyperbolicIBVP1::Main(argc, argv);
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

    std::cout << "Program finished successfully." << std::endl;
    return 0;
}
