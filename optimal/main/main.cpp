#include <float.h>
#include <time.h>
#include <cmethods.h>
#include <grid/pibvp.h>
#include <grid/hibvp.h>
#include <ode/nlode1o.h>

//#include "problem1/problem1L1.h"
//#include "problem1/problem1L2.h"
//#include "problem1/problem1L3.h"
//#include "problem1/art_problem1.h"
//#include "problem1/loadedheatequation.h"

//#include "problem2P/1d/problem2.h"
//#include "problem2P/2d/ex/problem22dex1.h"
//#include "problem2P/2d/ex/problem22dex2.h"
//#include "problem2P/2d/ex/problem22dex3.h"
//#include "problem2P/2d/ex/problem22dex4.h"
//#include "problem2P/2d/ex/problem22dex5.h"

//#include "../problem2P/2d/ex/expol.h"
//#include "../problem2P/2d/cproblem2forward2d.h"
//#include "../problem2P/2d/cproblem2backward2d.h"
//#include "../problem2P/2d/dirakdelta.h"
#include "../problem2P/2d/ex/p2_article.h"

#include "../problem2P/problem2p_solver.h"
//#include "../problem2H/problem2h_solver.h"
#include "ivp/nlode1oex1.h"

#include <utils/matrix.h>
#include <utils/random.h>

#include "nonlinearequationex1.h"
#include "loadedlinearode1order.h"

#include <grid/hpibvp.h>
#include "heatequationibvp1.h"

#include <QtGui>
#include <r1minimize.h>

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    QGuiApplication app(argc, argv);

    //Problem2Article::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    Problem2PNeumann::Main(argc, argv);
    //Problem2HNDirichlet::Main(argc, argv);
    return 0;

    //BorderHyperbolic2D::Main(argc, argv);
    //BorderHyperbolic2DN::Main(argc, argv);
    //srand(time(NULL));

    //Problem2Article::Main(argc, argv);
    //IProblem2H::IProblem2H2D::Main(argc, argv);
    //IPrinter::printSeperatorLine("+");
    //Problem2HDirichlet::Main(argc, argv);
    //Problem2HNDirichlet::Main(argc, argv);
    //Problem2HNM::Main(argc, argv);
    //return 0;

    //LinearEquation::func1(a, b, c, d, e, x, N);

    //IPrinter::printVector(x,N,NULL,N);
    //IPrinter::printSeperatorLine();

    //func2(a, b, c, d, e, x, N);

    //MatrixTest::Main(argc, argv);

    //Problem2::Main(argc, argv);

    //CProblem2Forward2D::Main(argc, argv);
    //CProblem2Backward2D::Main(argc, argv);

    //Problem22D::Main(argc, argv);
    //Problem22DEx1::Main(argc, argv);
    //Problem22DEx2::Main(argc, argv);
    //Problem22DEx3::Main(argc, argv);
    //Problem22DEx4::Main(argc, argv);
    //Problem22DEx5::Main(argc, argv);
    //ExpOptimalLetters::Main(argc, argv);

    //NLLIParabolicIBVP::Main(argc, argv);

    //NonLinearEquationEx1::Main(argc, argv);

    //LoadedLinearODE1Order::Main(argc, argv);

    //Example1::Main(argc, argv);
    //SystemLinearODENonLocalContionsV::Main(argc, argv);
    //SystemLinearODENonLocalContionsM::Main(argc, argv);
    //Problem4Ex2::Main(argc, argv);
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
