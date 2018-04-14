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

#include "problem2P/1d/problem2.h"
#include "problem2P/2d/ex/problem22dex1.h"
#include "problem2P/2d/ex/problem22dex2.h"
#include "problem2P/2d/ex/problem22dex3.h"
#include "problem2P/2d/ex/problem22dex4.h"
#include "problem2P/2d/ex/problem22dex5.h"
#include "problem2P/2d/ex/expol.h"
#include "problem2P/2d/cproblem2forward2d.h"
#include "problem2P/2d/cproblem2backward2d.h"
#include "problem2P/2d/dirakdelta.h"
#include "problem2P/2d/ex/p2_article.h"

#include "problem2H/iproblem2h2d.h"

#include "problem4/problem4ex1.h"
#include "problem4/problem4ex2.h"

#include "problem5/nllparabolic.h"

//#include <../border/borderparabolicd.h>
//#include <../border/borderparabolicn.h>
//#include <../border/borderparabolic2d.h>
#include <../border/borderhyperbolic2d.h>

#include <../border/grid/parabolicibvp1.h>
#include <../border/grid/parabolicibvp2.h>

//#include <../border/grid/hyperbolicibvp1.h>
//#include <../border/grid/newtonheatequationex1.h>

//#include <../parabolic/1d/heatcontrol.h>
//#include <../parabolic/1d/heatcontrol1.h>

#include <../rnfunction/rosenbrock.h>
#include <../rnfunction/quadraticfunction.h>

#include "load_sys/slodenlcsv.h"
#include "load_sys/slodenlcsm.h"
#include "load_sys/slodenlcsv2.h"
#include "load_sys/lode1oex1.h"

#include "ivp/nlode1oex1.h"

#include <utils/matrix.h>
#include <utils/random.h>

#include "nonlinearequationex1.h"
#include "loadedlinearode1order.h"

#include <QtGui>

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    QGuiApplication app(argc, argv);


//    DoubleVector a(11, -0.5); a[0] = 0.0;
//    DoubleVector b(11, +2.000015);
//    DoubleVector c(11, -0.5); c[10] = 0.0;
//    DoubleMatrix w(11, 11, 0.0);
//    DoubleVector d(11, 0.0);
//    DoubleVector x(11, 0.0);

//    w[2][7] = 1;
//    w[7][2] = -5;
//    w[9][2] = +7;
//    w[9][8] = 1.0;
//    w[4][7] = 10;

//    d[0] = 1.000015;
//    d[1] = 2.00003;
//    d[2] = 11.000045;
//    d[3] = 4.00006;
//    d[4] = 85.000075;
//    d[5] = 6.00009;
//    d[6] = 7.000105;
//    d[7] = -6.99988;
//    d[8] = 14.000135;
//    d[9] = 23;
//    d[10] = 10.000075;


//    LinearEquation::func1(a.data(), b.data(), c.data(), d.data(), w.data(), x.data(), 11);
//    IPrinter::print(x);
//    return 0.0;







    //BorderHyperbolic2D::Main(argc, argv);
    //BorderHyperbolic2DN::Main(argc, argv);
    //srand(time(NULL));

    //Problem2Article::Main(argc, argv);
    IProblem2H::IProblem2H2D::Main(argc, argv);

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
