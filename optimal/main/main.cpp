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
#include "../problem2H/problem2h_solver.h"
#include "../problem2H/problem2h_ibvp.h"
#include "ivp/nlode1oex1.h"

#include <utils/matrix.h>
#include <utils/random.h>

#include "nonlinearequationex1.h"
#include "loadedlinearode1order.h"

#include <grid/hpibvp.h>
#include <grid/hyperbolicibvp1.h>
#include <grid/parabolicibvp1.h>
#include <grid/newtonheatequationex1.h>
#include "heatequationibvp1.h"
#include <deltagrid.h>

#include <QtGui>
#include <r1minimize.h>

double funcU(double x, double y)
{
    return sin(M_PI*x)*sin(M_PI*y);
}

double funcUx(double x, double y)
{
    return M_PI*cos(M_PI*x)*sin(M_PI*y);
}

double funcUy(double x, double y)
{
    return M_PI*sin(M_PI*x)*cos(M_PI*y);
}

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
//    DeltaGrid2D delta;
//    unsigned int N = 100, M = 100;
//    double hx = 0.01, hy = 0.01;
//    delta.initGrid(N, hx, M, hy);
//    SpacePoint p; p.x = 0.5400; p.y = 0.3000;
//    delta.distributeGauss(p, 1, 1);
//    DoubleMatrix u(M+1, N+1);
//    for (unsigned int m=0; m<=M; m++)
//    {
//        for (unsigned int n=0; n<=N; n++)
//        {
//            u[m][n] = funcU(n*hx, m*hy);
//        }
//    }
//    double U = 0.0, Ux = 0.0, Uy = 0.0;

//    double a, b, c;
//    a = delta.consentrateInPoint(u, b, c);
//    printf("%.6f %.6f %.6f\n%.6f %.6f %.6f\n %4d %4d\n", a, b, c,
//           funcU(p.x, p.y), funcUx(p.x, p.y), funcUy(p.x, p.y), delta.rx(), delta.ry());

//    //U = u[delta.ry()][delta.rx()];
//    for (unsigned int m=delta.minY(); m<=delta.maxY(); m++)
//    {
//        for (unsigned int n=delta.minX(); n<=delta.maxX(); n++)
//        {
//            double k = 1.0;
//            //            if (m==delta.minY()) k *= 0.5;
//            //            if (m==delta.maxY()) k *= 0.5;
//            //            if (n==delta.minX()) k *= 0.5;
//            //            if (n==delta.maxX()) k *= 0.5;
//            U += k*u[m][n]*delta.weight(n,m)*(hx*hy);

//            if (delta.isCenter(n,m))
//            {
//                const unsigned int rx = static_cast<const unsigned int>(delta.rx());
//                const unsigned int ry = static_cast<const unsigned int>(delta.ry());
//                const unsigned int px = static_cast<const unsigned int>(delta.p().x);
//                const unsigned int py = static_cast<const unsigned int>(delta.p().y);

////                Ux = (u[ry][rx+1] - u[ry][rx-1])/(2.0*hx);
////                Uy = (u[ry+1][rx] - u[ry-1][rx])/(2.0*hy);
////                Ux += ((delta.p().x-rx*hx)/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
////                Uy += ((delta.p().y-ry*hy)/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);


//                Ux = (u[ry][rx-2] - 8.0*u[ry][rx-1] + 8.0*u[ry][rx+1] - u[ry][rx+2])/(12.0*hx);
//                Uy = (u[ry-2][rx] - 8.0*u[ry-1][rx] + 8.0*u[ry+1][rx] - u[ry+2][rx])/(12.0*hy);

//                Ux += ((delta.p().x-rx*hx))*((-2.0*u[ry][rx-2] + 32.0*u[ry][rx-1] - 60.0*u[ry][rx] + 32.0*u[ry][rx+1] - 2.0*u[ry][rx+2])/(24.0*hx*hx));
//                Uy += ((delta.p().y-ry*hy))*((-2.0*u[ry-2][rx] + 32.0*u[ry-1][rx] - 60.0*u[ry][rx] + 32.0*u[ry+1][rx] - 2.0*u[ry+2][rx])/(24.0*hy*hy));
//            }
//        }
//    }
//    printf("%.6f %.6f %.6f\n%.6f %.6f %.6f\n %4d %4d\n", U, Ux, Uy, funcU(p.x, p.y), funcUx(p.x, p.y), funcUy(p.x, p.y), delta.rx(), delta.ry());
//    return 0;



    //srand(time(0));
    //CcIHyperbolicIBVP1::Main(argc, argv);
    //CdIHyperbolicIBVP1::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //ConjugateCC1IHyperbolicIBVP1::Main(argc, argv);
    Problem2HDirichlet::Main(argc, argv);

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

    //Problem2Article::Main(argc, argv);
    //IProblem2H::IProblem2H2D::Main(argc, argv);
    //IPrinter::printSeperatorLine("+");
    //Problem2HDirichlet::Main(argc, argv);
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
