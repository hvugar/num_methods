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
    unsigned int N = 10;
    double hx = 0.1;
    double ht = 0.05;
    double alpha = -((ht*ht)/(hx*hx));
    double betta = 1.0+2.0*((ht*ht)/(hx*hx));
    DoubleVector u00(N+1, 0.0);
    DoubleVector u10(N+1, 0.0);
    DoubleVector u15(N+1, 0.0);
    DoubleVector u20(N+1, 0.0);
    u10[5] = 50.0;

    double *a = new double[N-1]; for (unsigned int i=0; i<N-1; i++) a[i] = alpha; a[0] = 0.0;
    double *b = new double[N-1]; for (unsigned int i=0; i<N-1; i++) b[i] = betta;
    double *c = new double[N-1]; for (unsigned int i=0; i<N-1; i++) c[i] = alpha; c[N-2] = 0.0;
    double *d = new double[N-1]; for (unsigned int i=0; i<N-1; i++) d[i] = +0.0;
    double *x = new double[N-1];
    for (unsigned int i=1; i<=N-1; i++)
    {
        d[i-1] = 2.0*u10[i] - u00[i];
    }
    tomasAlgorithm(a, b, c, d, x, N-1);
    for (unsigned int i=1; i<=N-1; i++) u15[i] = x[i-1];
    IPrinter::print(u00, 10, 10);
    IPrinter::print(u10, 10, 10);
    IPrinter::print(u15, 10, 10);


    //    double hx = 0.1;
    //    double hy = 0.1;
    //    double ht = 0.05;
    //    double alpha = -((ht*ht)/(hx*hx));
    //    double betta = 2.0+2.0*((ht*ht)/(hx*hx));
    //    double gamma = +((ht*ht)/(hy*hy));

    //    unsigned int N = 10;
    //    DoubleMatrix u00((N+1), (N+1), 0.0);
    //    DoubleMatrix u10((N+1), (N+1), 0.0);
    //    DoubleMatrix u15((N+1), (N+1), 0.0);
    //    DoubleMatrix u20((N+1), (N+1), 0.0);
    //    u10[4][5] = 50.0;
    //    u10[6][5] = 50.0;
    //    u10[5][5] = 50.0;
    //    u10[4][5] = 50.0;
    //    u10[6][5] = 50.0;


    //    double *a = new double[N-1]; for (unsigned int i=0; i<N-1; i++) a[i] = alpha; a[0] = 0.0;
    //    double *b = new double[N-1]; for (unsigned int i=0; i<N-1; i++) b[i] = betta;
    //    double *c = new double[N-1]; for (unsigned int i=0; i<N-1; i++) c[i] = alpha; c[N-2] = 0.0;
    //    double *d = new double[N-1]; for (unsigned int i=0; i<N-1; i++) d[i] = +0.0;
    //    double *x = new double[N-1];

    //    IPrinter::print(u10, 10, 10);
    //    IPrinter::printSeperatorLine();
    //    for (unsigned int j=1; j<=N-1; j++)
    //    {
    //        for (unsigned int i=1; i<=N-1; i++)
    //        {
    //            d[i-1] = gamma*(u10[j][i-1] - 2.0*u10[j][i] + u10[j][i+1]) + (2.0*u10[j][i]+u10[j][i]-u00[j][i]);
    //        }
    //        tomasAlgorithm(a, b, c, d, x, N-1);
    //        for (unsigned int i=1; i<=N-1; i++) u15[j][i] = x[i-1];
    //    }
    //    IPrinter::print(u15, 10, 10);
    //    IPrinter::printSeperatorLine();
    //    for (unsigned int i=1; i<=N-1; i++)
    //    {
    //        for (unsigned int j=1; j<=N-1; j++)
    //        {
    //            d[j-1] = gamma*(u15[j][i-1] - 2.0*u15[j][i] + u15[j][i+1]) + (2.0*u15[j][i]+u10[j][i]-u00[j][i]);
    //        }
    //        tomasAlgorithm(a, b, c, d, x, N-1);
    //        for (unsigned int j=1; j<=N-1; j++) u20[j][i] = x[j-1];
    //    }
    //    IPrinter::print(u20, 10, 10);

    return 0;


    //    QGuiApplication app(argc, argv);

    //Problem2Article::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //Problem2PNeumann::Main(argc, argv);
    Problem2HNDirichlet::Main(argc, argv);
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
