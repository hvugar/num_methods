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
#include "../problem2H/problem2h_solver_delta.h"
#include "ivp/nlode1oex1.h"

#include "../problem1H/problem1h_example.h"

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
#include "linearode1storderex1.h"

class A : public NonLinearODE1stOrder
{
public:
    virtual double f(double x, double y, unsigned int k) const;
};

double A::f(double x, double y, unsigned int k) const
{
    double sigma = 0.01;
//    return  3.0*x*x;// + 4.0*(1.0/(sqrt(2.0*M_PI*sigma*sigma))) * exp(((x-0.5)*(x-0.5))/(-2.0*sigma*sigma));
//    return  4.0*x*x*x*exp(pow(x,4.0));// + 4.0*(1.0/(sqrt(2.0*M_PI*sigma*sigma))) * exp(((x-0.8)*(x-0.8))/(-2.0*sigma*sigma));

    double w = 0.0;
    if (k==80) w = 100.0;
    return  4.0*x*x*x*exp(pow(x,4.0)) + 4.0*w;
}


int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
//    A a;
//    const unsigned int N = 100;
//    DoubleVector ry; ry.resize(N+1);
//    a.setDimension(Dimension(0.01, 0, N));
//    a.cauchyProblem(0.0, 0.0, ry, OrdinaryDifferentialEquation::OdeSolverMethod::EULER);
//    IPrinter::print(ry,ry.length());

    LinearODE1stOrderEx1::Main(argc, argv);

    return 0;

    //srand(static_cast<unsigned int>(time(nullptr)));

    //NonLocalSystem::Main(argc, argv);

    //CcIHyperbolicIBVP1::Main(argc, argv);
    //CdIHyperbolicIBVP1::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //ConjugateCC1IHyperbolicIBVP1::Main(argc, argv);
    //Problem1HDirichlet1::Main(argc, argv);
    //Problem2HDirichlet::Main(argc, argv);
    Problem2HDirichletDelta::Main(argc, argv);

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
