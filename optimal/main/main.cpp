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
#include "problem4/problem4ex2.h"

//#include <../border/borderparabolicd.h>
//#include <../border/borderparabolicn.h>
//#include <../border/borderparabolic2d.h>
//#include <../border/borderhyperbolic2d.h>

#include <../border/grid/parabolicibvp1.h>
#include <../border/grid/parabolicibvp2.h>

//#include <../border/grid/hyperbolicibvp1.h>
//#include <../border/grid/newtonheatequationex1.h>

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

#include "nonlinearequationex1.h"
#include "load_sys/lode1oex1.h"

#include "matrixtest.h"
#include <utils/random.h>

#include "loadedlinearode1order.h"

#include "problem5/nllparabolic.h"

#include "problem2/1d/problem2.h"
#include "problem2/2d/problem22d.h"


int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    srand(time(NULL));

    //    unsigned int N = 11;
    //    double *a = (double*) malloc(sizeof(double)*N);
    //    double *b = (double*) malloc(sizeof(double)*N);
    //    double *c = (double*) malloc(sizeof(double)*N);
    //    double *d = (double*) malloc(sizeof(double)*N);
    //    double *x = (double*) malloc(sizeof(double)*N);
    //    double **e = (double**) malloc(sizeof(double*)*N);

    //    for (unsigned int i=0; i<N; i++)
    //    {
    //        e[i] = (double*) malloc(sizeof(double)*N);
    //        for (unsigned int j=0; j<N; j++) e[i][j] = 0.0;
    //    }

    //    for (unsigned int i=0; i<N; i++)
    //    {
    //        a[i] = -1.0;
    //        b[i] = +4.1;
    //        c[i] = -1.0;
    //        x[i] = i;

    //        e[i][0] = +1.0;
    //        e[i][5] = -2.0;
    //        e[i][7] = +3.0;
    //        e[i][10] = -1.0;
    //    }
    //    a[0]   = 0.0;
    //    c[N-1] = 0.0;
    
    //    d[0] = b[0]*x[0] + c[1]*x[1]
    //            + e[0][0]*x[0] + e[0][5]*x[5] + e[0][7]*x[7] + e[0][10]*x[10];
    //    for (unsigned int i=1; i<N-1; i++)
    //    {
    //        d[i] = a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1]
    //                + e[i][0]*x[0] + e[i][5]*x[5] + e[i][7]*x[7] + e[i][10]*x[10];
    //    }
    //    d[N-1] = a[N-1]*x[N-2] + b[N-1]*x[N-1]
    //            + e[N-1][0]*x[0] + e[N-1][5]*x[5] + e[N-1][7]*x[7] + e[N-1][10]*x[10];

    //    IPrinter::printVector(a,N,NULL,N);
    //    IPrinter::printVector(b,N,NULL,N);
    //    IPrinter::printVector(c,N,NULL,N);
    //    IPrinter::printVector(d,N,NULL,N);
    //    IPrinter::printVector(x,N,NULL,N);
    //    IPrinter::printSeperatorLine();
    //    for (unsigned int i=0; i<N; i++)
    //    {
    //        for (unsigned int j=0; j<N; j++)
    //        {
    //            printf("%14.10f ", e[i][j]);
    //        }
    //        puts("");
    //    }
    //    IPrinter::printSeperatorLine();
    //    for (unsigned int i=0; i<N; i++) x[i] = 0.0;
    
    //    LinearEquation::func1(a, b, c, d, e, x, N);

    //    IPrinter::printVector(x,N,NULL,N);

    //func2(a, b, c, d, e, x, N);


    //MatrixTest::Main(argc, argv);

    //Problem2::Main(argc, argv);
    Problem22D::Main(argc, argv);

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
