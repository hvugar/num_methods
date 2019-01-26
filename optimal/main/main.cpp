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
#include "nonlocal.h"

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    //srand(static_cast<unsigned int>(time(nullptr)));

    NonLocal nl;

    //DoubleVector a(101, 0.0);
    //a[0] = a[50] = a[100] = 1.0;
    //double b = a[0]*nl.x(0.0) + a[50]*nl.x(0.5) + a[100]*nl.x(1.0);
    //nl.solve(a, b, 0.01);

    unsigned int N = 100;
    unsigned int M = 3;
    double h = 0.01;
    std::vector<DoubleMatrix> C;
    C.resize(N+1);
    for (unsigned int n=0; n<=N; n++) C[n].resize(M,M,0.0);
    C[0].clear();   C[0].resize(M,M,1.0);
    C[50].clear();  C[50].resize(M,M,2.0);
    C[100].clear(); C[100].resize(M,M,3.0);

    for (unsigned int j=0; j<3; j++)
    {
        for (unsigned int i=0; i<3; i++)
        {
            C[0][j][i] = Random::value(0,1,4);
            C[50][j][i] = Random::value(0,1,4);
            C[100][j][i] = Random::value(0,1,4);
        }
    }

    DoubleVector d;
    d.resize(M);
    d[0] = (C[0])[0][0]*nl.x(0.0,1)   + (C[0])[0][1]*nl.x(0.0,2)   + (C[0])[0][2]*nl.x(0.0,3)+
            (C[50])[0][0]*nl.x(0.5,1)  + (C[50])[0][1]*nl.x(0.5,2)  + (C[50])[0][2]*nl.x(0.5,3)+
            (C[100])[0][0]*nl.x(1.0,1) + (C[100])[0][1]*nl.x(1.0,2) + (C[100])[0][2]*nl.x(1.0,3);

    d[1] = (C[0])[1][0]*nl.x(0.0,1)   + (C[0])[1][1]*nl.x(0.0,2)   + (C[0])[1][2]*nl.x(0.0,3)+
            (C[50])[1][0]*nl.x(0.5,1)  + (C[50])[1][1]*nl.x(0.5,2)  + (C[50])[1][2]*nl.x(0.5,3)+
            (C[100])[1][0]*nl.x(1.0,1) + (C[100])[1][1]*nl.x(1.0,2) + (C[100])[1][2]*nl.x(1.0,3);

    d[2] = (C[0])[2][0]*nl.x(0.0,1)   + (C[0])[2][1]*nl.x(0.0,2)   + (C[0])[2][2]*nl.x(0.0,3)+
            (C[50])[2][0]*nl.x(0.5,1)  + (C[50])[2][1]*nl.x(0.5,2)  + (C[50])[2][2]*nl.x(0.5,3)+
            (C[100])[2][0]*nl.x(1.0,1) + (C[100])[2][1]*nl.x(1.0,2) + (C[100])[2][2]*nl.x(1.0,3);

    printf("%f %f %f\n", d[0], d[1], d[2]);

    nl.solveSystem(C, d, h, M);

    return 0;


    //CcIHyperbolicIBVP1::Main(argc, argv);
    //CdIHyperbolicIBVP1::Main(argc, argv);
    //IPrinter::printSeperatorLine();
    //ConjugateCC1IHyperbolicIBVP1::Main(argc, argv);
    //Problem1HDirichlet1::Main(argc, argv);
    //Problem2HDirichlet::Main(argc, argv);

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
