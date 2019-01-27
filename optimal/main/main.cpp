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

    NonLocalSystem nl;

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
    //C[0].clear();   C[0].resize(M,M,1.0);
    //C[10].clear();  C[5].resize(M,M,2.0);
    //C[20].clear();  C[10].resize(M,M,3.0);

    for (unsigned int j=0; j<3; j++)
    {
        for (unsigned int i=0; i<3; i++)
        {
            C[0][j][i] = Random::value(0,1,4);
            C[N/2][j][i] = Random::value(0,1,4);
            C[N][j][i] = Random::value(0,1,4);
        }
    }

    DoubleVector d;
    d.resize(M);
    DoubleVector x00; x00 << nl.x(0.0,1) << nl.x(0.0,2) << nl.x(0.0,3);
    DoubleVector x05; x05 << nl.x(0.5,1) << nl.x(0.5,2) << nl.x(0.5,3);
    DoubleVector x10; x10 << nl.x(1.0,1) << nl.x(1.0,2) << nl.x(1.0,3);
    d = C[0]*x00 + C[N/2]*x05 + C[N]*x10;
    std::vector<DoubleVector> x;
    nl.solve(C, d, x, Dimension(h, 0, static_cast<int>(N)), M);

    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            if (n%(N/10)==0) printf("%14.10f ", x[n][m]);
        }
        puts("");
    }
    IPrinter::printSeperatorLine();

    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            if (n%(N/10)==0) printf("%14.10f ", nl.x(n*h, m+1));
        }
        puts("");
    }

    double norm[] = {0.0, 0.0, 0.0};
    for (unsigned int m=0; m<M; m++)
    {
        norm[m] = 0.0;
        for (unsigned int n=0; n<=N; n++)
        {
            norm[m] += (x[n][m]-nl.x(n*h, m+1))*(x[n][m]-nl.x(n*h, m+1));
        }
        norm[m] = sqrt(norm[m]);
    }
    printf("Norms: %.10f %.10f %.10f\n", norm[0], norm[1], norm[2]);

//    return 0;


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
