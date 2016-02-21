#include "headers.h"

#include <libxl.h>

using namespace libxl;

int main()
{
//    DoubleCube c;
//    BorderHyperbolic2D s;

//    s.qamma = 0.0;
//    s.calculateMVD(c, s.h1, s.h2, s.ht, s.N1, s.N2, s.M, s.a1, s.a2);
//    DoubleMatrix u1 = c[c.size()-1];
//    IPrinter::printMatrix(u1);
//    puts("---");

//    s.qamma = 0.2;
//    s.calculateU1(c, s.h1, s.h2, s.ht, s.N1, s.N2, s.M, s.a1, s.a2, s.qamma);
//    DoubleMatrix u2 = c[c.size()-1];
//    IPrinter::printMatrix(u2);

//    puts("---");
//    for (unsigned int j=0; j<=s.N2; j++)
//    {
//        for (unsigned int i=0; i<=s.N1; i++)
//        {
//            u1[j][i] = s.u(i, j, s.M);
//        }
//    }
//    IPrinter::printMatrix(u1);

//    HeatControl2DeltaF::main();
//    HeatControlDeltaX::main();
//    DiscreteHyperbolic1::main();
//    HyperbolicControlH::main();

    HyperbolicControl2DM::main();
    puts("--------------------------------------");
    //HyperbolicControl2DMX::main();
    puts("--------------------------------------");
    //HyperbolicControl2DMV::main();
    puts("--------------------------------------");
    return 0;
}
