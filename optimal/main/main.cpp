#include "headers.h"

int main()
{
    DoubleMatrix u;
    BorderHyperbolic2D s;
    s.calculateU(u, s.h1, s.h2, s.ht, s.N1, s.N2, s.M);
    IPrinter::printMatrix(u);
    puts("---");
    for (unsigned int j=0; j<=s.N2; j++)
    {
        for (unsigned int i=0; i<=s.N1; i++)
        {
            u[j][i] = s.u(i, j, s.M);
        }
    }
    IPrinter::printMatrix(u);

//    HeatControl2DeltaX::main();
//    HeatControlDeltaX::main();
//    DiscreteHyperbolic1::main();
//    HyperbolicControlH::main();

//    HyperbolicControl2DMX::main();

    return 0;
}
