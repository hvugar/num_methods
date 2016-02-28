#include "headers.h"

int main()
{
//    DoubleCube c;
//    BorderHyperbolic bh;
//    bh.calculateU1(c, bh.h1, bh.h2, bh.ht, bh.N1, bh.N2, bh.M, bh.a1, bh.a2, bh.qamma);
//    IPrinter::printMatrix(c[bh.M]);
//    puts("---");
//    DoubleMatrix m;
//    bh.calculateU1(m, bh.h1, bh.h2, bh.ht, bh.N1, bh.N2, bh.M, bh.a1, bh.a2, bh.qamma);
//    IPrinter::printMatrix(m);

//    for (unsigned int i=0; i<=bh.M; i++)
//    {
//        const DoubleMatrix &m = c[i];
//        puts("-----------------");
//        char buffer[20];
//        int n = sprintf(buffer, "data/file%d.txt", i);
//        buffer[n] = '\0';
//        FILE *file = fopen(buffer, "w");
//        IPrinter::printMatrix(m, bh.N2, bh.N1, NULL, file);
//        fclose(file);
//    }

    //HeatControl2DeltaF::main();
    //HeatControlDeltaX::main();
    //DiscreteHyperbolic1::main();
    //HyperbolicControlH::main();
    HyperbolicControl2D1::main();
    //puts("--------------------------------------");
    //HyperbolicControl2DMX::main();
    //puts("--------------------------------------");
    //HyperbolicControl2DMV::main();
    //puts("--------------------------------------");
    return 0;
}
