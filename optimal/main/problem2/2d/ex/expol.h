#ifndef EXPOPTIMALLETTERS_H
#define EXPOPTIMALLETTERS_H

#include "../abs/abstractproblem22d.h"
#include "../abs/jfunctional.h"

class ExpOptimalLetters
{
public:
    static void Main(int argc, char* argv[]);

    ExpOptimalLetters();

    static void Table1Y1();
    static void Table1Y2();

    void optimization(DoubleVector &prm0);
};

#endif // EXPOPTIMALLETTERS_H
