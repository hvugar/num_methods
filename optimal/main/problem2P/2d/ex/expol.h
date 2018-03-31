#ifndef EXPOPTIMALLETTERS_H
#define EXPOPTIMALLETTERS_H

#include "../abs/abstractproblem22d.h"
#include "../abs/jfunctional.h"
#include "../abs/pfunctional.h"

class ExpOptimalLetters
{
public:
    static void Main(int argc, char* argv[]);

    ExpOptimalLetters();

    static void Table1Y1();
    static void Table1Y2();

    static void Table2Y1();
    static void Table3Y2();
    static void Table2Y2();

    static void Table4();

    static void test();

    static void figure1();
};

#endif // EXPOPTIMALLETTERS_H
