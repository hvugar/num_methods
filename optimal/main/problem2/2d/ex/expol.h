#ifndef EXPOPTIMALLETTERS_H
#define EXPOPTIMALLETTERS_H

#include "../abs/abstractproblem22d.h"

class ExpOptimalLetters : public AbstactProblem22D
{
public:
    static void Main(int argc, char* argv[]);

    ExpOptimalLetters();

    void table1();
};

#endif // EXPOPTIMALLETTERS_H
