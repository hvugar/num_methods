#ifndef PROBLEM22DEX5_H
#define PROBLEM22DEX5_H

#include "../abs/abstractproblem22d.h"
#include <imaging.h>

class Problem22DEx5 : public AbstactProblem22D
{
public:
    static void Main(int argc, char* argv[]);

    Problem22DEx5();

    static void Table1Y1();
    static void experiment1();
    static void experiment2();

private:
};

#endif // PROBLEM22DEX5_H
