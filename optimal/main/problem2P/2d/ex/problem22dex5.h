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
    static void Table1Y2();

    static void experiment1();
    static void experiment2();

    virtual double fx(const DoubleVector &prms) const;
    virtual void gradient(const DoubleVector &prms, DoubleVector &g);

private:
    vector<double> fis;
    vector<double> thetas;
};

#endif // PROBLEM22DEX5_H
