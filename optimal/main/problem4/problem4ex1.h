#ifndef PROBLEM5EX1_H
#define PROBLEM5EX1_H

#include <load_sys/islodenlcsv.h>
#include <load_sys/islodenlcsm.h>
#include <function.h>

#include "zetta0.h"
#include "zetta1.h"
#include "zetta2.h"

#include "../nonlinearfunctionsystem.h"

#define SAMPLE_1

class Problem4Ex1 : public NonLinearFunctionSystem
{
public:
    static void Main(int agrc, char *argv[]);

    Problem4Ex1(const ODEGrid &grid);
    void initialize();

    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int k, unsigned int row, unsigned int col, unsigned int i) const;
    virtual double C(double t, unsigned int k, unsigned int row) const;

    virtual double g(unsigned int num, unsigned int row) const;

    virtual double g(const DoubleVector &x, unsigned int num, unsigned int row) const;

    double X(double t, unsigned int num) const;
    double dX(double t, unsigned int num) const;

    unsigned int L1 = 3;

    virtual double fx(const DoubleVector &x, unsigned int num) const;

    void printResult();

    std::vector<DoubleVector> zm0;
    std::vector<std::vector<DoubleVector>> zm1;
    std::vector<std::vector<DoubleVector>> zm2;

private:
    ODEGrid mgrid;
};

#endif // PROBLEM5EX1_H
