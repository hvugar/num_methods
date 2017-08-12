#ifndef PROBLEM5EX1_H
#define PROBLEM5EX1_H

#include <load_sys/islodenlcs.h>
#include <load_sys/islodenlcsm.h>

#include "zetta0.h"
#include "zetta1.h"
#include "zetta2.h"

#define SAMPLE_1

class Problem5Ex1
{
public:
    static void Main(int agrc, char *argv[]);

    Problem5Ex1(const ODEGrid &grid);
    void initialize();

    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int k, unsigned int row) const;
    virtual double C(double t, unsigned int k, unsigned int row, unsigned int col, unsigned int i) const;
    virtual double g(double t, unsigned int k, unsigned int row, unsigned int i) const;

    double X(double t, unsigned int k, unsigned int i) const;
    double dX(double t, unsigned int k, unsigned int i) const;

    // Load points count
    unsigned int L0 = 2;
    // Conditions count
    unsigned int L1 = 3;

private:
    ODEGrid mgrid;
};

#endif // PROBLEM5EX1_H
