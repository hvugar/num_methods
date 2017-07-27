#ifndef PROBLEM4_H
#define PROBLEM4_H

#include <grid/grid.h>
#include <vector2d.h>

class Problem4
{
public:
    Problem4(const Dimension &time);

    double X(unsigned int k) const;

    void calculate(DoubleVector &x);

    Dimension mtime;

    double a(const TimeNode &tn) const;
    double b(const TimeNode &tn) const;
    double b1(const TimeNode &tn) const;
    double b2(const TimeNode &tn) const;
};

#endif // PROBLEM4_H
