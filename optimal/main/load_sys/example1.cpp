#include "example1.h"

void Example1::Main(int argc, char *argv[])
{
    Example1 ex1;
    ex1.setGrid(ODEGrid(Dimension(0.1, 10, 0)));

    std::vector<Condition> cs;

    Condition c0;
    c0.time = 0.1;
    c0.mtrx.resize(ex1.equationsNumber(), ex1.equationsNumber());
    c0.mtrx.at(0,0) = 10.2;
    c0.index = 0;
    cs.push_back(c0);

    Condition c1;
    c1.time = 0.403;
    c1.mtrx.resize(ex1.equationsNumber(), ex1.equationsNumber());
    c1.mtrx.at(0,0) = 4.2;
    c1.index = 0;
    cs.push_back(c1);

    Condition c2;
    c2.time = 1.0;
    c2.mtrx.resize(ex1.equationsNumber(), ex1.equationsNumber());
    c2.mtrx.at(0,0) = 3.5;
    c2.index = 0;
    cs.push_back(c2);

    DoubleVector bt(ex1.equationsNumber());

    ex1.highOder2Accuracy(cs, bt);
}

double Example1::A(double x, unsigned int i, unsigned int row, unsigned int col) const
{
    return 0.0;
}

double Example1::B(double x, unsigned int i, unsigned int row) const
{
    return 0.0;
}

unsigned int Example1::equationsNumber() const
{
    return 1.0;
}
