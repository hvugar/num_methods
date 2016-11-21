#ifndef BORDERTEST1_H
#define BORDERTEST1_H

#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <cmethods.h>

//#define SAMPLE_1
#define SAMPLE_2

class BorderTest1
{
public:
    BorderTest1();
    virtual ~BorderTest1() {}

    double a(unsigned int i) const;
    double b(unsigned int i) const;

    void calculate1(DoubleVector &x);
    void calculate2(DoubleVector &x);
    double X(unsigned int i) const;

    double h = 0.0001;
    unsigned int N = 10000;

    double alpha;
    double betta;

    static void Main(int argc, char* argv[]);
};

#endif // BORDERTEST1_H
