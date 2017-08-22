#ifndef HEATEXAMPLE1_H
#define HEATEXAMPLE1_H

#include <pde/parabolicequation.h>
#include <printer.h>

class HeatExample1 : public IParabolicEquation
{
public:
    HeatExample1();
    virtual ~HeatExample1() {}

    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    static void main(int argc, char ** argv);

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double a;

};

#endif // HEATEXAMPLE1_H
