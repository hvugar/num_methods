#ifndef ROSENBROCK_H
#define ROSENBROCK_H

#include "function.h"

class Rosenbrock : public RnFunction
{
public:
    virtual double fx(const std::vector<double>& x);

    static void Main();
};

#endif // ROSENBROCK_H
