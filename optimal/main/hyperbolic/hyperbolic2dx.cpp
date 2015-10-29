#include "hyperbolic2dx.h"

Hyperbolic2DX::Hyperbolic2DX()
{
}

Hyperbolic2DX::~Hyperbolic2DX()
{}

double Hyperbolic2DX::fx(const DoubleVector& x)
{
    return 0.0;
}

void Hyperbolic2DX::gradient(const DoubleVector& x, DoubleVector& g, double gradient_step)
{}

void Hyperbolic2DX::project(DoubleVector &x, int index)
{}

void Hyperbolic2DX::print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn)
{}
