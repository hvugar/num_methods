#include "hyperbolicequation.h"

HyperbolicEquation::HyperbolicEquation(double t0, double t1, double x0, double x1, double a, unsigned int M, unsigned int N)
    : t0(t0), t1(t1), x0(x0), x1(x1 ), a(a), M(M), N(N)
{
    ht = (t1 - t0) / M;
    hx = (x1 - x0) / N;
}

HyperbolicEquation::~HyperbolicEquation()
{

}

void HyperbolicEquation::setLengthInterval(double x0, double x1)
{
    this->x0 = x0;
    this->x1 = x1;
    hx = (x1 - x0) / N;
}

void HyperbolicEquation::setTimeInterval(double t0, double t1)
{
    this->t0 = t0;
    this->t1 = t1;
    ht = (t1 - t0) / M;
}

void HyperbolicEquation::setPartNumber(unsigned int M, unsigned int N)
{
    this->M = M;
    this->N = N;
    this->ht = (t1-t0)/M;
    this->hx = (x1-x0)/N;
}

void HyperbolicEquation::calculate(DoubleVector &u)
{

}
