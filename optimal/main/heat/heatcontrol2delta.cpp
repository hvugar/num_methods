#include "heatcontrol2delta.h"

HeatControl2Delta::HeatControl2Delta() : RnFunction()
{

}

HeatControl2Delta::~HeatControl2Delta()
{

}

double HeatControl2Delta::fx(const DoubleVector &x)
{
    return 0.0;
}

void HeatControl2Delta::gradient(double step, const DoubleVector &x, DoubleVector &g)
{

}

void HeatControl2Delta::calculateU(const DoubleVector &f)
{}

void HeatControl2Delta::calculateP(const DoubleVector &f, DoubleVector &g)
{}

