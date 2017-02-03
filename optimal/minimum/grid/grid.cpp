#include "grid.h"

Dimension::Dimension(double step, unsigned int maxN, unsigned int minN)
    : mstep(step), mmaxN(maxN), mminN(minN)
{
//    mmin = mminN*mstep;
//    mmax = mmaxN*mstep;
}

double Dimension::step() const { return mstep; }

//double Dimension::min() const { return mmin; }

//double Dimension::max() const { return mmax; }

unsigned int Dimension::minN() const { return mminN; }

unsigned int Dimension::maxN() const { return mmaxN; }

unsigned int Dimension::sizeN() const { return mmaxN-mminN; }
