#include "grid.h"
#include <cstdio>

Dimension::Dimension(double step, unsigned int maxN, unsigned int minN)
    : mstep(step), mmaxN(maxN), mminN(minN)
{

    union {
       double d;
       unsigned int i;
    } a;

    a.d = mstep;
    printf("1***************** %u %f\n", a.i, a.d);
//    mmin = mminN*mstep;
//    mmax = mmaxN*mstep;
}

double Dimension::step() const { return mstep; }

//double Dimension::min() const { return mmin; }

//double Dimension::max() const { return mmax; }

unsigned int Dimension::minN() const { return mminN; }

unsigned int Dimension::maxN() const { return mmaxN; }

unsigned int Dimension::sizeN() const { return mmaxN-mminN; }
