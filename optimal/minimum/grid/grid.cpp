#include "grid.h"
#include <cstdio>

Dimension::Dimension(double step, unsigned int maxN, unsigned int minN)
    : mstep(step), mmaxN(maxN), mminN(minN)
{}

double Dimension::step() const { return mstep; }

unsigned int Dimension::minN() const { return mminN; }

unsigned int Dimension::maxN() const { return mmaxN; }

unsigned int Dimension::sizeN() const { return mmaxN-mminN; }
