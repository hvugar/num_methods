#ifndef R1MINIMIZE_H
#define R1MINIMIZE_H

#include "function.h"

void stranghLineSearch(double x, double step, double &a, double &b, R1Function *fn);
double goldenSectionSearch(double &a, double &b, double &x, R1Function *f, double epsilon);

#endif // R1MINIMIZE_H
