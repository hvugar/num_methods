#ifndef _PENALTY_H_
#define _PENALTY_H_

#include "minimum.h"
#include "methods.h"

void penalty_method(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double R);

#endif // _PENALTY_H_