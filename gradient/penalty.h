#ifndef _PENALTY_H_
#define _PENALTY_H_

#include "minimum.h"
#include "methods.h"

/**
 * @brief
 * @param
 * @param
 * @param
 * @param
 * @return
 */
void penalty_method(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int k, double r1, double r2, double epsilon);

#endif // _PENALTY_H_