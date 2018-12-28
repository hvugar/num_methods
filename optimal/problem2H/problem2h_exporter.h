#ifndef PROBLEM2H_EXPORTER_H
#define PROBLEM2H_EXPORTER_H

#include "problem2h_global.h"

#ifdef __cplusplus
extern "C" {
#endif

void PROBLEM2HSHARED_EXPORT init_pr();

double PROBLEM2HSHARED_EXPORT call_fx(double *x);

void PROBLEM2HSHARED_EXPORT call_gr(double *x, double *g, unsigned int size);

void PROBLEM2HSHARED_EXPORT setPenaltyR(double r);

#ifdef __cplusplus
}
#endif

#endif // PROBLEM2H_EXPORTER_H
