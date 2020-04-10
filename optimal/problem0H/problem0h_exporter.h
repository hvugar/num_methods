#ifndef PROBLEM2H_EXPORTER_H
#define PROBLEM2H_EXPORTER_H

//#include "problem0h_global.h"

#define PROBLEM0H_LIBRARY_ML

#ifdef _WIN32
#if defined(PROBLEM0H_LIBRARY_ML)
#   define PROBLEM0HSHARED_EXPORT_ML __declspec(dllexport)
#else
#   define PROBLEM0HSHARED_EXPORT_ML __declspec(dllimport)
#endif
#else
#   define PROBLEM0HSHARED_EXPORT_ML
#endif


#ifdef __cplusplus
extern "C" {
#endif

void PROBLEM0HSHARED_EXPORT_ML init_problem();

void PROBLEM0HSHARED_EXPORT_ML get_vector_size(double size);

void PROBLEM0HSHARED_EXPORT_ML init_strt_vector(double* x);

double PROBLEM0HSHARED_EXPORT_ML call_fx(double *x);

void PROBLEM0HSHARED_EXPORT_ML call_gr(double *x, double *g, unsigned int size);

void PROBLEM0HSHARED_EXPORT_ML setPenaltyR(double r);

#ifdef __cplusplus
}
#endif

#endif // PROBLEM2H_EXPORTER_H
