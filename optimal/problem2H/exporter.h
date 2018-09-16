#ifndef EXPORTER_H
#define EXPORTER_H

#if defined(PROBLEM2H_LIBRARY)
#  define EXPORTERSHARED_EXPORT __declspec(dllexport)
#else
#  define EXPORTERSHARED_EXPORT __declspec(dllimport)
#endif

#ifdef __cplusplus
extern "C" {
#endif

void EXPORTERSHARED_EXPORT init_pr();

double EXPORTERSHARED_EXPORT call_fx(double *x);

void EXPORTERSHARED_EXPORT call_gr(double *x, double *g, unsigned int size);

void EXPORTERSHARED_EXPORT setPenaltyR(double r);

#ifdef __cplusplus
}
#endif


#endif // EXPORTER_H
