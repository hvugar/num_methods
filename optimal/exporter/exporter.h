#ifndef EXPORTER_H
#define EXPORTER_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(EXPORTER_LIBRARY)
#  define EXPORTERSHARED_EXPORT __declspec(dllexport)
#else
#  define EXPORTERSHARED_EXPORT __declspec(dllimport)
#endif

double EXPORTERSHARED_EXPORT init_pr(double *x);
double EXPORTERSHARED_EXPORT call_fx(double *x);
double EXPORTERSHARED_EXPORT call_gr(double *x);

#ifdef __cplusplus
}
#endif


#endif // EXPORTER_H
