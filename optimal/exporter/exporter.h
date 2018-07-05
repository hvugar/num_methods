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

EXPORTERSHARED_EXPORT double init_pr(double *x);
EXPORTERSHARED_EXPORT double call_fx(double *x);
EXPORTERSHARED_EXPORT void call_gr(double *x, double *g);

#ifdef __cplusplus
}
#endif


#endif // EXPORTER_H
