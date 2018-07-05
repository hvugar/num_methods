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

void EXPORTERSHARED_EXPORT init_pr();

double EXPORTERSHARED_EXPORT call_fx(double *x);

void EXPORTERSHARED_EXPORT call_gr(double *x, double *g);

double EXPORTERSHARED_EXPORT add_fx(double x, double y);

#ifdef __cplusplus
}
#endif


#endif // EXPORTER_H
