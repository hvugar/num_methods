#ifndef __MINIMUM_H
#define __MINIMUM_H

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*R1Function)(double);
typedef double (*RnFunction)(double*, int);

double golden_section_search_min(R1Function fx, double a, double b, double epsilon);
double straight_line_search_metod(R1Function fx, double x0, double dx, double *a, double *b);

#ifdef __cplusplus
}
#endif


#endif // __MINIMUM_H