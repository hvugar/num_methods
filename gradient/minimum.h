#ifndef __MINIMUM_H
#define __MINIMUM_H

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*R1Function)(double);
typedef double (*RnFunction)(double*, int);

/**
 * @brief
 * @param fx
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double straight_line_search_metod(R1Function fx, double x0, double dx, double *a, double *b, int *count);

/**
 * @brief
 * @param fx
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double golden_section_search_min(R1Function fx, double a, double b, double epsilon, int *count);

/**
 * @brief
 * @param fx
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double search_method_dck(R1Function f, double x0, double dx, double *a, double *b);

/**
 * @brief
 * @param fx
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double search_method_pauella(R1Function f, double x0, double dx, double *a, double *b);

/**
 * @brief Этап установления границ интервала. Метод Свенна
 * @param f     function
 * @param x0    initial point
 * @param dx
 * @param a
 * @param b
 * @return
 */
void search_interval_svenn(R1Function f, double x0, double dx, double *a, double *b);

#ifdef __cplusplus
}
#endif


#endif // __MINIMUM_H