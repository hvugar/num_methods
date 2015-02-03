#ifndef __MINIMUM_H
#define __MINIMUM_H

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*R1Function)(double);
typedef double (*RnFunction)(double*, int);

/**
 * @brief Метод золотого сечения
 * @param f
 * @param a
 * @param b
 * @param epsilon
 * @return
 */
double straight_line_search_metod(R1Function fx, double x0, double dx, double *a, double *b);

/**
 * @brief
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double golden_section_search_min(R1Function fx, double a, double b, double epsilon);

/**
 * @brief
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double search_method_dck(R1Function f, double x0, double dx, double *a, double *b);

/**
 * @brief
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double search_method_pauella(R1Function f, double x0, double dx, double epsilon, double *a, double *b);

/**
 * @brief Этап установления границ интервала. Метод Свенна
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
void search_interval_svenn(R1Function f, double x0, double dx, double *a, double *b);

/**
 * @brief Метод деления интервала пополам
 * @param f
 * @param epsilon
 * @param a
 * @param b
 * @return
 */
void halph_interval_method(R1Function f, double epsilon, double *a, double *b);

/**
 * @brief Метод Ньютона - Рафсона
 * @param f
 * @param x0
 * @param epsilon
 * @return
 */
double newton_raphson(R1Function f, double x0, double epsilon);

double derivative_1(R1Function f, double x, double h);

double derivative_2(R1Function f, double x, double h);

#ifdef __cplusplus
}
#endif


#endif // __MINIMUM_H