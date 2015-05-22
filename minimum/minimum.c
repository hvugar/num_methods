#include "methods.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

/**
 * @brief         Методы прямого поиска
 * @param f       Целевая функция
 * @param x0      Произвольно выбранная начальная точка
 * @param dx      Величина шага
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @return
 */
double straight_line_search_metod(R1Function f, double x0, double dx, double *a, double *b)
{
    if ( dx == 0.0 )
    {
        fputs("dx cannot be zero", stderr);
        exit(-1);
    }

    double y0 = f(x0);
    double y1 = f(x0 - dx);
    double y2 = f(x0 + dx);

    // if y1 and y2 are both greater than y0
    // minimum point is inside x1 and x2
    if (y1 >= y0 && y0 <= y2)
    {
        *a = x0 - dx;
        *b = x0 + dx;
        double c = (*a + *b) / 2.0;
        return c;
    }

    // if y1 and y2 are both lesser than y0
    // in this case there is not minimum
    // function is not unimodal
    if (y1 <= y0 && y0 >= y2)
    {
        fputs("Function is not unimodal", stderr);
        exit(-1);
    }

    {
        if ( y1 >= y0 && y0 >= y2 )
        {
            while ( y2 < y0 )
            {
                x0 = x0 + dx;
                y1 = y0;
                y0 = y2;
                y2 = f(x0 + fabs(dx));
            }
            *a = x0 - fabs(dx);
            *b = x0 + fabs(dx);
            double c = (*a + *b) / 2.0;
            return c;
        }
		
		if ( y1 <= y0 && y0 <= y2 )
        {
            while ( y1 < y0 )
            {
                x0 = x0 - fabs(dx);
                y2 = y0;
                y0 = y1;
                y1 = f(x0 - fabs(dx));
            }
            *a = x0 - fabs(dx);
            *b = x0 + fabs(dx);
            double c = (*a + *b) / 2.0;
            return c;
        }
    }

    /*
    if ( y1 <= y0 && y0 <= y2 )
    {
        dx = -1.0 * dx;
        y0 = f(x0);
        y1 = f(x0 - dx);
        y2 = f(x0 + dx);
    }

    while ( y2 < y0 )
    {
        x0 = x0 + dx;
        y1 = y0;
        y0 = y2;
        y2 = f(x0 + dx);
    }
*/

    *a = x0 - fabs(dx);
    *b = x0 + fabs(dx);
    double c = (*a + *b) / 2.0;
    return c;
}

/**
 * @brief           Метод Свенна. Установления границ интервала.
 * @param f         Целевая функция
 * @param x0        Произвольно выбранная начальная точка
 * @param dx        Величина шага
 * @param a         Начальная точка отрезка
 * @param b         Конечнная точка отрезка
 * @return
 */
void search_interval_svenn(R1Function f, double x0, double dx, double *a, double *b)
{
    if ( dx == 0.0 )
    {
        fputs("dx cannot be zero", stderr);
        exit(-1);
    }
	
    double y0 = f(x0);
    double y1 = f(x0 - dx);
    double y2 = f(x0 + dx);

    // if y1 and y2 are both greater than y0
    // minimum point is inside x1 and x2
    if (y1 >= y0 && y0 <= y2)
    {
        *a = x0 - dx;
        *b = x0 + dx;
        return;
    }
	
    // if y1 and y2 are both lesser than y0
    // in this case there is not minimum
    // function is not unimodal
    if (y1 <= y0 && y0 >= y2)
    {
        fputs("Function is not unimodal", stderr);
		exit(-1);
    }

	{
		if ( y1 >= y0 && y0 >= y2 )
		{
		}
	}
	
    /* if in point x1 = x0 + dx value of function is greater than value of in point x0
       then dx must negativ
    */
    if ( y1 < y0 && y0 < y2 )
    {
        dx *= -1;
    }

    int k = 0;
    x0 = x0 + dx;
    double y = f(x0);
    printf("k=%d %8.2f %8.2f %8.2f\n", k, x0, y0, y);

    while ( y <= y0 )
    {
        y0 = y;

        k = k + 1;
        x0 = x0 + pow(2, k)*dx;
        y = f(x0);

        printf("k=%d %8.2f %8.2f %8.2f\n", k, x0, y0, y);
    }

    *a = x0 - pow(2, k)*dx - pow(2, k-1)*dx;
    *b = x0;
}

/**
 * @brief         Метод равномерного поиска
 * @param f       Целевая функция
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param n       Количество вычислений функции
 */
void uniform_line_search_method(R1Function f, double *a, double *b, int n)
{
    double h = ((*a) - (*b)) / (n+1);
    int i;
    double x[n];
    double y[n];
    for (i=0; i<n; i++)
    {
        x[i] = (*a) + (i+1) * h;
        y[i] = f(x[i]);
    }
    double min = y[0];
    int k=0;
    for (i=1; i<n; i++)
    {
        if (y[i] < min) { k = i; }
    }
    if (i==0) { *b = x[1]; } else
        if (i==(k-1)) { *a = x[k-2]; } else
        { *a = x[k-1]; *b = x[k+1]; }

}

/**
 * @brief         Метод перебора
 * @param f       Целевая функция
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param n       Количество вычислений функции
 */
void bruteforce_line_search_method1(R1Function f, double *a, double *b, int n)
{
    double h = ((*a) - (*b)) / (n+1);
    int i;
    double x[n+2];
    double y[n+2];
    for (i=0; i<n+2; i++)
    {
        x[i] = (*a) + i * h;
        y[i] = f(x[i]);
    }
    double min = y[0];
    int k=0;
    for (i=1; i<n+2; i++)
    {
        if (y[i] < min) { k = i; }
    }
    if (i==0) { *b = x[1]; } else
        if (i==(k-1)) { *a = x[k-2]; } else
        { *a = x[k-1]; *b = x[k+1]; }
}

/**
 * @brief         Метод золотого сечения
 * @param f       Целевая функция
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param epsilon
 * @return
 */
double golden_section_search_min(R1Function f, double a, double b, double epsilon)
{
    //double sqrt_5 = 2.2360679774997896964091736687313
    //double phi = (sqrt(5) + 1.0) / 2.0;
    double phi = 1.6180339887498948482045868343656;
    //double phi = (sqrt(5) - 1) / 2.0;
    //double phi = 1.6180339887498948482045868343656;

    double x1 = NAN;
    double x2 = NAN;

    double y1 = 0.0;
    double y2 = 0.0;

    // Lazimi epsilon deqiqliyini alana qeder
    // iterasiyalari davam edirik
    while ( fabs(b-a) > epsilon )
    {
        if (isnan(x1))
        {
            x1 = b - fabs(b-a)/phi;
            y1 = f(x1);
        }

        if (isnan(x2))
        {
            x2 = a + fabs(b-a)/phi;
            y2 = f(x2);
        }

        if (y1 >= y2)
        {
            a = x1;
            x1 = x2;    // Tapilmish x2 noqtesi ve bu noqtede funkisiyanin qiymeti
            y1 = y2;    // sonraki iterasiyada x1 qiymeti kimi istifade olunacaq.
            x2 = NAN;   // x2 novbeti iterasiyada axtarilacaq
        }
        else
        {
            b = x2;
            x2 = x1;    // Tapilmish x1 noqtesi ve bu noqtede funkisiyanin qiymeti
            y2 = y1;    // sonraki iterasiyada x2 qiymeti kimi istifade olunacaq.
            x1 = NAN;   // x1 novbeti iterasiyada axtarilacaq
        }
    }

    double c = (a+b)/2.0;

    if (f(a)<f(b)) c = a;
    if (f(a)>f(b)) c = b;

    return c;
}

/**
 * @brief         Метод деления интервала пополам
 * @param f       Целевая функция
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param epsilon
 * @return
 */
double halph_interval_method(R1Function f, double *a, double *b, double epsilon)
{
    double L = *b - *a;

    while ( L > epsilon )
    {
        double xm = (*a+*b)/2.0;
        double x1 = *a + L/4.0;
        double x2 = *b - L/4.0;

        double f_xm = f(xm);
        double f_x1 = f(x1);
        double f_x2 = f(x2);

        if (f_x1 < f_xm)
        {
            *b = xm;
            xm = x1;
        }
        else
        {
            if ( f_x2 < f_xm )
            {
                *a = xm;
                xm = x2;
            }
            else
            {
                *a = x1;
                *b = x2;
            }
        }
        L = *b - *a;
    }
	
	return L;
}

/**
 *
 *
 *
**/
double R1Minimize(R1Function f, double line_step, double gold_epsilon)
{
    double a,b;
    double alpha0 = 0.0;
    straight_line_search_metod(f, alpha0, line_step, &a, &b);
    double alpha = golden_section_search_min(f, a, b, gold_epsilon);
    if ( f(alpha) > f(alpha0) ) alpha = alpha0;
    return alpha;
}

/**
 * @brief
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double search_method_dck(R1Function f, double x0, double dx, double *a, double *b)
{ 
    return 0;
}

/**
 * @brief
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double search_method_pauella(R1Function f, double x0, double dx, double epsilon, double *a, double *b) 
{
    /*
    double x1, x2, x3, xm, xt;
    double y1, y2, y3, ym, yt;

    x1 = x0;
    y1 = f(x1);
    x2 = x1 + dx;
    y2 = f(x2);
    if (y1 > y2) {
        x3 = x1 + dx*2.0;
    } else {
        x3 = x1 - dx;
    }

    do {
        if ( y1 < y2 ) {
            if ( y1 < y3 ) {
                xm = x1;
                ym = y1;
            } else {
                xm = x3;
                ym = y3;
            }
        } else {
            if ( y2 < y3 ) {
                xm = x2;
                ym = y2;
            } else {
                xm = x3;
                ym = y3;
            }
        }

        double a1 = (y2-y1)/(x2-x1);
        double a2 = ((y3-y1)/(x3-x1) - (y2-y1)/(x2-x1))/(x3-x2);
        double xt = (x2+x1)/2.0 - (a1/(2.0*a2));
        yt = f(xt);

        double dist_1 = xt - x1;
        double dist_2 = xt - x2;
        double dist_3 = xt - x3;

        if (fabs(dist_1) < fabs(dist_2) && fabs(dist_1) < fabs(dist_3))
            if (dist_1 < 0)
            {

            }

        printf("%f %f %f %f %f %f %f %f %f %f %f %f %f\n", x1, x2, x3, y1, y2, y3, xm, ym, a1, a2, xt, yt, (ym - yt)/yt);
        //break;
    } while ( fabs((ym - yt)/yt) > epsilon && fabs((xm - xt)/xt) > epsilon ) ;
*/
    return 0.0;
}

double derivative_1(R1Function f, double x, double h)
{
    return (f(x+h)-f(x-h)) / (2*h);
}

double derivative_2(R1Function f, double x, double h)
{
    return (f(x+2*h)-2*f(x)+f(x-2*h)) / (4*h*h);
}

/**
 * @brief Метод Ньютона - Рафсона
 * @param f
 * @param x0
 * @param epsilon
 * @return
 */
double newton_raphson(R1Function f, double x0, double epsilon)
{
    double dx = 0.00001;
    double f_1 = derivative_1(f, x0, dx);
    double f_2 = derivative_2(f, x0, dx);

    double x = x0 - f_1/f_2;
    printf("%f %f %f\n", f_1, f_2, x);

    while ( fabs(f_1) > epsilon ) {
        f_1 = derivative_1(f, x, dx);
        f_2 = derivative_2(f, x, dx);
        x = x - f_1/f_2;
        printf("%f %f %f\n", f_1, f_2, x);
    }

    return 0;
}
