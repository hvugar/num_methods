#include "minimum.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/**
 * @brief Метод золотого сечения
 * @param f
 * @param a
 * @param b
 * @param epsilon
 * @return
 */
double golden_section_search_min(R1Function f, double a, double b, double epsilon)
{
    double phi = (1 + sqrt(5)) / 2;

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
            x1 = b - (b-a)/phi;
            y1 = f(x1);
        }

        if (isnan(x2))
        {
            x2 = a + (b-a)/phi;
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

    return (a+b)/2;
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

    if (y1 >= y0 && y2 >= y0)
    {
        *a = x0 - dx;
        *b = x0 + dx;
    }

    if ( y2 > y0 )
    {
        dx = -dx;
        double ty = y2;
        y2 = y1;
        y1 = ty;
    }

    while ( y2 < y0 )
    {
        x0 = x0 + dx;
        y1 = y0;
        y0 = y2;
        y2 = f(x0 + dx);
    }

    if ( dx > 0 )
    {
        *a = x0 - dx;
        *b = x0 + dx;
    }
    else
    {
        *a = x0 + dx;
        *b = x0 - dx;
    }

    return (*a+*b)/2.0;
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
    double x1 = x0;
    double x2;
    double x3;
	double xt;
    double xm;
	
    double y1;
	double y2;
	double y3;
    double ym;
	double yt;

    do
    {
        x2 = x1 + dx;
        y1 = f(x1);
        y2 = f(x2);

        if ( y1 > y2 )
        {
            x3 = x1 + 2*dx;
            y3 = f(x3);
        }
        else
        {
            x3 = x1 - dx;
            y3 = f(x3);
        }

        if ( y1 < y2 )
        {
            if ( y1 < y3 )
            {
                xm = x1;
                ym = y1;
            }
            else
            {
                xm = x3;
                ym = y3;
            }
        }
        else
        {
            if ( y2 < y3 )
            {
                xm = x2;
                ym = y2;
            }
            else
            {
                xm = x3;
                ym = y3;
            }
        }
		
		double a1 = (y2-y1)/(x2-x1);
		double a2 = ((y3-y1)/(x3-x1) - (y2-y1)/(x2-x1))/(x3-x2);
		double xt = (x2+x1)/2.0 - (a1/(2.0*a2));
		yt = f(xt);
		
		printf("%f %f %f %f %f %f %f %f %f %f %f %f %f\n", x1, x2, x3, y1, y2, y3, xm, ym, a1, a2, xt, yt, (ym - yt)/yt);
		//break;
    } while ( (ym - yt)/yt > epsilon && (xm - xt)/xt > epsilon );
}

/**
 * @brief Этап установления границ интервала. Метод Свенна
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
void search_interval_svenn(R1Function f, double x0, double dx, double *a, double *b)
{
    double y1 = f(x0 - dx);
    double y0 = f(x0);
    double y2 = f(x0 + dx);

    /* if y1 and y2 are both lesser than y0 then function is not unimodal */
    if (y0 >= y1 && y0 >= y2)
    {
        *a = 0.0;
        *b = 0.0;
        fputs("Function is not unimodal.\n", stderr);
        return;
    }

    /* if y1 and y2 are both greater than y0 then x0 is minimum */
    if (y0 <= y1 && y0 <= y2)
    {
        *a = x0 - dx;
        *b = x0 + dx;
        return;
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
 * @brief Метод деления интервала пополам
 * @param f
 * @param epsilon
 * @param a
 * @param b
 * @return
 */
void halph_interval_method(R1Function f, double epsilon, double *a, double *b)
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
}
