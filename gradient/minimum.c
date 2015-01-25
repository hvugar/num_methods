#include "minimum.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double golden_section_search_min(R1Function fx, double a, double b, double epsilon, int *count)
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
             y1 = fx(x1);
			 (*count)++;
        }

        if (isnan(x2))
        {
            x2 = a + (b-a)/phi;
            y2 = fx(x2);
			(*count)++;
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

double straight_line_search_metod(R1Function f, double x0, double dx, double *a, double *b, int *count)
{
	if ( dx == 0.0 )
	{
		fputs("dx cannot be zero", stderr);
		exit(-1);
	}

	double y0 = f(x0);
	double y1 = f(x0 - dx);
	double y2 = f(x0 + dx);
	(*count) += 2;
	
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
		(*count) += 1;
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
