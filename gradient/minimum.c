#include "minimum.h"
#include <math.h>

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

double straight_line_search_metod1(R1Function f, double x0, double dx, double *a, double *b, int *count)
{
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
		y1 = y2;
	}
	
	while ( y2 < y0 )
	{
		x0 = x0 + dx;
		y1 = y0;
		y0 = y2;
		y2 = f(x0 + h);
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
	
	return (a+b)/2.0;
}

double straight_line_search_metod1(R1Function fx, double x0, double dx, double *a, double *b, int *count)
{
	double y0 = 0.0;
	double y1 = 0.0;
	double y2 = 0.0;
	
	// if at next and last point of x0 function is greater
	// then decrease the dx to half
	y0 = fx( x0 );
//	(*count)++;
	y1 = fx( x0 - dx );
//	(*count)++;
	y2 = fx( x0 + dx );
//	(*count)++;
	
	while (y1 > y0 && y2 > y0)
	{
		dx = dx / 2.0;
		y1 = fx( x0 - dx );
		y2 = fx( x0 + dx );
	}
	
	if ( y2 > y0 )
	{
		dx = -1 * dx;
		y2 = fx( x0 + dx );
	}

    double x1 = x0;
    double x2 = x1 + dx;
	
    y1 = fx( x1 );
	//(*count)++;
    y2 = fx( x2 );
	//(*count)++;
	
    int i = 1;
    while ( y2 <= y1 )
    {
		x1 = x2;
        y1 = y2;
		
		i++;
        x2 = x0 + i * dx;
        y2 = fx(x2);

		(*count)++;
		
		dx = dx*2;
    }
	
	if ( dx < 0 )
	{
		*a = x0 + (i) * dx;
		*b = x0 + (i-2) * dx;	
	}
	
	
	if ( dx > 0 )
	{
	    *a = x0 + (i-2) * dx;
		*b = x0 + (i+0) * dx;
	}
	
	return (*a+*b)/2;
}
