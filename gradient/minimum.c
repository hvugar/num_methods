#include "minimum.h"
#include <math.h>

double golden_section_search_min(R1Function fx, double a, double b, double epsilon)
{
    double phi = (1 + sqrt(5)) / 2;

    double x1 = NAN;
    double x2 = NAN;

    double y1 = 0;
    double y2 = 0;

    // Lazimi epsilon deqiqliyini alana qeder
    // iterasiyalari davam edirik
    while ( fabs(b-a) > epsilon )
    {
        if (isnan(x1))
        {
             x1 = b - (b-a)/phi;
             y1 = fx(x1);
        }

        if (isnan(x2))
        {
            x2 = a + (b-a)/phi;
            y2 = fx(x2);
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

void straight_line_search_metod(R1Function fx, double x0, double dx, double *a, double *b)
{
    int i = 1;
    double x1 = x0;
    double x2 = x0 + dx;
    double y1 = fx( x1 );
    double y2 = fx( x2 );
	
//    if ( y2 > y1 )
//	{
//        dx = -1 * dx;
//	}
	
//	x2 = x0 + dx;
//	y2 = fx(x2);
	
    while ( y2 <= y1 )
    {
        i++;
        x1 = x2;
        x2 = x0 + i * dx;
        y1 = y2;
        y2 = fx(x2);
    }

    *a = x1 - dx;
    *b = x2 - dx;
}
