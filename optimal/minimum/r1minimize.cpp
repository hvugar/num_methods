#include "r1minimize.h"
#include <stdio.h>
#include <math.h>

R1Minimize::R1Minimize()
{
}

R1Minimize::~R1Minimize()
{
}

double R1Minimize::straightLineSearch(double x0, double dx, double &a, double &b)
{
    if ( dx == 0.0 )
    {
        fputs("dx cannot be zero", stderr);
        return NAN;
    }

    double y0 = fx(x0);
    double y1 = fx(x0 - dx);
    double y2 = fx(x0 + dx);

    // if y1 and y2 are both greater than y0 then minimum point is inside x1 and x2
    if (y1 >= y0 && y0 <= y2)
    {
        a = x0 - dx;
        b = x0 + dx;
        double c = (a + b) / 2.0;
        return c;
    }

    // if y1 and y2 are both lesser than y0 then there is not minimum. function is not unimodal
    if (y1 <= y0 && y0 >= y2)
    {
        fputs("Function is not unimodal", stderr);
        return NAN;
    }

    {
        if ( y1 >= y0 && y0 >= y2 )
        {
            while ( y2 < y0 )
            {
                x0 = x0 + fabs(dx);
                y1 = y0;
                y0 = y2;
                y2 = fx(x0 + fabs(dx));
            }
            a = x0 - fabs(dx);
            b = x0 + fabs(dx);
            double c = (a + b) / 2.0;
            return c;
        }

        if ( y1 <= y0 && y0 <= y2 )
        {
            while ( y1 < y0 )
            {
                x0 = x0 - fabs(dx);
                y2 = y0;
                y0 = y1;
                y1 = fx(x0 - fabs(dx));
            }
            a = x0 - fabs(dx);
            b = x0 + fabs(dx);
            double c = (a + b) / 2.0;
            return c;
        }
    }

    a = x0 - fabs(dx);
    b = x0 + fabs(dx);
    double c = (a + b) / 2.0;
    return c;
}
