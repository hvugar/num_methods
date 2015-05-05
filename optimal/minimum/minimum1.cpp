#include "methods.h"

double straight_line_search_metod(R1Function *f, double x0, double dx, double &a, double &b)
{
    if (f == 0) {
        fputs("function is null", stderr);
        exit(-1);
    }

    if ( dx == 0.0 )
    {
        fputs("dx cannot be zero", stderr);
        exit(-1);
    }

    double y0 = f->fx(x0);
    double y1 = f->fx(x0 - dx);
    double y2 = f->fx(x0 + dx);

    // if y1 and y2 are both greater than y0
    // minimum point is inside x1 and x2
    if (y1 >= y0 && y0 <= y2)
    {
        a = x0 - dx;
        b = x0 + dx;
        double c = (a + b) / 2.0;
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
                y2 = f->fx(x0 + fabs(dx));
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
                y1 = f->fx(x0 - fabs(dx));
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

double golden_section_search_min(R1Function *f, double a, double b, double epsilon)
{
    double phi = 1.6180339887498948482045868343656;
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
            y1 = f->fx(x1);
        }

        if (isnan(x2))
        {
            x2 = a + fabs(b-a)/phi;
            y2 = f->fx(x2);
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

    if (f->fx(a)<f->fx(b)) c = a;
    if (f->fx(a)>f->fx(b)) c = b;

    return c;
}
