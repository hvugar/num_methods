#include "r1minimize.h"
#include <math.h>

#include <stdexcept>
#include <string>

void stranghLineSearch(double x, double step, double &a, double &b, R1Function *f)
{
    if ( f == NULL )
    {
        std::string msg = "in function \"stranghLineSearch\" function pointer is null.";
        throw std::runtime_error(msg);
    }

    if ( step <= 0.0 )
    {
        std::string msg = "in function \"stranghLineSearch\" step value is less than zero.";
        throw std::invalid_argument(msg);
    }

    double fstep = fabs(step);

    double y0 = f->fx(x);
    double y1 = f->fx(x - fstep);
    double y2 = f->fx(x + fstep);

    // if y1 and y2 are both greater than y0 then minimum point is inside x1 and x2
    if (y1 >= y0 && y0 <= y2)
    {
        a = x - fstep;
        b = x + fstep;
        return;
    }

    // if y1 and y2 are both lesser than y0 then there is not minimum. function is not unimodal
    if (y1 <= y0 && y0 >= y2)
    {
        a = b = NAN;
        fputs("Function is not unimodal\n", stderr);
        return;
    }

    if ( y1 >= y0 && y0 >= y2 )
    {
        while ( y2 < y0 )
        {
            x = x + fstep;
            y1 = y0;
            y0 = y2;
            y2 = f->fx(x + fstep);
        }
        a = x - fstep;
        b = x + fstep;
    }

    if ( y1 <= y0 && y0 <= y2 )
    {
        while ( y1 < y0 )
        {
            x = x - fstep;
            y2 = y0;
            y0 = y1;
            y1 = f->fx(x - fstep);
        }
        a = x - fstep;
        b = x + fstep;
    }

    a = x - fstep;
    b = x + fstep;
}

double goldenSectionSearch(double &a, double &b, double &x, R1Function *f, double epsilon)
{
    if ( f == NULL )
    {
        std::string msg = "in function \"goldenSectionSearch\" function pointer is null.";
        throw std::runtime_error(msg);
    }

    if ( epsilon <= 0.0 )
    {
        std::string msg = "in function \"goldenSectionSearch\" epsilon value is equal to or less than zero.";
        throw std::invalid_argument(msg);
    }

    double phi = 1.6180339887498948482045868343656;

    double x1 = NAN;
    double x2 = NAN;

    double y1 = 0.0;
    double y2 = 0.0;

    // Lazimi epsilon deqiqliyini alana qeder iterasiyalari davam edirik
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
            x1 = x2;    // Tapilmish x2 noqtesi ve bu noqtede funksiyanin qiymeti
            y1 = y2;    // sonraki iterasiyada x1 qiymeti kimi istifade olunacaq.
            x2 = NAN;   // x2 novbeti iterasiyada axtarilacaq
        }
        else
        {
            b = x2;
            x2 = x1;    // Tapilmish x1 noqtesi ve bu noqtede funksiyanin qiymeti
            y2 = y1;    // sonraki iterasiyada x2 qiymeti kimi istifade olunacaq.
            x1 = NAN;   // x1 novbeti iterasiyada axtarilacaq
        }
    }

    double c = (a+b)/2.0;

    double fa = f->fx(a);
    double fb = f->fx(b);

    //if (fa==fb) c = (a+b)/2.0;

    if (fa<fb)  c = a;
    if (fa>fb)  c = b;

    x = c;
    return c;
}
