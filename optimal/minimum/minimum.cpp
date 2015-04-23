#include "methods.h"

double straight_line_search_metod(R1Function f, double x0, double dx, double &a, double &b)
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
        a = x0 - dx;
        b = x0 + dx;
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
        a = x0 - dx;
        b = x0 + dx;
    }
    else
    {
        a = x0 + dx;
        b = x0 - dx;
    }

    if (a > b)
    {
        fprintf(stderr, "Value of a is grater than value of b. a=%.10f, b=%.10f\n", a, b);
        system("pause");
    }

    double c = (a + b) / 2.0;

    return c;
}

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

    double c = (a+b)/2.0;

    if (f(a)<f(b)) c = a;
    if (f(a)>f(b)) c = b;

    return c;
}

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

double R1Minimize(R1Function f, double line_step, double gold_epsilon)
{
    double a,b;
    double alpha0 = 0.0;
    straight_line_search_metod(f, alpha0, line_step, a, b);
    double alpha = golden_section_search_min(f, a, b, gold_epsilon);
    if ( f(alpha) > f(alpha0) ) alpha = alpha0;
    return alpha;
}
