#include "print.h"

void printer1(RnFunction f, double *x, int n, ...)
{
    va_list ap;
    va_start(ap, 0);

    int i = va_arg(ap, int);
    double* gr = va_arg(ap, double*);
    double nr = va_arg(ap, double);
	nr = nr;
    double alpha = va_arg(ap, double);

    if (i == 0)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|grad1      \t|grad2      \t|alpha      \t|count      \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+--------+\n");
    }

    double y = f(x,n);

    printf("%d\t", i);
    x[0]>=0.0 ? printf("|+%.8f\t", fabs(x[0])) : printf("|%.8f\t", x[0]);
    x[1]>=0.0 ? printf("|+%.8f\t", fabs(x[1])) : printf("|%.8f\t", x[1]);
    y>=0.0 ? printf("|%+0.8f\t", y) : printf("|%0.8f\t", y);
    gr[0]>=0.0 ? printf("|%+.8f\t", gr[0]) : printf("|%.8f\t", gr[0]);
    gr[1]>=0.0 ? printf("|%+.8f\t", gr[1]) : printf("|%.8f\t", gr[1]);
    alpha>=0.0 ? printf("|%+.8f\t", alpha) : printf("|%.8f\t", alpha);
    printf("|%d\t", n);
    printf("\n");

    va_end(ap);
}

void printer2(RnFunction f, double *x, int n, ...)
{
    va_list ap;
    va_start(ap, 0);

    int i = va_arg(ap, int);
    int count = va_arg(ap, int);
    double* g = va_arg(ap, double*);
    double* s = va_arg(ap, double*);
	double nr = va_arg(ap, double);
    double alpha = va_arg(ap, double);

    if (i == 0)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|grad1      \t|grad2      \t|s1         \t|s2       \t|norm         \t|alpha       \t|say      \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---\n");
    }
	
	double y = f(x,n);

    printf("%d\t", i);
    x[0]>=0 ? printf("|+%.8f\t", x[0]) : printf("|%.8f\t", x[0]);
    x[1]>=0 ? printf("|+%.8f\t", x[1]) : printf("|%.8f\t", x[1]);
    y>=0 ? printf("|+%.8f\t", y) : printf("|%.8f\t", y);
    g[0]>=0 ? printf("|+%.8f\t", g[0]) : printf("|%.8f\t", g[0]);
    g[1]>=0 ? printf("|+%.8f\t", g[1]) : printf("|%.8f\t", g[1]);
    s[0]>=0 ? printf("|+%.8f\t", s[0]) : printf("|%.8f\t", s[0]);
    s[1]>=0 ? printf("|+%.8f\t", s[1]) : printf("|%.8f\t", s[1]);
	nr>=0 ? printf("|+%.8f\t", nr) : printf("|%.8f\t", nr);
    alpha>=0 ? printf("|+%.8f\t", alpha) : printf("|%.8f\t", alpha);
    printf("|%d\t", count);
    printf("\n");

    va_end(ap);
}

void printer3(RnFunction f, double *x, int n, ...)
{
    //int n_args;
    va_list ap;
    va_start(ap, 0);

    int i = va_arg(ap, int);
    double* grad = va_arg(ap, double*);
    double grad_norm = va_arg(ap, double);
    double alpha = va_arg(ap, double);

    if (i == 0)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|grad1      \t|grad2      \t|grad_norm  \t|delta_x1  \t|delta_x2 \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+\n");
    }

    printf("%d\t", i);
    x[0]>=0 ? printf("|+%.6f\t", x[0]) : printf("|%.6f\t", x[0]);
    x[1]>=0 ? printf("|+%.6f\t", x[1]) : printf("|%.6f\t", x[1]);
    f(x,n)>=0 ? printf("|+%.6f\t", f(x,n)) : printf("|%.6f\t", f(x,n));
    grad[0]>=0 ? printf("|+%.6f\t", grad[0]) : printf("|%.6f\t", grad[0]);
    grad[1]>=0 ? printf("|+%.6f\t", grad[1]) : printf("|%.6f\t", grad[1]);
    grad_norm >=0 ? printf("|+%.6f\t", grad_norm) : printf("|%.6f\t", grad_norm);
    -grad[0]/grad_norm >=0 ? printf("|+%.6f\t", -grad[0]/grad_norm) : printf("|%.6f\t", -grad[0]/grad_norm);
    -grad[1]/grad_norm >=0 ? printf("|+%.6f\t", -grad[1]/grad_norm) : printf("|%.6f\t", -grad[1]/grad_norm);
    printf("|%f", alpha);
    printf("\n");

    va_end(ap);
}

void printer4(RnFunction P, double *x, int n, ...)
{
    //int n_args;
    va_list ap;
    va_start(ap, 0);

    int iter = va_arg(ap, int);
    //	int count = va_arg(ap, int);
    //	double* s = va_arg(ap, double*);
    //	double* s1 = va_arg(ap, double*);

    if (iter == 0)
    {
        printf("No\t|x1      \t|x2      \t|x3      \t|x4      \t|x5      \t|x6         \t|x7      \t|x8      \t|x9      \t|x10      \t|F(x)      \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+\n");
    }

    double y = P(x,n);

    printf("%d\t", iter);
    x[0]>=0 ? printf("|+%.10f\t", x[0]) : printf("|%.10f\t", x[0]);
    x[1]>=0 ? printf("|+%.10f\t", x[1]) : printf("|%.10f\t", x[1]);
    x[2]>=0 ? printf("|+%.10f\t", x[2]) : printf("|%.10f\t", x[2]);
    x[3]>=0 ? printf("|+%.10f\t", x[3]) : printf("|%.10f\t", x[3]);
    x[4]>=0 ? printf("|+%.10f\t", x[4]) : printf("|%.10f\t", x[4]);
    x[5]>=0 ? printf("|+%.10f\t", x[5]) : printf("|%.10f\t", x[5]);
    x[6]>=0 ? printf("|+%.10f\t", x[6]) : printf("|%.10f\t", x[6]);
    x[7]>=0 ? printf("|+%.10f\t", x[7]) : printf("|%.10f\t", x[7]);
    x[8]>=0 ? printf("|+%.10f\t", x[8]) : printf("|%.10f\t", x[8]);
    x[9]>=0 ? printf("|+%.10f\t", x[9]) : printf("|%.10f\t", x[9]);
    y>=0 ? printf("|+%.10f\t", y) : printf("|%.10f\t", y);
    //	s[0]>=0 ? printf("|+%.6f\t", s[0]) : printf("|%.6f\t", s[0]);
    //	s[1]>=0 ? printf("|+%.6f\t", s[1]) : printf("|%.6f\t", s[1]);
    //	s1[0]>=0 ? printf("|+%.6f\t", s1[0]) : printf("|%.6f\t", s1[0]);
    //	s1[1]>=0 ? printf("|+%.6f\t", s1[1]) : printf("|%.6f\t", s1[1]);
    //	printf("|%d\t", count);
    printf("\n");

    va_end(ap);
}

void printX(char *label, double *x, int n)
{
    int i;
    printf("double %s[] = \t{", label);
    for (i=0; i<n; i++)
    {
        printf("%12.8f", x[i]);
        if (i != n-1 )
            printf(", ");
    }
    printf("};");
    printf("\n");
}

void _print1(char *s, double *a, int n)
{
    int i;
    printf("double %s[] =\t{", s);
    //    printf("%s[] =\t{", s);
    for (i=0; i<n; i++)
    {
        if ( i%((n-1)/10) == 0 )
        {
            if (a[i] < 0)
            {
                printf("%12.8f", a[i]);
            }
            else
            {
                printf("%+12.8f", a[i]);
            }
        }
        if ( i%((n-1)/10) == 0 && i != n-1 )
        {
            printf(", ");
        }
    }
    printf("};");
    printf("\n");
}

void _seperator()
{
    puts("-------------------------------------------------------------------------------------------------------------------------------------------------------------------");
}

void _printM(double **x, int m, int n)
{
    int i,j;
    for (j=0; j<m; j++)
    {
        if (j%(m/10)==0)
        {
			printf("%6d|", j);
            for (i=0; i<n; i++)
            {
                if (i%(n/10)==0)
                    printf("%12.6f", x[j][i]);
            }
            puts("");
        }
    }
}
