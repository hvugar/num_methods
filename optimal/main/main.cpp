#include <stdio.h>
#include <methods.h>

double f1(double *x, int n)
{
    return 0.0;
}

int main()
{
    double x[] = {0.0, 0.0};
    double g[] = {0.0, 0.0};
    gradient(f1, x, 2, 0.1, g);
    return 0;
}

