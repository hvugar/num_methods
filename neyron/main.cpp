#include <stdio.h>
#include <math.h>

int main()
{
    int n;
    printf("n: ");
    scanf("%d", &n);

    double* w = new double[n];
    double* x = new double[n];

    double wi;
    double xi;
    for (int i=0; i<n; i++)
    {
        printf("w[%2d]: ", i);
        scanf("%lf", &wi);
        printf("x[%2d]: ", i);
        scanf("%lf", &xi);

        w[i] = wi;
        x[i] = xi;
    }

    double sum = 0.0;
    for (int i=0; i<n; i++)
    {
        sum += w[i]*x[i];
    }
    double y = 1.0/(1.0+exp(-sum));

    printf("y: %f\n", y);

    delete [] x;
    delete [] w;

    return 0;
}

