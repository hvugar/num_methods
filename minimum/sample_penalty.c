#include "methods.h"

extern double f_rosenbrock(double *x, int n); 

void penalty_sample1()
{
    double f(double *x, int n)  { return (x[0]-4.0)*(x[0]-4.0)+(x[1]-4.0)*(x[1]-4.0); }
    double h1(double *x, int n) { return x[0] + x[1] - 5.0; }
    double g1(double *x, int n) { return 5.0 - x[0] - x[1]; }

    int n = 2;
    double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = +10;
    x[1]    = +10;

    int m = 0;
    RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);
    h[0] = h1;

    int p = 1;
    RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
    g[0] = g1;

    double r1 = 100.00; //G
    double r2 = 1.00;   //H

    double epsilon = 0.01;
    penalty_method(f, x, n, h, m, g, p, r1, r2, epsilon);

    free(x);
    free(h);
    free(g);
}

void penalty_sample2()
{
    // x* = [0.0, 0.5]
    printf("Penalty sample 2. x* = [0.5, 0.0]\n");

    double f(double *x, int n) { return x[0]*x[0] + x[1]*x[1] + 2*x[1]; }

    double h1(double *x, int n) { return x[0]*x[0] + x[1]*x[1] - 1.0; }
    double g1(double *x, int n) { return x[0] + 2*x[1] - 0.5; }
    double g2(double *x, int n) { return x[0]; }
    double g3(double *x, int n) { return x[1]; }

    int n = 2;

    double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = +0.5;
    x[1]    = +0.5;

    int m = 1;
    RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);
    h[0] = h1;

    int p = 3;
    RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
    g[0] = g1;
    g[1] = g2;
    g[2] = g3;

    double r1 = 0.01;
    double r2 = 1.0;
    double epsilon = 0.001;
    penalty_method(f, x, n, h, m, g, p, r1, r2, epsilon);

    free(x);
    free(h);
    free(g);
}

void penalty_sample3()
{
    double f(double *x, int n) { return x[0]*x[0]; }

    double g1(double *x, int n) { return x[0]-1.0; }
    double g2(double *x, int n) { return 2.0-x[0]; }

    int n = 1;
    int m = 0;
    int p = 2;
    double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = +1.5;

    RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);

    RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
    g[0] = g1;
    g[1] = g2;

    double r1 = 1.00;
    double r2 = 1.00;
    double epsilon = 0.05;
    penalty_method(f, x, n, h, m, g, p, r1, r2, epsilon);

    free(x);
    free(h);
    free(g);

}

void penalty_sample4()
{
    // x* = [1.001282493, 4.89871752]

    printf("Penalty sample 4. x* = [1.001282493, 4.89871752]\n");

    int n = 2;
    double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = +2.0;
    x[1]    = +4.0;

    int m = 1;
    RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);
    double h1(double *x, int n) { return 25.0 - x[0]*x[0] - x[1]*x[1]; }
    h[0] = h1;

    int p = 3;
    RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
    double g1(double *x, int n) { return 10*x[0] - x[0]*x[0] + 10*x[1] - x[1]*x[1]- 34.0; }
    g[0] = g1;
    double g2(double *x, int n) { return x[0]; }
    g[1] = g2;
    double g3(double *x, int n) { return x[1]; }
    g[2] = g3;

    double r1 = 0.01;
    double r2 = 10.0;
    double epsilon = 0.0000001;
    double f(double *x, int n) { return 4*x[0] - x[1]*x[1] - 12.0; }
    penalty_method(f, x, n, h, m, g, p, r1, r2, epsilon);

    free(x);
    free(h);
    free(g);
}
