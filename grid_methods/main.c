#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double (*R1Function)(double);
typedef double (*RnFunction)(double);

double f(double, double);

//Initial condition
double fi(double x);
//Boundary conditions
double m1(double t);
double m2(double t);

void tomas_algorithm(double **a, double* b, double* x, int size);

int main(int argc, char** argv)
{
    double dx = 0.1;
    double dt = 0.1;
    double a  = 1.0;

    int M = 10; // x uzre
    int N = 10; // t uzre

    int j = 0;
    int i = 0;

    int size = 9;

    double** A = (double**) malloc(sizeof(double*)*size);
    double*  b = (double*)  malloc(sizeof(double)*size);
    double*  u = (double*)  malloc(sizeof(double)*size);

    for (i=0; i<size; i++)
    {
        A[i] = (double*)malloc(sizeof(double)*size);
        for (j=0; j<size; j++)
            A[i][j] = 0.0;
        b[i] = 0.0;
        u[i] = 0.0;
    }

    for (j=0; j<N; j++)
    {
        printf("Layer: %4d ", j+1);

        if ( j == 0 )
        {
            for (i=1; i<=M-1; i++)
            {
                if ( i == 1 )
                    b[i-1] = fi(dx*i) + dt*f(dx*i,dt*(j+1)) + (dt/(dx*dx))*m1(dt*(j+1));
                else if ( i == (M-1) )
                    b[i-1] = fi(dx*i) + dt*f(dx*i,dt*(j+1)) + (dt/(dx*dx))*m2(dt*(j+1));
                else
                    b[i-1] = fi(dx*i) + dt*f(dx*i,dt*(j+1));
            }
        }
        else
        {
            for (i=1; i<=M-1; i++)
            {
                if ( i == 1 )
                    b[i-1] = u[i] + dt*f(dx*i,dt*(j+1)) + (dt/(dx*dx))*m1(dt*(j+1));
                else if ( i == (M-1) )
                    b[i-1] = u[i] + dt*f(dx*i,dt*(j+1)) + (dt/(dx*dx))*m2(dt*(j+1));
                else
                    b[i-1] = u[i] + dt*f(dx*i,dt*(j+1));
            }
        }

        for (i=1; i<=M-1; i++)
        {
            if ( i == 1 )
            {
                A[i-1][0] = 1 + (2*dt)/(dx*dx);
                A[i-1][1] = -dt/(dx*dx);
            }
            else if ( i == (M-1) )
            {
                A[i-1][i-2] = -dt/(dx*dx);
                A[i-1][i-1] = 1 + (2*dt)/(dx*dx);
            }
            else
            {
                A[i-1][i-2] = -dt/(dx*dx);
                A[i-1][i-1] = 1 + (2*dt)/(dx*dx);
                A[i-1][i-0] = -dt/(dx*dx);
            }
        }

        tomas_algorithm(A, b, u, size);

        for (i=0; i<size; i++)
            printf("%8.4f ", u[i]);

		puts("");
    }

    free(u);
    free(b);
    free(A);

    return 0;
}

double f(double x, double t) { return -1.0; }

double fi(double x) { return x*x; }

double m1(double t) { return t; }

double m2(double t) { return 1.0+t; }


