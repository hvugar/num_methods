#include "diff_method.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * @brief Implicit difference scheme/Неявная разностная схема
 * @param f
 * @param fi		Initial condition function at t=0
 * @param m1		Boundary condition function at x=0
 * @param m2		Boundary condition function at x=l
 * @oaram a
 * @param delta_x	step for x
 * @param delta_t	step for t
 * @param l
 * @param t
 */
void implicit_difference_scheme(FxtFunction f, fiFunction fi, m1Function m1, m2Function m2, double a, double dx, double dt, double l, double t) 
{
    double di = l / dx;
    double dj = t / dt;

    int M = (int) di; // x uzre
    int N = (int) dj; // t uzre

    int j = 0;
    int i = 0;

    int size = M-1;

    printf("%d %d\n", M, N);

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
}

//
/**
 * @brief Explicit difference scheme/Явная разностная схема
 * @param f
 * @param fi	Initial condition function at t=0
 * @param m1	Boundary condition function at x=0
 * @param m2	Boundary condition function at x=l
 * @param a
 * @param dx	step for x
 * @param dt	step for t
 * @param l
 * @param t
 */
void explicit_difference_scheme(FxtFunction f, fiFunction fi, m1Function m1, m2Function m2, double a, double dx, double dt, double l, double t) 
{
    double di = l / dx;
    double dj = t / dt;

    int M = (int) di + 1; // x uzre
    int N = (int) dj + 1; // t uzre

    int j = 0;
    int i = 0;
	
	double** u = (double**) malloc(sizeof(double*) * N);
	for (j=0; j<N; j++)
	{
		u[j] = (double*) malloc(sizeof(double) * M);
		for (i=0; i<M; i++)
			u[j][i] = 0.0;
	}
	
	for (i=0; i<M; i++)
	{
		u[0][i] = fi(i*dx);
	}

	for (j=0; j<N; j++)
	{
		u[j][0]   = m1(j*dt);
		u[j][M-1] = m2(j*dt);
	}
	
	
	for (j=0; j<N-1; j++)
	{
		for (i=1; i<M-1; i++)
		{
			double r = 2.0;
			double k = a*(dt/(dx*dx));
			double u1 = (u[j][i-1] - 2.0 * u[j][i] + u[j][i+1]);
			double u2 = -dt;
			double u3 = k * u1;
			u[j+1][i] = u3 + u[j][i] + u2;
		}
	}
	
	for (j=0; j<N; j++)
	{
		for (i=0; i<M; i++)
		{
			printf("%.14f ", u[j][i]);
		}
		printf("\n");
		for (i=0; i<M; i++)
		{
			printf("%.14f ", (i*dx)*(i*dx)+(j*dt));
		}
		printf("\n----------------\n");
	}
	
	for (j=0; j<N; j++)
	{
		free(u[j]);
		u[j] = NULL;
	}
	free(u);
	u = NULL;
}
