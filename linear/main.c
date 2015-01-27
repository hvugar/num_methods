#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

void simple_iteration(double** A, double* b, double* x0, double* x, int n, double epsilon);

int main(int argc, char** argv)
{
	int i,j;
	int n = 4;

	double** a = (double**) malloc( sizeof(double*) * n );
	for (i=0; i<n; i++)
		a[i] = (double*) malloc( sizeof(double) * n );
	double* b = (double*) malloc( sizeof(double) * n );
	double* x0 = (double*) malloc( sizeof(double) * n );
	double* x = (double*) malloc( sizeof(double) * n );
	double epsilon = 0.000001;
	
	for (i=0; i<n; i++)
	{
		x[i] = 0.0;
		for (j=0; j<n; j++)
			a[i][j] = 0.0;
	}
	
	a[0][0] = +1.0;
	a[0][1] = +2.0;
	a[0][2] = -5.0;
	a[0][3] = -4.0;
	
	a[1][0] = +2.0;
	a[1][1] = -1.0;
	a[1][2] = +3.0;
	a[1][3] = +1.0;

	a[2][0] = +5.0;
	a[2][1] = +1.0;
	a[2][2] = -1.0;
	a[2][3] = +4.0;

	a[3][0] = +0.0;
	a[3][1] = +1.0;
	a[3][2] = +5.0;
	a[3][3] = +4.0;
	
	b[0] = -22.0;
	b[1] = +13.0;
	b[2] = +20.0;
	b[3] = +33.0;
	
	simple_iteration(a, b, x0, x, n, epsilon);
	
	printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);

	return 0;
}

void simple_iteration(double** a, double* b, double* x0, double* x, int n, double epsilon)
{
	int i,j;

	// create  copy of matrix and vector
	double*  b1 = (double*)  malloc( sizeof(double)  * n );
	double** a1 = (double**) malloc( sizeof(double*) * n );
	
	memcpy( b1, b, sizeof(double)*n );
	for (i=0; i<n; i++)
	{
		a1[i] = (double*) malloc( sizeof(double)*n );
		memcpy( a1[i], a[i], sizeof(double)*n );
	}
	
	// checking
	for (i=0; i<n; i++) { if (a[i][i] == 0.0) { } }

	// initializing
	for (i=0; i<n; i++)
	{
		b[i] = b[i] / a[i][i];
		
		for (j=0; j<n; j++) { if (i != j) { a[i][j] = -(a[i][j] / a[i][i]); } }
		a[i][i] = 0.0;
	}
	
	printf("b[0] = %2.4f\n", b[0]);
	printf("b[1] = %2.4f\n", b[1]);
	printf("b[2] = %2.4f\n", b[2]);
	printf("b[3] = %2.4f\n", b[3]);
	puts("---");
	printf("%2.4f %2.4f %2.4f %2.4f\n", a[0][0], a[0][1], a[0][2], a[0][3]);
	printf("%2.4f %2.4f %2.4f %2.4f\n", a[1][0], a[1][1], a[1][2], a[1][3]);
	printf("%2.4f %2.4f %2.4f %2.4f\n", a[2][0], a[2][1], a[2][2], a[2][3]);
	printf("%2.4f %2.4f %2.4f %2.4f\n", a[3][0], a[3][1], a[3][2], a[3][3]);
	puts("---");
	printf("x0[0] = %2.4f\n", x0[0]);
	printf("x0[1] = %2.4f\n", x0[1]);
	printf("x0[2] = %2.4f\n", x0[2]);
	printf("x0[3] = %2.4f\n", x0[3]);
	
	for (i=0; i<n; i++) { x[i] = b[i]; }

	int k=0;
	while ( 1 )
	{
		/* saving last values */
		for (i=0; i<n; i++)
		{
			x0[i] = x[i];
		}		
		
		for (i=0; i<n; i++)
		{
			x[i] = b[i];
			
			for (j=0; j<n; j++)
			{
				x[i] = x[i] + a[i][j] * x0[j];
			}
		}
		
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
		
		double max = 0.0;
		for (i=0; i<n; i++)
		{
			double diff = fabs(x[i] - x0[i]);
			if (diff > max) max = diff;
		}
		printf("%f\n", max);
		
		if (max < epsilon)
			break;
			
//		k++;
//		if (k>10)
//		break;			
	}
	//printf("%f %f %f %f\n", x0[0], x0[1], x0[2], x0[3]);
}

