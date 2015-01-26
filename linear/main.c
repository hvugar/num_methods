#include <stdio.h>
#include <stdlib.h>
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
	
	for (i=0; i<n; i++)
	{
		if (a[i][i] == 0.0)
		{
			//assert(a[i][i] != 0.0);
			//fputs("Dioqanal element is equal to zero\n", stderr);	
		}
	}
	
	for (i=0; i<n; i++)
	{
		b[i] = b[i] / a[i][i];
		
		for (j=0; j<n; j++)
			if (i != j) 
				a[i][j] = -(a[i][j] / a[i][i]);
		a[i][i] = 0.0;
		x0[i] = b[i];
		
		printf("b[%d] = %f\n", i, b[i]);
	}
	
	int k=0;
	while ( 1 )
	{
		for (i=0; i<n; i++)
		{
			x[i] = x0[i];
		}

		for (i=0; i<n; i++)
		{
			x0[i] = b[i];
			for (j=0; j<n; j++)
			{
				x0[i] += a[i][j] * x[j];
			}
		}
		
		double max = 0.0;
		for (i=0; i<n; i++)
		{
			double diff = fabs(x[i] - x0[i]);
			if (diff > max) max = diff;
		}
		printf("%f\n", max);
		
		if (max < epsilon)
			break;
			
		printf("%f %f %f %f\n", x0[0], x0[1], x0[2], x0[3]);
		k++;
		if (k>10)
		break;
				
	}
}

