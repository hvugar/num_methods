#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void simple_iteration(double** A, double* b, double* x0, double* x, int n, double epsilon);

int main(int argc, char** argv)
{
	int i,j;
	int n = 3;

	double** a = (double**) malloc( sizeof(double*) * n );
	for (i=0; i<n; i++)
		a[i] = (double*) malloc( sizeof(double) * n );
	double* b = (double*) malloc( sizeof(double) * n );
	double* x0 = (double*) malloc( sizeof(double) * n );
	double* x = (double*) malloc( sizeof(double) * n );
	double epsilon = 0.0001;
	
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
			a[i][j] = 0.0;
	}
	
	simple_iteration(a, b, x0, x, n, epsilon);


	return 0;
}

void simple_iteration(double** a, double* b, double* x0, double* x, int n, double epsilon)
{
	int i,j;
	
	for (i=0; i<n; i++)
	{
		if (a[i][i] == 0.0)
		{
			assert(a[i][i] != 0.0);
			fputs("Dioqanal element is equal to zero\n", stderr);	
		}
	}
}

