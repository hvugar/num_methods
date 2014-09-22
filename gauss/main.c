#include <stdio.h>
#include <stdlib.h>

void print_matrix(float** a, float* b, int n)
{
	int i,j;
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			printf("%10.3f", a[i][j]);
		}
		printf("%10.3f", b[i]);
		printf("\n");
	}
}

void check_matrix(float** a, float* b, int n)
{
	int i,j,k;
	for (i=0; i<n; i++)
	{
		int c = 0;
		for (j=0; j<n; j++)
		{
			if (a[i][j] != 0) c++;
		}
		if (c==0)
		{
			if (b[i] == 0) 
				fprintf(stderr, "Equation has infinity solutions");
			else
				fprintf(stderr, "Equation has not any solution");
		}			
	}
}

void gaussian_elimination(float** a, float* b, float* x, int n)
{
	int i,j,k;
	for (k=0; k<n-1; k++)
	{
		for (i=(k+1); i<n; i++)
		{
			double f = a[i][k] / a[k][k];
			for (j=k; j<n; j++)
			{
				a[i][j] = a[i][j] - a[k][j] * f;
			}
			b[i] = b[i] - b[k] * f;
		}
		check_matrix(a, b, n);
	}
	
	for (i=(n-1); i>=0; i--)
	{
		for (j=(n-1); j>i; j--)
			b[i] -= (a[i][j] * x[j]);
		
		x[i] = b[i] / a[i][i];
	}
}

int main(int argc, char** argv)
{
	/********************************************
	
	a11*x1 + a12*x2 + a13*x3 + ... + a1n*xn = b1
	a21*x1 + a22*x2 + a23*x3 + ... + a2n*xn = b2
	a31*x1 + a22*x2 + a33*x3 + ... + a3n*xn = b3
	............................................
	an1*x1 + an2*x2 + an3*x3 + ... + ann*xn = bn
	
	********************************************/

	int i,j,k;

	int N = 4;

	/* */
	float** a = (float**) malloc( sizeof(float*) * N );
	a[0] = (float*) malloc ( sizeof(float) * N);
	a[1] = (float*) malloc ( sizeof(float) * N);
	a[2] = (float*) malloc ( sizeof(float) * N);
	a[3] = (float*) malloc ( sizeof(float) * N);
	
	a[0][0] = 3.0;  a[0][1] = 3.0;  a[0][2] = -5.0; a[0][3] = 1.0;
	a[1][0] = 11.0; a[1][1] = -2.0; a[1][2] = -4.0; a[1][3] = 3.0;
	a[2][0] = 1.0;  a[2][1] = 1.0;  a[2][2] = -3.0; a[2][3] = 2.0;
	a[3][0] = 2.0;  a[3][1] = 3.0;  a[3][2] = 2.0;  a[3][3] = -1.0;
	
	/* */
	float* b = (float*) malloc ( sizeof(float) * N);
	b[0] = 14.0;
	b[1] = 61.0;
	b[2] = 17.0;
	b[3] = 6.0;
	
	float* x = (float*) malloc ( sizeof(float) * N);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	
	print_matrix(a, b, N);
	
	gaussian_elimination(a, b, x, N);
	
	puts("\n");
	
	for (j=0; j<N; j++) 
		printf("%10.3f", x[j]);

	return 0;
}
