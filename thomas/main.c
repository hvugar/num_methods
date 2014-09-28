#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
	/************************************************
	* b1*x1+c1*x2                                                   = d1,	a1 = 0
	* a2*x1+b2*x2+c2*x3                                             = d2
	*       a3*x2+b3*x3+c3*x3                                       = d3
 	*             a4*x3+b4*x4+c4*x5                                 = d4
	*      ---------------------------------------------
							a(n-1)*x(n-2)+b(n-1)*x(n-1)+c(n-1)*xn = d(n-1)
										      an*x(n-1)+    bn*xn = dn,     cn = 0	
	************************************************/
	
	int N = 0;
	
	printf("Enter number of linear equation of system: ");
	scanf("%d", &N);
	
	if (N==1)
	{
		puts("Number of equations must be at least 2");
		return 0;
	}
	
	int i = N;
	for (i=1; i<=N; i++)
	{
		if (i == 1)
		{
			printf("b1*x1 + c1*x2 = d1\n");
		}
		else if (i == N)
		{
			printf("a%d*x%d + b%d*x%d = d%d\n", i, i-1, i, i, i);
		}
		else
		{
			printf("a%d*x%d + b%d*x%d + c%d*x%d = d%d\n", i, i-1, i, i, i, i+1, i);
		}
		
	}
	
	float *a = (float*) malloc ( sizeof(float) * N );
	float *b = (float*) malloc ( sizeof(float) * N );
	float *c = (float*) malloc ( sizeof(float) * N );
	float *d = (float*) malloc ( sizeof(float) * N );
	
	for (i=0; i<N; i++)
	{
		if (i == 0)
		{
			printf("b1 = ");
			scanf("%f", &b[1]);
			printf("c1 = ");
			scanf("%f", &c[1]);
			printf("d1 = ");
			scanf("%f", &d[1]);
		}
		else if (i == N-1)
		{
			printf("a%d = ", i+1);
			scanf("%f", &a[i]);
			printf("b%d = ", i+1);
			scanf("%f", &b[i]);
			printf("c%d= ", i+1);
			scanf("%f", &c[i]);
			printf("d%d = ", i+1);
			scanf("%f", &d[i]);
		}
		else
		{
			printf("a%d*x%d + b%d*x%d + c%d*x%d = d%d\n", i, i-1, i, i, i, i+1, i);
		}
	}
	
	printf("Number of linear equation is %d\n", N);
	
	return 0;
}
