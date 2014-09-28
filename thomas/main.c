#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
	/************************************************
	*
	*
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
	for (i=0; i<N; i++)
	{
		if (i == 0)
		printf("b1*x1 + c1*x2 = d1\n");
		
		if (i == N-1)
		printf("a%d*x%d + b%d*x%d = d%d\n", i+1, i, i+1, i+1, i+1);
	}
	
	printf("Number of linear equation is %d\n", N);
	
	return 0;
}
