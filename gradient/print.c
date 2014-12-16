#include "print.h"

void print(int iter, double *x, double *s, double *s1, int n, RnFunction f, int count)
{
	if (iter == 0)
	{
		printf("No\t|x1      \t|x2      \t|f(x)      \t|grad1      \t|grad2      \t|grad1/norm \t|grad2/norm\t|say      \t");
		printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---\n");
	}
	
	printf("%d\t", iter);
	x[0]>=0 ? printf("|+%.6f\t", x[0]) : printf("|%.6f\t", x[0]);
	x[1]>=0 ? printf("|+%.6f\t", x[1]) : printf("|%.6f\t", x[1]);
	f(x,n)>=0 ? printf("|+%.6f\t", f(x,n)) : printf("|%.6f\t", f(x,n));
	s[0]>=0 ? printf("|+%.6f\t", s[0]) : printf("|%.6f\t", s[0]);
	s[1]>=0 ? printf("|+%.6f\t", s[1]) : printf("|%.6f\t", s[1]);
	s1[0]>=0 ? printf("|+%.6f\t", s1[0]) : printf("|%.6f\t", s1[0]);
	s1[1]>=0 ? printf("|+%.6f\t", s1[1]) : printf("|%.6f\t", s1[1]);
	printf("|%d\t", count);
	printf("\n");
}