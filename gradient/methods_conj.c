#include "methods.h"

void print(int iter, double *x, double *s, double *s1, int n, RnFunction f, int count);
double minimize(RnFunction f, double *x, double *grad, int n, double alpha0, double step, double epsilon, int *count);


void conjugate_gradient_method(RnFunction f, double *x, int n, double line_eps, double gold_eps, double grad_eps, double epsilon)
{
    int iter = 0;
    int i = 0;
    double *s = (double*) malloc(sizeof(double) * n);
	double *s1 = (double*) malloc(sizeof(double) * n);
	int count = 0;
    
	for (i=0; i<n; i++) 
	{ 
		s[i] = 0.0;
		s1[i] = 0.0; 
	}

    do
    {
        // First iteration
        double* gr1 = (double*) malloc(sizeof(double) * n);
        gradient(f, x, n, grad_eps, gr1);
		
        for (i=0; i<n; i++) 
		{
			s[i] = -gr1[i];
		}
		
		double ss = grad_module(s, n);
		
		for (i=0; i<n; i++) 
		{
			s1[i] = s[i] / ss;
		}

        double gr1_mod = 0.0;
        for (i=0; i<n; i++) 
		{
			gr1_mod += gr1[i]*gr1[i];
		}
		
        int k = 0;
		count += 2*n;
        do
        {
			print(iter, x, s, s1, n, f, count);
			count = 0;

			double alpha0 = 0.0;
			double alpha = minimize(f, x, s1, n, alpha0, line_eps, gold_eps, &count);

            for (i=0; i<n; i++) 
			{
				x[i] = x[i] + alpha * s1[i];
			}

			double* gr2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, grad_eps, gr2);
			count += 2*n+1;

            double gr2_mod = 0.0;
            for (i=0; i<n; i++) 
			{
				gr2_mod += gr2[i]*gr2[i];
			}

            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;

			for (i=0; i<n; i++) 
			{
				s[i] = -gr2[i] + s[i] * w;
			}
			
			ss = grad_module(s, n);
			for (i=0; i<n; i++) 
			{
				s1[i] = s[i] / ss;
			}

			iter++;

            free(gr2);
            gr2 = NULL;
			
            k++;
	    } while ( k < n );
		//puts("***");

        
		free(gr1);
        gr1 = NULL;
        
		double mod_s = s[0]*s[0] + s[1]*s[1];
        
		if ( mod_s < epsilon ) break;
		
    } while ( 1 );
    
	free(s1);
    s1 = NULL;

	free(s);
    s = NULL;
}

double minimize(RnFunction f, double *x, double *grad, int n, double alpha0, double line_eps, double gold_eps, int *count)
{
	double argmin(double alpha)
	{
		(*count)++;
		int j;
		for (j=0; j<n; j++) x[j] = x[j] + alpha * grad[j];
		double result = f(x, n);
		for (j=0; j<n; j++) x[j] = x[j] - alpha * grad[j];
		return result;
	}
	
	double a,b;
	straight_line_search_metod(argmin, alpha0, line_eps, &a, &b, count);
	double min = golden_section_search_min(argmin, a, b, gold_eps, count);
	
	return min; 
}

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
