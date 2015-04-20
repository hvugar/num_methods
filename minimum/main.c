#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"
#include "optimal.h"

void print1(char* s, double* a, int n);
void seperator();
double fx1(double t, double x1, double x2, double u);
double fx2(double t, double x1, double x2, double u);
double fp1(double t, double x1, double x2, double p1, double p2, double u);
double fp2(double t, double x1, double x2, double p1, double p2, double u);
double dIdu(double t, double x1, double x2, double p1, double p2, double u);

void runga_kutta_system(RmFunction *f, double t0, double *x0, double t, double *x, int n, double h);
void runga_kutta_system1(double t0, double x01, double x02, double t1, double x11, double x12, int n, double h, double u);

double JSum(double *t, double *x1, double *x2, int n, double *u, int N)
{
    double sum = 0.0;
    int i=0;
    for (i=0; i<(N-1); i++)
    {
        int j=i+1;
        double fj = (x1[j]-t[j]*t[j]*t[j])*(x1[j]-t[j]*t[j]*t[j]) + x2[j]*x2[j] - t[j]*t[j] + (2*u[j] - 1.0)*(2*u[j] - t[j]);
        double fi = (x1[i]-t[i]*t[i]*t[i])*(x1[i]-t[i]*t[i]*t[i]) + x2[i]*x2[i] - t[i]*t[i] + (2*u[i] - 1.0)*(2*u[i] - t[i]);
        sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
    }
    //double x[] = { x1[N-1], x2[N-1] };
    return sum;
}

int main(int argc, char** argv)
{
    int i;
    int N = 1000001;
    double h = 0.000001;

    //double t0 = 0.0;
    //double t1 = 1.0;

    //double x10 = 0.0;
    //double x20 = 0.0;

    double *t  = (double*) malloc( sizeof(double) * N );
    double *u  = (double*) malloc( sizeof(double) * N );
    double *x1 = (double*) malloc( sizeof(double) * N );
    double *x2 = (double*) malloc( sizeof(double) * N );
    double *p1 = (double*) malloc( sizeof(double) * N );
    double *p2 = (double*) malloc( sizeof(double) * N );
    double *gr = (double*) malloc( sizeof(double) * N );

    for (i=0; i<N; i++)
    {
        t[i] = i*h;
        u[i] = t[i]/2.0;
        x1[i] = 0.0;
        x2[i] = 0.0;
    }
    print1("t", t, N);
    print1("u", u, N);
    seperator();
    print1("x1", x1, N);
    print1("x2", x2, N);

    double k11;
    double k12;
    double k21;
    double k22;
    double k31;
    double k32;
    double k41;
    double k42;

    //if (t1<t0) h = -fabs(h);

	double j1, j2;
    //while (1)
    {
    //print1("t", t, N);
    print1("u", u, N);
	//seperator();
        i=0;
        x1[i] = 0.0;
        x2[i] = 0.0;
        h = +0.000001;
		for (i=0; i<N-1; i++)
//        while (i<N-1)
        {
            k11 = fx1(t[i],       x1[i],                 x2[i],                 u[i]);
            k12 = fx2(t[i],       x1[i],                 x2[i],                 u[i]);
            k21 = fx1(t[i]+h/2.0, x1[i] + (h/2.0) * k11, x2[i] + (h/2.0) * k12, u[i]);
            k22 = fx2(t[i]+h/2.0, x1[i] + (h/2.0) * k11, x2[i] + (h/2.0) * k12, u[i]);
            k31 = fx1(t[i]+h/2.0, x1[i] + (h/2.0) * k21, x2[i] + (h/2.0) * k22, u[i]);
            k32 = fx2(t[i]+h/2.0, x1[i] + (h/2.0) * k21, x2[i] + (h/2.0) * k22, u[i]);
            k41 = fx1(t[i]+h,     x1[i] + (h) * k31,     x2[i] + (h) * k32,     u[i]);
            k42 = fx2(t[i]+h,     x1[i] + (h) * k31,     x2[i] + (h) * k32,     u[i]);

            x1[i+1] = x1[i] + (h/6.0) * (k11 + 2*k21 + 2*k31 + k41);
            x2[i+1] = x2[i] + (h/6.0) * (k12 + 2*k22 + 2*k32 + k42);

            //t0 = t0 + h;
            //i++;
            //		x1[i] = x10;
            //		x2[i] = x20;
        }

        seperator();
        print1("x1", x1, N);
        print1("x2", x2, N);
        seperator();

        //if (t1<t0) h = -fabs(h);

        //double p10 = 0.0;
        //double p20 = -2.0 * (x2[i] - 1.0);
        //t0 = 1.0;
        //t1 = 0.0;
		
        i=N-1;
        p1[i] = 0.0;
        p2[i] = -2.0 * (x2[i] - 1.0);
        h = -0.000001;
        //while (fabs(t1-t0) >= fabs(h/2))
        while (i>0)
        {
            k11 = fp1(t[i],       x1[i], x2[i], p1[i],                 p2[i],                 u[i]);
            k12 = fp2(t[i],       x1[i], x2[i], p1[i],                 p2[i],                 u[i]);
			
            k21 = fp1(t[i]+h/2.0, x1[i], x2[i], p1[i] + (h/2.0) * k11, p2[i] + (h/2.0) * k12, u[i]);
            k22 = fp2(t[i]+h/2.0, x1[i], x2[i], p1[i] + (h/2.0) * k11, p2[i] + (h/2.0) * k12, u[i]);
            
			k31 = fp1(t[i]+h/2.0, x1[i], x2[i], p1[i] + (h/2.0) * k21, p2[i] + (h/2.0) * k22, u[i]);
            k32 = fp2(t[i]+h/2.0, x1[i], x2[i], p1[i] + (h/2.0) * k21, p2[i] + (h/2.0) * k22, u[i]);
            
			k41 = fp1(t[i]+h,     x1[i], x2[i], p1[i] + (h) * k31,     p2[i] + (h) * k32,     u[i]);
            k42 = fp2(t[i]+h,     x1[i], x2[i], p1[i] + (h) * k31,     p2[i] + (h) * k32,     u[i]);

            p1[i-1] = p1[i] + (h/6.0) * (k11 + 2*k21 + 2*k31 + k41);
            p2[i-1] = p2[i] + (h/6.0) * (k12 + 2*k22 + 2*k32 + k42);

            //t0 = t0 + h;
            i--;
            //p1[i] = p10;
            //p2[i] = p20;
        }
        print1("p1", p1, N);
        print1("p2", p2, N);
        seperator();
		printf("%.10f\n", JSum(t, x1, x2, 2, u, N));

        for (i=0; i<N; i++)
        {
            gr[i] = dIdu(t[i], x1[i], x2[i], p1[i], p2[i], u[i]);
        }
        print1("gr", gr, N);

        double argmin1(double alpha)
        {
            int i;
            //double u1[N];
			double *u1  = (double*) malloc( sizeof(double) * N );
            for (i=0; i<N; i++) u1[i] = u[i] - alpha * gr[i];
			double J = JSum(t, x1, x2, 2, u1, N);
			free(u1);
            return J;
        }

        double a,b;
        double alpha0 = 0.0;
        straight_line_search_metod(argmin1, alpha0, 0.001, &a, &b);
        double alpha = golden_section_search_min(argmin1, a, b, 0.00000001);
        if ( argmin1(alpha) > argmin1(alpha0) ) alpha = alpha0;
        printf("alpha\t%14.10f\n", alpha);
		
		j1 = JSum(t, x1, x2, 2, u, N);
        for (i=0; i<N; i++)
        {
            u[i] = u[i] - alpha*gr[i];
        }
        j2 = JSum(t, x1, x2, 2, u, N);

        double s=0.0;
        for (i=0; i<N; i++)
            s += (t[i]-1.0)*(t[i]-1.0)*(t[i]-1.0)/3.0;
        printf("%.10f %.10f\n", j1, j2);


    }

    free(p1);
    free(p2);
    free(x1);
    free(x2);
    free(u);
    free(t);

    return 0;
}

double fx1(double t, double x1, double x2, double u)
{
    return 3.0*x2*x2;
}

double fx2(double t, double x1, double x2, double u)
{
    return x1 + x2 - 2.0*u - t*t*t + 1.0;
}

double fp1(double t, double x1, double x2, double p1, double p2, double u)
{
    return 2.0 * (x1 - t*t*t) - p2;
}

double fp2(double t, double x1, double x2, double p1, double p2, double u)
{
    return 2.0 * (x2 - t) - 6.0 * x2 * p1 - p2;
}

double dIdu(double t, double x1, double x2, double p1, double p2, double u)
{
    return -4.0 * ( 2.0 * u - t ) - 2.0 *p2;
}

void print1(char *s, double *a, int n)
{
    int i;
    printf("%s =\t{", s);
    for (i=0; i<n; i++)
    {
        if ( i%((n-1)/10) == 0 )
            printf("%12.8f", a[i]);
        //if (i != n-1 )
        //	printf(", ");
    }
    printf("};");
    printf("\n");
}

void seperator()
{
    puts("---------------------------------------------------------------------------------------------------------------------------------------------");
}
