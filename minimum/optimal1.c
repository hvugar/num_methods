#include "optimal.h"

void smp1_control()
{
    int i;
    int n = 2;
    double h1 = 0.0;		//step for time
    double h2 = 0.000001;	//step for runga_kutta or euler method
    double t0 = 0.0;
    double t1 = 1.0;
    double x10 = 0.0;
    double x20 = 0.0;
    int M = 10;
	h1 = (t1-t0)/M;
    int N = M + 1;
	printf("%.10f\n", h1);

    double *t  = (double*) malloc( sizeof(double) * N );
    double *u  = (double*) malloc( sizeof(double) * N );
    double *x1 = (double*) malloc( sizeof(double) * N );
    double *x2 = (double*) malloc( sizeof(double) * N );
    double *p1 = (double*) malloc( sizeof(double) * N );
    double *p2 = (double*) malloc( sizeof(double) * N );
    double *gr = (double*) malloc( sizeof(double) * N );
	
	memset(t, 0, sizeof(double) * N);
	memset(u, 0, sizeof(double) * N);
	memset(x1, 0, sizeof(double) * N);
	memset(x2, 0, sizeof(double) * N);
	memset(p1, 0, sizeof(double) * N);
	memset(p2, 0, sizeof(double) * N);
	memset(gr, 0, sizeof(double) * N);
	
    for (i=0; i<N; i++)
    {
        t[i] = i*h1;
        u[i] = sin(t[i]);
    }

    double j1,j2;
    do
    {
        printX("u", u, N);
        x1[0] = x10;
        x2[0] = x20;
        for (i=0; i<N-1; i++)
        {
            double x0[] = { x1[i], x2[i] };
            double _x[] = { 0.0, 0.0 };
            smp1_RungaKuttaSystem1(t[i], x0, t[i+1], _x, n, h2, u[i]);
            x1[i+1] = _x[0];
            x2[i+1] = _x[1];
        }
        printX("x1", x1, N);
        printX("x2", x2, N);

        p1[N-1] = 0.0;
        p2[N-1] = 2 * (x2[N-1] - 1.0);        
        for (i=(N-1); i>0; i--)
        {
            double p0[] = { p1[i], p2[i] };
            double _x[] = { x1[i], x2[i] };
            double _p[] = { 0.0, 0.0 };
            smp1_RungaKuttaSystem2(t[i], p0, t[i-1], _p, n, h2, _x, u[i]);
            p1[i-1] = _p[0];
            p2[i-1] = _p[1];
        }
        printX("p1", p1, N);
        printX("p2", p2, N);
		
        for (i=0; i<N; i++)
        {
            double _x[] = { x1[i], x2[i] };
			double _p[] = { p1[i], p2[i] };
			gr[i] = smp1_Du(t[i], _x, _p, n, u[i]);
         }
        printX("gr", gr, N);

        double argmin1(double alpha)
        {
            int i;
			double u1[N];
            for (i=0; i<N; i++) u1[i] = u[i] - alpha * gr[i];
            return smp1_JSum(t, x1, x2, n, u1, N);
        }


        double a,b;
        double alpha0 = 0.0;
        straight_line_search_metod(argmin1, alpha0, 0.001, &a, &b);
        double alpha = golden_section_search_min(argmin1, a, b, 0.00000001);
        if ( argmin1(alpha) > argmin1(alpha0) ) alpha = alpha0;
        printf("alpha\t%14.10f\n", alpha);

		double U = 0.0;
		for (i=0;i<N;i++)
			U += u[i];
		printf("u\t%14.10f\n", U/N);

        //printf("J1\t%14.10f %14.10f\n", J1, JSum(u, N));
        j1 = smp1_JSum(t, x1, x2, n, u, N) - (x2[N-1]-1)*(x2[N-1]-1);
        for (i=0; i<N; i++)
        {
            u[i] = u[i] - alpha*gr[i];
        }
        j2 = smp1_JSum(t, x1, x2, n, u, N) - (x2[N-1]-1)*(x2[N-1]-1);
        printf("J1=%.16f\nJ2=%.16f\n", j1, j2);
		
		
		
        puts("--------------------------------");
		//if (fabs( j1 - j2 ) < 0.00001) break;
    } while (1);

    free(gr);
    free(p2);
    free(p1);
    free(x2);
    free(x1);
    free(u);
}

double smp1_f0(double t, double *x, int n, double u)
{
    return (x[0]-t*t*t)*(x[0]-t*t*t) + x[1]*x[1] - t*t + (2*u-1.0)*(2*u-1.0);
}

double smp1_f1(double t, double *x, int n, double u)
{
    return 3.0*x[1]*x[1];// - t;
}

double smp1_f2(double t, double *x, int n, double u)
{
    return x[0] + x[1] - 2.0*u - t*t*t + 1.0;
}

double smp1_Hamilton(double t, double *x, double *psi, int n, double u)
{
    return -1.0 * smp1_f0(t, x, n, u) + psi[0] * smp1_f1(t, x, n, u) + psi[1] * smp1_f2(t, x, n, u);
}

double smp1_Dpsi1(double t, double *x, double *psi, int n, double u)
{
    double _x1[n];
    double _x2[n];
    double h = 0.000001;
    _x1[0] = x[0] + h;
    _x1[1] = x[1];
    _x2[0] = x[0] - h;
    _x2[1] = x[1];
    double dH = -1.0*(smp1_Hamilton(t, _x1, psi, n, u) - smp1_Hamilton(t, _x2, psi, n, u)) / (2.0 * h);
    return dH;
//	return 2.0 * (x[0] - t*t*t) - psi[1];
}

double smp1_Dpsi2(double t, double *x, double *psi, int n, double u)
{
    double _x1[n];
    double _x2[n];
    double h = 0.000001;
    _x1[0] = x[0];
    _x1[1] = x[1] + h;
    _x2[0] = x[0];
    _x2[1] = x[1] - h;
    double dH = -1.0*(smp1_Hamilton(t, _x1, psi, n, u) - smp1_Hamilton(t, _x2, psi, n, u)) / (2.0 * h);
    return dH;
//	return 2.0 * x[1] - 6.0*x[1]*psi[0] - psi[1];
}

double smp1_Du(double t, double *x, double *psi, int n, double u)
{
    double h = 0.000001;
    double u1 = u + h;
    double u2 = u - h;
    double dH = +1.0*(smp1_Hamilton(t, x, psi, n, u1) - smp1_Hamilton(t, x, psi, n, u2)) / (2.0 * h);
    return dH;
//	return -2.0 * ( 2.0 * u + psi[1] - 1.0 );
}

double smp1_JSum(double *t, double *x1, double *x2, int n, double *u, int N)
{
	double sum = 0.0;
	int i=0;
	for (i=0; i<(N-1); i++)
	{
		int j=i+1;
	//	double x[n];
	//	x[0] = x1[j];
	//	x[1] = x2[j];
	//	double fj = smp1_f0(t[j], x, n, u[j]);
	//	x[0] = x1[i];
	//	x[1] = x2[i];
	//	double fi = smp1_f0(t[i], x, n, u[i]);
		double fj = (x1[j]-t[j]*t[j]*t[j])*(x1[j]-t[j]*t[j]*t[j]) + x2[j]*x2[j] - t[j]*t[j] + (2*u[j] - 1.0)*(2*u[j] - 1.0);
		double fi = (x1[i]-t[i]*t[i]*t[i])*(x1[i]-t[i]*t[i]*t[i]) + x2[i]*x2[i] - t[i]*t[i] + (2*u[i] - 1.0)*(2*u[i] - 1.0);
		sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
	}
	sum = sum + (x2[N-1]-1.0) * (x2[N-1]-1.0);
	return sum;
}

void smp1_RungaKuttaSystem1(double t0, double *x0, double t, double *x, int n, double h, double u)
{
    double xf1(double _t, double *_x, int _n) { return smp1_f1(_t, _x, _n, u); }
    double xf2(double _t, double *_x, int _n) { return smp1_f2(_t, _x, _n, u); }
    RmFunction xf[n];
    xf[0] = xf1;
    xf[1] = xf2;
    return RungaKuttaSystem(xf, t0, x0, t, x, n, h);
}

void smp1_RungaKuttaSystem2(double t0, double *p0, double t, double *p, int n, double h, double *x, double u)
{
    double __psi1(double _t, double *_p, int _n) { return smp1_Dpsi1(_t, x, _p, _n, u); }
    double __psi2(double _t, double *_p, int _n) { return smp1_Dpsi2(_t, x, _p, _n, u); }
    RmFunction fp[n];
    fp[0] = __psi1;
    fp[1] = __psi2;
    return RungaKuttaSystem(fp, t0, p0, t, p, n, h);
}
