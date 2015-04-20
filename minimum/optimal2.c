#include "optimal.h"

void printX1(char *label, double *x, int n);

void smp1_control()
{
    int i;
    int n = 2;
    double h1 = 0.001;		//step for time
    double h2 = 0.001;	//step for runga_kutta or euler method
    double t0 = 0.0;
    double t1 = 1.0;
    double x10 = 0.0;
    double x20 = 0.0;
    int M = 1000;
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
        u[i] = 0.5*t[i];
    }

    double j1,j2;
    do
    {
        x1[0] = x10;
        x2[0] = x20;
        {
            double *k1 = (double*) malloc( sizeof(double) * n );
            double *k2 = (double*) malloc( sizeof(double) * n );
            double *k3 = (double*) malloc( sizeof(double) * n );
            double *k4 = (double*) malloc( sizeof(double) * n );
            double *xc = (double*) malloc( sizeof(double) * n );

			double _t0 = t0;
			double _t1 = t1;
            int i=0;
            while ((_t1-_t0) >= h1)
            {
				double x[] = {x1[i], x2[i]};
			
                // Calculating k1 vector
                k1[0] = smp1_f1(_t0, x, n, u[i]);
                k1[1] = smp1_f2(_t0, x, n, u[i]);

                // Calculating k2 vector
                xc[0] = x[0] + (h1/2.0) * k1[0];
                xc[1] = x[1] + (h1/2.0) * k1[1];

                k2[0] = smp1_f1(_t0+h1/2.0, xc, n, u[i]);
                k2[1] = smp1_f2(_t0+h1/2.0, xc, n, u[i]);

                // Calculating k3 vector
                xc[0] = x[0] + (h1/2.0) * k2[0];
                xc[1] = x[1] + (h1/2.0) * k2[1];

                k3[0] = smp1_f1(_t0+h1/2.0, xc, n, u[i]);
                k3[1] = smp1_f2(_t0+h1/2.0, xc, n, u[i]);

                // Calculating k1 vector
                xc[0] = x[0] + h1 * k3[0];
                xc[1] = x[1] + h1 * k3[1];

                k4[0] = smp1_f1(_t0+h1, xc, n, u[i]);
                k4[1] = smp1_f2(_t0+h1, xc, n, u[i]);

                // Calculating y
                x1[i+1] = x1[i] + (h1/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
                x2[i+1] = x2[i] + (h1/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);

                _t0 = _t0 + h1;
				i++;
            }

            free(xc);
            free(k1);
            free(k2);
            free(k3);
            free(k4);
        }

        //        for (i=0; i<N-1; i++)
        //        {
        //            double _x0[] = { x1[i], x2[i] };
        //            double _x1[] = { 0.000, 0.000 };
        //            smp1_RungaKuttaSystem1(t[i], _x0, t[i+1], _x1, n, h2, u[i]);
        //            x1[i+1] = _x1[0];
        //            x2[i+1] = _x1[1];
        //        }
		
        p1[N-1] = 0.0;
        p2[N-1] = -2*(x2[N-1] - 1.0);
        {
            double *k1 = (double*) malloc( sizeof(double) * n );
            double *k2 = (double*) malloc( sizeof(double) * n );
            double *k3 = (double*) malloc( sizeof(double) * n );
            double *k4 = (double*) malloc( sizeof(double) * n );
            double *pc = (double*) malloc( sizeof(double) * n );

			double _t0 = t0;
			double _t1 = t1;
            int i=N-1;
            while ((_t1-_t0) >= h1)
            {
				double x[] = {x1[i], x2[i]};
				double p[] = {p1[i], p2[i]};
			
                // Calculating k1 vector
                k1[0] = smp1_Dpsi1(_t1, x, p, n, u[i]);
                k1[1] = smp1_Dpsi2(_t1, x, p, n, u[i]);

                // Calculating k2 vector
                pc[0] = p[0] + (h1/2.0) * k1[0];
                pc[1] = p[1] + (h1/2.0) * k1[1];

                k2[0] = smp1_Dpsi1(_t1+h1/2.0, pc, x, n, u[i]);
                k2[1] = smp1_Dpsi2(_t1+h1/2.0, pc, x, n, u[i]);

                // Calculating k3 vector
                pc[0] = p[0] + (h1/2.0) * k2[0];
                pc[1] = p[1] + (h1/2.0) * k2[1];

                k3[0] = smp1_Dpsi1(_t1+h1/2.0, pc, x, n, u[i]);
                k3[1] = smp1_Dpsi2(_t1+h1/2.0, pc, x, n, u[i]);

                // Calculating k1 vector
                pc[0] = p[0] + h1 * k3[0];
                pc[1] = p[1] + h1 * k3[1];

                k4[0] = smp1_Dpsi1(_t1+h1, pc, x, n, u[i]);
                k4[1] = smp1_Dpsi2(_t1+h1, pc, x, n, u[i]);

                // Calculating y
                p1[i-1] = p1[i] + (h1/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
                p2[i-1] = p2[i] + (h1/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);

                _t1 = _t1 - h1;
				i--;
            }

            free(pc);
            free(k1);
            free(k2);
            free(k3);
            free(k4);
        }
		printX1("u", u, N);
        printX1("x1", x1, N);
        printX1("x2", x2, N);

//		double _x[] = { x1[N-1], x2[N-1] };
//		double xg[] = { 0.00000, 0.00000 };
//		gradient(smp1_F, _x, n, 0.000001, xg);		
//		p1[N-1] = -xg[0];
//        p2[N-1] = -xg[1];
//        for (i=(N-1); i>0; i--)
//        {
//            double p0[] = { p1[i], p2[i] };
//            double _x[] = { x1[i], x2[i] };
//            double _p[] = { 0.0, 0.0 };
//            smp1_RungaKuttaSystem2(t[i], p0, t[i-1], _p, n, h2, _x, u[i]);
//            p1[i-1] = _p[0];
//            p2[i-1] = _p[1];
//        }
        printX1("p1", p1, N);
        printX1("p2", p2, N);

        for (i=0; i<N; i++)
        {
            double _x[] = { x1[i], x2[i] };
            double _p[] = { p1[i], p2[i] };
            gr[i] = smp1_Du(t[i], _x, _p, n, u[i]);
        }
        printX1("gr", gr, N);

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
		double _x[] = { x1[N-1], x2[N-1] };
        j1 = smp1_JSum(t, x1, x2, n, u, N);// - smp1_F(_x, n);
        for (i=0; i<N; i++)
        {
            u[i] = u[i] - alpha*gr[i];
        }
        j2 = smp1_JSum(t, x1, x2, n, u, N);// - smp1_F(_x, n);
        //printX1("u", u, N);
        printf("J1=%.16f\nJ2=%.16f %d\n", j1, j2, (j1-j2) > 0.0);

        puts("--------------------------------");
        //if (fabs( j1 - j2 ) < 0.00001) break;
        //break;
    } while ((j1-j2)>0.0);
	//} while(1);

    free(gr);
    free(p2);
    free(p1);
    free(x2);
    free(x1);
    free(u);
}

double smp1_F(double *x, int n)
{
    return (x[1] - 1.0) * (x[1] - 1.0);
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
    //double h = 0.0000001;
    //double _x1[] = { x[0] + h, x[1] };
    //double _x2[] = { x[0] - h, x[1] };
    //return -1.0*((smp1_Hamilton(t, _x1, psi, n, u) - smp1_Hamilton(t, _x2, psi, n, u)) / (2.0 * h));
    return 2.0 * (x[0] - t*t*t) - psi[1];
}

double smp1_Dpsi2(double t, double *x, double *psi, int n, double u)
{
    //double h = 0.0000001;
    //double _x1[] = { x[0], x[1] + h };
    //double _x2[] = { x[0], x[1] - h };
    //return -1.0*((smp1_Hamilton(t, _x1, psi, n, u) - smp1_Hamilton(t, _x2, psi, n, u)) / (2.0 * h));
    return 2.0 * x[1] - 6.0*x[1]*psi[0] - psi[1];
}

double smp1_Du(double t, double *x, double *psi, int n, double u)
{
    //double h = 0.0000001;
    //double u1 = u + h;
    //double u2 = u - h;
    //return +1.0*(smp1_Hamilton(t, x, psi, n, u1) - smp1_Hamilton(t, x, psi, n, u2)) / (2.0 * h);
    return +2.0 * ( 2.0 * u + psi[1] - 1.0 );
}

double smp1_dFdx1(double *x, int n)
{
    double h = 0.0000001;
    double _x1[] = { x[0] + h, x[1] };
    double _x2[] = { x[0] - h, x[1] };
    return -1.0*(smp1_F(_x1, n) - smp1_F(_x2, n)) / (2.0 * h);
}

double smp1_dFdx2(double *x, int n)
{
    double h = 0.0000001;
    double _x1[] = { x[0], x[1] + h };
    double _x2[] = { x[0], x[1] - h };
    return -1.0*(smp1_F(_x1, n) - smp1_F(_x2, n)) / (2.0 * h);
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
    //sum = sum + (x2[N-1]-1.0) * (x2[N-1]-1.0);
    double x[] = { x1[N-1], x2[N-1] };
    sum = sum + smp1_F(x, n);
    return sum;
}

void smp1_RungaKuttaSystem1(double t0, double *x0, double t, double *x, int n, double h, double u)
{
    //    double xf1(double _t, double *_x, int _n) { return smp1_f1(_t, _x, _n, u); }
    //    double xf2(double _t, double *_x, int _n) { return smp1_f2(_t, _x, _n, u); }
    //    RmFunction xf[] = { xf1, xf2 };
    //    RungaKuttaSystem(xf, t0, x0, t, x, n, h);

    memcpy(x, x0, sizeof(double)*n);
    if (fabs(t-t0) < fabs(h)) return;

    double *k1 = (double*) malloc( sizeof(double) * n );
    double *k2 = (double*) malloc( sizeof(double) * n );
    double *k3 = (double*) malloc( sizeof(double) * n );
    double *k4 = (double*) malloc( sizeof(double) * n );
    double *xc = (double*) malloc( sizeof(double) * n );

    if (t0 > t) h = -fabs(h);

    //while (fabs(t-t0)>=(fabs(h)))
    {
        int i=0;
        // Calculating k1 vector
        k1[0] = smp1_f1(t0, x, n, u);
        k1[1] = smp1_f2(t0, x, n, u);

        // Calculating k2 vector
        xc[0] = x[0] + (h/2.0) * k1[i];
        xc[1] = x[1] + (h/2.0) * k1[1];

        k2[0] = smp1_f1(t0+h/2.0, xc, n, u);
        k2[1] = smp1_f2(t0+h/2.0, xc, n, u);

        // Calculating k3 vector
        xc[0] = x[0] + (h/2.0) * k2[0];
        xc[1] = x[1] + (h/2.0) * k2[1];
        
        k3[0] = smp1_f1(t0+h/2.0, xc, n, u);
        k3[1] = smp1_f2(t0+h/2.0, xc, n, u);

        // Calculating k1 vector
        xc[0] = x[0] + h * k3[0];
        xc[1] = x[1] + h * k3[1];
        
        k4[0] = smp1_f1(t0+h, xc, n, u);
        k4[1] = smp1_f2(t0+h, xc, n, u);

        // Calculating y
        x[0] = x[0] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        x[1] = x[1] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);

        t0 = t0 + h;
    }

    free(xc);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
}

void smp1_RungaKuttaSystem2(double t0, double *p0, double t, double *p, int n, double h, double *x, double u)
{
    //    double __psi1(double _t, double *_p, int _n) { return smp1_Dpsi1(_t, x, _p, _n, u); }
    //    double __psi2(double _t, double *_p, int _n) { return smp1_Dpsi2(_t, x, _p, _n, u); }
    //    RmFunction fp[] = { __psi1, __psi2 };
    //    RungaKuttaSystem(fp, t0, p0, t, p, n, h);

    memcpy(p, p0, sizeof(double)*n);
    if (fabs(t-t0) < fabs(h)) return;

    double *k1 = (double*) malloc( sizeof(double) * n );
    double *k2 = (double*) malloc( sizeof(double) * n );
    double *k3 = (double*) malloc( sizeof(double) * n );
    double *k4 = (double*) malloc( sizeof(double) * n );
    double *pc = (double*) malloc( sizeof(double) * n );

    if (t0 > t) h = -fabs(h);

    //while (fabs(t-t0)>=(fabs(h)))
    {
        int i=0;
        // Calculating k1 vector
        k1[0] = smp1_Dpsi1(t0, x, p, n, u);
        k1[1] = smp1_Dpsi2(t0, x, p, n, u);

        // Calculating k2 vector
        pc[0] = p[0] + (h/2.0) * k1[i];
        pc[1] = p[1] + (h/2.0) * k1[1];
        k2[0] = smp1_Dpsi1(t0+h/2.0, x, pc, n, u);
        k2[1] = smp1_Dpsi2(t0+h/2.0, x, pc, n, u);

        // Calculating k3 vector
        pc[0] = p[0] + (h/2.0) * k2[0];
        pc[1] = p[1] + (h/2.0) * k2[1];
        k3[0] = smp1_Dpsi1(t0+h/2.0, x, pc, n, u);
        k3[1] = smp1_Dpsi2(t0+h/2.0, x, pc, n, u);

        // Calculating k1 vector
        pc[0] = p[0] + h * k3[0];
        pc[1] = p[1] + h * k3[1];
        k4[0] = smp1_Dpsi1(t0+h, x, pc, n, u);
        k4[1] = smp1_Dpsi2(t0+h, x, pc, n, u);

        // Calculating y
        p[0] = p[0] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        p[1] = p[1] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);

        t0 = t0 + h;
    }

    free(pc);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
}

void printX1(char *label, double *x, int n)
{
	int m = n / 10;
	int i;
	printf("double %s[] = \t{", label);
	for (i=0; i<n; i++)
	{
		//if (i%10000==0)
			printf("%12.8f", x[i]);
		if (i != n-1 )
			printf(", ");
	}
	printf("};");
	printf("\n");
}
