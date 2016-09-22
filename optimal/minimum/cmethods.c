#include "cmethods.h"
#include <float.h>

double derivative(CR1Function fx)
{
    return fx(5.0);
}

double derivative1(CR1Function fx, double x, double h)
{
    return (fx(x+h) - fx(x)) / h;
}

double derivative2(CR1Function fx, double x, double h)
{
    return (fx(x) - fx(x-h)) / h;
}

double derivative3(CR1Function fx, double x, double h)
{
    return (fx(x+h) - fx(x-h)) / (2.0*h);
}

void gradient1(CRnFunction fx, double *x, double *g, unsigned int n, double h)
{
    unsigned int i;
    double tx;
    double f0 = fx(x, n);
    for (i=0; i<n; i++)
    {
        tx = x[i];
        x[i] += h;
        double f = fx(x, n);
        g[i] = (f-f0)/h;
        x[i] = tx;
    }
}

void gradient2(CRnFunction fx, double *x, double *g, unsigned int n, double h)
{
    unsigned int i;
    double tx;
    double f0 = fx(x, n);
    for (i=0; i<n; i++)
    {
        tx = x[i];
        x[i] -= h;
        double f = fx(x, n);
        g[i] = (f-f0)/h;
        x[i] = tx;
    }
}

void gradient3(CRnFunction fx, double *x, double *g, unsigned int n, double h)
{
    unsigned int i;
    double tx;
    for (i=0; i<n; i++)
    {
        tx = x[i];
        x[i] = tx + h;
        double f1 = fx(x, n);
        x[i] = tx - h;
        double f2 = fx(x, n);
        g[i] = (f2-f1)/(2.0*h);
        x[i] = tx;
    }
}

double trapesium1(CR1Function fx, unsigned int n, double a, double b)
{
    double sum = 0.0;
    double h = (a-b)/n;
    unsigned int i;
    for (i=0; i<=n-1; i++)
    {
        double x = a + i*h;
        double f1 = fx(x);
        double f2 = fx(x+h);
        sum = sum + (f1+f2);
    }
    sum = (h/2.0)*sum;
    return sum;
}

double trapesium2(CR1Function fx, double h, double a, double b)
{
    unsigned int n = (unsigned int)round((b - a)/h);
    double sum = 0.0;
    unsigned int i;
    for (i=0; i<=(n-1); i++)
    {
        double x = a + i*h;
        double f1 = fx(x);
        double f2 = fx(x+h);
        sum = sum + (f1+f2);
    }
    sum = (h/2.0)*sum;
    return sum;
}

void tomasAlgorithm(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int n)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);

    const unsigned int j = (unsigned)0-1;
    unsigned int i;
    for (i=0; i<n; i++)
    {
        if (i==0)
        {
            p[0] = +d[0]/b[0];
            q[0] = -c[0]/b[0];
        }
        else
        {
            if(i==(n-1))
            {
                p[n-1] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[n-1] = 0.0;
            }
            else
            {
                p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[i] = -c[i]/(b[i]+a[i]*q[i-1]);
            }
        }
    }

    for (i=n-1; i != j; i--)
    {
        if (i==(n-1))
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];
        }
    }

    free(p);
    free(q);
}

void qovmaE(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e, unsigned int *E, unsigned int L)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);
    double **k = (double**) malloc(sizeof(double*)*L);
    for (unsigned int s=0; s<L; s++) k[s] = (double*)malloc(sizeof(double)*n);

    for (unsigned int i=0; i<n; i++)
    {
        if (i == 0)
        {
            p[0] = +d[0]/b[0];
            q[0] = -c[0]/b[0];

            for (unsigned int s=0; s<L; s++)
            {
                k[s][0] = -e[s]/b[0];
            }
        }
        else if (i == (n-1))
        {
            p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
            q[i] = 0.0;

            for (unsigned int s=0; s<L; s++) k[s][i] = 0.0;
        }
        else
        {
            p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
            q[i] = -c[i]/(b[i]+a[i]*q[i-1]);

            for (unsigned int s=0; s<L; s++)
            {
                if (i<(E[s]-1))
                    k[s][i] = -(a[i]*k[s][i-1])/(b[i]+a[i]*q[i-1]);
                else
                    k[s][i] = 0.0;
            }

            for (unsigned int s=0; s<L; s++) if (i==E[s]-1) q[i] += -(a[i]*k[s][i-1])/(b[i]+a[i]*q[i-1]);
        }
    }

    const unsigned int j = (unsigned)0-1;
    for (unsigned int i=n-1; i != j; i--)
    {
        if (i==(n-1))
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];

            for (unsigned int s=0; s<L; s++)
            {
                if (i<=E[s]-1)
                {
                    x[i] = x[i] + k[s][i]*x[E[s]];
                }
            }
        }
    }

    free(q);
    free(q);
    for (unsigned int s=0; s<L; s++) free(k[s]);
    free(k);
}

void qovma2(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);
    double *k = (double*)malloc(sizeof(double)*n);

    const unsigned int j = (unsigned)0-1;

    for (unsigned int i1=1; i1 <= n; i1++)
    {
        unsigned int i = n - i1;

        if (i == n-1)
        {
            p[i] = +d[i]/b[i];
            q[i] = -a[i]/b[i];
            k[i] = -e[i]/b[i];
        }
        else if (i == 1)
        {
            double m = b[i]+c[i]*q[i+1];
            p[i] = +(d[i]-c[i]*p[i+1])/m;
            q[i] = -(a[i]+c[i]*k[i+1])/m;
            k[i] = 0.0;
        }
        else if (i == 0)
        {
            double m = b[i]+c[i]*q[i+1];
            p[i] = +(d[i]-c[i]*p[i+1])/m;
            q[i] = 0.0;
            k[i] = 0.0;
        }
        else
        {
            double m = b[i]+c[i]*q[i+1];
            p[i] = +(d[i]-c[i]*p[i+1])/m;
            q[i] = -a[i]/m;
            k[i] = -(e[i]+c[i]*k[i+1])/m;
        }
    }
//    printf("%10.6f %10.6f %10.6f\n", p[6], q[6], k[6]);
//    printf("%10.6f %10.6f %10.6f\n", p[5], q[5], k[5]);
//    printf("%10.6f %10.6f %10.6f\n", p[4], q[4], k[4]);
//    printf("%10.6f %10.6f %10.6f\n", p[3], q[3], k[3]);
//    printf("%10.6f %10.6f %10.6f\n", p[2], q[2], k[2]);
//    printf("%10.6f %10.6f %10.6f\n", p[1], q[1], k[1]);
//    printf("%10.6f %10.6f %10.6f\n", p[0], q[0], k[0]);

    for (unsigned int i=0; i < n; i++)
    {
        if (i==0)
        {
            x[i] = p[i];
        }
        else if (i==1)
        {
            x[i] = p[i] + q[i]*x[i-1];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i-1] + k[i]*x[0];
        }
    }

    free(k);
    free(p);
    free(q);
}

void printMat(double **a, double *b, unsigned int n)
{
    for (unsigned int j=0; j<n; j++)
    {
        for (unsigned int i=0; i<n; i++)
            printf("%10.4f ", a[j][i]);
        printf("| %10.4f\n", b[j]);
    }
}

//void changeMatrixRow(double **a, double *b, unsigned int n)
//{
//}

void gaussianElimination(double **a, double *b, double *x, unsigned int n)
{
    C_UNUSED(x);

    unsigned int k;
    unsigned int j;
    unsigned int i;
    const unsigned int ui = (unsigned)0-1;


    for (k=0; k<n-1; k++)
    {
        //a[k][k] = 0.0;
        //        printf("\nIteration: %d matrix status\n", k);
        //        printMat(a, b, n);
        //        puts("matrix printed.");

        if (fabs(a[k][k]) <= DBL_EPSILON)
        {
            for (unsigned int k1 = k+1; k1 < n; k1++)
            {
                if (a[k][k1] != 0.0)
                {
                    //                    printf("Changing row: %d to row %d\n", k, k1);
                    for (unsigned int k2=k; k2 < n; k2++)
                    {
                        double _a = a[k1][k2];
                        a[k1][k2] = a[k][k2];
                        a[k][k2] = _a;
                    }
                    double _b = b[k1];
                    b[k1] = b[k];
                    b[k] = _b;
                    //                    printMat(a, b, n);
                    //                    printf("Row: %d changed to row %d\n", k, k1);
                    break;
                }
            }
        }

        //        printf("Relaxing\n");
        for (j=(k+1); j<n; j++)
        {
            double c = a[j][k]/a[k][k];
            for (i=k; i<n; i++) a[j][i] = a[j][i] - a[k][i] * c;
            b[j] = b[j] - b[k] *c;
        }
        //        printMat(a, b, n);
        //        printf("Relaxed\n");
    }

    for (i=(n-1); i!=ui; i--)
    {
        for (j=(n-1); j>i; j--) b[i] -= (a[i][j] * x[j]);
        x[i] = b[i] / a[i][i];
    }
    //    for (unsigned int i=0; i<n; i++) x[i] = 0.0;
}

void gaussJordanElimination(double **a, double *b, double *x, unsigned int n)
{
    C_UNUSED(a); C_UNUSED(b); C_UNUSED(x); C_UNUSED(n);
}

void euler(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq)
{
    C_UNUSED(x0); C_UNUSED(y0); C_UNUSED(xN); C_UNUSED(yN); C_UNUSED(N); C_UNUSED(x); C_UNUSED(y); C_UNUSED(eq);
}

void eulerMod(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq)
{
    C_UNUSED(x0); C_UNUSED(y0); C_UNUSED(xN); C_UNUSED(yN); C_UNUSED(N); C_UNUSED(x); C_UNUSED(y); C_UNUSED(eq);
}

void runge_kutta_rk3(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq)
{
    C_UNUSED(x0); C_UNUSED(y0); C_UNUSED(xN); C_UNUSED(yN); C_UNUSED(N); C_UNUSED(x); C_UNUSED(y); C_UNUSED(eq);
}

void runge_kutta_rk4(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq)
{
    C_UNUSED(yN);

    if (N == 0) return;

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    double h = (xN - x0) / N;
    //    x = (double*)malloc(sizeof(double)*(N+1));
    //    y = (double*)malloc(sizeof(double)*(N+1));

    if (h > 0)
    {
        x[0] = x0;
        y[0] = y0;
        unsigned int i = 1;
        for (i=1; i<=N; i++)
        {
            k1 = eq(x0, y0);
            k2 = eq(x0+h/2.0, y0+(h/2.0)*k1);
            k3 = eq(x0+h/2.0, y0+(h/2.0)*k2);
            k4 = eq(x0+h, y0+h*k3);
            y0 = y0 + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + h;
            x[i] = x0;
            y[i] = y0;
        }
        yN = y[N];
    }

    if (h < 0)
    {
        x[N] = x0;
        y[N] = y0;
        int i = N-1;
        for (i=N-1; i>=0; i--)
        {
            k1 = eq(x0,       y0);
            k2 = eq(x0+h/2.0, y0+(h/2.0)*k1);
            k3 = eq(x0+h/2.0, y0+(h/2.0)*k2);
            k4 = eq(x0+h,     y0+h*k3);
            y0 = y0 + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + h;
            x[i] = x0;
            y[i] = y0;
        }
        yN = y[0];
    }
}

void runge_kutta_rk4_system(double x0, double x1, double *y0, double **y, size_t n, unsigned int N, double h, ODE1stOrderEquationN *eq)
{
    C_UNUSED(x1);

    double *k1 = (double *)malloc(sizeof(double) * n);
    double *k2 = (double *)malloc(sizeof(double) * n);
    double *k3 = (double *)malloc(sizeof(double) * n);
    double *k4 = (double *)malloc(sizeof(double) * n);

    unsigned int i, j;
    for (i=0; i<n; i++) y[i][0] = y0[i];

    for (i=0; i<N; i++)
    {
        double *arg = (double*) malloc(sizeof(double) * n);

        // Calculating k1 vector
        for (j = 0; j<n; j++) arg[j] = y[j][i];
        for (j = 0; j<n; j++) k1[j] = eq[j](x0, arg, n);

        // Calculating k2 vector
        for (j = 0; j<n; j++) arg[j] = y[j][i] + (h/2.0) * k1[j];
        for (j = 0; j<n; j++) k2[j] = eq[j](x0+h/2.0, arg, n);

        // Calculating k3 vector
        for (j = 0; j<n; j++) arg[j] = y[j][i] + (h/2.0) * k2[j];
        for (j = 0; j<n; j++) k3[j] = eq[j](x0+h/2.0, arg, n);

        // Calculating k4 vector
        for (j = 0; j<n; j++) arg[j] = y[j][i] + h * k3[j];
        for (j = 0; j<n; j++) k4[j] = eq[j](x0+h, arg, n);

        // Calculating y
        for (j = 0; j<n; j++) y[j][i+1] = y[j][i] + (h/6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

        x0 = x0 + h;

        free(arg);
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
}


//void RungaKuttaSystem(RmFunction *f, double x0, const double *y0, double x, double *y, const int n, double h)
//{
//    //h = 0.000001;
//    memcpy(y, y0, sizeof(double)*n);
//    if (fabs(x-x0) < fabs(h)) return;

//    double *k1 = (double*) malloc( sizeof(double) * n );
//    double *k2 = (double*) malloc( sizeof(double) * n );
//    double *k3 = (double*) malloc( sizeof(double) * n );
//    double *k4 = (double*) malloc( sizeof(double) * n );
//    double *yc = (double*) malloc( sizeof(double) * n );

//    if (x0 > x) h = -fabs(h);

//    while (fabs(x-x0)>=(fabs(h)))
//    {
//        int i=0;
//        // Calculating k1 vector
//        for (i = 0; i<n; i++) k1[i] = f[i](x0, y, n);

//        // Calculating k2 vector
//        for (i = 0; i<n; i++) yc[i] = y[i] + (h/2.0) * k1[i];
//        for (i = 0; i<n; i++) k2[i] = f[i](x0+h/2.0, yc, n);

//        // Calculating k3 vector
//        for (i = 0; i<n; i++) yc[i] = y[i] + (h/2.0) * k2[i];
//        for (i = 0; i<n; i++) k3[i] = f[i](x0+h/2.0, yc, n);

//        // Calculating k1 vector
//        for (i = 0; i<n; i++) yc[i] = y[i] + h * k3[i];
//        for (i = 0; i<n; i++) k4[i] = f[i](x0+h, yc, n);

//        // Calculating y
//        for (i = 0; i<n; i++) y[i] = y[i] + (h/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

//        x0 = x0 + h;
//    }

//    free(yc);
//    free(k1);
//    free(k2);
//    free(k3);
//    free(k4);
//}
