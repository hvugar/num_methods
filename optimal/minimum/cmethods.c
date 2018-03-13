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
                //p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                //q[i] = -c[i]/(b[i]+a[i]*q[i-1]);

                double m = (b[i]+a[i]*q[i-1]);
                p[i] = +(d[i]-a[i]*p[i-1])/m;
                q[i] = -c[i]/m;
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

/**
 * @brief Метод прогонки с трехдиагональной матрицей
 * @param a
 * @param b
 * @param c
 * @param d
 * @param x
 * @param N
 */
void tomasAlgorithmL2R(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int N)
{
    // Корректность и устойчивость метода. ----------------------------------------------------------------------
//    int isCorrect1 = 0;
//    int isCorrect2 = 0;
//    if (fabs(c[0]) >= 0.0 && fabs(a[N-1]) >= 0.0 && fabs(b[0]) > 0.0 && fabs(b[N-1]) > 0.0)
//    {
//        isCorrect1 = 1;
//        isCorrect1 &= fabs(b[0]) >= fabs(c[0]);
//        isCorrect1 &= fabs(b[N-1]) >= fabs(a[N-1]);

//        for (unsigned int i=1; i<N-1; i++)
//        {
//            isCorrect1 &= (fabs(b[i]) >= fabs(a[i]) + fabs(c[i]));
//            isCorrect2 |= (fabs(b[i]) >  fabs(a[i]) + fabs(c[i]));
//        }
//    }
//    if (isCorrect1 == 0 || isCorrect2 == 0)
//    {
//        printf("Not correct %d %d!\n", isCorrect1, isCorrect2);
//        return;
//    }

    double *alpha = (double*)malloc(sizeof(double)*N);
    double *betta = (double*)malloc(sizeof(double)*N);

    unsigned int i;
    //unsigned int N0 = N-1;
    unsigned int N1 = N-2;

    // Прямой ход метода прогонки. Определение прогоночных коэффициентов.----------------------------------------
    alpha[0] = -c[0]/b[0];
    betta[0] = +d[0]/b[0];
    for (i=1; i<=N1; i++)
    {
        double m = b[i] + a[i]*alpha[i-1];
        alpha[i] = -c[i]/m;
        betta[i] = +(d[i]-a[i]*betta[i-1])/m;
    }
    double m = b[N-1] + a[N-1]*alpha[N-2];
    alpha[N-1] = 0.0;
    betta[N-1] = (d[N-1]-a[N-1]*betta[N-2])/m;
    //-----------------------------------------------------------------------------------------------------------

    const unsigned int U_INT32_MAX = (unsigned int)0-1;
    // Обратный ход метода прогонки.-----------------------------------------------------------------------------
    // Обратный ход метода прогонки начинается с вычисления хn.--------------------------------------------------
    x[N-1] = betta[N-1];
    // Остальные значения неизвестных находятся рекуррентно.-----------------------------------------------------
    for (i=N1; i != U_INT32_MAX; i--)
    {
        x[i] = alpha[i]*x[i+1] + betta[i];
    }

    free(betta); betta=NULL;
    free(alpha); alpha=NULL;
}

/**
 * @brief Метод прогонки с трехдиагональной матрицей
 * @param a
 * @param b
 * @param c
 * @param d
 * @param x
 * @param N
 */
void tomasAlgorithmR2L(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int N)
{
//    // Корректность и устойчивость метода. ----------------------------------------------------------------------
//    int isCorrect1 = 0;
//    int isCorrect2 = 0;
//    if (fabs(c[0]) >= 0.0 && fabs(a[N-1]) >= 0.0 && fabs(b[0]) > 0.0 && fabs(b[N-1]) > 0.0)
//    {
//        isCorrect1 = 1;
//        isCorrect1 &= fabs(b[0]) >= fabs(c[0]);
//        isCorrect1 &= fabs(b[N-1]) >= fabs(a[N-1]);

//        for (unsigned int i=1; i<N-1; i++)
//        {
//            isCorrect1 &= (fabs(b[i]) >= fabs(a[i]) + fabs(c[i]));
//            isCorrect2 |= (fabs(b[i]) >  fabs(a[i]) + fabs(c[i]));
//        }
//    }
//    if (isCorrect1 == 0 || isCorrect2 == 0)
//    {
//        printf("Not correct %d %d!\n", isCorrect1, isCorrect2);
//        return;
//    }

    double *alpha = (double*)malloc(sizeof(double)*N);
    double *betta = (double*)malloc(sizeof(double)*N);

    unsigned int i;
    //unsigned int N0 = N-1;
    //unsigned int N1 = N-2;

    // Прямой ход метода прогонки.
    // Определение прогоночных коэффициентов.
    alpha[N-1] = -a[N-1]/b[N-1];
    betta[N-1] = +d[N-1]/b[N-1];
    for (i=N-1; i>=1; i--)
    {
        double m = b[i] + c[i]*alpha[i+1];
        alpha[i] = -a[i]/m;
        betta[i] = +(d[i]-c[i]*betta[i+1])/m;
    }
    double m = b[0] + c[0]*alpha[1];
    alpha[0] = 0.0;
    betta[0] = (d[0]-c[0]*betta[1])/m;

    // Обратный ход метода прогонки.
    // Обратный ход метода прогонки начинается с вычисления хn.
    x[0] = betta[0];
    // Остальные значения неизвестных находятся рекуррентно.
    for (i=1; i<N; i++)
    {
        x[i] = alpha[i]*x[i-1] +  betta[i];
    }

    free(betta); betta=NULL;
    free(alpha); alpha=NULL;
}

void qovmaE(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e, unsigned int *E, unsigned int L)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);
    double **k = (double**) malloc(sizeof(double*)*L);
    unsigned int s;
    for (s=0; s<L; s++) k[s] = (double*)malloc(sizeof(double)*n);

    unsigned int i;
    for (i=0; i<n; i++)
    {
        if (i == 0)
        {
            p[0] = +d[0]/b[0];
            q[0] = -c[0]/b[0];

            unsigned int s;
            for (s=0; s<L; s++)
            {
                k[s][0] = -e[s]/b[0];
            }
        }
        else if (i == (n-1))
        {
            p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
            q[i] = 0.0;

            unsigned int s;
            for (s=0; s<L; s++) k[s][i] = 0.0;
        }
        else
        {
            p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
            q[i] = -c[i]/(b[i]+a[i]*q[i-1]);

            unsigned int s;
            for (s=0; s<L; s++)
            {
                if (i<(E[s]-1))
                    k[s][i] = -(a[i]*k[s][i-1])/(b[i]+a[i]*q[i-1]);
                else
                    k[s][i] = 0.0;
            }

            for (s=0; s<L; s++) if (i==E[s]-1) q[i] += -(a[i]*k[s][i-1])/(b[i]+a[i]*q[i-1]);
        }
    }

    const unsigned int ui = (unsigned)0-1;
    for (i=n-1; i != ui; i--)
    //for (i=n-1; i != UINT_MAX; i--)
    {
        if (i==(n-1))
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];

            unsigned int s;
            for (s=0; s<L; s++)
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

    for (s=0; s<L; s++) free(k[s]);
    free(k);
}

void qovma2(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);
    double *k = (double*)malloc(sizeof(double)*n);

    unsigned int i1;
    for (i1=1; i1 <= n; i1++)
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

    unsigned int i;
    for (i=0; i < n; i++)
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
    unsigned int j;
    for (j=0; j<n; j++)
    {
        unsigned int i;
        for (i=0; i<n; i++)
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
            unsigned int k1;
            for (k1 = k+1; k1 < n; k1++)
            {
                if (a[k][k1] != 0.0)
                {
                    unsigned int k2;
                    for (k2=k; k2 < n; k2++)
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
