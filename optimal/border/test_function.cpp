#include "test_function.h"

#define FUNCTION_T1
#define FUNCTION_X1
#define FUNCTION_Y1

double TestFunction::u(const TimeNodePDE &tn, const SpaceNodePDE &sn, Derivative derivative, const Dimension &dimX, const Dimension &dimY)
{
    double res = 0.0;
    switch (derivative)
    {
    case FunctionValue:
    {
        res = 0.0;
#ifdef FUNCTION_T1
        res += tn.t;
#endif
#ifdef FUNCTION_T2
        res += tn.t*tn.t;
#endif
#ifdef FUNCTION_T3
        res += tn.t*tn.t*tn.t;
#endif
/******************************************/
#ifdef FUNCTION_X1
        res += sn.x;
#endif
#ifdef FUNCTION_X2
        res += sn.x*sn.x;
#endif
#ifdef FUNCTION_X3
        res += sn.x*sn.x*sn.x;
#endif
/******************************************/
#ifdef FUNCTION_Y1
        res += sn.y;
#endif
#ifdef FUNCTION_Y2
        res += sn.y*sn.y;
#endif
#ifdef FUNCTION_Y3
        res += sn.y*sn.y*sn.y;
#endif

    } break;
    case TimeFirstDerivative:
    {
#ifdef FUNCTION_T1
        res = 1.0;
#endif
#ifdef FUNCTION_T2
        res = 2.0*tn.t;
#endif
#ifdef FUNCTION_T3
        res = 3.0*tn.t*tn.t;
#endif
    } break;
    case TimeSecondDerivative:
    {
#ifdef FUNCTION_T1
        res = 0.0;
#endif
#ifdef FUNCTION_T2
        res = 2.0;
#endif
#ifdef FUNCTION_T3
        res = 6.0*tn.t;
#endif
    } break;
    case SpaceFirstDerivativeX:
    {
#ifdef FUNCTION_X1
        res = 1.0;
#endif
#ifdef FUNCTION_X2
        res = 2.0*sn.x;
#endif
#ifdef FUNCTION_X3
        res = 3.0*sn.x*sn.x;
#endif
    } break;
    case SpaceSecondDerivativeX:
    {
#ifdef FUNCTION_X1
        res = 0.0;
#endif
#ifdef FUNCTION_X2
        res = 2.0;
#endif
#ifdef FUNCTION_X3
        res = 6.0*sn.x;
#endif
    } break;
    case SpaceFirstDerivativeY:
    {
#ifdef FUNCTION_Y1
        res = 1.0;
#endif
#ifdef FUNCTION_Y2
        res = 2.0*sn.y;
#endif
#ifdef FUNCTION_Y3
        res = 3.0*sn.y*sn.y;
#endif
    } break;
    case SpaceSecondDerivativeY:
    {
#ifdef FUNCTION_Y1
        res = 0.0;
#endif
#ifdef FUNCTION_Y2
        res = 2.0;
#endif
#ifdef FUNCTION_y3
        res = 6.0*sn.y;
#endif
    } break;
    case SpaceNorm:
    {
#ifdef FUNCTION_X1
        if (sn.i == dimX.min()) res = -1.0;
        if (sn.i == dimX.max()) res = +1.0;
#endif
#ifdef FUNCTION_X2
        if (sn.i == dimX.min()) res = -2.0*sn.x;
        if (sn.i == dimX.max()) res = +2.0*sn.x;
#endif
#ifdef FUNCTION_X3
        if (sn.i == dimX.min()) res = -3.0*sn.x*sn.x;
        if (sn.i == dimX.max()) res = +3.0*sn.x*sn.x;
#endif
#ifdef FUNCTION_Y1
        if (sn.j == dimY.min()) res = -1.0;
        if (sn.j == dimY.max()) res = +1.0;
#endif
#ifdef FUNCTION_Y2
        if (sn.j == dimY.min()) res = -2.0*sn.y;
        if (sn.j == dimY.max()) res = +2.0*sn.y;
#endif
#ifdef FUNCTION_Y3
        if (sn.j == dimY.min()) res = -3.0*sn.y*sn.y;
        if (sn.j == dimY.max()) res = +3.0*sn.y*sn.y;
#endif
    } break;
    }

    return res;
}

double test1()
{
    srand(time(nullptr));
    const size_t N = 1000;
    const double ht = 0.001;
    const double hx = 0.001;

    double *a = new double[N+1];
    double *b = new double[N+1];
    double *c = new double[N+1];
    double *d = new double[N+1];
    double *x = new double[N+1];
    double *u = new double[N+1];
    double **w = new double*[N+1];

    for (size_t i=0; i<=N; i++)
    {
        a[i] = -ht/(hx*hx);
        b[i] = +1.0+2.0*ht/(hx*hx);
        c[i] = -ht/(hx*hx);

        x[i] = (double)(rand() % 100) * 0.01;

        w[i] = new double[N+1];
        for (size_t j=0; j<=N; j++) w[i][j] = 0.0;
    }

    //w[0][30] = +325.0; w[0][60] = -128.0; w[0][90] = +229.0;
    //w[N][20] = -127.0; w[N][40] = +285.0; w[N][80] = +568.0;

    for (size_t i=0; i<=N; i++) { w[i][30] = rand()%100; w[i][50] = rand()%100; w[i][70] = rand()%100; }

    a[0] = 0.0;
    c[N] = 0.0;

    d[0] = b[0]*x[0] + c[0]*x[1];
    for (size_t i=1; i<N; i++)
    {
        d[i] = a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1];
    }
    d[N] = a[N]*x[N-1] + b[N]*x[N];

    for (size_t j=0; j<=N; j++)
    {
        for (size_t i=0; i<=N; i++) { d[j] += w[j][i] * x[i]; }
        //for (size_t i=0; i<=N; i++) { d[N] += w[N][i] * x[i]; }
    }

    FILE *file = fopen("d:/data.txt", "w");

    for (size_t i=0; i<=N; i++) { fprintf(file, "%10.2f", a[i]); } fputs("\n", file);
    for (size_t i=0; i<=N; i++) { fprintf(file, "%10.2f", b[i]); } fputs("\n", file);
    for (size_t i=0; i<=N; i++) { fprintf(file, "%10.2f", c[i]); } fputs("\n", file);
    fputs("---\n", file);
    for (size_t i=0; i<=N; i++) { fprintf(file, "%10.4f", d[i]); } fputs("\n", file);
    fputs("---\n", file);
    for (size_t i=0; i<=N; i++) { fprintf(file, "%10.6f", w[0][i]); } fputs("\n", file);
    for (size_t i=0; i<=N; i++) { fprintf(file, "%10.6f", w[N][i]); } fputs("\n", file);

    LinearEquation::func1(a, b, c, d, w, u, N+1);

    DoubleMatrix M(N+1, N+1, 0.0);
    DoubleVector B(N+1, 0.0);
    DoubleVector v(N+1);

    for (size_t i=0; i<=N; i++)
    {
        for (size_t j=0; j<=N; j++)
        {
            M[i][j] = w[i][j];
        }
        M[i][i] += b[i];
        B[i] = d[i];
    }
    for (size_t i=0; i<N; i++)  M[i][i+1] += c[i];
    for (size_t i=1; i<=N; i++) M[i][i-1] += a[i];

    IPrinter::printMatrix(12, 4, M);
    puts("---");
    IPrinter::printVector(12, 4, B);

    LinearEquation::GaussianElimination(M, B, v);

    fputs("---\n", file);
    for (size_t i=0; i<=N; i++) { fprintf(file, "%10.6f", x[i]); } fputs("\n", file);
    for (size_t i=0; i<=N; i++) { fprintf(file, "%10.6f", u[i]); } fputs("\n", file);
    for (size_t i=0; i<=N; i++) { fprintf(file, "%10.6f", v[i]); } fputs("\n", file);

    fclose(file);

    delete [] x;
    delete [] d;
    delete [] c;
    delete [] b;
    delete [] a;
}
