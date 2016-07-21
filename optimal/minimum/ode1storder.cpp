#include "ode1storder.h"

double max_abs(double *x, unsigned int N)
{
    double m = fabs(x[0]);
    for (unsigned int i=0; i<=N; i++) if (m<fabs(x[i])) m=fabs(x[i]);
    return m;
}

void LinearODE1stOrder::solveLinearBoundaryProblem(double t0, double x0, double tn, double xn, unsigned int N, double *x) const
{
    double h = (tn - t0)/N;

    double* da = (double*)malloc(sizeof(double)*(N-1));
    double* db = (double*)malloc(sizeof(double)*(N-1));
    double* dc = (double*)malloc(sizeof(double)*(N-1));
    double* dd = (double*)malloc(sizeof(double)*(N-1));
    double* rx = (double*)malloc(sizeof(double)*(N-1));

    for (unsigned int i=1; i<=N-1; i++)
    {
        double t = i*h;

        da[i-1] = -1.0;
        db[i-1] = -2.0*h*a(t, i);
        dc[i-1] = +1.0;
        dd[i-1] = 2.0*h*b(t, i);
    }

    da[0] = 0.0;
    dc[N-2] = 0.0;

    dd[0]   += x0;
    dd[N-2] -= xn;

    tomasAlgorithm(da, db, dc, dd, rx, N-1);

    x[0] = x0;
    for (unsigned int i=1; i<=N-1; i++)
    {
        x[i] = rx[i-1];
    }
    x[N] = xn;

    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

void NonLinearODE1stOrder::solveNonLinearBoundaryProblem(double t0, double a, double tn, double b, unsigned int N, double *x0) const
{
    class LinearODE1stOrderExt : public LinearODE1stOrder
    {
    public:

        virtual double a(double t, unsigned int i) const
        {
            double f1 = nl->f(i*h, x[i+1], i);
            double f2 = nl->f(i*h, x[i-1], i);
            //printf("%d %.10f %.10f %.10f\n", i,
            //       (f1 - f2) / (x[i+1]-x[i-1]),
            //       2.0 + sin(x[i]) + x[i]*cos(x[i]), h);
            return (f1 - f2) / (x[i+1]-x[i-1]);
            //return 2.0 + sin(x[i]) + x[i]*cos(x[i]);
        }

        virtual double b(double t, unsigned int i) const
        {
            C_UNUSED(t);
            //printf("%d %.10f %.10f\n", i,
            //       nl->f(i*h, x[i], i) - (x[i+1] - x[i-1])/(2.0*h),
            //       2.0*x[i] + x[i]*sin(x[i]) + 2.0*(i*h) - 2.0*(i*h)*(i*h) - (i*h)*(i*h)*sin((i*h)*(i*h)) - (x[i+1] - x[i-1])/(2.0*h));

            return nl->f(i*h, x[i], i) - (x[i+1] - x[i-1])/(2.0*h);
            //return 2.0*x[i] + x[i]*sin(x[i]) + 2.0*i*h - 2.0*i*h*i*h - i*h*i*h*sin(i*h*i*h) - (x[i+1] - x[i-1])/(2.0*h);
        }

        virtual double initial(double t, unsigned int i) const { return 0.0; }

        virtual double boundary(double t, unsigned int i) const
        {
            if (i == 0) return _a - x[0];
            if (i == N) return _b - x[N];
            return 0.0;
        }

        double *x;
        double h;
        const NonLinearODE1stOrder *nl;
        double _a;
        double _b;
        unsigned int N;
    };

    LinearODE1stOrderExt line;
    line.nl = this;
    line.x = x0;
    line.h = (tn - t0)/N;
    line._a = a;
    line._b = b;
    line.N = N;

    double *dx = (double*)malloc(sizeof(double)*(N+1));
    for (unsigned int i=0; i<=N; i++) { dx[i] = 0.0; }

    do {
        double dx0 = a - x0[0];
        double dx1 = b - x0[N];
        for (unsigned int i=0; i<=N; i++) if (i%(N/10)==0) printf("%14.8f", x0[i]);
        printf("\n");
        line.solveLinearBoundaryProblem(t0, dx0, b, dx1, N, dx);
        //for (unsigned int i=0; i<=N; i++) fprintf(file, "%14.8f", x0[i]);
        //fprintf(file, "\n");
        for (unsigned int i=0; i<=N; i++) if (i%(N/10)==0) printf("%14.8f", x0[i]);
        printf("\n");

        for (unsigned int i=0; i<=N; i++) { x0[i] += dx[i]; }
    } while ( max_abs(dx, N) > 0.0001);
}

//void ODE1stOrder::solveLinearBoundaryProblem(double t0, double x0, double tn, double xn, unsigned int N, double *x) const
//{
//    double h = (tn - t0)/N;

//    double* da = (double*)malloc(sizeof(double)*(N-1));
//    double* db = (double*)malloc(sizeof(double)*(N-1));
//    double* dc = (double*)malloc(sizeof(double)*(N-1));
//    double* dd = (double*)malloc(sizeof(double)*(N-1));
//    double* rx = (double*)malloc(sizeof(double)*(N-1));

//    for (unsigned int i=1; i<=N-1; i++)
//    {
//        double t = i*h;

//        da[i-1] = -1.0;
//        db[i-1] = -2.0*h*a(t);
//        dc[i-1] = +1.0;
//        dd[i-1] = 2.0*h*b(t);
//    }

//    da[0] = 0.0;
//    dc[N-2] = 0.0;

//    dd[0]   += x0;
//    dd[N-2] -= xn;

//    tomasAlgorithm(da, db, dc, dd, rx, N-1);

//    x[0] = x0;
//    for (unsigned int i=1; i<=N-1; i++)
//    {
//        x[i] = rx[i-1];
//    }
//    x[N] = xn;

//    free(rx);
//    free(dd);
//    free(dc);
//    free(db);
//    free(da);
//}

//void ODE1stOrder::solveNonLinearBoundaryProblem(double t0, double a, double tn, double b, unsigned int N, double *x) const
//{
//    double *x0 = (double*)malloc(sizeof(double)*(N+1));
//    for (unsigned int i=0; i<=N; i++) { x0[i] = 0.001; }
//    double *dx = (double*)malloc(sizeof(double)*(N+1));
//    for (unsigned int i=0; i<=N; i++) { dx[i] = 0.001; }

//    do {
//        double dx0 = a - x0[0];
//        double dx1 = b - x0[N];
//        solveLinearBoundaryProblem(t0, dx0, b, dx1, N, dx);

//        //for (unsigned int i=0; i<=N; i++) fprintf(file, "%14.6f", x0[i]);
//        //fprintf(file, "\n");

//        //for (unsigned int i=0; i<=N; i++) if (i%(N/10)==0) printf("%14.6f", x0[i]);
//        //printf("\n");

//        for (unsigned int i=0; i<=N; i++) { x0[i] += dx[i]; }
//    } while ( max_abs(dx, N) > 0.0001);
//}
