#include "headers.h"
#include "ode1storder.h"

//double fx1(double t, double x)
//{
//    return 3.0*x - 3.0*t*t + 2.0*t;
//}

//double a(double t) { return t; }
//double b(double t) { return 3.0*t*t - t*t*t*t; }
//double fx2(double t, double x) { return a(t)*x + b(t); }

double fxt(double t, double x)
{
    return 2.0*x + x*sin(x) + (2.0*t-2.0*t*t-t*t*sin(t*t));
}

double* x0 = 0;
double* dx = 0;
//double a(double t, unsigned int i) { return 6.0*x0[i]*x0[i]; }
//double b(double t, unsigned int i, double h) { return 2.0*x0[i]*x0[i]*x0[i] + 3.0*t*t - 2.0*t*t*t*t*t*t*t*t*t - (x0[i+1] - x0[i-1])/(2.0*h); }
double a(double t, unsigned int i) { return 2.0 + sin(x0[i]) + x0[i]*cos(x0[i]); }
double b(double t, unsigned int i, double h) { return fxt(t, x0[i]) - (x0[i+1] - x0[i-1])/(2.0*h); }
double max(double *x, unsigned int N)
{
    double m = fabs(x[0]);
    for (unsigned int i=0; i<=N; i++) if (m<fabs(x[i])) m=fabs(x[i]);
    return m;
}

void ode_1st_order_border(double t0, double x0, double tn, double xn, unsigned int N, double *x)
{
    double dt = (tn - t0)/N;

    double* da = (double*)malloc(sizeof(double)*(N-1));
    double* db = (double*)malloc(sizeof(double)*(N-1));
    double* dc = (double*)malloc(sizeof(double)*(N-1));
    double* dd = (double*)malloc(sizeof(double)*(N-1));
    double* rx = (double*)malloc(sizeof(double)*(N-1));

    for (unsigned int i=1; i<=N-1; i++)
    {
        double t = i*dt;

        da[i-1] = -1.0;
        db[i-1] = -2.0*dt*a(t, i);
        dc[i-1] = +1.0;
        dd[i-1] = 2.0*dt*b(t, i, dt);
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

class NonLinearODE1stOrder1 : public NonLinearODE1stOrder
{
public:
    virtual double f(double t, double x, unsigned i) const
    {
        return 2.0*x + x*sin(x) + (2.0*t - 2.0*t*t - t*t*sin(t*t));
    }
    virtual double boundary(double t, unsigned i) const
    {
        if (i==0) return x1;
        if (i==N) return x2;
        return 0.0;
    }
    virtual double initial(double t, unsigned i) const { return 0.0; }

    double t1;
    double t2;
    double x1;
    double x2;
    unsigned int N;
};

class LinearODE1stOrder1 : public LinearODE1stOrder
{
public:
    virtual double a(double t, unsigned int i) const { return t; }
    virtual double b(double t, unsigned int i) const { return 2.0*t - t*t*t; }
    virtual double boundary(double t, unsigned int i) const { if (i==0) return x0; if (i==N) return x1; return 0.0; }
    virtual double initial(double t, unsigned int i) const { return 0.0; }

    double t0;
    double x0;
    double t1;
    double x1;
    unsigned int N;
};

double f(double t, double x)
{
    return x*x + t*t;
}

int main(int argc, char ** argv)
{
    SampleBorderHyperBolic::main(argc, argv);
    return 0;
    unsigned int N = 1000;

//    double *x = (double*)malloc(sizeof(double)*(N+1));
//    double *t = (double*)malloc(sizeof(double)*(N+1));
//    for (unsigned int i=0; i<=N; i++)
//    {
//        t[i] = i*0.001;
//        x[i] = t[i]*t[i];
//    }
//    double *d = (double*)malloc(sizeof(double)*(N+1));
//    for (unsigned int i=1; i<=N-1; i++)
//    {
//        printf("%f %f %f %f %f %f\n",
//               t[i],
//               x[i],
//               2.0*x[i],
//               2.0*t[i],
//               (f(x[i+1], t[i])-f(x[i-1], t[i]))/(x[i+1]-x[i-1]),
//               (f(x[i], t[i+1])-f(x[i], t[i-1]))/(t[i+1]-t[i-1]));
//    }

//    return 0;

    //    double t0 = 0.0;
    //    double t1 = 1.0;
    //    double a = 0.0;
    //    double b = 1.0;

//    LinearODE1stOrder1 loo1;
//    loo1.t0 = 0.0;
//    loo1.t1 = 1.0;
//    loo1.x0 = 0.0;
//    loo1.x1 = 1.0;
//    loo1.N = 1000;
//    double *x1 = (double*)malloc(sizeof(double)*(N+1));
//    loo1.solveLinearBoundaryProblem(loo1.t0, loo1.x0, loo1.t1, loo1.x1, loo1.N, x1);
//    for (unsigned int i=0; i<=loo1.N; i++) if (i%100==0) printf("%14.6f", x1[i]);
//    printf("\n");

//    return 0.0;


    NonLinearODE1stOrder1 nla;
    nla.t1 = 0.0;
    nla.x1 = 0.0;
    nla.t2 = 1.0;
    nla.x2 = 1.0;
    nla.N = 1000;
    double *x0 = (double*)malloc(sizeof(double)*(N+1));
    for (unsigned int i=0; i<=nla.N; i++) x0[i] = sin(0.1*i);
    nla.solveNonLinearBoundaryProblem(nla.t1, nla.x1, nla.t2, nla.x2, nla.N, x0);

    //    FILE* file = fopen("data.txt", "w");

    //    x0 = (double*)malloc(sizeof(double)*(N+1));
    //    for (unsigned int i=0; i<=N; i++) { x0[i] = 0.1; }
    //    dx = (double*)malloc(sizeof(double)*(N+1));
    //    for (unsigned int i=0; i<=N; i++) { dx[i] = 0.0; }
    //    do {
    //        double dx0 = a - x0[0];
    //        double dx1 = b - x0[N];
    //        ode_1st_order_border(a, dx0, b, dx1, N, dx);

    //        for (unsigned int i=0; i<=N; i++) fprintf(file, "%14.6f", x0[i]);
    //        fprintf(file, "\n");

    //        for (unsigned int i=0; i<=N; i++) if (i%(N/10)==0) printf("%14.6f", x0[i]);
    //        printf("\n");

    //        for (unsigned int i=0; i<=N; i++) { x0[i] += dx[i]; }
    //    } while ( max(dx,N) > 0.0001);

    //    fclose(file);

    for (unsigned int i=0; i<=N; i++) if (i%(N/10)==0) printf("%14.8f", x0[i]);
    printf("\n");

    return 0;
}
