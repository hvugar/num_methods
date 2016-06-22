#include "rungekutta.h"
#include <math.h>

void RungeKutta::calculate(R2Function *f, double x0, double x1, double y0, double &y1, double dx)
{
    C_UNUSED(f);
    C_UNUSED(x0);
    C_UNUSED(x1);
    C_UNUSED(y0);
    C_UNUSED(y1);
    C_UNUSED(dx);
}

unsigned int RungeKutta::calculate(R2Function *f, double x0, double x1, double y0, DoubleVector &y, double dx)
{
    if (dx == 0.0) return 0;

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    int n = (int)(ceil(fabs(x1-x0)/fabs(dx)))+1;
    y.clear();
    y.resize(n);

    if (dx > 0.0)
    {
        y[0] = y0;
        for (int i=1; i<n; i++)
        {
            k1 = f->fx(x0,        y0);
            k2 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f->fx(x0+dx,     y0+dx*k3);

            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }

    if (dx < 0.0)
    {
        y[n-1] = y0;
        for (int i=n-2; i>=0; i--)
        {
            k1 = f->fx(x0,        y0);
            k2 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f->fx(x0+dx,     y0+dx*k3);

            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }

    return n;
}

void RungeKutta::calculate(R2Function *f, double x0, double y0, DoubleVector &y, double dx)
{
    if (dx == 0.0) return;

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    unsigned int n = y.size();
    if (dx > 0.0) y[0]   = y0;
    if (dx < 0.0) y[n-1] = y0;

    if (dx > 0.0)
    {
        y[0] = y0;
        for (unsigned int i=1; i<n; i++)
        {
            k1 = f->fx(x0,        y0);
            k2 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f->fx(x0+dx,     y0+dx*k3);
            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }

    if (dx < 0.0)
    {
        y[n-1] = y0;
        for (int i=n-2; i>=0; i--)
        {
            k1 = f->fx(x0,        y0);
            k2 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f->fx(x0+dx,     y0+dx*k3);
            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }
}

void RungeKutta::calculate(R2FunctionX f, double x0, double y0, DoubleVector &y, double dx)
{
    if (dx == 0.0) return;

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    unsigned int n = y.size();
    if (dx > 0.0) y[0]   = y0;
    if (dx < 0.0) y[n-1] = y0;

    if (dx > 0.0)
    {
        y[0] = y0;
        for (unsigned int i=1; i<n; i++)
        {
            k1 = f(x0,        y0);
            k2 = f(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f(x0+dx,     y0+dx*k3);
            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }

    if (dx < 0.0)
    {
        y[n-1] = y0;
        for (int i=n-2; i>=0; i--)
        {
            k1 = f(x0,        y0);
            k2 = f(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f(x0+dx,     y0+dx*k3);
            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }
}

void CauchyProblem::initial(double x0, double y0)
{
    this->mx0 = x0;
    this->my0 = y0;
}

void CauchyProblemSystem(std::vector<RnFunction*> fs, double x0, const DoubleVector &y0, DoubleMatrix &y, double x, double h, unsigned int N)
{
    unsigned int n = fs.size();

    DoubleVector k1(n);
    DoubleVector k2(n);
    DoubleVector k3(n);
    DoubleVector k4(n);

    for (unsigned int j=0; j<y.size(); j++) y[j].clear();

    y.resize(n);
    for (unsigned int j=0; j<n; j++)
    {
        y[j].resize(N+1);
        y[j][0] = y0[j];
    }

    for (unsigned int j=0; j<N; j++)
    {
        DoubleVector arg(n+1);

        // Calculating k1 vector
        arg[0] = x0;
        for (unsigned int i = 0; i<n; i++) arg[i+1] = y[j][i];
        for (unsigned int i = 0; i<n; i++) k1[i] = fs[i]->fx(arg);

        // Calculating k2 vector
        arg[0] = x0+h/2.0;
        for (unsigned int i = 0; i<n; i++) arg[i+1] = y[j][i] + (h/2.0) * k1[i];
        for (unsigned int i = 0; i<n; i++) k2[i] = fs[i]->fx(arg);

        // Calculating k3 vector
        arg[0] = x0+h/2.0;
        for (unsigned int i = 0; i<n; i++) arg[i+1] = y[j][i] + (h/2.0) * k2[i];
        for (unsigned int i = 0; i<n; i++) k2[i] = fs[i]->fx(arg);

        // Calculating k4 vector
        arg[0] = x0+h;
        for (unsigned int i = 0; i<n; i++) arg[i+1] = y[j][i] + h * k3[i];
        for (unsigned int i = 0; i<n; i++) k2[i] = fs[i]->fx(arg);

        // Calculating y
        for (unsigned int i = 0; i<n; i++) y[j+1][i] = y[j][i] + (h/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

        x0 = x0 + h;
    }

    k1.clear();
    k2.clear();
    k3.clear();
    k4.clear();




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


}
