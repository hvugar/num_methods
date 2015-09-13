#include "gridmethod.h"
#include <math.h>

GridMethod::GridMethod()
{
    this->alpha = 1.0;
}

void GridMethod::setLengthInterval(double x0, double x1)
{
    this->x0 = x0;
    this->x1 = x1;
}

void GridMethod::setTimeInterval(double t0, double t1)
{
    this->t0 = t0;
    this->t1 = t1;
}

void GridMethod::setLengthTimeStep(double dx, double dt)
{
    this->dx = dx;
    this->dt = dt;
    this->n = (unsigned int)(ceil(x1-x0)/dx) + 1;
    this->m = (unsigned int)(ceil(t1-t0)/dt) + 1;
}

void GridMethod::setLengthTimeStepCount(unsigned int n, unsigned int m)
{
    this->n = n;
    this->m = m;
    this->dx = (x1 - x0) / (n - 1);
    this->dt = (t1 - t0) / (m - 1);
}

void GridMethod::setF(R2Function* f)
{
    this->f = f;
}

void GridMethod::setM1(R1Function* m1)
{
    this->m1 = m1;
}

void GridMethod::setM2(R1Function* m2)
{
    this->m2 = m2;
}

void GridMethod::setFi(R1Function* fi)
{
    this->fi = fi;
}

void GridMethod::implicitDifferenceScheme()
{
    unsigned int k = n - 2;

    // initilizing

    for (unsigned int i=0; i<m; i++) u.push_back(DoubleVector(m));


    for (unsigned int j=0; j<m; j++)
    {
        u[j][0]   = m1->fx(dt * j);
        u[j][n-1] = m2->fx(dt * j);
    }

    for (unsigned int i=0; i<n; i++)
    {
        u[0][i] = fi->fx( dx * i );
    }

    /////////////////////////////////////////////////

    DoubleVector b(k);
    std::vector<DoubleVector> a;
    for (unsigned int i=0; i<k; i++)
    {
        a.push_back(DoubleVector(k));
    }

    ////////////////////////////////////////////////
    double c1 = -alpha * (dt / (dx*dx));
    double c2 = 1.0 + (2.0*alpha) * (dt / (dx*dx));

    for (unsigned int j=1; j<m; j++)
    {
        for (unsigned int i=1; i<=k; i++)
        {
            if ( i == 1 )
            {
                a[i-1][0] = c2;
                a[i-1][1] = c1;
                b[i-1] = -c1*u[j][0] + u[j-1][1] + dt * f->fx(i*dx, j*dt);
            }
            else if ( i == k )
            {
                a[i-1][i-2] = c1;
                a[i-1][i-1] = c2;
                b[i-1] = -c1*u[j][k+1] + u[j-1][k] + dt * f->fx(i*dx, j*dt);
            }
            else
            {
                a[i-1][i+0] = c1;
                a[i-1][i-1] = c2;
                a[i-1][i-2] = c1;
                b[i-1] = u[j-1][i] + dt * f->fx(i*dx, j*dt);
            }
        }

        DoubleVector x(k);
        tomas_algorithm(a, b, x);
        for (unsigned int i=0; i<k; i++) u[j][i+1] = x[i];
    }

    for (unsigned int i=0; i<k; i++)
    {
        a[i].clear();
    }
    a.clear();
    b.clear();

    for (unsigned int i=0; i<m; i++)
    {
        if (i%100==0 || i==999) {
            for (unsigned int j=0; j<n; j++)
            {
                if (j%100==0 || j==999) printf("%f ", u[i][j]);
            }
            puts("");
        }
    }
    printf("%d %d\n", m, n);
}

void GridMethod::tomas_algorithm(std::vector<DoubleVector> &a, const DoubleVector& b, DoubleVector& x)
{
    unsigned int size = x.size();

    DoubleVector p(size);
    DoubleVector q(size);

    for (unsigned int i=0; i<size; i++)
    {
        if (i==0)
        {
            p[i] = -a[i][i+1] / a[i][i];
            q[i] = b[i] / a[i][i];
        }
        else if (i==size-1)
        {
            p[i] = 0.0;
            q[i] = (b[i] - a[i][i-1]*q[i-1])/(a[i][i]+a[i][i-1]*p[i-1]);
        }
        else
        {
            p[i] = -(a[i][i+1])/(a[i][i]+a[i][i-1]*p[i-1]);
            q[i] = (b[i] - a[i][i-1]*q[i-1])/(a[i][i]+a[i][i-1]*p[i-1]);
        }
    }

    for (unsigned int i=size-1; i>=0; i--)
    {
        if (i==(size-1))
            x[i] = q[i];
        else
            x[i] = p[i]*x[i+1] + q[i];
    }

    p.clear();
    q.clear();
}

void printLayer(const std::vector<DoubleVector>& u)
{
    unsigned int m = u.size()/10;
    for (unsigned int j=0; j<u.size(); j++)
    {
        if (j%m==0)
        {
            for (unsigned int i=0; i<u[j].size(); i++)
            {
                if (i%m==0) printf("%14.10f ", u[j][i]);
            }
            puts("");
        }
    }
}

double _u(double x, double y, double t) { return x*x*x+y*y*y*y+t*t; }
double _f(double x, double y, double t) { return 2.0*t-6.0*x-12.0*y*y; }
double _fi(double x, double y) { return _u(x, y ,0.0); }
double _m1(double y, double t) { return _u(0.0, y, t); }
double _m2(double y, double t) { return _u(1.0, y, t); }
double _m3(double x, double t) { return _u(x, 0.0, t); }
double _m4(double x, double t) { return _u(x, 1.0, t); }

void GridMethod::printResult(unsigned int k, unsigned int N, unsigned int M)
{
    double l1 = 1.0;
    double l2 = 1.0;

    double dx = l1 / N;
    double dy = l2 / M;
    double dt = 0.001;

    DoubleMatrix u;

    u.resize(M+1);
    for (unsigned int j=0; j<M+1; j++)
    {
        u[j].resize(N+1);
    }

    for (unsigned int j=0; j<M+1; j++)
    {
        for (unsigned int i=0; i<N+1; i++)
        {
            u[j][i] = _u((dx*i), (dy*j), dt*(k+0.5));
        }
    }
    printf("k=%d\n", 2*k+1);
    printLayer(u);

    for (unsigned int j=0; j<M+1; j++)
    {
        for (unsigned int i=0; i<N+1; i++)
        {
            u[j][i] = _u((dx*i), (dy*j), dt*(k+1.0));
        }
    }

    printf("k=%d\n", 2*k+2);
    printLayer(u);
    puts("--------------------------------------------------");
}

void GridMethod::VariableDirectionsMethod(R2Function *fi, R2Function *m1, R2Function *m2, R2Function *m3, R2Function *m4, R3Function *f)
{
    unsigned int N = 1000;
    unsigned int M = 1000;
    unsigned int K = 4;
    double l1 = 1.0;
    double l2 = 1.0;
    double dx = l1 / N;
    double dy = l2 / M;
    double dt = 0.001;
    double a1 = 1.0;
    double a2 = 1.0;

    std::vector<DoubleMatrix> u;
    u.resize(2*K+1);

    for (unsigned int k=0; k<2*K+1; k++)
    {
        u[k].resize(M+1);
        for (unsigned int j=0; j<M+1; j++)
        {
            u[k][j].resize(N+1);
        }
    }

    for (unsigned int j=0; j<M+1; j++)
    {
        for (unsigned int i=0; i<N+1; i++)
        {
            u[0][j][i] = _fi((dx*i), (dy*j));
        }
    }
    printf("\nk = %d\n", 0);
    printLayer(u[0]);

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    double alpha1 = 0.0;
    double beta1  = 0.0;
    double alpha2 = 0.0;
    double beta2  = 0.0;

    for (unsigned int k=0; k<K; k++)
    {
        // Approximation on x direction
        alpha1 = -(a1*dt)/(2.0*dx*dx);
        beta1  = 1.0 + (a1*dt)/(dx*dx);
        alpha2 = (a2*dt)/(2.0*dy*dy);
        beta2  = 1.0 - (a2*dt)/(dy*dy);
        a.resize(N-1);
        b.resize(N-1);
        c.resize(N-1);
        d.resize(N-1);
        x.resize(N-1);
        for (unsigned int j=1; j<M; j++)
        {
            for (unsigned int i=1; i<N; i++)
            {
                a[i-1] = alpha1;
                b[i-1] = beta1;
                c[i-1] = alpha1;
                d[i-1] = alpha2*u[2*k][j-1][i] + beta2*u[2*k][j][i] + alpha2*u[2*k][j+1][i] + (dt/2.0)*_f(i*dx, j*dy, (k+0.5)*dt);
            }
            a[0]   = 0.0;
            c[N-2] = 0.0;
            d[0]   -= alpha1 * _m1(dy*j, dt*(k+0.5));
            d[N-2] -= alpha1 * _m2(dy*j, dt*(k+0.5));
            TomasAlgorithm(a, b, c, d, x);
            for (unsigned int i=1; i<N; i++)
            {
                u[2*k+1][j][i] = x[i-1];
            }
            u[2*k+1][j][0] = _m1(dy*j, dt*(k+0.5));
            u[2*k+1][j][N] = _m2(dy*j, dt*(k+0.5));
        }

        for (unsigned int i=0; i<=N; i++)
        {
            u[2*k+1][0][i] = _m3(dx*i, dt*(k+0.5));
            u[2*k+1][M][i] = _m4(dx*i, dt*(k+0.5));
        }
        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        printf("\nk = %d\n", 2*k+1);
        printLayer(u[2*k+1]);

        ///////////////////////////////////////////////////////////////////////
        // Approximation on y direktion
        alpha1 = -(a2*dt)/(2.0*dy*dy);
        beta1  = 1.0 + (a2*dt)/(dy*dy);
        alpha2 = (a1*dt)/(2.0*dx*dx);
        beta2  = 1.0 - (a1*dt)/(dx*dx);

        a.resize(M-1);
        b.resize(M-1);
        c.resize(M-1);
        d.resize(M-1);
        x.resize(M-1);
        for (unsigned int i=1; i<N; i++)
        {
            for (unsigned int j=1; j<M; j++)
            {
                a[j-1] = alpha1;
                b[j-1] = beta1;
                c[j-1] = alpha1;
                d[j-1] = alpha2*u[2*k+1][j][i-1] + beta2*u[2*k+1][j][i] + alpha2*u[2*k+1][j][i+1] + (dt/2.0)*_f(i*dx, j*dy, (k+0.5)*dt);
            }
            a[0]   = 0.0;
            c[M-2] = 0.0;
            d[0]   -= alpha1 * _m3(dx*i, dt*(k+1.0));
            d[M-2] -= alpha1 * _m4(dx*i, dt*(k+1.0));
            TomasAlgorithm(a, b, c, d, x);
            for (unsigned int j=1; j<M; j++)
            {
                u[2*k+2][j][i] = x[j-1];
            }
            u[2*k+2][0][i] = _m3(dx*i, dt*(k+1));
            u[2*k+2][M][i] = _m4(dx*i, dt*(k+1));
        }
        for (unsigned int j=0; j<=M; j++)
        {
            u[2*k+2][j][0] = _m1(dy*j, dt*(k+1));
            u[2*k+2][j][N] = _m2(dy*j, dt*(k+1));
        }
        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        printf("k = %d\n", 2*k+2);
        printLayer(u[2*k+2]);
    }
    printf("end\n");
}

void GridMethod::TomasAlgorithm(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x)
{
    if (x.size() != a.size() || x.size() != b.size() || x.size() != c.size() || x.size() != d.size())
        return;

    unsigned int n = x.size();
    DoubleVector p(n);
    DoubleVector q(n);

    for (unsigned int i=0; i<n; i++)
    {
        if (i==0)
        {
            p[0] = d[0]/b[0];
            q[0] = -c[0]/b[0];
        } else
            if(i==n-1)
            {
                p[n-1] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[n-1] = 0.0;
            }
            else
            {
                p[i] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[i] = -c[i]/(b[i]+a[i]*q[i-1]);
            }
    }

    for (int i=n-1; i>=0; i--)
    {
        if (i==n-1)
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];
        }
    }
}
