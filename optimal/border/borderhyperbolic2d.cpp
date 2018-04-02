#include "borderhyperbolic2d.h"

#define SAMPLE1

void BorderHyperbolic2D::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    BorderHyperbolic2D bp;
    bp.h1 = 0.01;
    bp.h2 = 0.01;
    bp.ht = 0.005;
    bp.N1 = 100;
    bp.N2 = 100;
    bp.M  = 200;
    bp.a1 = 1.0;
    bp.a2 = 1.0;
    bp.lambda = 0.1;

    //    DoubleMatrix m1;
    //    bp.calculate(m1, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, bp.a1, bp.a2);
    //    IPrinter::printMatrix(14,10,m1);
    //    IPrinter::printSeperatorLine();

    //    DoubleMatrix m2;
    //    bp.calculateMVD(m2, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, bp.a1, bp.a2);
    //    IPrinter::printMatrix(14,10,m2);
    //    IPrinter::printSeperatorLine();

    DoubleMatrix m3;
    bp.ht = 0.005;
    bp.M  = 200;
    //bp.calculateMVD2(m3, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, bp.a1, bp.a2);
    bp.calculateMVD3(m3, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, bp.a1, bp.a2, bp.lambda);
    IPrinter::printMatrix(14,10,m3);
    IPrinter::printSeperatorLine();
}

double BorderHyperbolic2D::u(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    double x1 UNUSED_PARAM = i*h1;
    double x2 UNUSED_PARAM = j*h2;
    double t  UNUSED_PARAM = k*ht;
#ifdef SAMPLE1
    return x1*x1 + x2*x2 + t*t;
#endif
#ifdef SAMPLE2
    return x1*x1 + x1*x2 + sin(x1) + cos(x2) + exp(t);
#endif
#ifdef SAMPLE3
    return sin(x1)*sin(x1) + cos(x2) + x1*x2 + t*t*t;
#endif
#ifdef SAMPLE4
    return sin(t) + exp(x1) + x2*x2;
#endif
#ifdef SAMPLE5
    //return x1*x1 + x1*x2 + 5.0*sin(10.0*x1) + 3.0*cos(5.0*x2) + 4.0*cos(5.0*x1*t);
    return x1*x1 + x1*x2 + k1*sin(k2*x1) + k3*cos(k4*x2) + k5*cos(k6*x1*t);
#endif
#ifdef SAMPLE6
    return 2.0*sin(3.0*x1) + 4.0*cos(2.0*x2) + sin(2.0*t);
#endif
}

double BorderHyperbolic2D::initial1(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    return u(i, j, 0);
}

double BorderHyperbolic2D::initial2(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    C_UNUSED(i);
    C_UNUSED(j);
#ifdef SAMPLE1
    return 0.0;
#endif
#ifdef SAMPLE2
    return 1.0;
#endif
#ifdef SAMPLE3
    return 0.0;
#endif
#ifdef SAMPLE4
    return 1.0;
#endif
#ifdef SAMPLE5
    return 0.0;
#endif
#ifdef SAMPLE6
    return 2.0;
#endif
}

double BorderHyperbolic2D::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    return u(i, j, k);
}

double BorderHyperbolic2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 UNUSED_PARAM = i*h1;
    double x2 UNUSED_PARAM = j*h2;
    double t  UNUSED_PARAM = k*ht;

#ifdef SAMPLE1
    return 2.0 - 2.0*(a1*a1) - 2.0*(a2*a2);// + 2.0*t*lambda;
#endif
#ifdef SAMPLE2
    return exp(t) + (a1*a1)*(sin(x1) - 2.0) + (a2*a2)*cos(x2);// + qamma*exp(t);
#endif
#ifdef SAMPLE3
    return 6*t - (a1*a1)*(2.0*cos(2.0*x1)) + (a2*a2)*cos(x2) ;
#endif
#ifdef SAMPLE4
    return -sin(t) - (a1*a1)*exp(x1) - (a2*a2)*2.0;
#endif
#ifdef SAMPLE5
    //return -100.0*x1*x1*cos(5.0*x1*t) + a1*a1*(500.0*sin(10.0*x1) + 100.0*t*t*cos(5.0*x1*t) - 2.0) + a2*a2*(75.0*cos(5.0*x2));
    return -k5*k6*k6*x1*x1*cos(k6*x1*t)
            - a1*a1*(2.0-k1*k2*k2*sin(k2*x1) - k5*k6*k6*t*t*cos(k6*x1*t))
            - a2*a2*(-k3*k4*k4*cos(k4*x2));
#endif
#ifdef SAMPLE6
    return -4.0*sin(2.0*t) + a1*a1*(18.0*sin(3.0*x1)) + a2*a2*(16.0*cos(2.0*x2)) + qamma*2.0*cos(2.0*t);
#endif
}

void BorderHyperbolic2D::calculate(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    u.resize(N2+1, N1+1);

    DoubleMatrix u0(N2+1, N1+1);
    DoubleMatrix u1(N2+1, N1+1);

    DoubleVector da(N1-1);
    DoubleVector db(N1-1);
    DoubleVector dc(N1-1);
    DoubleVector dd(N1-1);
    DoubleVector rx(N1-1);

    double x1_a = -(a1*a1*ht*ht)/(h1*h1);
    double x1_b  = 1.0 + (2.0*a1*a1*ht*ht)/(h1*h1);
    double x1_c = (a2*a2*ht*ht)/(h2*h2);

    double x2_a = -(a2*a2*ht*ht)/(h2*h2);
    double x2_b  = 1.0 + (2.0*a2*a2*ht*ht)/(h2*h2);
    double x2_c = (a1*a1*ht*ht)/(h1*h1);

    for (unsigned int k=1; k<=M; k++)
    {
        if (k==1)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = initial1(i, j);
                    u1[j][i] = u0[j][i] + ht*initial2(i, j);
                }
            }

            //            IPrinter::printMatrix(14,10,u0);
            //            IPrinter::printSeperatorLine();
            //            IPrinter::printMatrix(14,10,u1);
            //            IPrinter::printSeperatorLine();
        }
        else
        {
            if ((k % 2) == 0)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[0][i]  = boundary(i, 0, k);
                    u[N2][i] = boundary(i, N2, k);
                }

                // Approximation to x1 direction
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da[i-1] = x1_a;
                        db[i-1] = x1_b;
                        dc[i-1] = x1_a;
                        dd[i-1] = x1_c*(u1[j-1][i] - 2.0*u1[j][i] + u1[j+1][i]) + 2.0*u1[j][i] - u0[j][i] + (ht*ht) * f(i, j, k);
                    }

                    da[0]     = 0.0;
                    dc[N1-2]  = 0.0;

                    u[j][0]  = boundary(0, j, k);
                    u[j][N1] = boundary(N1, j, k);

                    dd[0]    -= x1_a * u[j][0];
                    dd[N1-2] -= x1_a * u[j][N1];

                    tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.length());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        u[j][i] = rx[i-1];
                    }
                }

            }
            else
            {
                for (unsigned int j=0; j<=N2; j++)
                {
                    u[j][0]  = boundary(0, j, k);
                    u[j][N1] = boundary(N1, j, k);
                }

                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da[j-1] = x2_a;
                        db[j-1] = x2_b;
                        dc[j-1] = x2_a;
                        dd[j-1] = x2_c*(u1[j][i-1] - 2.0*u1[j][i] + u1[j][i+1]) + 2.0*u1[j][i] - u0[j][i] + (ht*ht) * f(i, j, k);
                    }
                    da[0]     = 0.0;
                    dc[N2-2]  = 0.0;

                    u[0][i]  = boundary(i, 0, k);
                    u[N2][i] = boundary(i, N2, k);

                    dd[0]    -= x2_a * u[0][i];
                    dd[N2-2] -= x2_a * u[N2][i];

                    tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.length());

                    for (unsigned int j=1; j<N2; j++)
                    {
                        u[j][i] = rx[j-1];
                    }
                }
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = u1[j][i];
                    u1[j][i] = u[j][i];
                }
            }

            //            IPrinter::printMatrix(14,10,u1);
            //            IPrinter::printSeperatorLine();
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();

    //    da2.clear();
    //    db2.clear();
    //    dc2.clear();
    //    dd2.clear();
    //    rx2.clear();
}


void BorderHyperbolic2DN::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    BorderHyperbolic2DN bpn;
    bpn.hx = 0.01;
    bpn.hy = 0.01;
    bpn.ht = 0.005;
    bpn.Nx = 100;
    bpn.Ny = 100;
    bpn.M  = 200;
    bpn.a1 = 1.0;
    bpn.a2 = 1.0;
    bpn.lambda0 = 0.00;
    bpn.lambda1 = 0.01;

    DoubleMatrix m3;
    bpn.ht = 0.001;
    bpn.M  = 1000;
    //bp.calculateMVD2(m3, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, bp.a1, bp.a2);
    bpn.calculateMVD3(m3, bpn.hx, bpn.hy, bpn.ht, bpn.Nx, bpn.Ny, bpn.M, bpn.a1, bpn.a2);
    IPrinter::printMatrix(14,10,m3);
    IPrinter::printSeperatorLine();
}

void BorderHyperbolic2DN::calculateMVD3(DoubleMatrix &u, double hx, double hy, double ht, unsigned int Nx, unsigned Ny, unsigned int M, double a1, double a2) const
{
    u.resize(Ny+1, Nx+1);

    DoubleMatrix u0(Ny+1, Nx+1);
    DoubleMatrix u1(Ny+1, Nx+1);
    DoubleMatrix u05(Ny+1, Nx+1);
    DoubleMatrix u15(Ny+1, Nx+1);

    DoubleVector da1(Nx+1);
    DoubleVector db1(Nx+1);
    DoubleVector dc1(Nx+1);
    DoubleVector dd1(Nx+1);
    DoubleVector rx1(Nx+1);

    DoubleVector da2(Ny+1);
    DoubleVector db2(Ny+1);
    DoubleVector dc2(Ny+1);
    DoubleVector dd2(Ny+1);
    DoubleVector rx2(Ny+1);

    double x1_a = -0.5*(a1*a1*ht*ht)/(hx*hx);
    double x1_b = +1.0 + (a1*a1*ht*ht)/(hx*hx) + 1.5*lambda0*ht;
    double x1_c = +0.5*(a2*a2*ht*ht)/(hy*hy);
    double x1_d = +(a1*a1*ht*ht*lambda1)/(hx);

    double x2_a = -0.5*(a2*a2*ht*ht)/(hy*hy);
    double x2_b = +1.0 + (a2*a2*ht*ht)/(hy*hy) + 1.5*lambda0*ht;
    double x2_c = +0.5*(a1*a1*ht*ht)/(hx*hx);
    double x2_d = +(a2*a2*ht*ht*lambda1)/(hy);

    double hh = 0.5*ht;

    // initial conditions
    for (unsigned int j=0; j<=Ny; j++)
    {
        for (unsigned int i=0; i<=Nx; i++)
        {
            u0[j][i] = initial1(i, j);
        }
    }

    for (unsigned int j=0; j<=Ny; j++)
    {
        for (unsigned int i=0; i<=Nx; i++)
        {
            u1[j][i] = u0[j][i] + ht*initial2(i, j);
            u15[j][i] = u0[j][i] + hh*initial2(i, j);

            double sum = f(i,j,0);

            if (i==0)       sum += a1*a1*(u0[j][i]-2.0*u0[j][i+1]+u0[j][i+2])/(hx*hx);
            else if (i==Nx) sum += a1*a1*(u0[j][i-2]-2.0*u0[j][i-1]+u0[j][i])/(hx*hx);
            else            sum += a1*a1*(u0[j][i-1]-2.0*u0[j][i]+u0[j][i+1])/(hx*hx);

            if (j==0)       sum += a2*a2*(u0[j][i]-2.0*u0[j+1][i]+u0[j+2][i])/(hy*hy);
            else if (j==Ny) sum += a2*a2*(u0[j-2][i]-2.0*u0[j-1][i]+u0[j][i])/(hy*hy);
            else            sum += a2*a2*(u0[j-1][i]-2.0*u0[j][i]+u0[j+1][i])/(hy*hy);

            sum -= lambda0*initial2(i, j);

            u1[j][i] += 0.5*ht*ht*sum;

            u15[j][i] += 0.5*hh*hh*sum;
        }
    }

//    IPrinter::printMatrix(14,10,u0);
//    IPrinter::printSeperatorLine();
//    IPrinter::printMatrix(14,10,u1);
//    IPrinter::printSeperatorLine();
//    IPrinter::printMatrix(14,10,u15);
//    IPrinter::printSeperatorLine();
//    return;

    for (unsigned int k=2; k<=M; k++)
    {
        // Approximation to x direction

        double t = k*ht-0.5*ht;
        //double t2 = t*t;

        for (unsigned int j=0; j<=Ny; j++)
        {
            for (unsigned int i=0; i<=Nx; i++)
            {
                if (j==0)        dd1[i] = x1_c*(u1[0][i]    - 2.0*u1[1][i]    + u1[2][i]);
                if (j>0 && j<Ny) dd1[i] = x1_c*(u1[j-1][i]  - 2.0*u1[j][i]    + u1[j+1][i]);
                if (j==Ny)       dd1[i] = x1_c*(u1[Ny-2][i] - 2.0*u1[Ny-1][i] + u1[Ny][i]);

                if (i==0)
                {
                    da1[0] = 0.0;
                    db1[0] = x1_b - x1_d;
                    dc1[0] = 2.0*x1_a;
                    dd1[0] += 0.5*(u1[j][i] - u0[j][i]) + u1[j][i]
                            + 0.5*lambda0*ht*(4.0*u1[j][i]-u15[j][i])
                            + 0.5*ht*ht*(f(i, j, 0) + 2.0*t*lambda0)
                            - x1_d*mu(k, t)
                            + (a1*a1*ht*ht*h(0, j, k, t, 1))/hx;
                }
                else if (i==Nx)
                {
                    da1[Nx] = 2.0*x1_a;
                    db1[Nx] = x1_b - x1_d;
                    dc1[Nx] = 0.0;
                    dd1[Nx] += 0.5*(u1[j][i] - u0[j][i]) + u1[j][i]
                            + 0.5*lambda0*ht*(4.0*u1[j][i]-u15[j][i])
                            + 0.5*ht*ht*(f(i, j, 0) + 2.0*t*lambda0)
                            - x1_d*mu(k, t)
                            + (a1*a1*ht*ht*h(Nx, j, k, t, 2))/hx;
                }
                else
                {
                    da1[i] = x1_a;
                    db1[i] = x1_b;
                    dc1[i] = x1_a;
                    dd1[i] += 0.5*(u1[j][i] - u0[j][i]) + u1[j][i]
                            + 0.5*lambda0*ht*(4.0*u1[j][i]-u15[j][i])
                            + 0.5*ht*ht*(f(i, j, 0) + 2.0*t*lambda0);
                }
            }

            tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.length());

            for (unsigned int i=0; i<=Nx; i++)
            {
                u05[j][i] = rx1[i];
            }
        }

        //IPrinter::printSeperatorLine();
        //IPrinter::printMatrix(14,10,u05);
        //IPrinter::printSeperatorLine();
        //return;

        // Approximation to y direction

        t = k*ht;
        //t2 = t*t;

        for (unsigned int i=0; i<=Nx; i++)
        {
            for (unsigned int j=0; j<=Ny; j++)
            {
                if (i==0)        dd2[j] = x2_c*(u05[j][0]    - 2.0*u05[j][1]    + u05[j][2]);
                if (i>0 && i<Nx) dd2[j] = x2_c*(u05[j][i-1]  - 2.0*u05[j][i]    + u05[j][i+1]);
                if (i==Nx)       dd2[j] = x2_c*(u05[j][Nx-2] - 2.0*u05[j][Nx-1] + u05[j][Nx]);

                if (j==0)
                {
                    da2[0] = 0.0;
                    db2[0] = x2_b - x2_d;
                    dc2[0] = 2.0*x2_a;
                    dd2[0] += 0.5*(u1[j][i] - u0[j][i]) + u05[j][i]
                            + 0.5*lambda0*ht*(4.0*u05[j][i]-u1[j][i])
                            + 0.5*ht*ht*(f(i, j, 0) + 2.0*t*lambda0)
                            - x2_d*mu(k, t)
                            + (a2*a2*ht*ht*h(i, 0, k, t, 3))/hy;
                }
                else if (j==Ny)
                {
                    da2[Ny] = 2.0*x2_a;
                    db2[Ny] = x2_b - x2_d;
                    dc2[Ny] = 0.0;
                    dd2[Ny] += 0.5*(u1[j][i] - u0[j][i]) + u05[j][i]
                            + 0.5*lambda0*ht*(4.0*u05[j][i]-u1[j][i])
                            + 0.5*ht*ht*(f(i, j, 0) + 2.0*t*lambda0)
                            - x2_d*mu(k, t)
                            + (a2*a2*ht*ht*h(i, Ny, k, t, 4))/hy;
                }
                else
                {
                    da2[j] = x2_a;
                    db2[j] = x2_b;
                    dc2[j] = x2_a;
                    dd2[j] += 0.5*(u1[j][i] - u0[j][i]) + u05[j][i]
                            + 0.5*lambda0*ht*(4.0*u05[j][i]-u1[j][i])
                            + 0.5*ht*ht*(f(i, j, 0) + 2.0*t*lambda0);
                }
            }

            tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.length());

            for (unsigned int j=0; j<=Ny; j++)
            {
                u[j][i] = rx2[j];
            }
        }

        //IPrinter::printSeperatorLine();
        //IPrinter::printMatrix(14,10,u);
        //IPrinter::printSeperatorLine();
        //break;

        for (unsigned int j=0; j<=Ny; j++)
        {
            for (unsigned int i=0; i<=Nx; i++)
            {
                u0[j][i] = u1[j][i];
                u1[j][i] = u[j][i];
                u15[j][i] = u05[j][i];
            }
        }
    }

    da1.clear();
    db1.clear();
    dc1.clear();
    dd1.clear();
    rx1.clear();

    da2.clear();
    db2.clear();
    dc2.clear();
    dd2.clear();
    rx2.clear();
}

double BorderHyperbolic2DN::u(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    double x UNUSED_PARAM = i*hx;
    double y UNUSED_PARAM = j*hy;
    double t  UNUSED_PARAM = k*ht;
    return x*x + y*y + t*t;
}

double BorderHyperbolic2DN::initial1(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    return u(i, j, 0);
}

double BorderHyperbolic2DN::initial2(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    return 0.0;
}

double BorderHyperbolic2DN::boundary(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    return NAN;
}

double BorderHyperbolic2DN::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    //double x1 UNUSED_PARAM = i*h1;
    //double x2 UNUSED_PARAM = j*h2;
    //double t  UNUSED_PARAM = k*ht;
    return 2.0 - 2.0*(a1*a1) - 2.0*(a2*a2);// + 2.0*t*lambda;
}

double BorderHyperbolic2DN::mu(unsigned int k, double t) const
{
    return 5.0;
}

double BorderHyperbolic2DN::h(unsigned int i, unsigned int j, unsigned int k, double t, unsigned int n) const
{
    double x = i*hx;
    double y = j*hy;

    if (n==1) return 0.0 - lambda1*(x*x+y*y+t*t-mu(k,t));
    if (n==2) return 2.0 - lambda1*(x*x+y*y+t*t-mu(k,t));
    if (n==3) return 0.0 - lambda1*(x*x+y*y+t*t-mu(k,t));;
    if (n==4) return 2.0 - lambda1*(x*x+y*y+t*t-mu(k,t));;
    return NAN;
}


