#include "examplem.h"

void ExampleM::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    DoubleMatrix a(2,2);
    a[0][0] = 1.0; a[0][1] = 1.2;
    a[1][0] = 0.8; a[0][1] = 0.3;
    a.inverse();

    ExampleM e;

//    DoubleMatrix rx(3, e.N+1);
//    for (unsigned int k=0; k<=e.N; k++)
//    {
//        rx[0][k] = e.f(0,k);
//        rx[1][k] = e.f(1,k);
//        rx[2][k] = e.f(2,k);
//    }
//    IPrinter::printVector(e.w,e.p,rx.row(0));
//    IPrinter::printVector(e.w,e.p,rx.row(1));
//    IPrinter::printVector(e.w,e.p,rx.row(2));

    DoubleMatrix rx;
    e.calculateRX(rx);

    e.calculate2R2LV1(rx);
//    e.calculate4R2LV1(rx);
//    e.calculate6R2LV1(rx);
}

ExampleM::ExampleM()
{
    h = 0.01;
    N = 100;
    w = 14;
    p = 10;
}

void ExampleM::calculate2R2LV1(const DoubleMatrix &rx)
{
    DoubleMatrix* betta = new DoubleMatrix[N+1];
    for (unsigned int i=0; i<=N; i++) betta[i].resize(3,3,0.0);

    DoubleMatrix alpha(N+1, 3);


    betta[0][0][0] = +1.0; betta[0][0][1] = +1.0; betta[0][0][2] = +1.0;
    betta[0][1][0] = +1.0; betta[0][1][1] = +1.0; betta[0][1][2] = +1.0;
    betta[0][2][0] = +1.0; betta[0][2][1] = +1.0; betta[0][2][2] = +1.0;

    betta[N][0][0] = +1.5; betta[N][0][1] = +1.5; betta[N][0][2] = +1.5;
    betta[N][1][0] = +1.5; betta[N][1][1] = +1.5; betta[N][1][2] = +1.5;
    betta[N][2][0] = +1.5; betta[N][2][1] = +1.5; betta[N][2][2] = +1.5;

    DoubleVector eta = DoubleVector(betta[0]*rx.col(0)) + DoubleVector(betta[N]*rx.col(N));

    for (unsigned int k=0; k<=N-2; k++)
    {
//        double m = +3.0 + 2.0*h*a(k);
//        alpha[k][1] = +4.0/m;
//        alpha[k][2] = -1.0/m;
//        alpha[k][0] = -2.0*h*b(k)/m;
    }

//    IPrinter::printSeperatorLine();

//    for (unsigned int k=0; k<=N-2; k++)
//    {
//        betta[k+1] = betta[k+1] + alpha[k][1]*betta[k];
//        betta[k+2] = betta[k+2] + alpha[k][2]*betta[k];
//        eta        = eta        - alpha[k][0]*betta[k];
//    }

//    DoubleMatrix M(3,3);
//    DoubleVector A(3);
//    DoubleVector x(3);

//    M[0][0] = 0.0;
//    M[0][1] = betta[N-1];
//    M[0][2] = betta[N-0];
//    A[0] = eta;

//    M[1][0] = +1.0;
//    M[1][1] = +2.0*h*a(N-1);
//    M[1][2] = -1.0;
//    A[1] = -2.0*h*b(N-1);

//    M[2][0] = +1.0;
//    M[2][1] = -4.0;
//    M[2][2] = +3.0-2.0*h*a(N);
//    A[2] = +2.0*h*b(N);

//    //printf("det: %.10f\n", M.determinant());
//    //printf("%18.10f %18.10f\n", M[0][0]*rx[N-2] + M[0][1]*rx[N-1] + M[0][2]*rx[N-0], A[0]);
//    //printf("%18.10f %18.10f\n", M[1][0]*rx[N-2] + M[1][1]*rx[N-1] + M[1][2]*rx[N-0], A[1]);
//    //printf("%18.10f %18.10f\n", M[2][0]*rx[N-2] + M[2][1]*rx[N-1] + M[2][2]*rx[N-0], A[2]);

//    //IPrinter::printSeperatorLine();
//    //IPrinter::print(M,3,3);
//    //IPrinter::printSeperatorLine();

//    GaussianElimination(M, A, x);

//    //printf("%.10f %.10f %.10f\n", x[0], x[1], x[2]);

//    //IPrinter::printSeperatorLine();
//    //printf("%18.10f %18.10f\n", M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2], A[0]);
//    //printf("%18.10f %18.10f\n", M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2], A[1]);
//    //printf("%18.10f %18.10f\n", M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2], A[2]);

//    DoubleVector nx(N+1);
//    nx[N-0] = x[2];
//    nx[N-1] = x[1];
//    nx[N-2] = x[0];

//    //IPrinter::printSeperatorLine();
//    //printf("%d %.10f %.10f\n", N-0, rx[N-0], nx[N-0]);
//    //printf("%d %.10f %.10f\n", N-1, rx[N-1], nx[N-1]);
//    //printf("%d %.10f %.10f\n", N-2, rx[N-2], nx[N-2]);

//    betta[N-1] -= alpha[N-2][1]*betta[N-2];
//    betta[N-0] -= alpha[N-2][2]*betta[N-2];
//    eta        += alpha[N-2][0]*betta[N-2];

//    for (unsigned int k=N-2; k>0; k--)
//    {
//        betta[k+0] -= alpha[k-1][1]*betta[k-1];
//        betta[k+1] -= alpha[k-1][2]*betta[k-1];
//        eta        += alpha[k-1][0]*betta[k-1];

//        nx[k-1] = eta;
//        for (unsigned int i=k; i<=N; i++)
//        {
//            nx[k-1] -= betta[i]*nx[i];
//        }
//        nx[k-1] /= betta[k-1];
//    }

//    IPrinter::printVector(w,p,nx);

//    FILE *file =fopen("data_rx.txt", "a");
//    IPrinter::printVector(14,10,nx,"nv2",nx.size(),0,0,file);
//    fclose(file);
}

void ExampleM::calculateRX(DoubleMatrix &rx)
{
    rx.clear();
    rx.resize(3, N+1, 0.0);

    for (unsigned int k=0; k<=N; k++)
    {
        rx[0][k] = f(0,k);
        rx[1][k] = f(1,k);
        rx[2][k] = f(2,k);
    }
    IPrinter::printVector(w,p,rx.row(0));
    IPrinter::printVector(w,p,rx.row(1));
    IPrinter::printVector(w,p,rx.row(2));

//    //    FILE *file1 =fopen("data_rx.txt", "w");
//    //    IPrinter::printVector(14,10,rx,"rx0",rx.size(),0,0,file1);
//    //    fclose(file1);
}

double ExampleM::a(unsigned int i, unsigned int j, unsigned int k) const
{
    if (i==0 && j==0) return +2.0;
    if (i==0 && j==1) return +3.0;
    if (i==0 && j==2) return -1.0;

    if (i==1 && j==0) return +4.0;
    if (i==1 && j==1) return +6.0;
    if (i==1 && j==2) return -2.0;

    if (i==2 && j==0) return -1.0;
    if (i==2 && j==1) return +1.0;
    if (i==2 && j==2) return -1.0;

    return NAN;
}

double ExampleM::b(unsigned int i, unsigned int k) const
{
    double t = k*h;

#ifdef SAMPLE_1
    if (i==0) return -(+2.0*sin(20.0*t*t) + 3.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+40.0*t*cos(20.0*t*t));
    if (i==1) return -(+4.0*sin(20.0*t*t) + 6.0*(cos(10.0*t) - sin(20.0*t)) - 2.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (-10.0*sin(10.0*t) - 20.0*cos(20.0*t));
    if (i==2) return -(-1.0*sin(20.0*t*t) + 1.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+3.0*t*t - 16.0*cos(8.0*t)*sin(8.0*t));
#endif
#ifdef SAMPLE_2
    if (i==0) return -(+2.0*sin(t*t) + 3.0*(cos(t) - sin(t)) - 1.0*(t*t*t - sin(t)*sin(t))) + (+2.0*t*cos(t*t));
    if (i==1) return -(+4.0*sin(t*t) + 6.0*(cos(t) - sin(t)) - 2.0*(t*t*t - sin(t)*sin(t))) + (-sin(t) - cos(t));
    if (i==2) return -(-1.0*sin(t*t) + 1.0*(cos(t) - sin(t)) - 1.0*(t*t*t - sin(t)*sin(t))) + (+3.0*t*t - 2.0*cos(t)*sin(t));
#endif
    return NAN;
}

double ExampleM::f(unsigned int i, unsigned int k) const
{
    double t = k*h;
#ifdef SAMPLE_1
    if (i==0) return sin(20.0*t*t);
    if (i==1) return cos(10.0*t) - sin(20.0*t);
    if (i==2) return t*t*t - sin(8.0*t)*sin(8.0*t);
#endif
#ifdef SAMPLE_2
    if (i==0) return sin(t*t);
    if (i==1) return cos(t) - sin(t);
    if (i==2) return t*t*t - sin(t)*sin(t);
#endif
    return NAN;
}
