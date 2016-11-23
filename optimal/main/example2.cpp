#include "example2.h"

double my_rand()
{
    return ((rand() % 100) + 1) / 100.0;
}

void Example2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Example2 e2;
//    e2.sample_n3();
    e2.sample_n4();
//    e2.sample_n6();
    e2.sample_n4_test1();
    e2.sample_n4_test2();
}

Example2::Example2() {}

void Example2::calculateLeft2RightSample()
{
    unsigned int N = 1000;
    unsigned int K = 2;
    double h = 0.001;

    FILE *file1 = fopen("data_x.txt", "w");

    DoubleVector x0(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        double t = i*h;
        x0.at(i) = 5.0*sin(t) + 3.0*t;
    }
    IPrinter::printVector(x0,"x0:",x0.size(),0,0,file1);

    DoubleMatrix a(N-K+1, K+1);

    DoubleVector x1(N+1);
    x1.at(N)   = 5.0*sin(N*h) + 3.0*N*h;
    x1.at(N-1) = 5.0*sin((N-1)*h) + 3.0*(N-1)*h;
    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        double t = (i+1)*h;
        double a1 = (2.0-h*h)/(h+1.0);
        double a2 = (h-1.0)/(h+1.0);
        double a0 = ((h*h)/(h+1.0))*(3.0*t-10.0*cos(t)-6.0);

        a.at(i,0) = a0;
        a.at(i,1) = a1;
        a.at(i,2) = a2;

        x1.at(i) = a.at(i,1)*x1.at(i+1) + a.at(i,2)*x1.at(i+2) + a.at(i,0);
    }
    IPrinter::printVector(x1,"x1:",N+1,0,0,file1);

    unsigned int n11 = 5*N/10;
    unsigned int n12 = 8*N/10;
    unsigned int n21 = 2*N/10;
    unsigned int n22 = 4*N/10;

    DoubleMatrix beta(K, N+1);
    beta.at(0,0) = 0.5; beta.at(0,1) = 1.5;
    beta.at(1,0) = +1.0; beta.at(1,n21) = +0.1; beta.at(1,n22) = +1.3; beta.at(1,N-1) = -1.2; beta.at(1,N) = +1.1;

    DoubleVector qamma(K);
    qamma.at(0) = 0.0;
    qamma.at(1) = 0.0;
    for (unsigned int i=0; i<=N; i++)
    {
        qamma.at(0) += beta.at(0,i)*x1.at(i);
        qamma.at(1) += beta.at(1,i)*x1.at(i);
    }
    printf("%.10f %.10f\n", qamma.at(0), qamma.at(1));

    DoubleVector x2(N+1, 0.0);
    calculateLeft2Right(N,K,a,beta,qamma,x2);
    IPrinter::printVector(x2,"x2:",x2.size(),0,0,file1);
    fclose(file1);
}

void Example2::calculateLeft2Right(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x)
{
    for (unsigned int eq=0; eq<K; eq++)
    {
        for (unsigned int i=0; i<=N-K; i++)
        {
            for (unsigned int j=1; j<=K; j++)
            {
                beta.at(eq,i+j) = beta.at(eq,i+j) + beta.at(eq,i)*a.at(i,j);
            }
            qamma.at(eq) = qamma.at(eq) - beta.at(eq,i)*a.at(i,0);
        }
    }

    DoubleMatrix m(K,K);

    for (unsigned int j=0; j<K; j++)
    {
        for (unsigned int i=0; i<K; i++) m.at(j,i) = beta.at(j,N-K+1+i);
    }

//    printf("Determinant: %.14f\n", m.determinant());
//    IPrinter::printMatrix(18, 14, m, m.rows(), m.cols());

    DoubleVector b(K);
    for (unsigned int i=0; i<K; i++) b.at(i) = qamma.at(i);

    DoubleVector x1(K);
    GaussianElimination(m, b, x1);
    for (unsigned int i=0; i<K; i++) x.at(N-i) = x1.at((K-1)-i);
    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x.at(i) = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x.at(i) += a.at(i,j)*x.at(i+j);
    }
}

void Example2::calculateRight2LeftSample()
{
    unsigned int N = 1000;
    unsigned int K = 2;
    double h = 0.001;

    FILE *file1 = fopen("data_x.txt", "w");

    DoubleVector x0(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        double t = i*h;
        x0.at(i) = 5.0*sin(t) + 3.0*t;
    }
    IPrinter::printVector(x0,"x0:",N+1,0,0,file1);

    DoubleMatrix a(N-K+1, K+1);

    DoubleVector x1(N+1);
    x1.at(0) = 5.0*sin(0*h) + 3.0*0*h;
    x1.at(1) = 5.0*sin(1*h) + 3.0*1*h;
    for (unsigned int i=2; i<=N; i++)
    {
        double t = (i-1)*h;
        double a1 = (h*h-2.0)/(h-1.0);
        double a2 = (h+1.0)/(h-1.0);
        double a0 = (-(h*h)/(h-1.0))*(3.0*t-10.0*cos(t)-6.0);

        a.at(i-2,0) = a0;
        a.at(i-2,1) = a1;
        a.at(i-2,2) = a2;

        x1.at(i) = a.at(i-2,1)*x1.at(i-1) + a.at(i-2,2)*x1.at(i-2) + a.at(i-2,0);
    }
    IPrinter::printVector(x1,"x1:",N+1,0,0,file1);

    unsigned int n11 = 5*N/10;
    unsigned int n12 = 8*N/10;
    unsigned int n21 = 2*N/10;
    unsigned int n22 = 4*N/10;

    DoubleMatrix beta(K, N+1);
    beta.at(0,0) = 0.5; beta.at(0,N) = 1.5;
    beta.at(1,0) = +1.0; beta.at(1,n21) = +0.1; beta.at(1,n22) = +1.3; beta.at(1,N-1) = -1.2; beta.at(1,N) = +1.1;

    DoubleVector qamma(K);
    qamma.at(0) = 0.0;
    qamma.at(1) = 0.0;
    for (unsigned int i=0; i<=N; i++)
    {
        qamma.at(0) += beta.at(0,i)*x1.at(i);
        qamma.at(1) += beta.at(1,i)*x1.at(i);
    }
    printf("%.10f %.10f\n", qamma.at(0), qamma.at(1));

    DoubleVector x2(N+1, 0.0);
    calculateRight2Left(N,K,a,beta,qamma,x2, x1);
    IPrinter::printVector(x2,"x2:",x2.size(),0,0,file1);
    fclose(file1);
}

void Example2::calculateRight2Left(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x, const DoubleVector &x0 UNUSED_PARAM)
{
    for (unsigned int eq=0; eq<K; eq++)
    {
        for (unsigned int i=N; i>=K; i--)
        {
            for (unsigned int j=1; j<=K; j++)
            {
                beta.at(eq,i-j) = beta.at(eq,i-j) + beta.at(eq,i)*a.at(i-K,j);
            }
            qamma.at(eq) = qamma.at(eq) - beta.at(eq,i)*a.at(i-K,0);
        }
    }

    DoubleMatrix m(K,K);
    for (unsigned int j=0; j<K; j++)
    {
        for (unsigned int i=0; i<K; i++) m.at(j,i) = beta.at(j,i);
    }
    DoubleVector b(K);
    for (unsigned int i=0; i<K; i++) b.at(i) = qamma.at(i);

    DoubleVector x1(K);
    GaussianElimination(m, b, x1);
    for (unsigned int i=0; i<K; i++) x.at(i) = x1.at(i);
    for (unsigned int i=K; i <= N; i++)
    {
        x[i] = a.at(i-K,0);
        for (unsigned int j=1; j<=K; j++)
        {
            x.at(i) += a.at(i-K,j)*x.at(i-j);
        }
    }
}

double Example2::X(double t UNUSED_PARAM) const
{
#ifdef SAMPLE_1
    return t*t;
#endif
#ifdef SAMPLE_2
    return 3.0*sin(t);
#endif
#ifdef SAMPLE_3
    return (t*t-t)*sin(50.0*t);
#endif
#ifdef SAMPLE_4
    return exp(4.0*t);
#endif
}

double Example2::A(double t UNUSED_PARAM) const
{
#ifdef SAMPLE_1
    return 3.0;
#endif
#ifdef SAMPLE_2
    return t;
#endif
#ifdef SAMPLE_3
    return t;
#endif
#ifdef SAMPLE_4
    return (4.0+t);
#endif
}

double Example2::B(double t UNUSED_PARAM) const
{
#ifdef SAMPLE_1
    return 2.0*t - 3.0*t*t;
#endif
#ifdef SAMPLE_2
    return 3.0*cos(t) - 3.0*t*sin(t);
#endif
#ifdef SAMPLE_3
    return 50.0*(t*t-t)*cos(50.0*t) + (2.0*t-1)*sin(50.0*t) - t*(t*t-t)*sin(50.0*t);
#endif
#ifdef SAMPLE_4
    return -t*exp(4.0*t);
#endif
}

void Example2::sample1()
{
    unsigned int N = 1000;
    unsigned int K = 2;
    double h = 0.001;

    FILE *file1 = fopen("data_sample.txt", "w");
    DoubleVector x0(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        double t = i*h;
        x0.at(i) = X(t);
    }
    IPrinter::printVector(x0,"x0:",N+1,0,0,file1);

    DoubleMatrix a(N-K+1, K+1);

    DoubleVector x1(N+1);
    x1.at(N)   = X(N*h);
    x1.at(N-1) = X((N-1)*h);
    //x1.at(N-2) = -2.0*h*A((N-1)*h)*x1.at(N-1) + x1.at(N) - 2.0*h*B((N-1)*h);
    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        double t = (i+1)*h;
        a.at(i,0) = -2.0*h*B(t);
        a.at(i,1) = -2.0*h*A(t);
        a.at(i,2) = 1.0;

        x1.at(i) = a.at(i,1)*x1.at(i+1) + a.at(i,2)*x1.at(i+2) + a.at(i,0);
    }
    IPrinter::printVector(x1,"x1:",N+1,0,0,file1);

    //    unsigned int n11 = 5*N/10;
    //    unsigned int n12 = 8*N/10;
    unsigned int n21 = 2*N/10;
    unsigned int n22 = 4*N/10;

    DoubleMatrix beta(K, N+1);
    beta.at(0,0) = -3.0-2.0*h*A(0.0); beta.at(0,1) = +4.0; beta.at(0,2) = -1.0;
    beta.at(1,0) = +1.0; beta.at(1,N/2) = -2.1; beta.at(1,N) = +1.1;

    DoubleVector qamma(K);
    qamma.at(0) = 2.0*h*B(0.0);
    qamma.at(1) = 0.0;
    for (unsigned int i=0; i<=N; i++)
    {
        //qamma.at(0) += beta.at(0,i)*x0.at(i);
        qamma.at(1) += beta.at(1,i)*x1.at(i);
    }
    printf("%.10f %.10f\n", qamma.at(0), qamma.at(1));

    DoubleVector x2(N+1, 0.0);
    calculateLeft2Right(N,K,a,beta,qamma,x2);
    IPrinter::printVector(x2,"x2:",x2.size(),0,0,file1);
    fclose(file1);
}

void Example2::sample2()
{
    unsigned int N = 100;
    unsigned int K = 4;
    double h = 1.0/N;

    FILE *file1 = fopen("data_sample_n4.txt", "w");
    DoubleVector x0(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        double t = i*h;
        x0.at(i) = X(t);
    }
    IPrinter::printVector(14, 10, x0, "x0:", x0.size(), 0, 0, file1);

    DoubleMatrix a(N-K+1, K+1);
    DoubleVector x1(N+1);
    x1.at(N-0) = X((N-0)*h);
    x1.at(N-1) = X((N-1)*h);
    x1.at(N-2) = X((N-2)*h);
    x1.at(N-3) = X((N-3)*h);
    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        //        double t  = (i+2)*h;
        //        a.at(i,0) = +12.0*h*B(t);
        //        a.at(i,1) = +8.0;
        //        a.at(i,2) = +12.0*h*A(t);
        //        a.at(i,3) = -8.0;
        //        a.at(i,4) = +1.0;

        double t  = i*h;
        double m = -25.0 - 12.0*h*A(t);
        a.at(i,0) = +(12.0*h*B(t))/m;
        a.at(i,1) = -48.0/m;
        a.at(i,2) = +36.0/m;
        a.at(i,3) = -16.0/m;
        a.at(i,4) = +3.0/m;

        x1.at(i) = a.at(i,1)*x1.at(i+1) + a.at(i,2)*x1.at(i+2) + a.at(i,3)*x1.at(i+3) + a.at(i,4)*x1.at(i+4) + a.at(i,0);
    }
    IPrinter::printVector(18, 14, x1, "x1:", x1.size(), 0, 0, file1);

    FILE *file2 = fopen("alpha.txt", "w");
    for (unsigned int i=0; i<a.rows(); i++)
    {
        fprintf(file2, "%4d %.10f %.10f %.10f %.10f %.10f\n", i, a.at(i,1), a.at(i,2), a.at(i,3), a.at(i,4), a.at(i,0));
    }
    fclose(file2);

    DoubleMatrix beta(K, N+1);

    //    beta.at(0,0) = -25.0-12.0*h*A(0.0);
    //    beta.at(0,1) = +48.0;
    //    beta.at(0,2) = -36.0;
    //    beta.at(0,3) = +16.0;
    //    beta.at(0,4) = -3.0;
    beta.at(0,0)   = +2.5;
    beta.at(0,N)   = +0.8;

    //    beta.at(1,0) = -3.0;
    //    beta.at(1,1) = -10.0-12.0*h*A(h);
    //    beta.at(1,2) = +18.0;
    //    beta.at(1,3) = -6.0;
    //    beta.at(1,4) = +1.0;
    beta.at(1,0)   = +3.5;
    beta.at(1,N)   = +1.1;

    //    beta.at(2,N-4) = -1.0;
    //    beta.at(2,N-3) = +6.0;
    //    beta.at(2,N-2) = -18.0;
    //    beta.at(2,N-1) = +10.0-12.0*h*A((N-1)*h);
    //    beta.at(2,N-0) = +3.0;
    beta.at(2,0)   = +5.5;
    beta.at(2,N)   = +2.1;

    //    beta.at(3,0) = +25.0-12.0*h*A(N*h);
    //    beta.at(3,1) = +3.0;
    //    beta.at(3,2) = -16.0;
    //    beta.at(3,3) = +36.0;
    //    beta.at(3,4) = -48.0;

    beta.at(3,0)   = -0.5;
    beta.at(3,N)   = +1.1;

    DoubleVector qamma(K);
    qamma.at(0) = 0.0;//+12.0*h*B(0.0);
    qamma.at(1) = 0.0;//+12.0*h*B(h);
    qamma.at(2) = 0.0;//+12.0*h*B((N-1)*h);
    qamma.at(3) = 0.0;//beta.at(3,0)*x1.at(0) + beta.at(3,N/2)*x1.at(N/2) + beta.at(3,N)*x1.at(N);
    for (unsigned int i=0; i<=N; i++)
    {
        qamma.at(0) += beta.at(0,i)*x1.at(i);
        qamma.at(1) += beta.at(1,i)*x1.at(i);
        qamma.at(2) += beta.at(2,i)*x1.at(i);
        qamma.at(3) += beta.at(3,i)*x1.at(i);
    }

    DoubleVector x(N+1, 0.0);
    calculateLeft2Right(N,K,a,beta,qamma, x);
    IPrinter::printVector(18, 14, x, "x2:", x.size(), 0, 0, file1);

    fclose(file1);
}

void Example2::sample_n3()
{
    FILE *file1 = fopen("sample_n3_x.txt", "w");

    unsigned int N = 1000;
    unsigned int K = 3;
    double h = 1.0/N;

    DoubleVector x0(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        double t = i*h;
        x0.at(i) = X(t);
    }
    IPrinter::printVector(14, 10, x0, "x0:", x0.size(), 0, 0, file1);

    DoubleMatrix a(N-K+1, K+1);

    DoubleVector x1(N+1);

    x1.at(N-0) = X((N-0)*h);
    x1.at(N-1) = X((N-1)*h);
    x1.at(N-2) = X((N-2)*h);
    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        double t  = i*h;
        double m  = 11.0+6.0*h*A(t);
        a.at(i,0) = -6.0*h*B(t)/m;
        a.at(i,1) = +18.0/m;
        a.at(i,2) = -9.0/m;
        a.at(i,3) = +2.0/m;

        x1.at(i) = a.at(i,1)*x1.at(i+1) + a.at(i,2)*x1.at(i+2) + a.at(i,3)*x1.at(i+3) + a.at(i,0);
    }
    IPrinter::printVector(14, 10, x1, "x1:", x1.size(), 0, 0, file1);

    DoubleMatrix beta(K, N+1, 0.0);

    beta.at(0,N-3) = -2.0;
    beta.at(0,N-2) = -3.0 - 6.0*h*A((N-2)*h);
    beta.at(0,N-1) = +6.0;
    beta.at(0,N-0) = -1.0;

    beta.at(1,N-3) = +1.0;
    beta.at(1,N-2) = -6.0;
    beta.at(1,N-1) = +3.0-6.0*h*A((N-1)*h);
    beta.at(1,N-0) = +2.0;

    beta.at(2,0) = +3.0;
    beta.at(2,N) = -2.0;

    DoubleVector qamma(K, 0.0);
    for (unsigned int i=0; i<=N; i++)
    {
        qamma.at(0) += beta.at(0,i)*x1.at(i);
        qamma.at(1) += beta.at(1,i)*x1.at(i);
        qamma.at(2) += beta.at(2,i)*x1.at(i);
    }

    DoubleVector x2(N+1, 0.0);
    calculateLeft2Right(N, K, a, beta, qamma, x2);
    IPrinter::printVector(14, 10, x2, "x2:", x2.size(), 0, 0, file1);
    IPrinter::printVector(14, 10, x0, "x0:", 10);
    IPrinter::printVector(14, 10, x2, "x2:", 10);
    fclose(file1);
}

void Example2::sample_n4()
{
    FILE *file1 = fopen("sample_n4_x.txt", "w");

    unsigned int N = 1000;
    unsigned int K = 4;
    double h = 1.0/N;

    DoubleVector x0(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        double t = i*h;
        x0.at(i) = X(t);
    }
    IPrinter::printVector(14, 10, x0, "x0:", x0.size(), 0, 0, file1);

    DoubleMatrix a(N-K+1, K+1);

    DoubleVector x1(N+1);

    x1.at(N-0) = X((N-0)*h);
    x1.at(N-1) = X((N-1)*h);
    x1.at(N-2) = X((N-2)*h);
    x1.at(N-3) = X((N-3)*h);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        double t  = i*h;
        double m  = 25.0+12.0*h*A(t);
        a.at(i,0) = -12.0*h*B(t)/m;
        a.at(i,1) = +48.0/m;
        a.at(i,2) = -36.0/m;
        a.at(i,3) = +16.0/m;
        a.at(i,4) = -3.0/m;

        x1.at(i) = a.at(i,1)*x1.at(i+1) + a.at(i,2)*x1.at(i+2) + a.at(i,3)*x1.at(i+3) + a.at(i,4)*x1.at(i+4) + a.at(i,0);
    }
    IPrinter::printVector(14, 10, x1, "x1:", x1.size(), 0, 0, file1);

    DoubleMatrix beta(K, N+1, 0.0);

    beta.at(0,N-4) = -3.0;
    beta.at(0,N-3) = -10.0-12.0*h*A((N-3)*h);
    beta.at(0,N-2) = +18.0;
    beta.at(0,N-1) = -6.0;
    beta.at(0,N-0) = +1.0;

    beta.at(1,N-4) = +1.0;
    beta.at(1,N-3) = -8.0;
    beta.at(1,N-2) = -12.0*h*A((N-2)*h);
    beta.at(1,N-1) = +8.0;
    beta.at(1,N-0) = -1.0;

    beta.at(2,N-4) = -1.0;
    beta.at(2,N-3) = +6.0;
    beta.at(2,N-2) = -18.0;
    beta.at(2,N-1) = +10.0-12.0*h*A((N-1)*h);
    beta.at(2,N-0) = +3.0;

    beta.at(3,0) = +3.0;
    beta.at(3,N) = -2.0;

    DoubleVector qamma(K);
    qamma.at(0) = 0.0;
    qamma.at(1) = 0.0;
    qamma.at(2) = 0.0;
    qamma.at(3) = 0.0;
    for (unsigned int i=0; i<=N; i++)
    {
        qamma.at(0) += beta.at(0,i)*x1.at(i);
        qamma.at(1) += beta.at(1,i)*x1.at(i);
        qamma.at(2) += beta.at(2,i)*x1.at(i);
        qamma.at(3) += beta.at(3,i)*x1.at(i);
    }

    DoubleVector x2(N+1, 0.0);
    calculateLeft2Right(N, K, a, beta, qamma, x2);
    IPrinter::printVector(14, 10, x2, "x2:", x2.size(), 0, 0, file1);
    IPrinter::printVector(14, 10, x0, "x0:", 10);
    IPrinter::printVector(14, 10, x2, "x2:", 10);
    fclose(file1);
}

void Example2::sample_n6()
{
    FILE *file1 = fopen("sample_n6_x.txt", "w");

    unsigned int N = 1000;
    unsigned int K = 6;
    double h = 1.0/N;

    DoubleVector x0(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        double t = i*h;
        x0.at(i) = X(t);
    }
    IPrinter::printVector(14, 10, x0, "x0:", x0.size(), 0, 0, file1);

    DoubleMatrix a(N-K+1, K+1);

    DoubleVector x1(N+1);

    x1.at(N-0) = X((N-0)*h);
    x1.at(N-1) = X((N-1)*h);
    x1.at(N-2) = X((N-2)*h);
    x1.at(N-3) = X((N-3)*h);
    x1.at(N-4) = X((N-4)*h);
    x1.at(N-5) = X((N-5)*h);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        double t  = i*h;
        double m  = 147.0+60.0*h*A(t);

        a.at(i,0) = -60.0*h*B(t)/m;
        a.at(i,1) = +360.0/m;
        a.at(i,2) = -450.0/m;
        a.at(i,3) = +400.0/m;
        a.at(i,4) = -225.0/m;
        a.at(i,5) = +72.0/m;
        a.at(i,6) = -10.0/m;

        x1.at(i) = a.at(i,1)*x1.at(i+1) + a.at(i,2)*x1.at(i+2) + a.at(i,3)*x1.at(i+3) + a.at(i,4)*x1.at(i+4) + a.at(i,5)*x1.at(i+5) + a.at(i,6)*x1.at(i+6) + a.at(i,0);
    }
    IPrinter::printVector(14, 10, x1, "x1:", x1.size(), 0, 0, file1);

    DoubleMatrix beta(K, N+1, 0.0);

    beta.at(0,N-6) = -10.0;
    beta.at(0,N-5) = -77.0-60.0*h*A((N-5)*h);
    beta.at(0,N-4) = +150.0;
    beta.at(0,N-3) = -100.0;
    beta.at(0,N-2) = +50.0;
    beta.at(0,N-1) = -15.0;
    beta.at(0,N-0) = +2.0;

    beta.at(1,N-6) = +2.0;
    beta.at(1,N-5) = -24.0;
    beta.at(1,N-4) = -35.0-60.0*h*A((N-4)*h);
    beta.at(1,N-3) = +80.0;
    beta.at(1,N-2) = -30.0;
    beta.at(1,N-1) = +8.0;
    beta.at(1,N-0) = -1.0;

    beta.at(2,N-6) = -1.0;
    beta.at(2,N-5) = +9.0;
    beta.at(2,N-4) = -45.0;
    beta.at(2,N-3) = -60.0*h*A((N-3)*h);
    beta.at(2,N-2) = +45.0;
    beta.at(2,N-1) = -9.0;
    beta.at(2,N-0) = +1.0;

    beta.at(3,N-6) = +1.0;
    beta.at(3,N-5) = -8.0;
    beta.at(3,N-4) = +30.0;
    beta.at(3,N-3) = -80.0;
    beta.at(3,N-2) = +35.0-60.0*h*A((N-2)*h);
    beta.at(3,N-1) = +24.0;
    beta.at(3,N-0) = -2.0;

    beta.at(4,N-6) = -2.0;
    beta.at(4,N-5) = +15.0;
    beta.at(4,N-4) = -50.0;
    beta.at(4,N-3) = +100.0;
    beta.at(4,N-2) = -150.0;
    beta.at(4,N-1) = +77.0-60.0*h*A((N-1)*h);
    beta.at(4,N-0) = +10.0;

    //    beta.at(5,N-6) = +10.0;
    //    beta.at(5,N-5) = -72.0;
    //    beta.at(5,N-4) = +225.0;
    //    beta.at(5,N-3) = -400.0;
    //    beta.at(5,N-2) = +450.0;
    //    beta.at(5,N-1) = -360.0;
    //    beta.at(5,N-0) = +147.0-60.0*h*A((N-0)*h);

    beta.at(5,0) = +3.0;
    //beta.at(3,N/2) = -1.0;
    beta.at(5,N) = -2.0;

    DoubleVector qamma(K);
    qamma.at(0) = 0.0;
    qamma.at(1) = 0.0;
    qamma.at(2) = 0.0;
    qamma.at(3) = 0.0;
    qamma.at(4) = 0.0;
    qamma.at(5) = 0.0;
    for (unsigned int i=0; i<=N; i++)
    {
        qamma.at(0) += beta.at(0,i)*x1.at(i);
        qamma.at(1) += beta.at(1,i)*x1.at(i);
        qamma.at(2) += beta.at(2,i)*x1.at(i);
        qamma.at(3) += beta.at(3,i)*x1.at(i);
        qamma.at(4) += beta.at(4,i)*x1.at(i);
        qamma.at(5) += beta.at(5,i)*x1.at(i);
    }

    DoubleVector x2(N+1, 0.0);
    calculateLeft2Right(N, K, a, beta, qamma, x2);
    IPrinter::printVector(14, 10, x2, "x2:", x2.size(), 0, 0, file1);
    IPrinter::printVector(14, 10, x0, "x0:", 10);
    IPrinter::printVector(14, 10, x2, "x2:", 10);
    fclose(file1);
}

void Example2::sample_n4_test1()
{
    unsigned int N = 1000;
    double h = 0.001;
    unsigned int K = 4;

    DoubleVector x(N+1);
    DoubleVector p3(N+1);
    DoubleVector p2(N+1);
    DoubleVector p1(N+1);
    DoubleVector p0(N+1);
    DoubleVector q(N+1);

    x.at(0) = X(0*h);
    x.at(1) = X(1*h);
    x.at(2) = X(2*h);
    x.at(3) = X(3*h);

    for (unsigned int i=K; i<=N; i++)
    {
        double t = i*h;
        double m = 12.0*h*A(t)-25.0;

        if (i==K)
        {
            p3.at(i) = -48.0/m;
            p2.at(i) = +36.0/m;
            p1.at(i) = -16.0/m;
            p0.at(i) = +3.0/m;
            q.at(i)  = -12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(3) + p2.at(i)*x.at(2) + p1.at(i)*x.at(1) + p0.at(i)*x.at(0) + q.at(i);
            //printf("%4d %.10f\n", i, x.at(i));
        }
        if (i == K+1)
        {
            p3.at(i) = (-48.0/m)*p3.at(i-1) + (+36.0/m);
            p2.at(i) = (-48.0/m)*p2.at(i-1) + (-16.0/m);
            p1.at(i) = (-48.0/m)*p1.at(i-1) + (+3.0/m);
            p0.at(i) = (-48.0/m)*p0.at(i-1);
            q.at(i)  = (-48.0/m)*q.at(i-1) - 12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(3) + p2.at(i)*x.at(2) + p1.at(i)*x.at(1) + p0.at(i)*x.at(0) + q.at(i);
            //printf("%4d %.10f\n", i, x.at(i));
        }
        if (i == K+2)
        {
            p3.at(i) = (-48.0/m)*p3.at(i-1) + (+36.0/m)*p3.at(i-2) + (-16.0/m);
            p2.at(i) = (-48.0/m)*p2.at(i-1) + (+36.0/m)*p2.at(i-2) + (+3.0/m);
            p1.at(i) = (-48.0/m)*p1.at(i-1) + (+36.0/m)*p1.at(i-2);
            p0.at(i) = (-48.0/m)*p0.at(i-1) + (+36.0/m)*p0.at(i-2);
            q.at(i)  = (-48.0/m)*q.at(i-1)  + (+36.0/m)*q.at(i-2) - 12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(3) + p2.at(i)*x.at(2) + p1.at(i)*x.at(1) + p0.at(i)*x.at(0) + q.at(i);
            //printf("%4d %.10f\n", i, x.at(i));
        }
        if (i == K+3)
        {
            p3.at(i) = (-48.0/m)*p3.at(i-1) + (+36.0/m)*p3.at(i-2) + (-16.0/m)*p3.at(i-3) + (+3.0/m);
            p2.at(i) = (-48.0/m)*p2.at(i-1) + (+36.0/m)*p2.at(i-2) + (-16.0/m)*p2.at(i-3);
            p1.at(i) = (-48.0/m)*p1.at(i-1) + (+36.0/m)*p1.at(i-2) + (-16.0/m)*p1.at(i-3);
            p0.at(i) = (-48.0/m)*p0.at(i-1) + (+36.0/m)*p0.at(i-2) + (-16.0/m)*p0.at(i-3);
            q.at(i)  = (-48.0/m)*q.at(i-1)  + (+36.0/m)*q.at(i-2)  + (-16.0/m)*q.at(i-3) - 12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(3) + p2.at(i)*x.at(2) + p1.at(i)*x.at(1) + p0.at(i)*x.at(0) + q.at(i);
            //printf("%4d %.10f\n", i, x.at(i));
        }
        if (i >= 2*K)
        {
            p3.at(i) = (-48.0/m)*p3.at(i-1) + (+36.0/m)*p3.at(i-2) + (-16.0/m)*p3.at(i-3) + (+3.0/m)*p3.at(i-4);
            p2.at(i) = (-48.0/m)*p2.at(i-1) + (+36.0/m)*p2.at(i-2) + (-16.0/m)*p2.at(i-3) + (+3.0/m)*p2.at(i-4);
            p1.at(i) = (-48.0/m)*p1.at(i-1) + (+36.0/m)*p1.at(i-2) + (-16.0/m)*p1.at(i-3) + (+3.0/m)*p1.at(i-4);
            p0.at(i) = (-48.0/m)*p0.at(i-1) + (+36.0/m)*p0.at(i-2) + (-16.0/m)*p0.at(i-3) + (+3.0/m)*p0.at(i-4);
            q.at(i)  = (-48.0/m)*q.at(i-1)  + (+36.0/m)*q.at(i-2)  + (-16.0/m)*q.at(i-3) +  (+3.0/m)*q.at(i-4) - 12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(3) + p2.at(i)*x.at(2) + p1.at(i)*x.at(1) + p0.at(i)*x.at(0) + q.at(i);
            //printf("%4d %.10f\n", i, x.at(i));
        }
    }

    IPrinter::printVector(14,10,x,"x5:");

    q.clear();
    p0.clear();
    p1.clear();
    p2.clear();
    p3.clear();
    x.clear();
}

void Example2::sample_n4_test2()
{
    unsigned int N = 1000;
    double h = 0.001;
    unsigned int K = 4;

    DoubleVector x(N+1);
    DoubleVector p3(N+1);
    DoubleVector p2(N+1);
    DoubleVector p1(N+1);
    DoubleVector p0(N+1);
    DoubleVector q(N+1);

    x.at(N-0) = X((N-0)*h);
    x.at(N-1) = X((N-1)*h);
    x.at(N-2) = X((N-2)*h);
    x.at(N-3) = X((N-3)*h);

    for (unsigned int i=N-K; i!=UINT32_MAX; i--)
    {
        double t = i*h;
        double m = 12.0*h*A(t)+25.0;

        if (i==N-K)
        {
            p3.at(i) = +48.0/m;
            p2.at(i) = -36.0/m;
            p1.at(i) = +16.0/m;
            p0.at(i) = -3.0/m;
            q.at(i)  = -12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(N-3) + p2.at(i)*x.at(N-2) + p1.at(i)*x.at(N-1) + p0.at(i)*x.at(N-0) + q.at(i);
            //printf("%4d %.10f %.10f\n", i, t, x.at(i));
        }
        if (i == N-K-1)
        {
            p3.at(i) = (48.0/m)*p3.at(i+1) + (-36.0/m);
            p2.at(i) = (48.0/m)*p2.at(i+1) + (+16.0/m);
            p1.at(i) = (48.0/m)*p1.at(i+1) + (-3.0/m);
            p0.at(i) = (48.0/m)*p0.at(i+1);
            q.at(i)  = (48.0/m)*q.at(i+1) - 12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(N-3) + p2.at(i)*x.at(N-2) + p1.at(i)*x.at(N-1) + p0.at(i)*x.at(N-0) + q.at(i);
            //printf("%4d %.10f %.10f\n", i, t, x.at(i));
        }
        if (i == N-K-2)
        {
            p3.at(i) = (48.0/m)*p3.at(i+1) + (-36.0/m)*p3.at(i+2) + (+16.0/m);
            p2.at(i) = (48.0/m)*p2.at(i+1) + (-36.0/m)*p2.at(i+2) + (-3.0/m);
            p1.at(i) = (48.0/m)*p1.at(i+1) + (-36.0/m)*p1.at(i+2);
            p0.at(i) = (48.0/m)*p0.at(i+1) + (-36.0/m)*p0.at(i+2);
            q.at(i)  = (48.0/m)*q.at(i+1)  + (-36.0/m)*q.at(i+2) - 12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(N-3) + p2.at(i)*x.at(N-2) + p1.at(i)*x.at(N-1) + p0.at(i)*x.at(N-0) + q.at(i);
            //printf("%4d %.10f %.10f\n", i, t, x.at(i));
        }
        if (i == N-K-3)
        {
            p3.at(i) = (+48.0/m)*p3.at(i+1) + (-36.0/m)*p3.at(i+2) + (+16.0/m)*p3.at(i+3) + (-3.0/m);
            p2.at(i) = (+48.0/m)*p2.at(i+1) + (-36.0/m)*p2.at(i+2) + (+16.0/m)*p2.at(i+3);
            p1.at(i) = (+48.0/m)*p1.at(i+1) + (-36.0/m)*p1.at(i+2) + (+16.0/m)*p1.at(i+3);
            p0.at(i) = (+48.0/m)*p0.at(i+1) + (-36.0/m)*p0.at(i+2) + (+16.0/m)*p0.at(i+3);
            q.at(i)  = (+48.0/m)*q.at(i+1)  + (-36.0/m)*q.at(i+2)  + (+16.0/m)*q.at(i+3) - 12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(N-3) + p2.at(i)*x.at(N-2) + p1.at(i)*x.at(N-1) + p0.at(i)*x.at(N-0) + q.at(i);
            //printf("%4d %.10f %.10f\n", i, t, x.at(i));
        }
        if (i <= N-2*K)
        {
            p3.at(i) = (+48.0/m)*p3.at(i+1) + (-36.0/m)*p3.at(i+2) + (+16.0/m)*p3.at(i+3) + (-3.0/m)*p3.at(i+4);
            p2.at(i) = (+48.0/m)*p2.at(i+1) + (-36.0/m)*p2.at(i+2) + (+16.0/m)*p2.at(i+3) + (-3.0/m)*p2.at(i+4);
            p1.at(i) = (+48.0/m)*p1.at(i+1) + (-36.0/m)*p1.at(i+2) + (+16.0/m)*p1.at(i+3) + (-3.0/m)*p1.at(i+4);
            p0.at(i) = (+48.0/m)*p0.at(i+1) + (-36.0/m)*p0.at(i+2) + (+16.0/m)*p0.at(i+3) + (-3.0/m)*p0.at(i+4);
            q.at(i)  = (+48.0/m)*q.at(i+1)  + (-36.0/m)*q.at(i+2)  + (+16.0/m)*q.at(i+3) +  (-3.0/m)*q.at(i+4) - 12.0*h*B(t)/m;
            x.at(i)  = p3.at(i)*x.at(N-3) + p2.at(i)*x.at(N-2) + p1.at(i)*x.at(N-1) + p0.at(i)*x.at(N-0) + q.at(i);
            //printf("%4d %.10f %.10f\n", i, t, x.at(i));
        }
    }

    IPrinter::printVector(14,10,x,"x5:");

    q.clear();
    p0.clear();
    p1.clear();
    p2.clear();
    p3.clear();
    x.clear();
}

