#include "example2.h"

double my_rand()
{
    return ((rand() % 100) + 1) / 100.0;
}

void Example2::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    Example2 e2;
    e2.sample_n4();
    //    e2.sample1();
}

Example2::Example2()
{
}

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
    calculateLeft2Right(N,K,a,beta,qamma,x2, x1);
    IPrinter::printVector(x2,"x2:",x2.size(),0,0,file1);
    fclose(file1);
}

void Example2::calculateLeft2Right(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x, const DoubleVector &x0)
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

    DoubleVector b(K);
    for (unsigned int i=0; i<K; i++) b.at(i) = qamma.at(i);

    DoubleVector x1(K);
    GaussianElimination(m, b, x1);
    //double a1 = m.at(1,1) - m.at(0,1)*m.at(1,0)/m.at(0,0);
    //double b1 = b.at(1) - b.at(0)*m.at(1,0)/m.at(0,0);
    //x1.at(1) = b1/a1;
    //x1.at(0) = -(m.at(0,1)/m.at(0,0))*x1.at(1) + b.at(0)/m.at(0,0);

    printf("x0 %18.14f %18.14f\n", x0.at(N-1), x0.at(N));
    printf("x1 %18.14f %18.14f\n", x1.at(0), x1.at(1));

    double c1 = m.at(0,0)*x0.at(N-1)+m.at(0,1)*x0.at(N);
    double e1 = m.at(0,0)*x1.at(0)+m.at(0,1)*x1.at(1);

    double c2 = m.at(1,0)*x0.at(N-1)+m.at(1,1)*x0.at(N);
    double e2 = m.at(1,0)*x1.at(0)+m.at(1,1)*x1.at(1);

    printf("a00: %18.14f a01: %18.14f b0 : %18.14f x0: %18.14f x1: %18.14f\n", m.at(0,0), m.at(0,1), b.at(0), c1, e1);
    printf("a10: %18.14f a11: %18.14f b1 : %18.14f x0: %18.14f x1: %18.14f\n", m.at(1,0), m.at(1,1), b.at(1), c2, e2);
    printf("k00: %18.14f k01: %18.14f k02: %18.14f\n", m.at(0,0)/m.at(1,0), m.at(0,1)/m.at(1,1), b.at(0)/b.at(1));
    printf("k01: %18.14f k11: %18.14f k12: %18.14f\n", m.at(1,0)/m.at(0,0), m.at(1,1)/m.at(0,1), b.at(1)/b.at(0));

    for (unsigned int i=0; i<K; i++) x.at(N-i) = x1.at((K-1)-i);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x[i] += a.at(i,j)*x.at(i+j);
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

void Example2::calculateRight2Left(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x, const DoubleVector &x0)
{
    for (unsigned int eq=0; eq<K; eq++)
    {
        for (unsigned int i=N; i>=K; i--)
        {
            for (unsigned int j=1; j<=K; j++)
            {
                beta.at(eq,i-j) = beta.at(eq,i-j) + beta.at(eq,i)*a.at(i-K,j);
            }
            qamma[eq] = qamma[eq] - beta.at(eq,i)*a.at(i-K,0);
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
    //double a1 = m.at(1,1) - m.at(0,1)*m.at(1,0)/m.at(0,0);
    //double b1 = b.at(1) - b.at(0)*m.at(1,0)/m.at(0,0);
    //x1.at(1) = b1/a1;
    //x1.at(0) = -(m.at(0,1)/m.at(0,0))*x1.at(1) + b.at(0)/m.at(0,0);

    printf("x0 %18.14f %18.14f\n", x0.at(0), x0.at(1));
    printf("x1 %18.14f %18.14f\n", x1.at(0), x1.at(1));

    double c1 = m.at(0,0)*x0.at(0)+m.at(0,1)*x0.at(1);
    double e1 = m.at(0,0)*x1.at(0)+m.at(0,1)*x1.at(1);

    double c2 = m.at(1,0)*x0.at(0)+m.at(1,1)*x0.at(1);
    double e2 = m.at(1,0)*x1.at(0)+m.at(1,1)*x1.at(1);

    printf("a00: %18.14f a01: %18.14f b0 : %18.14f x0: %18.14f x1: %18.14f\n", m.at(0,0), m.at(0,1), b.at(0), c1, e1);
    printf("a10: %18.14f a11: %18.14f b1 : %18.14f x0: %18.14f x1: %18.14f\n", m.at(1,0), m.at(1,1), b.at(1), c2, e2);
    printf("k00: %18.14f k01: %18.14f k02: %18.14f\n", m.at(0,0)/m.at(1,0), m.at(0,1)/m.at(1,1), b.at(0)/b.at(1));
    printf("k01: %18.14f k11: %18.14f k12: %18.14f\n", m.at(1,0)/m.at(0,0), m.at(1,1)/m.at(0,1), b.at(1)/b.at(0));

    for (unsigned int i=0; i<K; i++) x[i] = x1.at(i);

    for (unsigned int i=K; i <= N; i++)
    {
        x[i] = a.at(i-K,0);
        for (unsigned int j=1; j<=K; j++)
        {
            x[i] += a.at(i-K,j)*x[i-j];
        }
    }
}

double Example2::X(double t) const
{
    //    return 3.0*sin(t);
    return t*t;
    //    return (t*t-t)*sin(50.0*t);
}

double Example2::A(double t) const
{
    //    return t;
    return 3.0;
    //    return t;
}

double Example2::B(double t) const
{
    //    return 3.0*cos(t) - 3.0*t*sin(t);
    return 2.0*t - 3.0*t*t;
    //    return 50.0*(t*t-t)*cos(50.0*t) + (2.0*t-1)*sin(50.0*t) - t*(t*t-t)*sin(50.0*t);
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
    calculateLeft2Right(N,K,a,beta,qamma,x2, x1);
    IPrinter::printVector(x2,"x2:",x2.size(),0,0,file1);
    fclose(file1);
}

void Example2::sample_n4()
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
    IPrinter::printVector(18, 14, x0, "x0:", x0.size(), 0, 0, file1);

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
    calculateLeft2Right4(N,K,a,beta,qamma, x, x1);
    IPrinter::printVector(18, 14, x, "x2:", x.size(), 0, 0, file1);
    fclose(file1);
}

void Example2::calculateLeft2Right4(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x, const DoubleVector &x0)
{
    //    FILE *f = fopen("betta_1.txt", "w");
    //    IPrinter::printVector(beta.row(eq), NULL, beta.cols(), 0, 0, f);
    //    for (unsigned int i=0; i<=N-K; i++)
    //    {
    //        beta.at(eq,i+1) = beta.at(eq,i+1) + beta.at(eq,i)*a.at(i,1);
    //        beta.at(eq,i+2) = beta.at(eq,i+2) + beta.at(eq,i)*a.at(i,2);
    //        beta.at(eq,i+3) = beta.at(eq,i+3) + beta.at(eq,i)*a.at(i,3);
    //        beta.at(eq,i+4) = beta.at(eq,i+4) + beta.at(eq,i)*a.at(i,4);
    //        qamma.at(eq)    = qamma.at(eq)    - beta.at(eq,i)*a.at(i,0);
    //        IPrinter::printVector(beta.row(eq), NULL, beta.cols(), 0, 0, f);
    //    }
    //    fclose(f);

    for (unsigned int eq=0; eq<K; eq++)
    {
        FILE *f = NULL;
        if (eq==0) f = fopen("betta_1.txt", "w");
        if (eq==1) f = fopen("betta_2.txt", "w");
        if (eq==2) f = fopen("betta_3.txt", "w");
        if (eq==3) f = fopen("betta_4.txt", "w");
        IPrinter::printVector(beta.row(eq), NULL, beta.cols(), 0, 0, f);
        for (unsigned int i=0; i<=N-K; i++)
        {
            beta.at(eq,i+1) = beta.at(eq,i+1) + beta.at(eq,i)*a.at(i,1);
            beta.at(eq,i+2) = beta.at(eq,i+2) + beta.at(eq,i)*a.at(i,2);
            beta.at(eq,i+3) = beta.at(eq,i+3) + beta.at(eq,i)*a.at(i,3);
            beta.at(eq,i+4) = beta.at(eq,i+4) + beta.at(eq,i)*a.at(i,4);
            qamma.at(eq)    = qamma.at(eq)    - beta.at(eq,i)*a.at(i,0);
            IPrinter::printVector(beta.row(eq), NULL, beta.cols(), 0, 0, f);
            //            for (unsigned int j=1; j<=K; j++)
            //            {
            //                beta.at(eq,i+j) = beta.at(eq,i+j) + beta.at(eq,i)*a.at(i,j);
            //            }
        }
        fclose(f);
    }

    DoubleMatrix m(K,K);
    m.at(0,0) = beta.at(0,N-3);
    m.at(0,1) = beta.at(0,N-2);
    m.at(0,2) = beta.at(0,N-1);
    m.at(0,3) = beta.at(0,N-0);

    m.at(1,0) = beta.at(1,N-3);
    m.at(1,1) = beta.at(1,N-2);
    m.at(1,2) = beta.at(1,N-1);
    m.at(1,3) = beta.at(1,N-0);

    m.at(2,0) = beta.at(2,N-3);
    m.at(2,1) = beta.at(2,N-2);
    m.at(2,2) = beta.at(2,N-1);
    m.at(2,3) = beta.at(2,N-0);

    m.at(3,0) = beta.at(3,N-3);
    m.at(3,1) = beta.at(3,N-2);
    m.at(3,2) = beta.at(3,N-1);
    m.at(3,3) = beta.at(3,N-0);

    puts("------------------------------------");
    IPrinter::printMatrix(18, 14, m, m.rows(), m.cols());
    puts("------------------------------------");

    printf("Determinant: %.14f\n", m.determinant());

    //    for (unsigned int j=0; j<K; j++)
    //    {
    //        for (unsigned int i=0; i<K; i++) m.at(j,i) = beta.at(j,N-K+1+i);
    //    }

    DoubleVector b(K);
    //    for (unsigned int i=0; i<K; i++) b.at(i) = qamma.at(i);
    b.at(0) = qamma.at(0);
    b.at(1) = qamma.at(1);
    b.at(2) = qamma.at(2);
    b.at(3) = qamma.at(3);

    puts("------------------------------------------------------------------");
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | %18.14f\n", m.at(0,0), m.at(0,1), m.at(0,2), m.at(0,3), b.at(0),
           m.at(0,0)*x0.at(N-3)+m.at(0,1)*x0.at(N-2)+m.at(0,2)*x0.at(N-1)+m.at(0,3)*x0.at(N-0));
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | %18.14f\n", m.at(1,0), m.at(1,1), m.at(1,2), m.at(1,3), b.at(1),
           m.at(1,0)*x0.at(N-3)+m.at(1,1)*x0.at(N-2)+m.at(1,2)*x0.at(N-1)+m.at(1,3)*x0.at(N-0));
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | %18.14f\n", m.at(2,0), m.at(2,1), m.at(2,2), m.at(2,3), b.at(2),
           m.at(2,0)*x0.at(N-3)+m.at(2,1)*x0.at(N-2)+m.at(2,2)*x0.at(N-1)+m.at(2,3)*x0.at(N-0));
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | %18.14f\n", m.at(3,0), m.at(3,1), m.at(3,2), m.at(3,3), b.at(3),
           m.at(3,0)*x0.at(N-3)+m.at(3,1)*x0.at(N-2)+m.at(3,2)*x0.at(N-1)+m.at(3,3)*x0.at(N-0));
    puts("------------------------------------------------------------------");

    //DoubleMatrix m1 = m;
    //DoubleVector b1 = b;

    DoubleVector x1(K);
    //printf("%18.10f %18.10f %18.10f %18.10f\n", m.at(0,0), m.at(0,1), m.at(0,2), m.at(0,3));
    //printf("%18.10f %18.10f %18.10f %18.10f\n", m.at(1,0), m.at(1,1), m.at(1,2), m.at(1,3));
    //printf("%18.10f %18.10f %18.10f %18.10f\n", m.at(2,0), m.at(2,1), m.at(2,2), m.at(2,3));
    //printf("%18.10f %18.10f %18.10f %18.10f\n", m.at(3,0), m.at(3,1), m.at(3,2), m.at(3,3));
    GaussianElimination(m, b, x1);
    //printf("%18.10f %18.10f %18.10f %18.10f\n", m.at(0,0), m.at(0,1), m.at(0,2), m.at(0,3));
    //printf("%18.10f %18.10f %18.10f %18.10f\n", m.at(1,0), m.at(1,1), m.at(1,2), m.at(1,3));
    //printf("%18.10f %18.10f %18.10f %18.10f\n", m.at(2,0), m.at(2,1), m.at(2,2), m.at(2,3));
    //printf("%18.10f %18.10f %18.10f %18.10f\n", m.at(3,0), m.at(3,1), m.at(3,2), m.at(3,3));

    puts("------------------------------------------------------------------");
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | %18.14f\n", m.at(0,0), m.at(0,1), m.at(0,2), m.at(0,3), b.at(0),
           m.at(0,0)*x1.at(0)+m.at(0,1)*x1.at(1)+m.at(0,2)*x1.at(2)+m.at(0,3)*x1.at(3));
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | %18.14f\n", m.at(1,0), m.at(1,1), m.at(1,2), m.at(1,3), b.at(1),
           m.at(1,0)*x1.at(0)+m.at(1,1)*x1.at(1)+m.at(1,2)*x1.at(2)+m.at(1,3)*x1.at(3));
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | %18.14f\n", m.at(2,0), m.at(2,1), m.at(2,2), m.at(2,3), b.at(2),
           m.at(2,0)*x1.at(0)+m.at(2,1)*x1.at(1)+m.at(2,2)*x1.at(2)+m.at(2,3)*x1.at(3));
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | %18.14f\n", m.at(3,0), m.at(3,1), m.at(3,2), m.at(3,3), b.at(3),
           m.at(3,0)*x1.at(0)+m.at(3,1)*x1.at(1)+m.at(3,2)*x1.at(2)+m.at(3,3)*x1.at(3));
    puts("------------------------------------------------------------------");

    printf("x0 %18.14f %18.14f %18.14f %18.14f\n", x0.at(N-3), x0.at(N-2), x0.at(N-1), x0.at(N));
    printf("x1 %18.14f %18.14f %18.14f %18.14f\n", x1.at(0), x1.at(1), x1.at(2), x1.at(3));

    double c1 = m.at(0,0)*x0.at(N-3)+m.at(0,1)*x0.at(N-2)+m.at(0,2)*x0.at(N-1)+m.at(0,3)*x0.at(N);
    double e1 = m.at(0,0)*x1.at(0)+m.at(0,1)*x1.at(1)+m.at(0,2)*x1.at(2)+m.at(0,3)*x1.at(3);

    double c2 = m.at(1,0)*x0.at(N-3)+m.at(1,1)*x0.at(N-2)+m.at(1,2)*x0.at(N-1)+m.at(1,3)*x0.at(N);
    double e2 = m.at(1,0)*x1.at(0)+m.at(1,1)*x1.at(1)+m.at(1,2)*x1.at(2)+m.at(1,3)*x1.at(3);

    puts("------------------------------------------------------------------");
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | x1: %18.14f | x0: %18.14f\n", m.at(0,0), m.at(0,1), m.at(0,2), m.at(0,3), b.at(0), c1, e1);
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | x1: %18.14f | x0: %18.14f\n", m.at(1,0), m.at(1,1), m.at(1,2), m.at(1,3), b.at(1), c2, e2);
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | x1: %18.14f | x0: %18.14f\n", m.at(1,0), m.at(1,1), m.at(1,2), m.at(1,3), b.at(2), c2, e2);
    printf("%18.14f %18.14f %18.14f %18.14f | %18.14f | x1: %18.14f | x0: %18.14f\n", m.at(1,0), m.at(1,1), m.at(1,2), m.at(1,3), b.at(3), c2, e2);
    puts("------------------------------------------------------------------");

    //printf("k00: %18.14f k01: %18.14f k02: %18.14f\n", m.at(0,0)/m.at(1,0), m.at(0,1)/m.at(1,1), b.at(0)/b.at(1));
    //printf("k01: %18.14f k11: %18.14f k12: %18.14f\n", m.at(1,0)/m.at(0,0), m.at(1,1)/m.at(0,1), b.at(1)/b.at(0));

    for (unsigned int i=0; i<K; i++) x.at(N-i) = x1.at((K-1)-i);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x[i] += a.at(i,j)*x.at(i+j);
    }
}
