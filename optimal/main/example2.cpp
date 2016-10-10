#include "example2.h"
#include <math.h>

Example2::Example2()
{
    init();
}

double my_rand()
{
    return ((rand() % 100) + 1) / 100.0;
}

void Example2::init()
{
    unsigned int N = 50;
    unsigned int K = 2;

    DoubleMatrix a(N-K+1, K+1);
    DoubleVector x(N+1);

    for (unsigned int i=0; i<=N; i++)
    {
        x[i] = (5.0/(i+3.0))*cos((M_PI/50.0)*i);
    }
    IPrinter::printVector(x,"x:");
    printf("%12.8f %12.8f\n", x[N-1], x[N]);

    for (unsigned int i=0; i<(N+1)-K; i++)
    {
        a.at(i,1) = (i+1)/20.0;
        a.at(i,2) = ((i+1)*(i+2))/2000.0;
        a.at(i,0) = x[i] - a.at(i,1)*x.at(i+1) - a.at(i,2)*x.at(i+2);
    }

//    for (unsigned int i=0; i<(N+1)-K; i++)
//    {
//        a.at(i,1) = my_rand();//(i+1)/20.0;
//        a.at(i,2) = my_rand();//((i+1)*(i+2))/2000.0;
//        a.at(i,0) = my_rand();//x[i] - a.at(i,1)*x.at(i+1) - a.at(i,2)*x.at(i+2);
//    }


//    for (unsigned int i=0; i<K; i++) x.at(N-i) = my_rand();
//    for (unsigned int i=N-K; i != UINT_MAX; i--)
//    {
//        x[i] = a.at(i,0);
//        for (unsigned int j=1; j<=K; j++) x.at(i) += a.at(i,j)*x[i+j];
//    }
//    IPrinter::printVector(x,"x:");

    DoubleMatrix beta(K,N+1);
    beta.at(0,0) = 2.0; beta.at(0,1) = -1.0;
    beta.at(1,2) = 1.0; beta.at(1,4) = +3.0; beta.at(1,N-1) = -2.0; beta.at(1,N) = 1.0;

    DoubleVector qamma(K);
    qamma[0] = beta.at(0,0)*x[0] + beta.at(0,1)*x[1];
    qamma[1] = beta.at(1,2)*x[2] + beta.at(1,4)*x[4] + beta.at(1,N-1)*x[N-1] + beta.at(1,N)*x[N];

    DoubleVector x1(N+1);
    for (unsigned int i=0; i<=N; i++) x1[i] = 0.0;
    calculate1(N,K,a,beta,qamma,x1);
    IPrinter::printVector(x1,"x:");


//    puts("-----------------------------------");
//    DoubleMatrix a1(N-K+1, K+1);
//    for (unsigned int i=0; i<(N+1)-K; i++)
//    {
//        a1.at(i,1) = -a.at(i,1)/a.at(i,2);
//        a1.at(i,2) = 1.0/a.at(i,2);
//        a1.at(i,0) = -a.at(i,0)/a.at(i,2);
//    }
//    puts("---");
//    DoubleMatrix beta1(K,N+1);
//    beta1.at(0,0) = 2.0; beta1.at(0,1) = -1.0;
//    beta1.at(1,2) = 1.0; beta1.at(1,4) = +3.0; beta1.at(1,N-1) = -2.0; beta1.at(1,N) = 1.0;
//    DoubleVector qamma1(K);
//    qamma1[0] = beta1.at(0,0)*x[0] + beta1.at(0,1)*x[1];
//    qamma1[1] = beta1.at(1,2)*x[2] + beta1.at(1,4)*x[4] + beta1.at(1,N-1)*x[N-1] + beta1.at(1,N)*x[N];
//    DoubleVector x2(N+1);
//    for (unsigned int i=0; i<=N; i++) x2[i] = 0.0;
//    calculate2(N,K,a1,beta1,qamma1,x2);
//    IPrinter::printVector(x2,"x:");
}

void Example2::calculate()
{
    unsigned int N = 1000;
    unsigned int K = 2;

    DoubleMatrix a(N-K+1, K+1);

    for (unsigned int j=0; j<N-K+1; j++)
    {
        for (unsigned int i=0; i<K+1; i++)
        {
            a.at(j,i) = pow(-1,j)*my_rand();
        }
    }

    DoubleVector x(N+1);
    for (unsigned int i=0; i<K; i++) x.at(N-i) = 1.0;//my_rand();
    printf("%f %f\n", x[N-1], x[N]);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x.at(i) += a.at(i,j)*x[i+j];
    }
    IPrinter::printVector(x,"x:");

    DoubleMatrix beta(K,N+1);
    beta.at(0,0) = my_rand(); beta.at(0,1) = my_rand();
    beta.at(1,2) = my_rand(); beta.at(1,4) = my_rand(); beta.at(1,N-1) = my_rand(); beta.at(1,N) = my_rand();

    DoubleVector qamma(K);
    qamma[0] = beta.at(0,0)*x[0] + beta.at(0,1)*x[1];
    qamma[1] = beta.at(1,2)*x[2] + beta.at(1,4)*x[4] + beta.at(1,N-1)*x[N-1] + beta.at(1,N)*x[N];

    DoubleVector x1(N+1);
    for (unsigned int i=0; i<=N; i++) x1[i] = 0.0;
    calculate1(N,K,a,beta,qamma,x1);
    IPrinter::printVector(x1,"x:");
}

void Example2::calculate1(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x)
{
    for (unsigned int eq=0; eq<K; eq++)
    {
        for (unsigned int i=0; i<=N-K; i++)
        {
            for (unsigned int j=1; j<=K; j++)
            {
                beta.at(eq,i+j) = beta.at(eq,i+j) + beta.at(eq,i)*a.at(i,j);
            }
            qamma[eq] = qamma[eq] - beta.at(eq,i)*a.at(i,0);
        }
    }

    DoubleMatrix m(K,K);
    for (unsigned int j=0; j<K; j++)
    {
        for (unsigned int i=K-1; i != UINT_MAX; i--) m.at(j,(K-1)-i) = beta.at(j,N-i);
    }
    DoubleVector b(K);
    for (unsigned int i=0; i<K; i++) b.at(i) = qamma.at(i);

    DoubleVector x1(K);
    GaussianElimination(m, b, x1);

    printf("%12.8f %12.8f\n", x1[0], x1[1]);

    for (unsigned int i=0; i<K; i++) x[N-i] = x1.at((K-1)-i);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x[i] += a.at(i,j)*x[i+j];
    }
}

void Example2::calculate2(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x)
{
    for (unsigned int eq=0; eq<K; eq++)
    {
        for (unsigned int i=N; i>=K; i--)
        {
            for (unsigned int j=1; j<=K; j++)
            {
                //printf("%d %d %d %f %f %f\n", eq, i-K, j, a.at(i-K,0), a.at(i-K,1), a.at(i-K,2));
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

    printf("%12.8f %12.8f\n", x1[0], x1[1]);

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
