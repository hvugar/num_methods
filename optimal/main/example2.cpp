#include "example2.h"

Example2::Example2()
{}

double my_rand()
{
   return ((rand() % 400) + 1 - 200) / 1000.0;
}

void Example2::calculate()
{
    unsigned int N = 100;
    unsigned int K = 2;

    DoubleMatrix a(N-K+1, K+1);

    for (unsigned int j=0; j<N-K+1; j++)
    {
        for (unsigned int i=0; i<K+1; i++)
        {
            a.at(j,i) = my_rand();
        }
    }

    DoubleVector x(N+1);
    for (unsigned int i=0; i<K; i++) x.at(N-i) = my_rand();
    printf("%f %f %f\n", x[N-2], x[N-1], x[N]);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x.at(i) += a.at(i,j)*x[i+j];
    }
    IPrinter::printVector(x,"x:");

    DoubleMatrix beta(K,N+1);
    beta.at(0,0) = my_rand(); beta.at(0,1) = my_rand();
    beta.at(1,2) = my_rand(); beta.at(1,4) = my_rand(); beta.at(1,N-1) = my_rand(); beta.at(1,N) = my_rand();
    //beta.at(2,0) = my_rand(); beta.at(2,5) = my_rand(); beta.at(2,10) = my_rand();

    DoubleVector qamma(K);
    qamma[0] = beta.at(0,0)*x[0] + beta.at(0,1)*x[1];
    qamma[1] = beta.at(1,2)*x[2] + beta.at(1,4)*x[4] + beta.at(1,N-1)*x[N-1] + beta.at(1,N)*x[N];
    //qamma[2] = beta.at(2,0)*x[0] + beta.at(2,5)*x[5] + beta.at(2,10)*x[10];

    for (unsigned int i=0; i<=N; i++) x[i] = 0.0;
    calculate(N,K,a,beta,qamma,x);
    IPrinter::printVector(x,"x:");
}

void Example2::calculate(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x)
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

    for (unsigned int i=0; i<K; i++) x[N-i] = x1.at((K-1)-i);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x[i] += a.at(i,j)*x[i+j];
    }
}
