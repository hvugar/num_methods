#include "example2.h"
#include <math.h>

double my_rand()
{
    return ((rand() % 100) + 1) / 100.0;
}

void Example2::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    Example2 e2;
    e2.init2();
}

Example2::Example2()
{
}

void Example2::init1()
{
    unsigned int N = 1000;
    unsigned int K = 2;

    DoubleMatrix a(N-K+1, K+1);
    DoubleVector x(N+1);

    for (unsigned int i=0; i<=N; i++)
    {
        x[i] = (5.0/(i+3.0))*cos((M_PI/50.0)*i);
    }

    for (unsigned int i=0; i<(N+1)-K; i++)
    {
        a.at(i,1) = (i+1)/20.0;
        a.at(i,2) = ((i+1)*(i+2))/2000.0;
        a.at(i,0) = x[i] - a.at(i,1)*x.at(i+1) - a.at(i,2)*x.at(i+2);
    }

    FILE* file1 = fopen("data.txt", "w");
    IPrinter::printVector(x,NULL,x.size(),0,0,file1);
    fclose(file1);

    DoubleMatrix beta(K,N+1);
    beta.at(0,0) = 2.0; beta.at(0,1) = -1.0;
    beta.at(1,2) = 1.0; beta.at(1,4) = +3.0; beta.at(1,N-1) = -2.0; beta.at(1,N) = 1.0;

    DoubleVector qamma(K);
    qamma[0] = beta.at(0,0)*x[0] + beta.at(0,1)*x[1];
    qamma[1] = beta.at(1,2)*x[2] + beta.at(1,4)*x[4] + beta.at(1,N-1)*x[N-1] + beta.at(1,N)*x[N];

    DoubleVector x1(N+1);
    for (unsigned int i=0; i<=N; i++) x1[i] = 0.0;
    calculate1(N,K,a,beta,qamma,x1, x);
    FILE* file2 = fopen("data.txt", "a");
    IPrinter::printVector(x1,NULL,x1.size(),0,0,file2);
    fclose(file2);


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

void Example2::init2()
{
    unsigned int N = 1000;
    unsigned int K = 2;

    DoubleMatrix a(N-K+1, K+1);
    DoubleVector x(N+1, 0.0);

    for (unsigned int i=0; i<(N+1)-K; i++)
    {
//        a.at(i,1) = 0.9;
//        a.at(i,2) = 0.5;
//        a.at(i,0) = 0.2;

        a.at(i,1) = my_rand();
        a.at(i,2) = my_rand();
        a.at(i,0) = my_rand();
    }

    for (unsigned int i=N; i>N-K; i--)
    {
        x.at(i) = my_rand();
//        x.at(i) = 0.1*i;
    }

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x.at(i) += a.at(i,j)*x.at(i+j);
    }
    FILE* file1 = fopen("data.txt", "w");
    IPrinter::printVector(x,NULL,x.size(),0,0,file1);
    fclose(file1);

    DoubleMatrix beta(K,N+1);
    for (unsigned int i=0; i<beta.cols(); i++)
    {
        beta.at(0, i) = 0.0;
        beta.at(1, i) = 0.0;
    }

    unsigned int n11 = 5*N/10;
    unsigned int n12 = 8*N/10;
    unsigned int n21 = 2*N/10;
    unsigned int n22 = 4*N/10;

    //beta.at(0,0) = my_rand(); beta.at(0,N) = my_rand();
    //beta.at(0,N-1) = +1.2; beta.at(0,N) = -1.1;
    //beta.at(0,0) = +0.4; beta.at(0,N) = +0.5;
    beta.at(0,0) = +0.5; beta.at(0,n11) = +0.3; beta.at(0,n12) = -0.8; beta.at(0,N-1) = +1.4; beta.at(0,N) = +2.1;
    beta.at(1,0) = +1.0; beta.at(1,n21) = +0.1; beta.at(1,n22) = +1.3; beta.at(1,N-1) = -1.2; beta.at(1,N) = +1.1;

    DoubleVector qamma(K);
    //qamma.at(0) = beta.at(0,0)*x.at(0) + beta.at(0,1)*x.at(1);
    //qamma.at(1) = beta.at(1,0)*x.at(0) + beta.at(1,1)*x.at(1);

    qamma.at(0) = 0.0;
    for (unsigned int i=0; i<=N; i++)
    {
        qamma.at(0) += beta.at(0,i)*x.at(i);
        qamma.at(1) += beta.at(1,i)*x.at(i);
    }

    //qamma.at(0) = beta.at(0,n11)*x.at(n11) + beta.at(0,n12)*x.at(n12) + beta.at(0,N-1)*x.at(N-1) + beta.at(0,N)*x.at(N);
    //qamma.at(1) = beta.at(1,n21)*x.at(n21) + beta.at(1,n22)*x.at(n22) + beta.at(1,N-1)*x.at(N-1) + beta.at(1,N)*x.at(N);

    DoubleVector x1(N+1);
    for (unsigned int i=0; i<=N; i++) x1[i] = 0.0;
    calculate1(N,K,a,beta,qamma,x1, x);
    FILE* file2 = fopen("data.txt", "a");
    IPrinter::printVector(x1,NULL,x1.size(),0,0,file2);
    fclose(file2);

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

void Example2::calculate1(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x, const DoubleVector &x0)
{
    FILE *file1 = fopen("data1.txt", "w");
    FILE *file2 = fopen("data2.txt", "w");

    IPrinter::printVector(beta.row(0),NULL,beta.row(0).size(),0,0,file1);
    IPrinter::printVector(beta.row(1),NULL,beta.row(1).size(),0,0,file2);
    for (unsigned int eq=0; eq<K; eq++)
    {
        for (unsigned int i=0; i<=N-K; i++)
        {
//            for (unsigned int j=1; j<=K; j++)
//            {
//                beta.at(eq,i+j) = beta.at(eq,i+j) + beta.at(eq,i)*a.at(i,j);
//            }
            beta.at(eq,i+1) = beta.at(eq,i+1) + beta.at(eq,i)*a.at(i,1);
            beta.at(eq,i+2) = beta.at(eq,i+2) + beta.at(eq,i)*a.at(i,2);
            qamma.at(eq) = qamma.at(eq) - beta.at(eq,i)*a.at(i,0);

            if (eq == 0)
            {
                IPrinter::printVector(beta.row(0),NULL,beta.row(0).size(),0,0,file1);
            }

            if (eq == 1)
            {
                IPrinter::printVector(beta.row(1),NULL,beta.row(1).size(),0,0,file2);
            }
        }
    }

    fclose(file1);
    fclose(file2);

    DoubleMatrix m(K,K);
    m.at(0,0) = beta.at(0,N-1); m.at(0,1) = beta.at(0,N);
    m.at(1,0) = beta.at(1,N-1); m.at(1,1) = beta.at(1,N);

//    for (unsigned int j=0; j<K; j++)
//    {
//        for (unsigned int i=K-1; i != UINT_MAX; i--) m.at(j,(K-1)-i) = beta.at(j,N-i);
//    }
    DoubleVector b(K);
    b.at(0) = qamma.at(0);
    b.at(1) = qamma.at(1);

//    for (unsigned int i=0; i<K; i++) b.at(i) = qamma.at(i);

    DoubleVector x1(K);

    //GaussianElimination(m, b, x1);
    double a1 = m.at(1,1) - m.at(0,1)*m.at(1,0)/m.at(0,0);
    double b1 = b.at(1) - b.at(0)*m.at(1,0)/m.at(0,0);
    x1.at(1) = b1/a1;
    printf("%18.14f %18.14f\n", b1, a1);
    x1.at(0) = -(m.at(0,1)/m.at(0,0))*x1.at(1) + b.at(0)/m.at(0,0);

    printf("x0 %18.14f %18.14f\n", x0.at(N-1), x0.at(N));
    printf("x1 %18.14f %18.14f\n", x1.at(0), x1.at(1));
    double c1 = m.at(0,0)*x0.at(N-1)+m.at(0,1)*x0.at(N);
    double c2 = m.at(1,0)*x0.at(N-1)+m.at(1,1)*x0.at(N);
    double e1 = m.at(0,0)*x1.at(0)+m.at(0,1)*x1.at(1);
    double e2 = m.at(1,0)*x1.at(0)+m.at(1,1)*x1.at(1);

    printf("%18.14f %18.14f %18.14f %18.14f %18.14f\n", m.at(0,0), m.at(0,1), b.at(0), c1, e1);
    printf("%18.14f %18.14f %18.14f %18.14f %18.14f\n", m.at(1,0), m.at(1,1), b.at(1), c2, e2);
    printf("%18.14f %18.14f %18.14f %18.14f %18.14f\n", m.at(1,0)/m.at(0,0), m.at(1,1)/m.at(0,1), b.at(1)/b.at(0), 0.0, 0.0);

    for (unsigned int i=0; i<K; i++) x.at(N-i) = x1.at((K-1)-i);

    for (unsigned int i=N-K; i != UINT_MAX; i--)
    {
        x[i] = a.at(i,0);
        for (unsigned int j=1; j<=K; j++) x[i] += a.at(i,j)*x.at(i+j);
    }
}

void Example2::init3()
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
    qamma[0] = beta.at(0,0)*x.at(0) + beta.at(0,1)*x.at(1);
    qamma[1] = beta.at(1,2)*x.at(2) + beta.at(1,4)*x.at(4) + beta.at(1,N-1)*x.at(N-1) + beta.at(1,N)*x.at(N);

    DoubleVector x1(N+1);
    for (unsigned int i=0; i<=N; i++) x1[i] = 0.0;
    calculate1(N,K,a,beta,qamma,x1, x);
    IPrinter::printVector(x1,"x:");
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
