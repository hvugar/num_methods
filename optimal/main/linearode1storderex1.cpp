#include "linearode1storderex1.h"

#define EXAMPLE_1

void LinearODE1stOrderEx1::Main(int argc, char **argv)
{
    LinearODE1stOrderEx1 nl;

    unsigned int N = 100;
    unsigned int M = 3;
    double h = 0.01;

    std::vector<NonLocalCondition> C;
    C.resize(3);
    C[0].i = 0;
    C[0].n = PointNodeODE(0.0, 0);
    C[0].m.resize(M,M,0.0);

    C[1].i = 1;
    C[1].n = PointNodeODE(0.5, 50);
    C[1].m.resize(M,M,0.0);

    C[2].i = 2;
    C[2].n = PointNodeODE(1.0, 100);
    C[2].m.resize(M,M,0.0);

    for (unsigned int r=0; r<M; r++)
    {
        for (unsigned int c=0; c<M; c++)
        {
            C[0].m[r][c] = Random::value(0,1,4);
            C[1].m[r][c] = Random::value(0,1,4);
            C[2].m[r][c] = Random::value(0,1,4);
        }
    }

    DoubleVector x00; for (unsigned int m=1; m<=M; m++) x00 << nl.x(C[0].n,m);
    DoubleVector x05; for (unsigned int m=1; m<=M; m++) x05 << nl.x(C[1].n,m);
    DoubleVector x10; for (unsigned int m=1; m<=M; m++) x10 << nl.x(C[2].n,m);

    DoubleVector d = C[0].m*x00 + C[1].m*x05 + C[2].m*x10;

    std::vector<DoubleVector> x;
    nl.setDimension(Dimension(h, 0, static_cast<int>(N)));

    nl.transferOfCondition(C, d, x, 2, M);

    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            if (n%(N/10)==0) printf("%14.10f ", x[n][m]);
        }
        puts("");
    }
    IPrinter::printSeperatorLine();

    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            if (n%(N/10)==0) printf("%14.10f ", nl.x(n*h, m+1));
        }
        puts("");
    }
    IPrinter::printSeperatorLine();

    double norm[] = {0.0, 0.0, 0.0};
    for (unsigned int m=0; m<M; m++)
    {
        norm[m] = 0.0;
        for (unsigned int n=0; n<=N; n++)
        {
            norm[m] += (x[n][m]-nl.x(n*h, m+1))*(x[n][m]-nl.x(n*h, m+1));
        }
        norm[m] = sqrt(norm[m]);
    }
    printf("Norms: %.10f %.10f %.10f\n", norm[0], norm[1], norm[2]);
}

LinearODE1stOrderEx1::LinearODE1stOrderEx1()
{
}

LinearODE1stOrderEx1::~LinearODE1stOrderEx1()
{}

double LinearODE1stOrderEx1::A(const PointNodeODE &node, unsigned int r, unsigned int c) const
{
#ifdef EXAMPLE_1
    return 1.0;
#endif
#ifdef EXAMPLE_2
    return 1.0;
#endif
#ifdef EXAMPLE_3
    return 1.0;
#endif
#ifdef EXAMPLE_4
    return 1.0;
#endif
}

double LinearODE1stOrderEx1::B(const PointNodeODE &node, unsigned int r) const
{
    double t =  node.x;

#ifdef EXAMPLE_1
    if (r == 1) return -(A(node,1,1)*x(node,1)+A(node,1,2)*x(node,2)+A(node,1,3)*x(node,3)) + 1.0;
    if (r == 2) return -(A(node,2,1)*x(node,1)+A(node,2,2)*x(node,2)+A(node,2,3)*x(node,3)) + 2.0*t;
    if (r == 3) return -(A(node,3,1)*x(node,1)+A(node,3,2)*x(node,2)+A(node,3,3)*x(node,3)) + 2.0*t+1.0;
#endif
#ifdef EXAMPLE_2
    if (r == 1) return -(A(node,1,1)*x(node,1)+A(node,1,2)*x(node,2)+A(t,1,3)*x(node,3)) + 6.0*M_PI*cos(6.0*M_PI*t);
    if (r == 2) return -(A(node,2,1)*x(node,1)+A(node,2,2)*x(node,2)+A(t,2,3)*x(node,3)) - 8.0*M_PI*sin(8.0*M_PI*t);
    if (r == 3) return -(A(node,3,1)*x(node,1)+A(node,3,2)*x(node,2)+A(t,3,3)*x(node,3)) - 6.0*exp(-6.0*t);
#endif
#ifdef EXAMPLE_3
    if (r == 1) return -A(node,1,1)*x(node,1) + 2.0*t;
#endif
#ifdef EXAMPLE_4
    if (r == 1) return -A(node,1,1)*x(node,1) - 8.0*M_PI*sin(8.0*M_PI*t);
#endif

    throw std::exception();
    return NAN;
}

double LinearODE1stOrderEx1::x(const PointNodeODE &node, unsigned int r UNUSED_PARAM) const
{
    double t =  node.x;

#ifdef EXAMPLE_1
    if (r == 1) return t;
    if (r == 2) return t*t;
    if (r == 3) return t*t+t;
#endif

#ifdef EXAMPLE_2
    if (r == 1) return sin(6.0*M_PI*t);
    if (r == 2) return cos(8.0*M_PI*t);
    if (r == 3) return exp(-6.0*t);
#endif
#ifdef EXAMPLE_3
    if (r == 1) return t*t;
#endif
#ifdef EXAMPLE_4
    if (r == 1) return cos(8.0*M_PI*t);
#endif

    throw std::exception();
    return NAN;
}
