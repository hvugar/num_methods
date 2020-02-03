#include "first_order_linear_ode.h"

#define EXAMPLE_7

#define _N 1000
#define _H 0.01

void FirstOrderLinearODEEx1::Main(int argc UNUSED_PARAM, char **argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    //srand(time(NULL));

    FirstOrderLinearODEEx1 nl;

    const unsigned int M = nl.count();
    const unsigned int N = _N;
    const double h = _H;
    //const unsigned int order = 2;

    std::vector<NonLocalCondition> C;
    C.push_back(NonLocalCondition(0, PointNodeODE(0.00,   0), DoubleMatrix(M,M,+1.8)));
    C.push_back(NonLocalCondition(1, PointNodeODE(2.5,  250), DoubleMatrix(M,M,+1.5)));
    C.push_back(NonLocalCondition(2, PointNodeODE(5.0,  500), DoubleMatrix(M,M,-2.5)));
    C.push_back(NonLocalCondition(3, PointNodeODE(7.5,  750), DoubleMatrix(M,M,+4.3)));
    C.push_back(NonLocalCondition(4, PointNodeODE(10.0, 1000), DoubleMatrix(M,M,+3.4)));
    DoubleVector d(M, 0.0);

    for (unsigned int s=0; s<C.size(); s++)
    {
        DoubleVector x; for (unsigned int m=1; m<=M; m++) x << nl.x(C[s].n, m);
        d += C[s].m*x;
    }
    IPrinter::printVector(d, nullptr, d.length());
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x;
    nl.setDimension(Dimension(h, 0, static_cast<int>(N)));

    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            //if (n%(N/100)==0) printf("%14.10f ", nl.x(n*h, m+1));
            if (n%(N/10)==0) printf("%14.6f ", nl.x(n*h, m+1));
            //printf("%14.10f ", nl.x(n*h, m+1));
        }
        puts("");
    }
    IPrinter::printSeperatorLine();

    //nl.transferOfCondition2(C, d, x, 2);
    //for (unsigned int m=0; m<M; m++) for (unsigned int n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.6f ", x[n][m]); nl.printNorms(x);
    //IPrinter::printSeperatorLine();

    puts("===== transferOfCondition =====");
    nl.transferOfCondition3(C, d, x, 4);
    for (unsigned int m=0; m<M; m++) for (unsigned int n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.6f ", x[n][m]); nl.printNorms(x);
    IPrinter::printSeperatorLine();

//    nl.transferOfCondition(C, d, x, 6);
    //for (unsigned int m=0; m<M; m++) for (unsigned int n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.6f ", x[n][m]); nl.printNorms(x);
    //IPrinter::printSeperatorLine();
}

FirstOrderLinearODEEx1::FirstOrderLinearODEEx1()
{
}

FirstOrderLinearODEEx1::~FirstOrderLinearODEEx1()
{}

double FirstOrderLinearODEEx1::A(const PointNodeODE &node, unsigned int r, unsigned int c) const
{
    C_UNUSED(node);
    C_UNUSED(r);
    C_UNUSED(c);

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
#ifdef EXAMPLE_5
    if (r==1 && c == 1) return +2.0;
    if (r==1 && c == 2) return +1.0;//node.x;
    if (r==1 && c == 3) return +3.0;
    if (r==1 && c == 4) return +5.0;

    if (r==2 && c == 1) return +4.0;//*node.x*node.x;
    if (r==2 && c == 2) return +1.0;
    if (r==2 && c == 3) return -2.0;
    if (r==2 && c == 4) return +8.0;

    if (r==3 && c == 1) return +6.0;
    if (r==3 && c == 2) return -1.0;//node.x;
    if (r==3 && c == 3) return +5.0;
    if (r==3 && c == 4) return +3.0;

    if (r==4 && c == 1) return +4.0;
    if (r==4 && c == 2) return -2.0;//node.x;
    if (r==4 && c == 3) return +1.0;
    if (r==4 && c == 4) return +1.0;
#endif

#if defined(EXAMPLE_6) ||defined(EXAMPLE_7) || defined(EXAMPLE_8) || defined(EXAMPLE_9) || defined(EXAMPLE_10)
    return -1.0;
#endif

    throw std::exception();
}

double FirstOrderLinearODEEx1::B(const PointNodeODE &node, unsigned int r) const
{
    //const unsigned int _count = count();
    //double result = dt(node, r);
    //for (unsigned int c=1; c<=_count; c++) result -= A(node,r,c)*x(node,c);
    //return result;

    return dt(node, r) - A(node,r,1)*x(node,1);
}

unsigned int FirstOrderLinearODEEx1::count() const
{
#if defined(EXAMPLE_1) || defined(EXAMPLE_2)
    return 3;
#endif
#if defined(EXAMPLE_3) || defined(EXAMPLE_4)
    return 1;
#endif
#if defined(EXAMPLE_5)
    return 4;
#endif
#if defined(EXAMPLE_6) || defined(EXAMPLE_7) || defined(EXAMPLE_8) || defined(EXAMPLE_9) || defined(EXAMPLE_10)
    return 1;
#endif
}

double FirstOrderLinearODEEx1::x(const PointNodeODE &node, unsigned int r UNUSED_PARAM) const
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
#ifdef EXAMPLE_5
    if (r == 1) return 0.4*sin(12.0*t) + 0.2*t + cos(t*t) - 0.5;
    if (r == 2) return 0.6*sin(4.0*t) + exp(t) + 0.8*sin(10.0*t*t*t);
    if (r == 3) return 0.8*sin(6.0*t*t);
    if (r == 4) return 0.4*exp(t)*sin(12.0*t) + 0.4;
#endif
#ifdef EXAMPLE_6
    if (r == 1) return sin(4.0*M_PI*t*t*t*t) + 0.6*exp(t) - 1.4;
#endif
#ifdef EXAMPLE_7
    if (r == 1) return 6.0*t*t*t*t*t*t - 5.0*t*t*t + t*t + 2.0;
#endif
#ifdef EXAMPLE_8
    if (r == 1) return t*t*t*t*t*t;
#endif
#ifdef EXAMPLE_9
    if (r == 1) return t*t;
#endif
#ifdef EXAMPLE_10
    if (r == 1) return t;
#endif

    throw std::exception();
}

double FirstOrderLinearODEEx1::dt(const PointNodeODE &node, unsigned int r UNUSED_PARAM) const
{
    double t =  node.x;

#ifdef EXAMPLE_1
    if (r == 1) return 1.0;
    if (r == 2) return 2.0*t;
    if (r == 3) return 2.0*t+1.0;
#endif
#ifdef EXAMPLE_2
    if (r == 1) return +6.0*M_PI*cos(6.0*M_PI*t);
    if (r == 2) return -8.0*M_PI*sin(8.0*M_PI*t);
    if (r == 3) return -6.0*exp(-6.0*t);
#endif
#ifdef EXAMPLE_3
    if (r == 1) return 2.0*t;
#endif
#ifdef EXAMPLE_4
    if (r == 1) return -8.0*M_PI*sin(8.0*M_PI*t);
#endif
#ifdef EXAMPLE_5
    if (r == 1) return 4.8*cos(12.0*t) - 2.0*t*sin(t*t) + 0.2;
    if (r == 2) return 2.4*cos(4.0*t) + exp(t) + 24.0*t*t*cos(10.0*t*t*t);
    if (r == 3) return 9.6*t*cos(6.0*t*t);
    if (r == 4) return 0.4*exp(t)*sin(12.0*t) + 4.8*exp(t)*cos(12.0*t);
#endif
#ifdef EXAMPLE_6
    if (r == 1) return 16.0*M_PI*t*t*t*cos(4.0*M_PI*t*t*t*t) + 0.6*exp(t);
#endif
#ifdef EXAMPLE_7
    if (r == 1) return 36.0*t*t*t*t*t - 15.0*t*t + 2.0*t;
#endif
#ifdef EXAMPLE_8
    if (r == 1) return 6.0*t*t*t*t*t;
#endif
#ifdef EXAMPLE_9
    if (r == 1) return 2.0*t;
#endif
#ifdef EXAMPLE_10
    if (r == 1) return 1.0;
#endif
    throw std::exception();
}

auto FirstOrderLinearODEEx1::initial(InitialCondition, unsigned int r) const -> double
{
    PointNodeODE node; node.i = 0; node.x = 0.0;
    return x(node, r);
}

void FirstOrderLinearODEEx1::printNorms(std::vector<DoubleVector> &_x) const
{
    const unsigned int M = count();
    const unsigned int N = _N;
    const double h = _H;

    printf("Norms: ");
    for (unsigned int m=0; m<M; m++)
    {
        double norm = 0.0;
        for (unsigned int n=0; n<=N; n++)
        {
            norm += (_x[n][m]-x(n*h, m+1))*(_x[n][m]-x(n*h, m+1));
        }
        norm = sqrt(norm);
        printf("%14.10f ", norm);
    }
    puts("");
}
