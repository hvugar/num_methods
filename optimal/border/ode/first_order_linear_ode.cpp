#include "first_order_linear_ode.h"

#define EXAMPLE_6

void FirstOrderLinearODEEx1::Main(int argc UNUSED_PARAM, char **argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    //srand(time(NULL));

    FirstOrderLinearODEEx1 nl;

    const unsigned int M = nl.count();
    const unsigned int N = 100;
    const double h = 0.01;
    const unsigned int order = 2;

    std::vector<NonLocalCondition> C;

    C.push_back(NonLocalCondition(0, PointNodeODE(0.000, static_cast<int>(0.00*N)), DoubleMatrix(M,M,0.0)));
    C.push_back(NonLocalCondition(1, PointNodeODE(0.250, static_cast<int>(0.25*N)), DoubleMatrix(M,M,0.0)));
    C.push_back(NonLocalCondition(2, PointNodeODE(0.500, static_cast<int>(0.50*N)), DoubleMatrix(M,M,0.0)));
    C.push_back(NonLocalCondition(3, PointNodeODE(0.750, static_cast<int>(0.75*N)), DoubleMatrix(M,M,0.0)));
    C.push_back(NonLocalCondition(4, PointNodeODE(1.000, static_cast<int>(1.00*N)), DoubleMatrix(M,M,0.0)));

    printf("%f %f %f %f %f\n", C[0].n.x, C[1].n.x, C[2].n.x, C[3].n.x, C[4].n.x);
    printf("%d %d %d %d %d\n", C[0].n.i, C[1].n.i, C[2].n.i, C[3].n.i, C[4].n.i);

    //for (unsigned i=0; i<C.size(); i++)
    {
        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                for (unsigned i=0; i<C.size(); i++) C[i].m[r][c] = Random::value(0,1,2);
                //C[i].m[r][c] = 10.0*Random::value(0,1,1)-2.0;

            }
        }
        //IPrinter::printMatrix(C[i].m, C[i].m.rows(), C[i].m.cols());
        //IPrinter::printSeperatorLine();
    }

    for (unsigned i=0; i<C.size(); i++)
    {
        IPrinter::printMatrix(C[i].m, C[i].m.rows(), C[i].m.cols());
        IPrinter::printSeperatorLine();
    }

    DoubleVector x000; for (unsigned int m=1; m<=M; m++) x000 << nl.x(C[0].n,m);
    DoubleVector x025; for (unsigned int m=1; m<=M; m++) x025 << nl.x(C[1].n,m);
    DoubleVector x050; for (unsigned int m=1; m<=M; m++) x050 << nl.x(C[2].n,m);
    DoubleVector x075; for (unsigned int m=1; m<=M; m++) x075 << nl.x(C[3].n,m);
    DoubleVector x100; for (unsigned int m=1; m<=M; m++) x100 << nl.x(C[4].n,m);

    IPrinter::printVector(x000, nullptr, x000.length());
    IPrinter::printVector(x025, nullptr, x025.length());
    IPrinter::printVector(x050, nullptr, x050.length());
    IPrinter::printVector(x075, nullptr, x075.length());
    IPrinter::printVector(x100, nullptr, x100.length());
    IPrinter::printSeperatorLine();

    DoubleVector d = C[0].m*x000 + C[1].m*x025 + C[2].m*x050 + C[3].m*x075 + C[4].m*x100;
    IPrinter::printVector(d, nullptr, d.length());
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x;
    nl.setDimension(Dimension(h, 0, static_cast<int>(N)));

    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            //if (n%(N/100)==0) printf("%14.10f ", nl.x(n*h, m+1));
            if (n%(N/10)==0) printf("%14.10f ", nl.x(n*h, m+1));
            //printf("%14.10f ", nl.x(n*h, m+1));
        }
        puts("");
    }
    IPrinter::printSeperatorLine();

    nl.transferOfCondition(C, d, x, 2);
    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.10f ", x[n][m]);
        puts("");
    }

    nl.printNorms(x);
    IPrinter::printSeperatorLine();

    nl.transferOfCondition(C, d, x, 4);
    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.10f ", x[n][m]);
        puts("");
    }

    nl.printNorms(x);
    IPrinter::printSeperatorLine();

    nl.transferOfCondition(C, d, x, 6);
    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.10f ", x[n][m]);
        puts("");
    }

    nl.printNorms(x);
    IPrinter::printSeperatorLine();
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
#ifdef EXAMPLE_6
    return 1.0;
#endif

    throw std::exception();
}

double FirstOrderLinearODEEx1::B(const PointNodeODE &node, unsigned int r) const
{
    const unsigned int _count = count();
    double result = dt(node, r);
    for (unsigned int c=1; c<=_count; c++) result -= A(node,r,c)*x(node,c);
    return result;
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
#if defined(EXAMPLE_6)
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
    const unsigned int N = 100;
    const double h = 0.01;

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
