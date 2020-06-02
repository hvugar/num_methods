#include "first_order_linear_ode.h"

#define EXAMPLE_31

void FirstOrderLinearODEIVP::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    CauchyProblemExample();
}

void FirstOrderLinearODEIVP::CauchyProblemExample()
{
    FirstOrderLinearODEIVP fnl;
    fnl.mx.clear();
    fnl.mx.resize(fnl.dimension().size());

    std::vector<DoubleVector> x1;
    fnl.solveInitialValueProblem(x1, ODESolverMethod::EULER);
    IPrinter::print(x1, fnl.count());
    IPrinter::printSeperatorLine();

    fnl.solveInitialValueProblem(ODESolverMethod::EULER);
    IPrinter::print(fnl.mx, fnl.count());
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x2;
    fnl.solveInitialValueProblem(x2, ODESolverMethod::RUNGE_KUTTA_2);
    IPrinter::print(x2, fnl.count());
    IPrinter::printSeperatorLine();

    fnl.solveInitialValueProblem(ODESolverMethod::RUNGE_KUTTA_2);
    IPrinter::print(fnl.mx, fnl.count());
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x3;
    fnl.solveInitialValueProblem(x3, ODESolverMethod::RUNGE_KUTTA_4);
    IPrinter::print(x3, fnl.count());
    IPrinter::printSeperatorLine();

    fnl.solveInitialValueProblem(ODESolverMethod::RUNGE_KUTTA_4);
    IPrinter::print(fnl.mx, fnl.count());
    IPrinter::printSeperatorLine();

    DoubleVector x0(fnl.count()); PointNodeODE n0;
    DoubleVector x4(fnl.count()); PointNodeODE n4;

    fnl.start(x0, n0);
    fnl.mx[0] = x0;

    for (size_t i=1; i<fnl.dimension().size(); i++)
    {
        fnl.next(x0, n0, x4, n4);
        fnl.mx[i] = x4;
        x0 = x4;
        n0 = n4;
    }
    IPrinter::print(fnl.mx, fnl.count());
    IPrinter::printSeperatorLine();
}

void FirstOrderLinearODEIVP::NonLocalConditionExample()
{
    FirstOrderLinearODEIVP nl;

    const size_t M = nl.count();
    const size_t N = TIME_MAX;
    const double h = TIME_STEP;
    //const size_t order = 2;

    double p1 = N*h;
    int p2 = N/100;
    std::vector<NonLocalCondition> C;
    C.push_back(NonLocalCondition(0, PointNodeODE(0.00*p1,    0*p2), DoubleMatrix(M,M,+1.8)));
    //    C.push_back(NonLocalCondition(1, PointNodeODE(0.25*p1,   25*p2), DoubleMatrix(M,M,+1.5)));
    //    C.push_back(NonLocalCondition(2, PointNodeODE(0.50*p1,   50*p2), DoubleMatrix(M,M,-2.5)));
    //    C.push_back(NonLocalCondition(3, PointNodeODE(0.75*p1,   75*p2), DoubleMatrix(M,M,+4.3)));
    C.push_back(NonLocalCondition(1, PointNodeODE(1.00*p1,  100*p2), DoubleMatrix(M,M,+3.4)));
    DoubleVector d(M, 0.0);

    for (size_t s=0; s<C.size(); s++)
    {
        Random::fillMatrix(C[s].m, 0, 10, 2);

        DoubleVector x; for (size_t m=1; m<=M; m++) x << nl.x(C[s].n, m);
        d += C[s].m*x;
    }
    IPrinter::printVector(d, nullptr, d.length());
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x;
    //nl.setDimension(Dimension(h, 0, static_cast<int>(N)));

    for (size_t m=0; m<M; m++)
    {
        for (size_t n=0; n<=N; n++)
        {
            //if (n%(N/100)==0) printf("%14.10f ", nl.x(n*h, m+1));
            if (n%(N/10)==0) printf("%14.8f ", nl.x(n*h, m+1));
            //printf("%14.10f ", nl.x(n*h, m+1));
        }
        puts("");
    }
    IPrinter::printSeperatorLine();

    printf("***** %14.8f %14.8f %14.8f %14.8f %14.8f\n", nl.x((N-4)*h, 1), nl.x((N-3)*h, 1), nl.x((N-2)*h, 1), nl.x((N-1)*h, 1), nl.x((N)*h, 1));
    printf("***** %14.8f %14.8f %14.8f %14.8f %14.8f\n", nl.x((0)*h, 1),   nl.x((1)*h, 1),   nl.x((2)*h, 1),   nl.x((3)*h, 1),   nl.x((4)*h, 1));

    //    puts("===== transferOfCondition =====");
    //    x.clear();
    //    nl.transferOfCondition(C, d, x, 4);
    //    for (size_t m=0; m<M; m++) { for (size_t n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.8f ", x[n][m]); printf("\n"); /*nl.printNorms(x);*/ }
    //    IPrinter::printSeperatorLine();

    //    puts("===== transferOfConditionN 2 0 =====");
    //    x.clear();
    //    nl.transferOfConditionN(C, d, x, 2, 0);
    //    for (size_t m=0; m<M; m++) { for (size_t n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.8f ", x[n][m]); printf("\n"); /*nl.printNorms(x);*/ }
    //    IPrinter::printSeperatorLine();

    const double epsilon = 0.0000001;
    printf(">>>> %20.14f %20.14f %20.14f\n", nl.x(0.50-epsilon, 1), nl.x(0.50, 1), nl.x(0.50+epsilon, 1));
    printf(">>>> %20.14f %20.14f %20.14f\n", nl.dt(0.50-epsilon, 1), nl.dt(0.50, 1), nl.dt(0.50+epsilon, 1));
    printf(">>>> %20.14f %20.14f %20.14f\n", nl.d2t(0.50-epsilon, 1), nl.d2t(0.50, 1), nl.d2t(0.50+epsilon, 1));
    printf(">>>> %20.14f %20.14f %20.14f\n", nl.d3t(0.50-epsilon, 1), nl.d3t(0.50, 1), nl.d3t(0.50+epsilon, 1));
    printf(">>>> %20.14f %20.14f %20.14f\n", nl.d4t(0.50-epsilon, 1), nl.d4t(0.50, 1), nl.d4t(0.50+epsilon, 1));
    //    printf(">>>> %20.14f %20.14f %20.14f\n", nl.d5t(0.50-epsilon, 1), nl.d5t(0.50, 1), nl.d5t(0.50+epsilon, 1));

    puts("===== transferOfConditionN 4 0 =====");
    x.clear();
    nl.transferOfConditionP(C, d, x, 4, 1);
    for (size_t m=0; m<M; m++) { for (size_t n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.8f ", x[n][m]); printf("\n"); /*nl.printNorms(x);*/ }
    IPrinter::printSeperatorLine();

    //    puts("===== transferOfConditionN 4 1 =====");
    //    x.clear();
    //    nl.transferOfConditionN(C, d, x, 4, 1);
    //    for (size_t m=0; m<M; m++) { for (size_t n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.8f ", x[n][m]); printf("\n"); /*nl.printNorms(x);*/ }
    //    IPrinter::printSeperatorLine();

    //    puts("===== transferOfConditionM4 =====");
    //    x.clear();
    //    nl.transferOfConditionM(C, d, x, 4);
    //    for (size_t m=0; m<M; m++) { for (size_t n=0; n<=N; n++) if (n%(N/10)==0) printf("%14.8f ", x[n][m]); printf("\n"); /*nl.printNorms(x);*/ }
    //    IPrinter::printSeperatorLine();
}

/*****************************************************************************************************/

void FirstOrderLinearODEFVP::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    CauchyProblemExample();
}

void FirstOrderLinearODEFVP::CauchyProblemExample()
{
    FirstOrderLinearODEFVP bnl;
    bnl.mx.clear();
    bnl.mx.resize(bnl.dimension().size());

    std::vector<DoubleVector> x1;
    bnl.solveFinalValueProblem(x1, ODESolverMethod::EULER);
    IPrinter::print(x1, bnl.count());
    IPrinter::printSeperatorLine();

    bnl.solveFinalValueProblem(ODESolverMethod::EULER);
    IPrinter::print(bnl.mx, bnl.count());
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x2;
    bnl.solveFinalValueProblem(x2, ODESolverMethod::RUNGE_KUTTA_2);
    IPrinter::print(x2, bnl.count());
    IPrinter::printSeperatorLine();

    bnl.solveFinalValueProblem(ODESolverMethod::RUNGE_KUTTA_2);
    IPrinter::print(bnl.mx, bnl.count());
    IPrinter::printSeperatorLine();

    std::vector<DoubleVector> x3;
    bnl.solveFinalValueProblem(x3, ODESolverMethod::RUNGE_KUTTA_4);
    IPrinter::print(x3, bnl.count());
    IPrinter::printSeperatorLine();

    bnl.solveFinalValueProblem(ODESolverMethod::RUNGE_KUTTA_4);
    IPrinter::print(bnl.mx, bnl.count());
    IPrinter::printSeperatorLine();
}

/*****************************************************************************************************/

double FirstOrderLinearODEIVP::A(const PointNodeODE &node, size_t r, size_t c) const
{
    C_UNUSED(node);
    C_UNUSED(r);
    C_UNUSED(c);

#if defined(EXAMPLE_11) || defined (EXAMPLE_13) || defined (EXAMPLE_14) || defined (EXAMPLE_15) || defined (EXAMPLE_16) || defined (EXAMPLE_31)
    return -1.0;
    //    return 10.0*sin(4.0*M_PI*node.x*node.x);
#endif
#ifdef EXAMPLE_12
    return -1.0;
#endif
#ifdef EXAMPLE_32
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

double FirstOrderLinearODEIVP::B(const PointNodeODE &node, size_t r) const
{
    const size_t M = count();
    double result = dt(node, r);
    for (size_t c=1; c<=M; c++) result -= A(node,r,c)*x(node,c);
    return result;
}

size_t FirstOrderLinearODEIVP::count() const
{
#if defined(EXAMPLE_11) || defined(EXAMPLE_12) || defined(EXAMPLE_13) || defined(EXAMPLE_14) || defined(EXAMPLE_15) || defined(EXAMPLE_16)
    return 1;
#endif
#if defined(EXAMPLE_31) || defined(EXAMPLE_32)
    return 3;
#endif
#if defined(EXAMPLE_5)
    return 4;
#endif
#if defined(EXAMPLE_6) || defined(EXAMPLE_7) || defined(EXAMPLE_8) || defined(EXAMPLE_9) || defined(EXAMPLE_10)
    return 1;
#endif

    throw std::exception();
}

auto FirstOrderLinearODEIVP::dimension() const -> Dimension { return Dimension(TIME_STEP, 0, TIME_MAX); }

auto FirstOrderLinearODEFVP::dimension() const -> Dimension { return Dimension(TIME_STEP, 0, TIME_MAX); }

auto FirstOrderLinearODEIVP::initial(InitialCondition, size_t r) const -> double
{
    PointNodeODE node; node.i = 0; node.x = 0.0;
    return x(node, r);
}

void FirstOrderLinearODEIVP::iterationInfo(double y, const PointNodeODE &node) const
{
    //const_cast<FirstOrderLinearODEEx1*>(this)->mx[node.i] = v;
}

void FirstOrderLinearODEIVP::iterationInfo(const DoubleVector &v, const PointNodeODE &node) const
{
    const_cast<FirstOrderLinearODEIVP*>(this)->mx[node.i] = v;
}

/*****************************************************************************************************/

double FirstOrderLinearODEFVP::A(const PointNodeODE &node, size_t r, size_t c) const
{
    C_UNUSED(node);
    C_UNUSED(r);
    C_UNUSED(c);

#if defined(EXAMPLE_11) || defined (EXAMPLE_13) || defined (EXAMPLE_14) || defined (EXAMPLE_15) || defined (EXAMPLE_16) || defined (EXAMPLE_31)
    return +1.0;
    //    return 10.0*sin(4.0*M_PI*node.x*node.x);
#endif
#ifdef EXAMPLE_12
    return -1.0;
#endif
#ifdef EXAMPLE_32
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

double FirstOrderLinearODEFVP::B(const PointNodeODE &node, size_t r) const
{
    const size_t M = count();
    double result = dt(node, r);
    for (size_t c=1; c<=M; c++) result -= A(node,r,c)*x(node,c);
    return result;
}

size_t FirstOrderLinearODEFVP::count() const
{
#if defined(EXAMPLE_11) || defined(EXAMPLE_12) || defined(EXAMPLE_13) || defined(EXAMPLE_14) || defined(EXAMPLE_15) || defined(EXAMPLE_16)
    return 1;
#endif
#if defined(EXAMPLE_31) || defined(EXAMPLE_32)
    return 3;
#endif
#if defined(EXAMPLE_5)
    return 4;
#endif
#if defined(EXAMPLE_6) || defined(EXAMPLE_7) || defined(EXAMPLE_8) || defined(EXAMPLE_9) || defined(EXAMPLE_10)
    return 1;
#endif

    throw std::exception();
}

auto FirstOrderLinearODEFVP::final(FinalCondition, size_t r) const -> double
{
    PointNodeODE node; node.i = dimension().max(); node.x = dimension().max()*dimension().step();
    return x(node, r);
}

auto FirstOrderLinearODEFVP::iterationInfo(double y, const PointNodeODE &node) const -> void
{
    //const_cast<FirstOrderLinearODEEx1*>(this)->mx[node.i] = v;
}

auto FirstOrderLinearODEFVP::iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void
{
    const_cast<FirstOrderLinearODEFVP*>(this)->mx[node.i] = v;
}

/*****************************************************************************************************/

double FirstOrderLinearSample::x(const PointNodeODE &node, size_t r UNUSED_PARAM) const
{
    double t =  node.x;

#ifdef EXAMPLE_11
    if (r == 1) return t*t;
#endif
#ifdef EXAMPLE_12
    if (r == 1) return cos(8.0*M_PI*t);
#endif
#ifdef EXAMPLE_13
    if (r == 1) return pow(fabs(t-0.5), 5);
#endif
#ifdef EXAMPLE_14
    if (r == 1)
    {
        if (node.x <= 1.0) return -0.25*t*t*t*t + t*t*t - 0.5*t*t - t - 1.0;
        if (node.x  > 1.0) return t*t - 2.0*t - 0.75;
    }
#endif
#ifdef EXAMPLE_15
    if (r == 1)
    {
        if (node.x <= 0.5) return t*t*t*t*t - 2.5*t*t*t*t + 5.0*t*t*t + + 5.0*t*t - 6.0*t + 2.0;
        if (node.x  > 0.5) return 2.5*t*t*t + 6.25*t*t - 6.3125*t + 2.03125;
    }
#endif
#ifdef EXAMPLE_16
    if (r == 1)
    {
        if (node.x <= 0.5) return t*t*t*t + 3.0*t*t*t - 5.5*t*t + 0.5;
        if (node.x  > 0.5) return 2.0*t*t*t*t + t*t*t - 4.0*t*t - 0.5*t + 0.5625;
    }
#endif
#ifdef EXAMPLE_31
    if (r == 1) return t;
    if (r == 2) return t*t;
    if (r == 3) return t*t+t;
#endif
#ifdef EXAMPLE_32
    if (r == 1) return sin(6.0*M_PI*t);
    if (r == 2) return cos(8.0*M_PI*t);
    if (r == 3) return exp(-6.0*t);
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

    throw std::runtime_error("double FirstOrderLinearODEEx1::x");
}

double FirstOrderLinearSample::dt(const PointNodeODE &node, size_t r UNUSED_PARAM) const
{
    double t =  node.x;
#ifdef EXAMPLE_11
    if (r == 1) return 2.0*t;
#endif
#ifdef EXAMPLE_12
    if (r == 1) return -8.0*M_PI*sin(8.0*M_PI*t);
#endif
#ifdef EXAMPLE_13
    if (r == 1) return 5.0*fabs(t-0.5)*pow((t-0.5), 3.0);
#endif
#ifdef EXAMPLE_14
    if (r == 1)
    {
        if (node.x <= 1.0) return -t*t*t + 3.0*t*t - t - 1.0;
        if (node.x  > 1.0) return 2.0*t - 2.0;
    }
#endif
#ifdef EXAMPLE_15
    if (r == 1)
    {
        if (node.x <= 0.5) return 5.0*t*t*t*t - 10.0*t*t*t + 15.0*t*t + + 10.0*t - 6.0;
        if (node.x  > 0.5) return 7.5*t*t + 12.5*t - 6.3125;
    }
#endif
#ifdef EXAMPLE_16
    if (r == 1)
    {
        if (node.x <= 0.5) return 4.0*t*t*t + 9.0*t*t - 11.0*t;
        if (node.x  > 0.5) return 8.0*t*t*t + 3.0*t*t - 8.0*t - 0.5;
    }
#endif
#ifdef EXAMPLE_31
    if (r == 1) return 1.0;
    if (r == 2) return 2.0*t;
    if (r == 3) return 2.0*t+1.0;
#endif
#ifdef EXAMPLE_32
    if (r == 1) return +6.0*M_PI*cos(6.0*M_PI*t);
    if (r == 2) return -8.0*M_PI*sin(8.0*M_PI*t);
    if (r == 3) return -6.0*exp(-6.0*t);
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

double FirstOrderLinearSample::d2t(const PointNodeODE &node, size_t r UNUSED_PARAM) const
{
    double t =  node.x;
#ifdef EXAMPLE_11
    if (r == 1) return 2.0*t;
#endif
#ifdef EXAMPLE_12
    if (r == 1) return -8.0*M_PI*sin(8.0*M_PI*t);
#endif
#ifdef EXAMPLE_13
    if (r == 1) return 5.0*fabs(t-0.5)*pow((t-0.5), 3.0);
#endif
#ifdef EXAMPLE_14
    if (r == 1)
    {
        if (node.x <= 1.0) return -t*t*t + 3.0*t*t - t - 1.0;
        if (node.x  > 1.0) return 2.0*t - 2.0;
    }
#endif
#ifdef EXAMPLE_15
    if (r == 1)
    {
        if (node.x <= 0.5) return 20.0*t*t*t - 30.0*t*t + 30.0*t + 10.0;
        if (node.x  > 0.5) return 15.0*t + 12.5;
    }
#endif
#ifdef EXAMPLE_16
    if (r == 1)
    {
        if (node.x <= 0.5) return 12.0*t*t + 18.0*t - 11.0;
        if (node.x  > 0.5) return 24.0*t*t + 6.0*t - 8.0;
    }
#endif
#ifdef EXAMPLE_31
    if (r == 1) return 1.0;
    if (r == 2) return 2.0*t;
    if (r == 3) return 2.0*t+1.0;
#endif
#ifdef EXAMPLE_32
    if (r == 1) return +6.0*M_PI*cos(6.0*M_PI*t);
    if (r == 2) return -8.0*M_PI*sin(8.0*M_PI*t);
    if (r == 3) return -6.0*exp(-6.0*t);
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

    throw std::runtime_error("FirstOrderLinearODEEx1::d2t");
}

double FirstOrderLinearSample::d3t(const PointNodeODE &node, size_t r UNUSED_PARAM) const
{
    double t =  node.x;
#ifdef EXAMPLE_11
    if (r == 1) return 2.0*t;
#endif
#ifdef EXAMPLE_12
    if (r == 1) return -8.0*M_PI*sin(8.0*M_PI*t);
#endif
#ifdef EXAMPLE_13
    if (r == 1) return 5.0*fabs(t-0.5)*pow((t-0.5), 3.0);
#endif
#ifdef EXAMPLE_14
    if (r == 1)
    {
        if (node.x <= 1.0) return -t*t*t + 3.0*t*t - t - 1.0;
        if (node.x  > 1.0) return 2.0*t - 2.0;
    }
#endif
#ifdef EXAMPLE_15
    if (r == 1)
    {
        if (node.x <= 0.5) return 60.0*t*t - 60.0*t + 30.0;
        if (node.x  > 0.5) return 15.0;
    }
#endif
#ifdef EXAMPLE_16
    if (r == 1)
    {
        if (node.x <= 0.5) return 24.0*t + 18.0;
        if (node.x  > 0.5) return 48.0*t + 6.0;
    }
#endif
#ifdef EXAMPLE_31
    if (r == 1) return 1.0;
    if (r == 2) return 2.0*t;
    if (r == 3) return 2.0*t+1.0;
#endif
#ifdef EXAMPLE_32
    if (r == 1) return +6.0*M_PI*cos(6.0*M_PI*t);
    if (r == 2) return -8.0*M_PI*sin(8.0*M_PI*t);
    if (r == 3) return -6.0*exp(-6.0*t);
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

    throw std::runtime_error("FirstOrderLinearODEEx1::d3t");
}

double FirstOrderLinearSample::d4t(const PointNodeODE &node, size_t r UNUSED_PARAM) const
{
    double t =  node.x;
#ifdef EXAMPLE_11
    if (r == 1) return 2.0*t;
#endif
#ifdef EXAMPLE_12
    if (r == 1) return -8.0*M_PI*sin(8.0*M_PI*t);
#endif
#ifdef EXAMPLE_13
    if (r == 1) return 5.0*fabs(t-0.5)*pow((t-0.5), 3.0);
#endif
#ifdef EXAMPLE_14
    if (r == 1)
    {
        if (node.x <= 1.0) return -t*t*t + 3.0*t*t - t - 1.0;
        if (node.x  > 1.0) return 2.0*t - 2.0;
    }
#endif
#ifdef EXAMPLE_15
    if (r == 1)
    {
        if (node.x <= 0.5) return 120.0*t - 60.0;
        if (node.x  > 0.5) return 0.0;
    }
#endif
#ifdef EXAMPLE_16
    if (r == 1)
    {
        if (node.x <= 0.5) return 24.0;
        if (node.x  > 0.5) return 48.0;
    }
#endif
#ifdef EXAMPLE_31
    if (r == 1) return 1.0;
    if (r == 2) return 2.0*t;
    if (r == 3) return 2.0*t+1.0;
#endif
#ifdef EXAMPLE_32
    if (r == 1) return +6.0*M_PI*cos(6.0*M_PI*t);
    if (r == 2) return -8.0*M_PI*sin(8.0*M_PI*t);
    if (r == 3) return -6.0*exp(-6.0*t);
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

    throw std::runtime_error("FirstOrderLinearODEEx1::d4t");
}

void FirstOrderLinearSample::printNorms(std::vector<DoubleVector> &_x, size_t k) const
{
    const size_t M = k;
    const size_t N = TIME_MAX;
    const double h = TIME_STEP;

    printf("Norms: ");
    for (size_t m=0; m<M; m++)
    {
        double norm = 0.0;
        for (size_t n=0; n<=N; n++)
        {
            norm += (_x[n][m]-x(n*h, m+1))*(_x[n][m]-x(n*h, m+1));
        }
        norm = sqrt(norm);
        printf("%14.10f ", norm);
    }
    puts("");
}


