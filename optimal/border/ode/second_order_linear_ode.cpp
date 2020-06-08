#include "second_order_linear_ode.h"

#define EXAMPLE_31

void SecondOrderLinearODEIBVP::Main(int argc, char **argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    CauchyProblemExample();
}

void SecondOrderLinearODEIBVP::CauchyProblemExample()
{
    SecondOrderLinearODEIBVP snl;
    snl.mx.clear();
    snl.mx.resize(snl.dimension().size());

    {
        IPrinter::printSeperatorLine("EULER");
//        std::vector<DoubleVector> x1;
//        fnl.solveInitialValueProblem(x1, ODESolverMethod::EULER);
//        IPrinter::print(x1, snl.count());

//        IPrinter::printSeperatorLine();
        snl.solveInitialValueProblem(ODESolverMethod::EULER);
        IPrinter::print(snl.mx, snl.count());
    }

    {
        IPrinter::printSeperatorLine("RUNGE_KUTTA_2");
//        std::vector<DoubleVector> x2;
//        fnl.solveInitialValueProblem(x2, ODESolverMethod::RUNGE_KUTTA_2);
//        IPrinter::print(x2, snl.count());

//        IPrinter::printSeperatorLine();
        snl.solveInitialValueProblem(ODESolverMethod::RUNGE_KUTTA_2);
        IPrinter::print(snl.mx, snl.count());
    }

    {
        IPrinter::printSeperatorLine("RUNGE_KUTTA_4");
//        std::vector<DoubleVector> x3;
//        fnl.solveInitialValueProblem(x3, ODESolverMethod::RUNGE_KUTTA_4);
//        IPrinter::print(x3, snl.count());

//        IPrinter::printSeperatorLine();
        snl.solveInitialValueProblem(ODESolverMethod::RUNGE_KUTTA_4);
        IPrinter::print(snl.mx, snl.count());
    }

    {
        IPrinter::printSeperatorLine("STEP");
        DoubleVector x0(snl.count()); PointNodeODE n0;
        DoubleVector x4(snl.count()); PointNodeODE n4;

        snl.start(x0, n0);
        snl.mx[0] = x0;

        for (size_t i=1, sz=snl.dimension().size(); i<sz; i++)
        {
            snl.next(x0, n0, x4, n4, ODESolverMethod::RUNGE_KUTTA_2);
            snl.mx[i] = x4;
            x0 = x4;
            n0 = n4;
        }
        IPrinter::print(snl.mx, snl.count());
        IPrinter::printSeperatorLine();
    }
}

double SecondOrderLinearODEIBVP::A(const PointNodeODE & node, size_t r, size_t c) const
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

double SecondOrderLinearODEIBVP::B(const PointNodeODE &node, size_t r, size_t c) const
{
    C_UNUSED(node);
    C_UNUSED(r);
    C_UNUSED(c);

#if defined(EXAMPLE_11) || defined (EXAMPLE_13) || defined (EXAMPLE_14) || defined (EXAMPLE_15) || defined (EXAMPLE_16) || defined (EXAMPLE_31)
    return +0.5;
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

double SecondOrderLinearODEIBVP::C(const PointNodeODE &node, size_t r) const
{
    const size_t m = count();
    double result = d2t(node, r);
    for (size_t c=1; c<=m; c++) result -= (A(node,r,c)*dt(node,c) + B(node,r,c)*x(node,c));
    return result;
}

size_t SecondOrderLinearODEIBVP::count() const
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

auto SecondOrderLinearODEIBVP::dimension() const -> Dimension { return Dimension(TIME_STEP_2ND, TIME_MIN_2ND, TIME_MAX_2ND); }

auto SecondOrderLinearODEIBVP::initial(InitialCondition c, size_t r) const -> double
{
    PointNodeODE node; node.i = TIME_MIN_2ND; node.x = node.i*TIME_STEP_2ND;
    return (c == InitialCondition::InitialValue) ? x(node, r) : dt(node, r);
}

void SecondOrderLinearODEIBVP::iterationInfo(const DoubleVector &v, const PointNodeODE &node) const
{
    const_cast<SecondOrderLinearODEIBVP*>(this)->mx.at(node.i) = v;

   //if (node.i%((dimension().size()-1)/10)==0) { printf("%6d: ", node.i); IPrinter::print(v, v.length()); }
}

//auto SecondOrderLinearODEIBVP::boundary(const PointNodeODE &, BoundaryConditionPDE &, size_t) const -> double
//{
//    throw runtime_error("");
//}

/****************************************************************************************************************************/

void SecondOrderLinearODEFBVP::Main(int argc, char **argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    CauchyProblemExample();
}

void SecondOrderLinearODEFBVP::CauchyProblemExample()
{
    SecondOrderLinearODEFBVP snl;
    snl.mx.clear();
    snl.mx.resize(snl.dimension().size());

    {
        IPrinter::printSeperatorLine("EULER");
//        std::vector<DoubleVector> x1;
//        fnl.solveInitialValueProblem(x1, ODESolverMethod::EULER);
//        IPrinter::print(x1, snl.count());

        snl.mx.clear();
        snl.mx.resize(snl.dimension().size());
        IPrinter::printSeperatorLine();
        snl.solveInitialValueProblem(ODESolverMethod::EULER);
        IPrinter::print(snl.mx, snl.count());
    }

    {
        IPrinter::printSeperatorLine("RUNGE_KUTTA_2");
//        std::vector<DoubleVector> x2;
//        fnl.solveInitialValueProblem(x2, ODESolverMethod::RUNGE_KUTTA_2);
//        IPrinter::print(x2, snl.count());

        snl.mx.clear();
        snl.mx.resize(snl.dimension().size());
        IPrinter::printSeperatorLine();
        snl.solveInitialValueProblem(ODESolverMethod::RUNGE_KUTTA_2);
        IPrinter::print(snl.mx, snl.count());
    }

    {
        IPrinter::printSeperatorLine("RUNGE_KUTTA_4");
//        std::vector<DoubleVector> x3;
//        fnl.solveInitialValueProblem(x3, ODESolverMethod::RUNGE_KUTTA_4);
//        IPrinter::print(x3, snl.count());

        snl.mx.clear();
        snl.mx.resize(snl.dimension().size());
        IPrinter::printSeperatorLine();
        snl.solveInitialValueProblem(ODESolverMethod::RUNGE_KUTTA_4);
        IPrinter::print(snl.mx, snl.count());
    }

//    {
//        IPrinter::printSeperatorLine("STEP");
//        DoubleVector x0(snl.count()); PointNodeODE n0;
//        DoubleVector x4(snl.count()); PointNodeODE n4;

//        snl.start(x0, n0);
//        snl.mx[0] = x0;

//        for (size_t i=1, sz=snl.dimension().size(); i<sz; i++)
//        {
//            snl.next(x0, n0, x4, n4);
//            snl.mx[i] = x4;
//            x0 = x4;
//            n0 = n4;
//        }
//        IPrinter::print(snl.mx, snl.count());
//        IPrinter::printSeperatorLine();
//    }
}

double SecondOrderLinearODEFBVP::A(const PointNodeODE & node, size_t r, size_t c) const
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

double SecondOrderLinearODEFBVP::B(const PointNodeODE &node, size_t r, size_t c) const
{
    C_UNUSED(node);
    C_UNUSED(r);
    C_UNUSED(c);

#if defined(EXAMPLE_11) || defined (EXAMPLE_13) || defined (EXAMPLE_14) || defined (EXAMPLE_15) || defined (EXAMPLE_16) || defined (EXAMPLE_31)
    return +0.5;
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

double SecondOrderLinearODEFBVP::C(const PointNodeODE &node, size_t r) const
{
    const size_t m = count();
    double result = d2t(node, r);
    for (size_t c=1; c<=m; c++) result -= (A(node,r,c)*dt(node,c) + B(node,r,c)*x(node,c));
    return result;}

size_t SecondOrderLinearODEFBVP::count() const
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

auto SecondOrderLinearODEFBVP::dimension() const -> Dimension { return Dimension(TIME_STEP_2ND, TIME_MIN_2ND, TIME_MAX_2ND); }

auto SecondOrderLinearODEFBVP::final(FinalCondition c, size_t r) const -> double
{
    PointNodeODE node; node.i = TIME_MAX_2ND; node.x = node.i*TIME_STEP_2ND;
    return (c == FinalCondition::FinalValue) ? x(node, r) : dt(node, r);
}

void SecondOrderLinearODEFBVP::iterationInfo(const DoubleVector &v, const PointNodeODE &node) const
{
    const_cast<SecondOrderLinearODEFBVP*>(this)->mx[node.i] = v;
}

//auto SecondOrderLinearODEFBVP::boundary(const PointNodeODE &, BoundaryConditionPDE &, size_t) const -> double
//{
//    throw runtime_error("");
//}

/*****************************************************************************************************/

double SecondOrderLinearSample::x(const PointNodeODE &node, size_t r UNUSED_PARAM) const
{
    double t =  node.x;

#ifdef EXAMPLE_11
    if (r == 1) return t*t + 1.0;
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
    if (r == 1) return t+0.2;
    if (r == 2) return t*t-1.0;
    if (r == 3) return t*t+t-0.2;
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

double SecondOrderLinearSample::dt(const PointNodeODE &node, size_t r UNUSED_PARAM) const
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

double SecondOrderLinearSample::d2t(const PointNodeODE &node, size_t r UNUSED_PARAM) const
{
    double t =  node.x;
#ifdef EXAMPLE_11
    if (r == 1) return 2.0;
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
    if (r == 1) return 0.0;
    if (r == 2) return 2.0;
    if (r == 3) return 2.0;
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

double SecondOrderLinearSample::d3t(const PointNodeODE &node, size_t r UNUSED_PARAM) const
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

double SecondOrderLinearSample::d4t(const PointNodeODE &node, size_t r UNUSED_PARAM) const
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

void SecondOrderLinearSample::printNorms(std::vector<DoubleVector> &_x, size_t k) const
{
    const size_t M = k;
    const size_t N = TIME_MAX_2ND;
    const double h = TIME_STEP_2ND;

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


