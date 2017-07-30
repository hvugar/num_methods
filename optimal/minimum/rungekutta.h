#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "global.h"
#include "function.h"
#include "matrix2d.h"
#include <vector>

typedef double (*R2FunctionX)(double x, double y);

class MINIMUMSHARED_EXPORT RungeKutta
{
public:
    static void calculate(R2Function *f, double x0, double x1, double y0, double &y1, double dx);
    static unsigned int calculate(R2Function *f, double x0, double x1, double y0, DoubleVector& y, double dx);
    static void calculate(R2Function *f, double x0, double y0, DoubleVector& y, double dx);

    static void calculate(R2FunctionX f, double x0, double y0, DoubleVector& y, double dx);
};

class MINIMUMSHARED_EXPORT CauchyProblem
{
public:
    enum SolutionMethods
    {
        RungeKutta4rdOrder,
        Euler
    };

    /**
     * @brief f Right side of ordinary differensial equation first order. dy/dx=f(x,y);
     * @param x Independent variable;
     * @param y
     * @return
     */
    virtual double f(double x, const DoubleVector &y) const = 0;

    void calculate(SolutionMethods method);

    static void rungeKutta(CauchyProblem *cp, double x0, double y0, double h, unsigned int N, DoubleVector &y);
    static void rungeKutta(std::vector<CauchyProblem*> cps, double x0, double h, unsigned int N, DoubleMatrix &my);

    static void euler1(std::vector<CauchyProblem*> cps, double x0, double h, unsigned int N, DoubleMatrix &my);
    static void euler2(std::vector<CauchyProblem*> cps, double x0, double h, unsigned int N, DoubleMatrix &my);

    double x0;
    double y0;

private:
    void calculateRungeKutta4rdOrder();
};

#endif // RUNGEKUTTA_H
