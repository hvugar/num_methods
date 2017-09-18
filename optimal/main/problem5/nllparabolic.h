#ifndef NONLINEARLOADEDPARABOLIC_H
#define NONLINEARLOADEDPARABOLIC_H

#include <grid/pibvp.h>
#include <nonlinearequation.h>
#include <utils/random.h>

/**
 * @brief The NLLIParabolicIBVP class
 * u_t = a^2 u_xx + f + \sum_{i=1}^L {g_s(x,t) u^2(x,t)}
 */
class NLLIParabolicIBVP : public IParabolicIBVP, NonLinearEquation
{
public:
    static void Main(int argc, char *argv[]);

    NLLIParabolicIBVP();
    virtual ~NLLIParabolicIBVP();

protected:
    double g(const SpaceNodePDE &sn, const TimeNodePDE &tn, unsigned int s) const;

    void solveEquation(DoubleMatrix &u, double a);

    void solveEquationM1(DoubleMatrix &u, double a);
    void solveEquationM2(DoubleMatrix &u, double a);
    void solveEquationM4(DoubleMatrix &u, double a);

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

protected:
    virtual double fx(const DoubleVector &x, unsigned int num) const;

public:
    double U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    UniformPDEGrid grid;
    DoubleMatrix *pu;
    int cur_m;
    DoubleVector ts;
};

#endif // NONLINEARLOADEDPARABOLIC_H
