#ifndef ISYSTEMLINEARODENONLOCALCONTIONSV2_H
#define ISYSTEMLINEARODENONLOCALCONTIONSV2_H

#include <vector2d.h>
#include <matrix2d.h>
#include <vector>
#include <ode/cauchyp.h>

class ISystemLinearODENonLocalContionsV2
{
public:
    ISystemLinearODENonLocalContionsV2(const ODEGrid& grid);

    virtual double A(TimeNode node, unsigned int row = 0, unsigned int col = 0) const = 0;
    virtual double B(TimeNode node, unsigned int s, unsigned int row = 0, unsigned int col = 0) const = 0;
    virtual double C(TimeNode node, unsigned int row = 0) const = 0;

    void calculateForward();
    void calculateBackward();

    void addCondition(const DoubleMatrix &alpha, double time, unsigned int nmbr);
    void addLoadPoint(const DoubleMatrix &betta, double time, unsigned int nmbr);
    void setRightSize(const DoubleVector &gamma);

    void setGrid(const ODEGrid &grid);

private:
    std::vector<DoubleMatrix> alphas;
    std::vector<double> atimes;
    std::vector<unsigned int> anmbrs;

    std::vector<DoubleMatrix> bettas;
    std::vector<double> btimes;
    std::vector<unsigned int> bnmbrs;

    DoubleVector gamma;

    ODEGrid mgrid;
};

#endif // ISYSTEMLINEARODENONLOCALCONTIONSV2_H
