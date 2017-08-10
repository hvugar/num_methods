#ifndef ISYSTEMLINEARODENONLOCALCONTIONS_H
#define ISYSTEMLINEARODENONLOCALCONTIONS_H

#include <vector2d.h>
#include <matrix2d.h>
#include <grid/grid.h>
#include <grid/cauchyp.h>
#include <vector>

class MINIMUMSHARED_EXPORT ISystemLinearODENonLocalContionsM
{
public:
    enum ConditionType
    {
        SeparatedLeft = 0,
        SeparatedRight = 1,
        NonSeparated = 2
    };

    struct Condition
    {
        ConditionType type;
        double time;
        unsigned int nmbr;
        DoubleMatrix alpha;
    };

    ISystemLinearODENonLocalContionsM(const ODEGrid& grid);

    virtual void calculateForward(DoubleMatrix &x);
    virtual void calculateBackward(DoubleMatrix &x);

    void setLeftSeparatedCondition(const Condition &lscs);
    void setRightSeparatedCondition(const Condition &lscs);
    void addNonSeparatedCondition(const Condition &nsc);
    const std::vector<Condition>& nonSeparatedConditions() const;
    void setBetta(const DoubleMatrix &betta);
    void setSystemOrder(unsigned int n);
    unsigned int systemOrder() const;

    const ODEGrid& grid() const;

public:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;

//private:
    std::vector<Condition> nscs;
    Condition lscs;
    Condition rscs;
    DoubleMatrix betta;
    unsigned int n0;
    unsigned int n1;
    unsigned int n2;
    unsigned int n;

private:
    virtual void calculateIntervalF(unsigned int s, unsigned int row, unsigned int col);
    //virtual void calculateIntervalB(unsigned int s, unsigned int r);

protected:
    ODEGrid mgrid;
};

#endif // ISYSTEMLINEARODENONLOCALCONTIONS_H
