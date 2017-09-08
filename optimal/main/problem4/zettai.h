#ifndef ZETTA1_H
#define ZETTA1_H

#include <load_sys/islodenlcsm.h>
#include <ode/lode1o.h>
#include <matrix2d.h>

class Problem4Ex1;

class Zettai : public ISystemLinearODENonLocalContionsM
{
public:
    Zettai(const Problem4Ex1 &p4, unsigned int i);

    virtual double A(const GridNodeODE &node, unsigned int row, unsigned int col) const;
    virtual double B(const GridNodeODE &node, unsigned int row, unsigned int col) const;
private:
    const Problem4Ex1 &p4;
    unsigned int i;
};

class Zettai1 : public LinearODE1stOrder
{
public:
    Zettai1(const Problem4Ex1 &p4, unsigned int i);
    virtual unsigned int equationsNumber() const;

    void calculateM(const std::vector<LinearODE1stOrder::Condition> &cs, const DoubleMatrix &betta, std::vector<std::vector<DoubleVector>> &zmi);

protected:
    virtual double A(const GridNodeODE &node, unsigned int row, unsigned int col) const;
    virtual double B(const GridNodeODE &node, unsigned int row) const;
    virtual double C(const GridNodeODE &node, unsigned int row, unsigned int col) const;
private:
    const Problem4Ex1 &p;
    unsigned int i;
    unsigned int cur_col;
};

#endif // ZETTA1_H
