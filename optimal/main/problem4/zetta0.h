#ifndef ZETTA0_H
#define ZETTA0_H

#include <load_sys/islodenlcsv.h>
#include <ode/lode1o.h>

class Problem4Ex1;

//class Zetta0 : public ISystemLinearODENonLocalContionsV
//{
//public:
//    Zetta0(const Problem4Ex1 &p4);

//    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
//    virtual double B(double t, unsigned int k, unsigned int row) const;
//private:
//    const Problem4Ex1 &p4;
//};

class Zetta01 : public LinearODE1stOrder
{
public:
    Zetta01(const Problem4Ex1 &p);
    virtual unsigned int equationsNumber() const;
protected:
    virtual double A(const GridNodeODE &node, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(const GridNodeODE &node, unsigned int row = 0) const;
private:
    const Problem4Ex1 &p;
};

#endif // ZETTA0_H
