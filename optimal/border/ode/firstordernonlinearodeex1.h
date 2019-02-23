#ifndef FIRST_ORDER_NON_LINEAR_ODE_EX1_H
#define FIRST_ORDER_NON_LINEAR_ODE_EX1_H

#include <ode/nlode1o.h>
#include <cmath>

class FirstOrderNonLinearODEEx1 : public FirstOrderNonLinearODE
{
public:
    virtual double f(double x, double y, unsigned int k) const;
};

class B
{
public:
    B(int i=1) : a(i) { std::cout << "B constructor " << a << std::endl; }
    virtual ~B() { std::cout << "B destructor" << std::endl; }

    int a;
};

#endif // FIRST_ORDER_NON_LINEAR_ODE_EX1_H
