#ifndef R1MINIMIZE_H
#define R1MINIMIZE_H

class R1Minimize
{
public:
    R1Minimize();
    virtual ~R1Minimize();
    virtual double straightLineSearch(double x0, double dx, double &a, double &b);

protected:
    virtual double fx(double x) = 0;
};

#endif // R1MINIMIZE_H
