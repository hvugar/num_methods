#ifndef ZETTA1_H
#define ZETTA1_H

#include <load_sys/islodenlcsm.h>

class Problem4Ex1;

class Zettai : public ISystemLinearODENonLocalContionsM
{
public:
    Zettai(const Problem4Ex1 &p4, unsigned int i);

    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int k, unsigned int row, unsigned int col) const;
private:
    const Problem4Ex1 &p4;
    unsigned int i;
};

#endif // ZETTA1_H
