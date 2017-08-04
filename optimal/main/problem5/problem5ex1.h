#ifndef PROBLEM5EX1_H
#define PROBLEM5EX1_H

#include <islodenlcs.h>

#define SAMPLE_1

class Problem5Ex1
{
public:
    Problem5Ex1();

    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int k, unsigned int row) const;
    virtual double B(double t, unsigned int k, unsigned int row, unsigned int col, unsigned int i) const;
    virtual double C(double t, unsigned int k, unsigned int row) const;
    virtual double g(double t, unsigned int k, unsigned int row, unsigned int i) const;

    double X(double t, unsigned int k, unsigned int i) const;

    unsigned int L0 = 2;
    unsigned int L1 = 3;
};

class HelpZ0 : public ISystemLinearODENonLocalContions
{
public:
    HelpZ0(const ODEGrid &grid, const Problem5Ex1 &p5);

    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const;
private:
    const Problem5Ex1 &p5;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HelpZ0::HelpZ0(const ODEGrid &grid, const Problem5Ex1 &p5) : ISystemLinearODENonLocalContions(grid), p5(p5) {}

double HelpZ0::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p5.A(t,k,row,col);
}

double HelpZ0::B(double t, unsigned int k, unsigned int row) const
{
    return p5.B(t, k, row);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class HelpZ1 : public ISystemLinearODENonLocalContions
{
public:
    HelpZ1(const ODEGrid &grid, const Problem5Ex1 &p5);

    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const;
private:
    const Problem5Ex1 &p5;
};

HelpZ1::HelpZ1(const ODEGrid &grid, const Problem5Ex1 &p5) : ISystemLinearODENonLocalContions(grid), p5(p5) {}

double HelpZ1::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p5.A(t,k,row,col);
}

double HelpZ1::B(double t, unsigned int k, unsigned int row) const
{
    return p5.B(t, k, row);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class HelpZ2 : public ISystemLinearODENonLocalContions
{
public:
    HelpZ2(const ODEGrid &grid, const Problem5Ex1 &p5);

    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const;
private:
    const Problem5Ex1 &p5;
};

HelpZ2::HelpZ2(const ODEGrid &grid, const Problem5Ex1 &p5) : ISystemLinearODENonLocalContions(grid), p5(p5) {}

double HelpZ2::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p5.A(t,k,row,col);
}

double HelpZ2::B(double t, unsigned int k, unsigned int row) const
{
    return p5.B(t, k, row);
}


#endif // PROBLEM5EX1_H
