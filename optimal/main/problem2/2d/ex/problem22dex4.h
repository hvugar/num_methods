#ifndef PROBLEM22DEX4_H
#define PROBLEM22DEX4_H

#include "../abstractproblem22d.h"

//---------------------------------------------------------------------------------------------------------------//
// Məqalə üçün hazırlanmış test məsələ
// İdarə nöqtələrinin sayı 2.
// Ölçü cihazlarının sayı 3.
//---------------------------------------------------------------------------------------------------------------//

class Problem22DEx4;
class Problem2Forward2DEx4;
class Problem2Backward2DEx4;

//---------------------------------------------------------------------------------------------------------------//

class Problem2Forward2DEx4 : public IProblem2Forward2D
{
public:
    double fi;

protected:
    virtual double initial(const SpaceNodePDE &) const;
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const;
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;
};

//---------------------------------------------------------------------------------------------------------------//

class Problem2Backward2DEx4 : public IProblem2Backward2D
{
public:
    Problem22DEx4 *p22dEx4;
    DoubleMatrix *u;

protected:
    virtual double initial(const SpaceNodePDE &) const;
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const;
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual double h(const SpaceNodePDE &) const;
    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const;
};

//---------------------------------------------------------------------------------------------------------------//


class Problem22DEx4 : public AbstactProblem22D
{
public:
    static void Main(int argc, char* argv[]);

    Problem22DEx4();
    virtual ~Problem22DEx4();

    void compareNandAGradients();

private:
    Problem2Forward2DEx4 *forward;
    Problem2Backward2DEx4 *backward;
};

#endif // PROBLEM22DEX4_H
