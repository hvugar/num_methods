#ifndef PARABOLICIBVP1_H
#define PARABOLICIBVP1_H

#include <grid/pibvp.h>
#include <time.h>

#define SAMPLE_0

class MINIMUMSHARED_EXPORT ParabolicIBVP1 : public ParabolicIBVP
{
public:
    void static Main(int argc, char* argv[]);

public:
    double U(unsigned int i, unsigned int j) const;
    double U(const SpaceNodePDE &sn, const TimeNodePDE& tn) const;

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE& tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector &, unsigned int) const;
};

class MINIMUMSHARED_EXPORT ParabolicIBVP2 : public ParabolicIBVP
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE& tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    void static Main(int argc, char* argv[]);
    double U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
};

class MINIMUMSHARED_EXPORT CCParabolicIBVP1 : public IHeatEquationIBVP
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual double U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    static void Main(int argc, char* argv[]);
};

#endif // PARABOLICIBVP1_H
