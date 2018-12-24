#ifndef HYPERBOLICIBVP1_H
#define HYPERBOLICIBVP1_H

#include <grid/hibvp.h>
#include <time.h>
#include <QPixmap>
#include <QGuiApplication>
#include "../imaging/imaging.h"

class MINIMUMSHARED_EXPORT HyperbolicIBVP1 : public HyperbolicIBVP
{
public:
    HyperbolicIBVP1();
    double U(unsigned int i, unsigned int j) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
public:
    static void Main(int argc, char* argv[]);
};

class MINIMUMSHARED_EXPORT CCIHyperbolicIBVP1 : public CCIHyperbolicIBVP
{
public:
    static void Main(int argc, char* argv[]);

    virtual void layerInfo(const DoubleVector &, unsigned int) const;
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

    double integralU1(const DoubleVector &u) const;
    double integralU2(const DoubleVector &u) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    double a = 1.0;

    DoubleVector u2;
    DoubleVector u1;
    DoubleVector u0;
    DoubleVector ut;

    DoubleMatrix u10;
    bool calculateU = true;
};

class MINIMUMSHARED_EXPORT CC1IHyperbolicIBVP1 : public CC1IHyperbolicIBVP
{
public:
    static void Main(int argc, char* argv[]);

    virtual ~CC1IHyperbolicIBVP1();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    double a = 1.0;
    double alpha = 0.001;
    DoubleMatrix u10;
};

#endif // HYPERBOLICIBVP1_H
