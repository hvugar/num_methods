#ifndef HYPERBOLICIBVP1_H
#define HYPERBOLICIBVP1_H

#include <grid/hibvp.h>
#include <time.h>
#include <QPixmap>
#include <QGuiApplication>
#include "../imaging/imaging.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MINIMUMSHARED_EXPORT CcIHyperbolicIBVP1 : public CcIHyperbolicIBVP
{
public:
    static void Main(int argc, char* argv[]);

    virtual void layerInfo(const DoubleVector &, unsigned int) const;
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

    double integralU1(const DoubleVector &u) const;
    double integralU2(const DoubleVector &u) const;
    double integralU1(const DoubleMatrix &u) const;
    double integralU2(const DoubleMatrix &u) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    double a = 1.0;
    unsigned int method = 0;

    DoubleVector vu0;
    DoubleVector vu1;
    DoubleVector vu2;
    DoubleVector vut;

    DoubleMatrix mu0;
    DoubleMatrix mu1;
    DoubleMatrix mu2;
    DoubleMatrix mut;

    bool calculateU = true;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MINIMUMSHARED_EXPORT CdIHyperbolicIBVP1 : public CdIHyperbolicIBVP
{
public:
    static void Main(int argc, char* argv[]);

    virtual ~CdIHyperbolicIBVP1();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    double a;
    double alpha;
    DoubleMatrix u10;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MINIMUMSHARED_EXPORT ConjugateCdIHyperbolicIBVP1 : public ConjugateCdIHyperbolicIBVP
{
public:
    static void Main(int argc, char* argv[]);

    virtual ~ConjugateCdIHyperbolicIBVP1();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    double a = 1.0;
    double alpha = 0.1;
    DoubleMatrix u10;
};

#endif // HYPERBOLICIBVP1_H
