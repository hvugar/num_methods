#ifndef WAVEEQUATIONIBVP_H
#define WAVEEQUATIONIBVP_H

#include "border_global.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BORDERSHARED_EXPORT WaveEquationIBVP : public IWaveEquationIBVP
{
public:
    static void Main(int argc, char* argv[]);

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector&, const TimeNodePDE&) const;
    virtual void layerInfo(const DoubleMatrix&, const TimeNodePDE&) const;

private:
    void saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void saveToTextF(const DoubleMatrix &u, const TimeNodePDE &tn) const;

    DoubleMatrix u1;
    DoubleMatrix u2;

    double integralUP(const DoubleVector &u) const;
    double integralUK(const DoubleVector &u) const;
    double integralUP(const DoubleMatrix &u) const;
    double integralUK(const DoubleMatrix &u) const;

    DoubleVector vu0;
    DoubleVector vu1;
    DoubleVector vu2;
    DoubleVector vut;

    DoubleMatrix mu0;
    DoubleMatrix mu1;
    DoubleMatrix mu2;

    DoubleMatrix mut;
    DoubleMatrix muc;

    bool calculateU = true;
    unsigned int method = 0;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BORDERSHARED_EXPORT ConjugateCdIHyperbolicIBVP1 : public IConjugateWaveEquationIBVP
{
public:
    static void Main(int argc, char* argv[]);

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector&, const TimeNodePDE&) const {}
    virtual void layerInfo(const DoubleMatrix&, const TimeNodePDE&) const;

private:
    double a = 1.0;
    double alpha = 0.1;
    DoubleMatrix u10;
};

#endif // WAVEEQUATIONIBVP_H
