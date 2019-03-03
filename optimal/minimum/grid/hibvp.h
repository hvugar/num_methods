#ifndef HYPERBOLICIBVP_H
#define HYPERBOLICIBVP_H

#include "ibvp.h"

/**
 * @brief The IHyperbolicIBVP class
 * u_tt(x,t) = a^2u_xx(x,t) + f(x,t), t in (0,T], x in (0,l)
 */
class MINIMUMSHARED_EXPORT IHyperbolicIBVP : public InitialBoundaryValueProblemPDE
{
public:
    virtual ~IHyperbolicIBVP();

protected:
    virtual double initial1(const SpaceNodePDE &sn) const = 0;
    virtual double initial2(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

protected:
    bool isHalfLayerEnabled() const;
    void setHalfLayerEnabled(bool) const;

private:
    bool _isHalfLayerEnabled = false;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////s

class MINIMUMSHARED_EXPORT CcIHyperbolicIBVP : public IHyperbolicIBVP
{
public:
    virtual ~CcIHyperbolicIBVP();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    void explicit_calculate_D1V1(DoubleVector &u, double a) const;
    void implicit_calculate_D1V1(DoubleVector &u, double a, double lambda=0.25) const;

    void explicit_calculate_D2V1(DoubleMatrix &u, double a) const;
    void implicit_calculate_D2V1(DoubleMatrix &u, double a, double lambda = 0.25) const;

private:
    void explicit_calculate_D1V1_initial(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a) const;

    void explicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a) const;
    void explicit_calculate_D2V1_border(DoubleMatrix &u, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const;
    void implicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a) const;
    void implicit_calculate_D2V1_border(DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn05, const TimeNodePDE &tn10) const;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MINIMUMSHARED_EXPORT CdIHyperbolicIBVP : public IHyperbolicIBVP
{
public:
    virtual ~CdIHyperbolicIBVP();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    void explicit_calculate_D1V1(DoubleVector &u, double a, double alpha) const;
    void implicit_calculate_D1V1(DoubleVector &u, double a, double alpha, double lambda=0.25) const;

    void explicit_calculate_D2V1(DoubleMatrix &u, double a, double alpha) const;
    void implicit_calculate_D2V1(DoubleMatrix &u, double a, double alpha, double lambda=0.25) const;

private:
    void explicit_calculate_D1V1_initial(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a, double sigma) const;

    void explicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const;
    void explicit_calculate_D2V1_border(DoubleMatrix &u, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const;
    void implicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const;
    void implicit_calculate_D2V1_border(DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn05, const TimeNodePDE &tn10) const;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MINIMUMSHARED_EXPORT ConjugateCdIHyperbolicIBVP : public IHyperbolicIBVP
{
public:
    virtual ~ConjugateCdIHyperbolicIBVP();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    void explicit_calculate_D1V1(DoubleVector &p, double a, double alpha) const;
    void implicit_calculate_D1V1(DoubleVector &p, double a, double alpha, double lambda=0.25) const;

    void explicit_calculate_D2V1(DoubleMatrix &p, double a, double alpha) const;
    void implicit_calculate_D2V1(DoubleMatrix &p, double a, double alpha, double lambda=0.25) const;

private:
    void explicit_calculate_D1V1_initial(DoubleVector &p00, DoubleVector &u10, unsigned int N, double hx, double ht, double a, double sigma) const;

    void explicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha, unsigned int L) const;
    void explicit_calculate_D2V1_border(DoubleMatrix &u, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const;
    void implicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha, unsigned int L) const;
    void implicit_calculate_D2V1_border(DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn05, const TimeNodePDE &tn10) const;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // HYPERBOLICIBVP_H
