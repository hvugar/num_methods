#ifndef HYPERBOLIC_IBVP_H
#define HYPERBOLIC_IBVP_H

#include "ibvp.h"

/**
 * @brief The IHyperbolicIBVP class
 * u_tt(x,t) = a^2u_xx(x,t) + f(x,t), t in (0,T], x in (0,l)
 */
class MINIMUMSHARED_EXPORT IHyperbolicIBVP : public InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const {}
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}
};

//--------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IWaveEquationIBVP : public IHyperbolicIBVP
{
public:
    explicit IWaveEquationIBVP(double waveSpeed = 1.0, double waveDissipation = 0.0);

    virtual double waveSpeed() const;
    virtual double waveDissipation() const;

    void setWaveSpeed(double waveSpeed);
    void setWaveDissipation(double waveDissipation);

    void explicit_calculate_D1V1() const;
    void implicit_calculate_D1V1() const;

    void explicit_calculate_D2V1() const;
    void implicit_calculate_D2V1() const;

    virtual double lambda() const;

private:
    void explicit_calculate_D1V1_initial(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a, double alpha) const;
    void explicit_calculate_D1V1_border(DoubleVector &u, unsigned int N, double hx, double ht, const TimeNodePDE &tn) const;

    void implicit_calculate_D1V1_initial(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a, double alpha) const;
    void implicit_calculate_D1V1_border(DoubleVector &u, unsigned int N, double hx, const TimeNodePDE &tn) const;

    void explicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const;
    void explicit_calculate_D2V1_border(DoubleMatrix &u, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const;

    void implicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const;
    void implicit_calculate_D2V1_border(DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn05, const TimeNodePDE &tn10) const;

protected:
    double _waveSpeed;
    double _waveDissipation;
};

//--------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IConjugateWaveEquationIBVP : public IHyperbolicIBVP
{
public:
    explicit IConjugateWaveEquationIBVP(double waveSpeed = 1.0, double waveDissipation = 0.0);

    virtual double waveSpeed() const;
    virtual double waveDissipation() const;

    void setWaveSpeed(double waveSpeed);
    void setWaveDissipation(double waveDissipation);

    void explicit_calculate_D1V1() const;
    void implicit_calculate_D1V1() const;

    void explicit_calculate_D2V1() const;
    void implicit_calculate_D2V1() const;

    virtual double lambda() const;

private:
    void explicit_calculate_D1V1_initial(DoubleVector &p00, DoubleVector &p10, unsigned int N, double hx, double ht, double a, double alpha) const;
    void explicit_calculate_D1V1_border(DoubleVector &p, unsigned int N, double hx, double ht, const TimeNodePDE &tn) const;

    void implicit_calculate_D1V1_initial(DoubleVector &p00, DoubleVector &p10, unsigned int N, double hx, double ht, double a, double alpha) const;
    void implicit_calculate_D1V1_border(DoubleVector &p, unsigned int N, double hx, const TimeNodePDE &tn) const;

    void explicit_calculate_D2V1_initial(DoubleMatrix &p00, DoubleMatrix &p10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha, unsigned int L) const;
    void explicit_calculate_D2V1_border(DoubleMatrix &p, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const;

    void implicit_calculate_D2V1_initial(DoubleMatrix &p00, DoubleMatrix &p05, DoubleMatrix &p10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha, unsigned int L) const;
    void implicit_calculate_D2V1_border(DoubleMatrix &p05, DoubleMatrix &p10, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn05, const TimeNodePDE &tn10) const;

protected:
    double _waveSpeed;
    double _waveDissipation;
};

//--------------------------------------------------------------------------------------------------------------//

#endif // HYPERBOLIC_IBVP_H
