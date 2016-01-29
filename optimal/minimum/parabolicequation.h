#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include "global.h"
#include "tomasmethod.h"
#include "printer.h"

/**
 * @brief Параболического уравнения
 */
class MINIMUMSHARED_EXPORT IParabolicEquation
{
public:
    virtual double fi(unsigned int i) const = 0;
    virtual double m1(unsigned int j) const = 0;
    virtual double m2(unsigned int j) const = 0;
    virtual double f(unsigned int i, unsigned int j) const = 0;
    /**
     * @brief calculateU
     * @param u
     * @param hx
     * @param ht
     * @param N
     * @param M
     * @param a положительная константа. число a^2 является коэффициентом температуропроводности
     */
    virtual void calculateU(DoubleVector &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
    /**
     * @brief calculateU
     * @param u
     * @param hx
     * @param ht
     * @param N
     * @param M
     * @param a положительная константа. число a^2 является коэффициентом температуропроводности
     */
    virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
};

/**
 * @brief Обратная параболического уравнения
 *
 */
class MINIMUMSHARED_EXPORT IBackwardParabolicEquation
{
public:
    virtual double bfi(unsigned int i) const = 0;
    virtual double bm1(unsigned int j) const = 0;
    virtual double bm2(unsigned int j) const = 0;
    virtual double bf(unsigned int i, unsigned int j) const = 0;
    virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
};

/**
 * @brief The IParabolicEquation2D class
 */
class MINIMUMSHARED_EXPORT IParabolicEquation2D
{
public:
    virtual double fi(unsigned int i, unsigned int j) const = 0;
    virtual double m1(unsigned int j, unsigned int k) const = 0;
    virtual double m2(unsigned int j, unsigned int k) const = 0;
    virtual double m3(unsigned int i, unsigned int k) const = 0;
    virtual double m4(unsigned int i, unsigned int k) const = 0;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const = 0;

    void calculateU(DoubleMatrix &u, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
    void calculateU(DoubleCube &u, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
};

class MINIMUMSHARED_EXPORT IBackwardParabolicEquation2D
{
public:
    virtual double bfi(unsigned int i, unsigned int j) const = 0;
    virtual double bm1(unsigned int j, unsigned int k) const = 0;
    virtual double bm2(unsigned int j, unsigned int k) const = 0;
    virtual double bm3(unsigned int i, unsigned int k) const = 0;
    virtual double bm4(unsigned int i, unsigned int k) const = 0;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const = 0;

    void calculateU(DoubleCube &psi, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
};

#endif // HEATEQUATION_H
