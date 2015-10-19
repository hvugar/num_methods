#ifndef DOUBLEVECTOR_H
#define DOUBLEVECTOR_H

#include <vector>
#include "global.h"

using namespace std;

class MINIMUMSHARED_EXPORT DoubleVector : public std::vector<double>
{
public:
    explicit DoubleVector(unsigned int size = 0, double value = 0.0);
    virtual ~DoubleVector();

    double L2Norm() const;
    double L1Norm() const;
    double LInfNorm() const;
    double EuclideanNorm() const;
    double EuclideanDistance(const DoubleVector&) const;
    void L2Normalize();
    void L1Normalize();
    void EuclideanNormalize();

    /* static methods */
    static double L2Norm(const DoubleVector&);
    static double L1Norm(const DoubleVector&);
    static double LInfNorm(const DoubleVector&);
    static double EuclideanNorm(const DoubleVector&);
    static double EuclideanDistance(const DoubleVector&, const DoubleVector&);
    static void L2Normalize(DoubleVector&);
    static void L1Normalize(DoubleVector&);
    static void EuclideanNormalize(DoubleVector&);
};

class MINIMUMSHARED_EXPORT DoubleMatrix : public std::vector<DoubleVector>
{
public:
    DoubleMatrix();
    virtual ~DoubleMatrix();
    void Clear();
    void Resize(unsigned int Ny, unsigned Nx);
};

class MINIMUMSHARED_EXPORT DoubleCube : public std::vector<DoubleMatrix>
{
public:
    DoubleCube();
    virtual ~DoubleCube();

    void Resize(unsigned int Nz, unsigned int Ny, unsigned Nx);
    void Clear();
};


#endif // DOUBLEVECTOR_H
