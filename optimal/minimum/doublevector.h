#ifndef DOUBLEVECTOR_H
#define DOUBLEVECTOR_H

#include "global.h"
#include <vector>
#include <string.h>
#include <stdlib.h>

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

    double min() const;
    double max() const;

    DoubleVector* mid(unsigned int s, unsigned int e) const;
};

class MINIMUMSHARED_EXPORT DoubleMatrix : public std::vector<DoubleVector>
{
public:
    DoubleMatrix(unsigned int m = 0, unsigned int n = 0, double value = 0.0);
    virtual ~DoubleMatrix();
    void Clear();
    void Resize(unsigned int m, unsigned int n);

    double min() const;
    double max() const;

    double determinant() const;
    DoubleMatrix* transpose() const;
    DoubleMatrix* inverse() const;
    DoubleMatrix* minor(size_t row, size_t col) const;
    DoubleMatrix* multiply(const DoubleMatrix &m) const;
    DoubleMatrix* multiply(const DoubleVector &v) const;

    DoubleMatrix operator+(const DoubleMatrix &A) const;
    DoubleMatrix operator-(const DoubleMatrix &A) const;
    DoubleMatrix operator*(const DoubleMatrix &A) const;
    DoubleMatrix operator/(const DoubleMatrix &A) const;
    DoubleMatrix operator*(double c) const;
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
