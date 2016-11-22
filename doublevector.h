#ifndef DOUBLEVECTOR_H
#define DOUBLEVECTOR_H

#include <vector>
#include <string.h>
#include <stdlib.h>

using namespace std;

class DoubleVector : public std::vector<double>
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
};

class DoubleMatrix : public std::vector<DoubleVector>
{
public:
    DoubleMatrix(unsigned int m = 0, unsigned int n = 0, double value = 0.0);
    virtual ~DoubleMatrix();
    void Clear();
    void Resize(unsigned int m, unsigned int n);
};

class DoubleCube : public std::vector<DoubleMatrix>
{
public:
    DoubleCube();
    virtual ~DoubleCube();

    void Resize(unsigned int Nz, unsigned int Ny, unsigned Nx);
    void Clear();
};

#endif // DOUBLEVECTOR_H
