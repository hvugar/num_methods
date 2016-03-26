#ifndef DOUBLEVECTOR_H
#define DOUBLEVECTOR_H

#include <vector>
#include "global.h"

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

    DoubleVector mid(unsigned int s, unsigned int e) const;

    /* static methods */
    static double L2Norm(const DoubleVector&);
    static double L1Norm(const DoubleVector&);
    static double LInfNorm(const DoubleVector&);
    static double EuclideanNorm(const DoubleVector&);
    static double EuclideanDistance(const DoubleVector&, const DoubleVector&);
    static void L2Normalize(DoubleVector&);
    static void L1Normalize(DoubleVector&);
    static void EuclideanNormalize(DoubleVector&);
    static DoubleVector null;
};

class MINIMUMSHARED_EXPORT DoubleMatrix : public std::vector<DoubleVector>
{
public:
    DoubleMatrix(unsigned int m = 0, unsigned int n = 0, double value = 0.0);
    virtual ~DoubleMatrix();
    void Clear();
    void Resize(unsigned int m, unsigned n);
};

class MINIMUMSHARED_EXPORT DoubleCube : public std::vector<DoubleMatrix>
{
public:
    DoubleCube();
    virtual ~DoubleCube();

    void Resize(unsigned int Nz, unsigned int Ny, unsigned Nx);
    void Clear();
};

class MINIMUMSHARED_EXPORT DblVector
{
public:
    DblVector(unsigned int n=0);
    virtual ~DblVector();

    unsigned int size() const;
    double at(unsigned int i) const;
    double* data() const;
    void add(double d);
    void insert(unsigned int i, double d);
    void remove(unsigned int i);
    void clear();
private:
    double* pdata;
    unsigned int msize;
};

#endif // DOUBLEVECTOR_H
