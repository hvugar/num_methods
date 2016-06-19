#ifndef MATRIX_H
#define MATRIX_H

#include "global.h"
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <math.h>

class MINIMUMSHARED_EXPORT CVector
{
public:
    CVector(size_t size = 0);
    virtual ~CVector();

    size_t size() const;
    void resize(size_t size);
    void clear();

    double L2Norm() const;
    double L1Norm() const;
    double LInfNorm() const;
    double EuclideanNorm() const;
    double EuclideanDistance(const CVector&) const;
    void L2Normalize();
    void L1Normalize();
    void EuclideanNormalize();

    double min() const;
    double max() const;

    CVector mid(size_t s, size_t e) const;

    double& operator[](size_t i);
    const double& operator[](size_t i) const;

    //    double at(unsigned int i) const;
    //    double* data() const;
    //    void add(double d);
    //    void insert(unsigned int i, double d);
    //    void remove(unsigned int i);

private:
    size_t msize;
    double* pdata;
};

class MINIMUMSHARED_EXPORT CMatrix
{
public:
    explicit CMatrix(size_t rows, size_t cols);
    virtual ~CMatrix();

    size_t rows() const;
    size_t columns() const;

    void clear();
    void resize(size_t m, size_t n);

    double min() const;
    double max() const;

    CVector& operator[](size_t i);
    const CVector& operator[](size_t i) const;

private:
    size_t mrows;
    size_t mcols;
    CVector *pvector;
};

#endif // MATRIX_H
