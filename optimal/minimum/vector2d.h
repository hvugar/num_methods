#ifndef VECTOR2D_H
#define VECTOR2D_H

#include "global.h"
#include <stdio.h>

class DoubleMatrix;

class MINIMUMSHARED_EXPORT DoubleVector
{
public:
    explicit DoubleVector(unsigned int size=0, double val = 0.0);
    DoubleVector(const double* data, unsigned int size);
    DoubleVector(const DoubleVector &vector);
    DoubleVector(const DoubleMatrix &matrix);
    virtual ~DoubleVector();

    void clear();
    void resize(unsigned int size, double value = 0.0);
    bool empty() const;
    double& at (unsigned int n);
    const double& at (unsigned int n) const;
    unsigned int size() const;
    double* data() NOEXCEPT;

    void randomData();

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
    DoubleVector mid(unsigned int s, unsigned int e) const;

    ///////
    void assign (unsigned int n, const double& value);
    double& back();
    const double& back() const;
    double& begin();
    const double& begin() const;
    unsigned int capacity() const;
    double& cbegin() const NOEXCEPT;
    double& cend() const NOEXCEPT;
    const double* data() const NOEXCEPT;
    double& end();
    const double& end() const;
    double& erase (unsigned int position);
    double& erase (unsigned int first, unsigned int last);
    double& front();
    const double& front() const;
    double& insert (unsigned int position, const double& val);
    void insert (unsigned int position, unsigned int n, const double& val);
    void insert (unsigned int position, unsigned int first, unsigned int last);
    unsigned int max_size() const;
    DoubleVector& operator= (const DoubleVector& x);
    double& operator[] (unsigned int n);
    double operator[] (unsigned int n) const;
    DoubleVector& operator <<(double value);
    DoubleVector& operator +(const DoubleVector &other);

    void swap(DoubleVector& x);

    void print(unsigned int cols, char* label = NULL, unsigned int start=0, unsigned int end=0, FILE* file=stdout);
    void print();

    friend class DoubleMatrix;
    friend MINIMUMSHARED_EXPORT DoubleVector operator*(double scalar, const DoubleVector &);
//    friend MINIMUMSHARED_EXPORT DoubleVector operator +(const DoubleVector&, const DoubleVector &);
//    friend MINIMUMSHARED_EXPORT DoubleVector operator -(const DoubleVector&, const DoubleVector &);
//    friend MINIMUMSHARED_EXPORT DoubleVector operator *(const DoubleVector&, const DoubleVector &);
//    friend MINIMUMSHARED_EXPORT DoubleVector operator *(double, const DoubleVector &);
private:
    unsigned int mSize;
    double *mData;
};

#endif // VECTOR2D_H
