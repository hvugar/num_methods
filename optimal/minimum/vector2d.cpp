#include "vector2d.h"
#include "matrix2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <exception>
#include <stdexcept>

DoubleVector::DoubleVector(unsigned int size, double val) : mSize(0), mData(NULL)
{
    if (size <= 0) return;

    mSize = size;
    mData = (double*) malloc(sizeof(double) * size);
    for (unsigned int i=0; i<size; i++) mData[i] = val;
}

DoubleVector::DoubleVector(const double *data, unsigned int size) : mSize(0), mData(NULL)
{
    if (size <= 0) return;

    mSize = size;
    mData = (double*) malloc(sizeof(double) * size);
    memcpy(mData, data, sizeof(double)*size);
}

DoubleVector::DoubleVector(const DoubleVector &vector) : mSize(0), mData(NULL)
{
    if (vector.mSize == 0) return;

    mSize = vector.mSize;
    mData = (double*) (malloc(sizeof(double)*mSize));
    memcpy(mData, vector.mData, sizeof(double)*mSize);
}

DoubleVector::DoubleVector(const DoubleMatrix &matrix) : mSize(0), mData(NULL)
{
    if (matrix.cols() != 1) return;

    mSize = matrix.rows();
    mData = (double*) (malloc(sizeof(double)*mSize));
    for (unsigned int i=0; i<mSize; i++) mData[i] = matrix.at(i,0);
}

DoubleVector::~DoubleVector()
{
    clear();
}

void DoubleVector::clear()
{
    if (mData != NULL)
    {
        free(mData);
        mData = NULL;
        mSize = 0;
    }
}

void DoubleVector::resize(unsigned int size, double value)
{
    if (size == mSize) return;

    if (size == 0) clear();

    if (size > 0)
    {
        if (mData == NULL)
        {
            mData = (double*) malloc(sizeof(double)*size);
            for (unsigned int i=0; i<size; i++) mData[i] = value;
            mSize = size;
        }
        else if (size != mSize)
        {
            double *ptr = (double *) realloc(mData, sizeof(double) * size);
            for (unsigned int i=mSize; i<size; i++) ptr[i] = value;
            mData = ptr;
            mSize = size;
        }
    }
}

bool DoubleVector::empty() const
{
    return mSize == 0;
}

double& DoubleVector::at(unsigned int n)
{
    if (n >= mSize)
    {
        throw std::out_of_range("out of range");
    }
    return mData[n];
}

const double& DoubleVector::at(unsigned int n) const
{
    if (n >= mSize)
    {
        throw std::out_of_range("out of range");
    }
    return mData[n];
}

unsigned int DoubleVector::size() const
{
    return mSize;
}

double* DoubleVector::data() NOEXCEPT
{
    return mData;
}

const double* DoubleVector::data() const NOEXCEPT
{
    return mData;
}

double DoubleVector::L2Norm() const
{
    double norm = 0.0;
    if (!empty())
    {
        for (unsigned int i=0; i < mSize; i++)
        {
            double item = mData[i];
            norm += item*item;
        }
        return sqrt(norm);
    }
    return norm;
}

double DoubleVector::L1Norm() const
{
    double norm = 0.0;
    if (!empty())
    {
        for (unsigned int i=0; i < mSize; i++)
        {
            norm += fabs(mData[i]);
        }
    }
    return norm;
}

double DoubleVector::LInfNorm() const
{
    double norm = 0.0;
    if (!empty())
    {
        norm = fabs(mData[0]);
        for (unsigned int i=1; i < mSize; i++)
        {
            if (norm < fabs(mData[i])) norm = fabs(mData[i]);
        }
    }
    return norm;
}

double DoubleVector::EuclideanNorm() const
{
    return L2Norm();
}

double DoubleVector::EuclideanDistance(const DoubleVector &p) const
{
    if ( mSize != p.mSize ) return INFINITY;

    double distance = 0.0;
    for (unsigned int i=0; i<mSize; i++)
    {
        double dx = mData[i]-p.mData[i];
        distance += dx*dx;
    }
    return sqrt(distance);
}

void DoubleVector::L2Normalize()
{
    double norm = L2Norm();
    if (norm > DBL_EPSILON) for (unsigned int i=0; i<mSize; i++) mData[i] /= norm;
}

void DoubleVector::L1Normalize()
{
    double norm = L1Norm();
    if (norm > DBL_EPSILON) for (unsigned int i=0; i<mSize; i++) mData[i] /= norm;
}

void DoubleVector::EuclideanNormalize()
{
    L2Normalize();
}

double DoubleVector::min() const
{
    double minimum = NAN;
    if (!empty())
    {
        minimum = mData[0];
        for (unsigned int i=1; i<mSize; i++) if (minimum > mData[i]) minimum = mData[i];
    }
    return minimum;
}

double DoubleVector::max() const
{
    double maximum = NAN;
    if (!empty())
    {
        maximum = mData[0];
        for (unsigned int i=1; i<mSize; i++) if (maximum < mData[i]) maximum = mData[i];
    }
    return maximum;
}

DoubleVector DoubleVector::mid(unsigned int s, unsigned int e) const
{
    DoubleVector vector((e-s)+1);
    for (unsigned int i=s; i<=e; i++) vector.mData[i-s] = mData[i];
    return vector;
}

DoubleVector& DoubleVector::operator =(const DoubleVector& other)
{
    if (this != &other)
    {
        unsigned int size = other.size();
        resize(size, 0.0);
        for (unsigned int i=0; i<size; i++) mData[i] = other.mData[i];
    }
    return *this;
}

double& DoubleVector::operator [](unsigned int n)
{
    return mData[n];
}

double DoubleVector::operator [](unsigned int n) const
{
    return mData[n];
}

DoubleVector& DoubleVector::operator <<(double value)
{
    if (mData==NULL)
    {
        mSize = 1;
        mData = (double*) malloc(sizeof(double)*mSize);
        mData[0] = value;
    }
    else
    {
        mSize++;
        mData = (double*)realloc(mData, sizeof(double)*mSize);
        mData[mSize-1] = value;
    }
    return *this;
}

DoubleVector& DoubleVector::operator +(const DoubleVector &other)
{
    if (mSize != other.mSize)
    {
        //throw std::exception("");
    }
    else
    {
        for (unsigned int i=0; i<mSize; i++)
        {
            mData[i] += other.mData[i];
        }
    }
    return *this;
}

//DoubleVector& DoubleVector::operator -(const DoubleVector &other)
//{
//    if (mSize != other.mSize)
//    {
//        //throw std::exception("");
//    }
//    else
//    {
//        for (unsigned int i=0; i<mSize; i++)
//        {
//            mData[i] -= other.mData[i];
//        }
//    }
//    return *this;
//}

DoubleVector operator *(double scalar, const DoubleVector &v)
{
    DoubleVector rv(v.mSize);
    for (unsigned int i=0; i<v.mSize; i++) rv.mData[i] = v.mData[i]*scalar;
    return rv;
}

DoubleVector operator *(const DoubleVector &v, double scalar)
{
    return scalar*v;
}

bool operator ==(const DoubleVector& vector1, const DoubleVector& vector2)
{
    if (vector1.mSize != vector2.mSize) return false;

    unsigned int length = vector1.mSize;
    for (unsigned int i=0; i<length; i++) if (vector1.mData[i] != vector2.mData[i]) return false;

    return true;
}

bool operator !=(const DoubleVector& vector1, const DoubleVector& vector2)
{
    return !(vector1 == vector2);
}

void DoubleVector::randomData()
{
    for (unsigned int i=0; i<mSize; i++) mData[i] = rand() % 100;
}

void DoubleVector::print()
{
    for (unsigned int i=0; i<mSize; i++)
    {
        printf("%10.6f ", mData[i]);
    }
    puts("");
}
