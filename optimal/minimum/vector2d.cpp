#include "vector2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <exception>
#include <stdexcept>
#include "matrix2d.h"

DoubleVector::DoubleVector(unsigned int size, double val) : mLength(0), mData(NULL)
{
    if (size <= 0) return;

    mLength = size;
    mData = (double*) malloc(sizeof(double) * size);
    for (unsigned int i=0; i<size; i++) mData[i] = val;
}

DoubleVector::DoubleVector(const double* data, unsigned int size) : mLength(0), mData(NULL)
{
    if (size <= 0) return;

    mLength = size;
    mData = (double*) malloc(sizeof(double) * size);
    memcpy(mData, data, sizeof(double)*size);
}

DoubleVector::DoubleVector(const DoubleVector &vector) : mLength(0), mData(NULL)
{
    if (vector.mLength == 0) return;

    mLength = vector.mLength;
    mData = (double*) (malloc(sizeof(double)*mLength));
    memcpy(mData, vector.mData, sizeof(double)*mLength);
}

DoubleVector::DoubleVector(const DoubleMatrix &matrix) : mLength(0), mData(NULL)
{
    if (matrix.cols() != 1) return;

    mLength = matrix.rows();
    mData = (double*) (malloc(sizeof(double)*mLength));
    for (unsigned int i=0; i<mLength; i++) mData[i] = matrix.at(i,0);
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
        mLength = 0;
    }
}

void DoubleVector::resize(unsigned int length, double value)
{
    if (length == mLength) return;

    if (length == 0) clear();

    if (length > 0)
    {
        if (mData == NULL)
        {
            mData = (double*) malloc(sizeof(double)*length);
            for (unsigned int i=0; i<length; i++) mData[i] = value;
            mLength = length;
        }
        else if (length != mLength)
        {
            double *ptr = (double *) realloc(mData, sizeof(double) * length);
            for (unsigned int i=mLength; i<length; i++) ptr[i] = value;
            mData = ptr;
            mLength = length;
        }
    }
}

bool DoubleVector::empty() const
{
    return mLength == 0;
}

double& DoubleVector::at(unsigned int n)
{
    if (n >= mLength)
    {
        throw std::out_of_range("out of range");
    }
    return mData[n];
}

const double& DoubleVector::at(unsigned int n) const
{
    if (n >= mLength)
    {
        throw std::out_of_range("out of range");
    }
    return mData[n];
}

unsigned int DoubleVector::length() const
{
    return mLength;
}

/********************************************************************
 *                      NORM
 *******************************************************************/

double DoubleVector::L1Norm() const
{
    double norm = 0.0;
    if (!empty())
    {
        for (unsigned int i=0; i < mLength; i++)
        {
            norm += fabs(mData[i]);
        }
    }
    return norm;
}

double DoubleVector::L2Norm() const
{
    double norm = 0.0;
    if (!empty())
    {
        for (unsigned int i=0; i < mLength; i++)
        {
            double item = mData[i];
            norm += item*item;
        }
        return sqrt(norm);
    }
    return norm;
}

double DoubleVector::LpNorm(unsigned int p) const
{
    double norm = 0.0;
    if (!empty())
    {
        for (unsigned int i=0; i < mLength; i++)
        {
            double item = mData[i];
            norm += pow(item, (double)p);
        }
        return pow(norm, 1.0/(double)p);
    }
    return norm;
}

double DoubleVector::LInfNorm() const
{
    double norm = 0.0;
    if (!empty())
    {
        norm = fabs(mData[0]);
        for (unsigned int i=1; i < mLength; i++)
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

void DoubleVector::L1Normalize()
{
    double norm = L1Norm();
    if (norm > DBL_EPSILON) for (unsigned int i=0; i<mLength; i++) mData[i] /= norm;
}

void DoubleVector::L2Normalize()
{
    double norm = L2Norm();
    if (norm > DBL_EPSILON) for (unsigned int i=0; i<mLength; i++) mData[i] /= norm;
}

void DoubleVector::EuclideanNormalize()
{
    L2Normalize();
}

/********************************************************************
 *                      NORM
 *******************************************************************/

double DoubleVector::min() const
{
    double minimum = NAN;
    if (!empty())
    {
        minimum = mData[0];
        for (unsigned int i=1; i<mLength; i++) if (minimum > mData[i]) minimum = mData[i];
    }
    return minimum;
}

double DoubleVector::max() const
{
    double maximum = NAN;
    if (!empty())
    {
        maximum = mData[0];
        for (unsigned int i=1; i<mLength; i++) if (maximum < mData[i]) maximum = mData[i];
    }
    return maximum;
}

DoubleVector DoubleVector::mid(unsigned int s, unsigned int e) const
{
    DoubleVector vector((e-s)+1);
    for (unsigned int i=s; i<=e; i++) vector.mData[i-s] = mData[i];
    return vector;
}

double DoubleVector::EuclideanDistance(const DoubleVector &p) const
{
    if ( mLength != p.mLength ) return INFINITY;

    double distance = 0.0;
    for (unsigned int i=0; i<mLength; i++)
    {
        double dx = mData[i]-p.mData[i];
        distance += dx*dx;
    }
    return sqrt(distance);
}

double* DoubleVector::data() NOEXCEPT
{
    return mData;
}

const double* DoubleVector::data() const NOEXCEPT
{
    return mData;
}

double& DoubleVector::operator [](unsigned int n)
{
    return mData[n];
}

double DoubleVector::operator [](unsigned int n) const
{
    return mData[n];
}

DoubleVector& DoubleVector::operator =(const DoubleVector& other)
{
    // the vector object holds reusable storage, such as a heap-allocated buffer mData

    if (this != &other)  // self-assignment check expected
    {
        clear();
        mLength = other.mLength;
        resize(mLength, 0.0);
        memcpy(mData, other.mData, sizeof(double)*mLength);
    }
    return *this;
}

DoubleVector& DoubleVector::operator +=(const DoubleVector& v)
{
    if (mLength != v.mLength)
    {
        fprintf(stderr, "DoubleVector& DoubleVector::operator +=(const DoubleVector& v) this:%d v:%d\n",
                mLength, v.mLength);
        fflush(stderr);

        throw DoubleMatrixException(1);
    }
    for (unsigned int i=0; i < mLength; i++)
    {
        mData[i] += v.mData[i];
    }
    return *this;
}

DoubleVector& DoubleVector::operator -=(const DoubleVector& v)
{
    if (mLength != v.mLength)
    {
        fprintf(stderr, "DoubleVector& DoubleVector::operator +=(const DoubleVector& v) this:%d v:%d\n",
                mLength, v.mLength);
        fflush(stderr);

        throw DoubleMatrixException(1);
    }
    for (unsigned int i=0; i < mLength; i++)
    {
        mData[i] -= v.mData[i];
    }
    return *this;
}

DoubleVector& DoubleVector::operator *=(const DoubleVector& v)
{
    if (mLength != v.mLength)
    {
        fprintf(stderr, "DoubleVector& DoubleVector::operator *=(const DoubleVector& v) this:%d v:%d\n",
                mLength, v.mLength);
        fflush(stderr);
        throw DoubleMatrixException(3);
    }

    for (unsigned int i=0; i<v.mLength; i++)
    {
        mData[i] *= v.mData[i];
    }
    return *this;
}

DoubleVector& DoubleVector::operator *=(double scalar)
{
    for (unsigned int i=0; i<mLength; i++)
    {
        mData[i] *= scalar;
    }
    return *this;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DoubleVector operator +(DoubleVector v1, const DoubleVector& v2)
{
    if (v1.mLength != v2.mLength)
    {
        fprintf(stderr, "DoubleVector operator +(DoubleVector v1, const DoubleVector& v2) v1:%d v2:%d\n",
                v1.mLength, v2.mLength);
        fflush(stderr);

        throw DoubleMatrixException(1);
    }

    unsigned int length = v1.mLength;
    for (unsigned int i=0; i<length; i++)
    {
        v1.mData[i] += v2.mData[i];
    }
    return v1;
}

DoubleVector operator -(DoubleVector v1, const DoubleVector& v2)
{
    if (v1.mLength != v2.mLength)
    {
        fprintf(stderr, "DoubleVector operator -(DoubleVector v1, const DoubleVector& v2) v1:%d v2:%d\n",
                v1.mLength, v2.mLength);
        fflush(stderr);

        throw DoubleMatrixException(1);
    }

    unsigned int length = v1.mLength;
    for (unsigned int i=0; i<length; i++)
    {
        v1.mData[i] -= v2.mData[i];
    }
    return v1;
}

DoubleVector operator *(DoubleVector v1, const DoubleVector& v2)
{
    if (v1.mLength != v2.mLength)
    {
        fprintf(stderr, "DoubleVector operator *(DoubleVector v1, const DoubleVector& v2) v1:%d v2:%d\n",
                v1.mLength, v2.mLength);
        fflush(stderr);
        throw DoubleMatrixException(3);
    }

    for (unsigned int i=0; i<v1.mLength; i++)
    {
        v1[i] *= v2[i];
    }
    return v1;
}

DoubleVector operator *(double scalar, DoubleVector v)
{
    for (unsigned int i=0; i<v.mLength; i++)
        v.mData[i] *= scalar;
    return v;
}

DoubleVector operator *(const DoubleVector v, double scalar)
{
    for (unsigned int i=0; i<v.mLength; i++)
        v.mData[i] *= scalar;
    return v;
}

bool operator ==(const DoubleVector& v1, const DoubleVector& v2)
{
    if (v1.mLength != v2.mLength) return false;

    unsigned int length = v1.mLength;
    for (unsigned int i=0; i<length; i++) if (v1.mData[i] != v2.mData[i]) return false;

    return true;
}

bool operator !=(const DoubleVector& v1, const DoubleVector& v2)
{
    return !(v1 == v2);
}

DoubleVector& DoubleVector::operator <<(double value)
{
    if (mData == NULL)
    {
        mLength = 1;
        mData = (double*) malloc(sizeof(double)*mLength);
        mData[0] = value;
    }
    else
    {
        mLength++;
        mData = (double*)realloc(mData, sizeof(double)*mLength);
        mData[mLength-1] = value;
    }
    return *this;
}

//std::ostream& operator <<(std::ostream& os, const DoubleVector& v)
//{
//    for (unsigned int i=0; i<v.mLength; i++) os << v.mData[i] << " ";
//    return os;
//}
