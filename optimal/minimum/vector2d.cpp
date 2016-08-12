#include "vector2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <exception>

DoubleVector::DoubleVector(unsigned int size, double val) : mSize(size), mData(NULL)
{
    if (size > 0)
    {
        mSize = size;
        mData = (double*) malloc(sizeof(double) * size);
        for (unsigned int i=0; i<size; i++) mData[i] = val;
        //        memset(mData, val, size);
    }
}

DoubleVector::DoubleVector(const double *data, unsigned int size) : mSize(size), mData(NULL)
{
    if (size > 0)
    {
        mSize = size;
        mData = (double*) malloc(sizeof(double) * size);
        memcpy(mData, data, sizeof(double)*size);
    }
}

DoubleVector::DoubleVector(const DoubleVector &vector)
{
    if (vector.mSize == 0) return;

    mSize = vector.mSize;
    mData = (double*) (malloc(sizeof(double)*mSize));
    memcpy(mData, vector.mData, sizeof(double)*mSize);
}

DoubleVector::~DoubleVector()
{
    clear();
}

void DoubleVector::clear()
{
    if (mData!=NULL)
    {
        free(mData);
        mData = NULL;
        mSize = 0;
    }
}

void DoubleVector::resize(unsigned int size, double value)
{
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
    return mData[n];
}

const double& DoubleVector::at(unsigned int n) const
{
    return mData[n];
}

unsigned int DoubleVector::size() const
{
    return mSize;
}

double* DoubleVector::data() noexcept
{
    return mData;
}

const double* DoubleVector::data() const noexcept
{
    return mData;
}

double DoubleVector::L2Norm() const
{
    double norm = 0.0;
    if (!empty())
    {
        for (unsigned int i=0; i<mSize; i++)
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
        for (unsigned int i=0; i<mSize; i++)
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
        for (unsigned int i=0; i<mSize; i++)
        {
            if ( i == 0 )
                norm = fabs(mData[i]);
            else if (norm < fabs(mData[i]))
                norm = fabs(mData[i]);
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
        distance += (at(i)-p.at(i))*(at(i)-p.at(i));
    }
    return sqrt(distance);
}

void DoubleVector::L2Normalize()
{
    double norm = L2Norm();
    if (norm != 0.0) for (unsigned int i=0; i<mSize; i++) mData[i] /= norm;
}

void DoubleVector::L1Normalize()
{
    double norm = L1Norm();
    for (unsigned int i=0; i<mSize; i++) mData[i] /= norm;
}

void DoubleVector::EuclideanNormalize()
{
    L2Normalize();
}

double DoubleVector::min() const
{
    double _min = NAN;
    if (!empty())
    {
        _min = at(0);
        for (unsigned int i=1; i<mSize; i++)
            if (_min > mData[i]) _min = mData[i];
    }
    return _min;
}

double DoubleVector::max() const
{
    double _max = 0.0;
    if (!empty())
    {
        _max = mData[0];
        for (unsigned int i=1; i<mSize; i++)
            if (_max < mData[i]) _max = mData[i];
    }
    return _max;
}

DoubleVector DoubleVector::mid(unsigned int s, unsigned int e) const
{
    DoubleVector vector((e-s)+1);
    for (unsigned int i=s; i<=e; i++) vector[i-s] = mData[i];
    return vector;
}

DoubleVector& DoubleVector::operator= (const DoubleVector& other)
{
    if (this != &other)
    {
        unsigned int size = other.size();
        resize(size, 0.0);
        for (unsigned int i=0; i<size; i++) mData[i] = other.mData[i];
    }
    return *this;
}

double& DoubleVector::operator[] (unsigned int n)
{
    return mData[n];
}

double DoubleVector::operator[] (unsigned int n) const
{
    return mData[n];
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
