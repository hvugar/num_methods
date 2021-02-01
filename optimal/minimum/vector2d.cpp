#include "vector2d.h"
#include "matrix2d.h"
#include "exceptions.h"
#include <time.h>
#include <math.h>
#include <float.h>

DoubleVector::DoubleVector(size_t length, double val) : mLength(length), mData(nullptr)
{
    if (length == 0) return;

    mData = static_cast<double*>(malloc(sizeof(double)*length));
    for (size_t i=0; i<length; i++) mData[i] = val;
}

DoubleVector::DoubleVector(const double* data, size_t length) : mLength(length), mData(nullptr)
{
    if (length == 0 || data == nullptr) return;

    mData = static_cast<double*>(malloc(sizeof(double)*mLength));
    memcpy(mData, data, sizeof(double)*mLength);
}

DoubleVector::DoubleVector(const DoubleVector &vector) : mLength(vector.mLength), mData(nullptr)
{
    if (vector.mLength == 0) return;

    mData = static_cast<double*>(malloc(sizeof(double)*mLength));
    memcpy(mData, vector.mData, sizeof(double)*mLength);
}

DoubleVector::DoubleVector(const DoubleMatrix &matrix) : mLength(0), mData(nullptr)
{
    if (matrix.cols() != 1) return;

    mLength = matrix.rows();
    mData = static_cast<double*>(malloc(sizeof(double)*mLength));
    for (size_t i=0; i<mLength; i++) mData[i] = matrix.at(i,0);
}

DoubleVector::~DoubleVector()
{
    clear();
}

void DoubleVector::clear()
{
    if (mData != nullptr)
    {
        free(mData);
        mData = nullptr;
        mLength = 0;
    }
}

void DoubleVector::resize(size_t length, double value)
{
    if (length == mLength) return;

    if (length == 0) clear();

    if (length > 0)
    {
        if (mData == nullptr)
        {
            mData = static_cast<double*>(malloc(sizeof(double)*length));
            for (size_t i=0; i<length; i++) mData[i] = value;
            mLength = length;
        }
        else if (length != mLength)
        {
            double *ptr = static_cast<double*>(realloc(mData, sizeof(double) * length));
            for (size_t i=mLength; i<length; i++) ptr[i] = value;
            mData = ptr;
            mLength = length;
        }
    }
}

bool DoubleVector::empty() const
{
    return mLength == 0;
}

double& DoubleVector::at(size_t n)
{
    if (n >= mLength) { throw double_vector_exception(2); }
    return mData[n];
}

const double& DoubleVector::at(size_t n) const
{
    if (n >= mLength) { throw double_vector_exception(2); }
    return mData[n];
}

size_t DoubleVector::length() const
{
    return mLength;
}

DoubleVector& DoubleVector::append(const double *data, size_t length)
{
    if (length == 0 || data == nullptr) return *this;

    if (mLength == 0)
    {
        mLength = length;
        mData = static_cast<double*>(malloc(sizeof(double)*mLength));
        memcpy(mData, data, sizeof(double)*mLength);
    }
    else
    {
        mData = static_cast<double*>(realloc(mData, sizeof(double)*(mLength+length)));
        memcpy(mData+mLength, data, sizeof(double)*length);
        mLength += length;
    }

    return *this;
}

DoubleVector& DoubleVector::append(const DoubleVector &v)
{
    this->append(v.data(), v.length());
    return *this;
}

/********************************************************************
 *                      NORM
 *******************************************************************/

//double DoubleVector::norm(IVectorNormalizer::Norm norm) const
//{
//    switch (norm) {
//    case IVectorNormalizer::EUCLIDEAN_NORM: IVectorNormalizer::EuclideanNorm(this); break;
//    case IVectorNormalizer::L1_NORM: IVectorNormalizer::L1Norm(this); break;
//    case IVectorNormalizer::L2_NORM: IVectorNormalizer::L2Norm(this); break;
//    case IVectorNormalizer::LInf_NORM: IVectorNormalizer::LInfNorm(this); break;
//    default:
//        break;
//    }

//    };
//    return IVectorNormalizer::
//}

double DoubleVector::EuclideanNorm() const
{
    double norm = 0.0;
    for (size_t i=0; i<mLength; i++)
    {
        double item = mData[i];
        norm += item*item;
    }
    return sqrt(norm);
}

double DoubleVector::L1Norm() const
{
    double norm = 0.0;
    for (size_t i=0; i<mLength; i++)
    {
        norm += fabs(mData[i]);
    }
    return norm;
}

double DoubleVector::L2Norm() const
{
    double norm = 0.0;
    for (size_t i=0; i<mLength; i++)
    {
        double item = mData[i];
        norm += item*item;
    }
    return sqrt(norm);
}

double DoubleVector::LpNorm(size_t p) const
{
    double norm = 0.0;
    if (!empty())
    {
        for (size_t i=0; i < mLength; i++)
        {
            double item = mData[i];
            norm += pow(item, static_cast<double>(p));
        }
        return pow(norm, 1.0/static_cast<double>(p));
    }
    return norm;
}

double DoubleVector::LInfNorm() const
{
    double norm = 0.0;
    if (!empty())
    {
        norm = fabs(mData[0]);
        for (size_t i=1; i < mLength; i++)
        {
            if (norm < fabs(mData[i])) norm = fabs(mData[i]);
        }
    }
    return norm;
}

DoubleVector& DoubleVector::EuclideanNormalize()
{
    double norm = EuclideanNorm();
    if (norm > DBL_EPSILON) for (size_t i=0; i<mLength; i++) mData[i] /= norm;
    return *this;
}


DoubleVector& DoubleVector::L1Normalize()
{
    double norm = L1Norm();
    if (norm > DBL_EPSILON) for (size_t i=0; i<mLength; i++) mData[i] /= norm;
    return *this;
}

DoubleVector& DoubleVector::L2Normalize()
{
    double norm = L2Norm();
    if (norm > DBL_EPSILON) for (size_t i=0; i<mLength; i++) mData[i] /= norm;
    return *this;
}

/********************************************************************
 *                      NORM
 *******************************************************************/

double DoubleVector::min() const
{
    double minimum = static_cast<double>(NAN);
    if (!empty())
    {
        minimum = mData[0];
        for (size_t i=1; i<mLength; i++) if (minimum > mData[i]) minimum = mData[i];
    }
    return minimum;
}

double DoubleVector::max() const
{
    double maximum = static_cast<double>(NAN);
    if (!empty())
    {
        maximum = mData[0];
        for (size_t i=1; i<mLength; i++) if (maximum < mData[i]) maximum = mData[i];
    }
    return maximum;
}

DoubleVector DoubleVector::mid(size_t s, size_t e) const
{
    DoubleVector vector((e-s)+1);
    //unsigned char *src = reinterpret_cast<unsigned char*>(mData) + s;
    //unsigned char *dst = reinterpret_cast<unsigned char*>(vector.mData);
    //std::memcpy(dst, src, static_cast<size_t>((e-s)+1));
    for (size_t i=s; i<=e; i++) vector.mData[i-s] = mData[i];
    return vector;
}

double DoubleVector::EuclideanDistance(const DoubleVector &p) const
{
    if ( mLength != p.mLength ) return static_cast<double>(INFINITY);

    double distance = 0.0;
    for (size_t i=0; i<mLength; i++)
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

double& DoubleVector::operator [](size_t n)
{
    return mData[n];
}

double DoubleVector::operator [](size_t n) const
{
    return mData[n];
}

DoubleVector& DoubleVector::operator =(const DoubleVector& other)
{
    // the vector object holds reusable storage, such as a heap-allocated buffer mData

    if (this == &other) return *this; // self-assignment check expected

    clear();
    mLength = other.mLength;
    mData = static_cast<double*>(malloc(sizeof(double)*mLength));
    memcpy(mData, other.mData, sizeof(double)*mLength);
    return *this;
}

DoubleVector& DoubleVector::operator +=(const DoubleVector& v)
{
    if (mLength != v.mLength) { throw double_vector_exception(1); }

    for (size_t i=0; i < mLength; i++)
    {
        mData[i] += v.mData[i];
    }
    return *this;
}

DoubleVector& DoubleVector::operator -=(const DoubleVector& v)
{
    if (mLength != v.mLength) { throw double_vector_exception(1); }

    for (size_t i=0; i < mLength; i++)
    {
        mData[i] -= v.mData[i];
    }
    return *this;
}

DoubleVector& DoubleVector::operator *=(const DoubleVector& v)
{
    if (mLength != v.mLength) { throw double_vector_exception(1); }

    for (size_t i=0; i<v.mLength; i++)
    {
        mData[i] *= v.mData[i];
    }
    return *this;
}

DoubleVector& DoubleVector::operator *=(double scalar)
{
    for (size_t i=0; i<mLength; i++)
    {
        mData[i] *= scalar;
    }
    return *this;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DoubleVector operator +(DoubleVector v1, const DoubleVector& v2)
{
    if (v1.mLength != v2.mLength) { throw double_vector_exception(1); }

    size_t length = v1.mLength;
    for (size_t i=0; i<length; i++)
    {
        v1.mData[i] += v2.mData[i];
    }
    return v1;
}

DoubleVector operator -(DoubleVector v1, const DoubleVector& v2)
{
    if (v1.mLength != v2.mLength) { throw double_vector_exception(1); }

    size_t length = v1.mLength;
    for (size_t i=0; i<length; i++)
    {
        v1.mData[i] -= v2.mData[i];
    }
    return v1;
}

DoubleVector operator *(DoubleVector v1, const DoubleVector& v2)
{
    if (v1.mLength != v2.mLength) { throw double_vector_exception(1); }

    for (size_t i=0; i<v1.mLength; i++)
    {
        v1[i] *= v2[i];
    }
    return v1;
}

DoubleVector operator *(double scalar, DoubleVector v)
{
    for (size_t i=0; i<v.mLength; i++) v.mData[i] *= scalar;
    return v;
}

DoubleVector operator *(const DoubleVector v, double scalar)
{
    for (size_t i=0; i<v.mLength; i++) v.mData[i] *= scalar;
    return v;
}

bool operator ==(const DoubleVector& v1, const DoubleVector& v2)
{
    if (v1.mLength != v2.mLength) return false;

    size_t length = v1.mLength;
    for (size_t i=0; i<length; i++)
    {
        if (fabs(v1.mData[i] - v2.mData[i]) >= DBL_EPSILON) return false;
    }

    return true;
}

bool operator !=(const DoubleVector& v1, const DoubleVector& v2)
{
    return !(v1 == v2);
}

DoubleVector& DoubleVector::operator <<(double value)
{
    if (mData == nullptr)
    {
        mLength = 1L;
        mData = static_cast<double*>(malloc(sizeof(double)*mLength));
        mData[0] = value;
    }
    else
    {
        mLength++;
        mData = static_cast<double*>(realloc(mData, sizeof(double)*mLength));
        mData[mLength-1] = value;
    }
    return *this;
}
