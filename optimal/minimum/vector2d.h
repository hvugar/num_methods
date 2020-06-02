#ifndef VECTOR2D_H
#define VECTOR2D_H

#include "global.h"

class DoubleMatrix;
class IVectorNormalizer;

class MINIMUMSHARED_EXPORT DoubleVector
{
public:
    explicit DoubleVector(size_t length = 0, double value = 0.0);
    explicit DoubleVector(const double* data, size_t length);
    DoubleVector(const DoubleVector &vector);
    DoubleVector(const DoubleMatrix &matrix);
    virtual ~DoubleVector();

    void clear();
    void resize(size_t length, double value = 0.0);
    bool empty() const;
    double& at (size_t n);
    const double& at (size_t n) const;
    size_t length() const;
    DoubleVector& append(const double *data, size_t length);
    DoubleVector& append(const DoubleVector &v);

    /********************************************************************
     *                               NORM
     ********************************************************************/

    double EuclideanNorm() const;
    double L1Norm() const;
    double L2Norm() const;
    double LpNorm(size_t p) const;
    double LInfNorm() const;

    DoubleVector& EuclideanNormalize();
    DoubleVector& L1Normalize();
    DoubleVector& L2Normalize();

    /********************************************************************
     *                               NORM
     ********************************************************************/

    double min() const;
    double max() const;
    DoubleVector mid(size_t s, size_t e) const;

    double EuclideanDistance(const DoubleVector&) const;

    double* data() NOEXCEPT;
    const double* data() const NOEXCEPT;

    /********************************************************************
     *                             OPERATORS
     ********************************************************************/

    double& operator [](size_t n);
    double operator [](size_t n) const;

    DoubleVector& operator =(const DoubleVector& vector);
    DoubleVector& operator +=(const DoubleVector& vector);
    DoubleVector& operator -=(const DoubleVector& vector);
    DoubleVector& operator *=(const DoubleVector& vector);
    DoubleVector& operator *=(double scalar);

    ///////////////////////////////////////////////////////////////////////////////////////////////

    friend MINIMUMSHARED_EXPORT DoubleVector operator +(DoubleVector v1, const DoubleVector& v2);
    friend MINIMUMSHARED_EXPORT DoubleVector operator -(DoubleVector v1, const DoubleVector& v2);
    friend MINIMUMSHARED_EXPORT DoubleVector operator *(DoubleVector v1, const DoubleVector& v2);
    friend MINIMUMSHARED_EXPORT DoubleVector operator *(double scalar, DoubleVector v);
    friend MINIMUMSHARED_EXPORT DoubleVector operator *(DoubleVector v, double scalar);

    friend MINIMUMSHARED_EXPORT bool operator ==(const DoubleVector& v1, const DoubleVector& v2);
    friend MINIMUMSHARED_EXPORT bool operator !=(const DoubleVector& v1, const DoubleVector& v2);

    DoubleVector& operator <<(double value);
    //friend MINIMUMSHARED_EXPORT std::ostream& operator <<(std::ostream& os, const DoubleVector& v);

    friend class DoubleMatrix;
    friend class IVectorNormalizer;
private:
    size_t mLength;
    double *mData;
};

#endif // VECTOR2D_H
