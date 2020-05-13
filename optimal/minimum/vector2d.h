#ifndef VECTOR2D_H
#define VECTOR2D_H

#include "global.h"

class DoubleMatrix;
class IVectorNormalizer;

class MINIMUMSHARED_EXPORT DoubleVector
{
public:
    explicit DoubleVector(unsigned int length = 0, double value = 0.0);
    explicit DoubleVector(const double* data, unsigned int length);
    DoubleVector(const DoubleVector &vector);
    DoubleVector(const DoubleMatrix &matrix);
    virtual ~DoubleVector();

    void clear();
    void resize(unsigned int length, double value = 0.0);
    bool empty() const;
    double& at (unsigned int n);
    const double& at (unsigned int n) const;
    unsigned int length() const;
    DoubleVector& append(const double *data, unsigned int length);
    DoubleVector& append(const DoubleVector &v);

    /********************************************************************
     *                               NORM
     ********************************************************************/

    double EuclideanNorm() const;
    double L1Norm() const;
    double L2Norm() const;
    double LpNorm(unsigned int p) const;
    double LInfNorm() const;

    DoubleVector& EuclideanNormalize();
    DoubleVector& L1Normalize();
    DoubleVector& L2Normalize();

    /********************************************************************
     *                               NORM
     ********************************************************************/

    double min() const;
    double max() const;
    DoubleVector mid(unsigned int s, unsigned int e) const;

    double EuclideanDistance(const DoubleVector&) const;

    double* data() NOEXCEPT;
    const double* data() const NOEXCEPT;

    /********************************************************************
     *                             OPERATORS
     ********************************************************************/

    double& operator [](unsigned int n);
    double operator [](unsigned int n) const;

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
    unsigned int mLength;
    double *mData;
    static unsigned int instanceCount;
};

#endif // VECTOR2D_H
