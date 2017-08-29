#ifndef MATRIX2D_H
#define MATRIX2D_H

#include "global.h"
#include "vector2d.h"
#include <exception>
#include <string>

class MINIMUMSHARED_EXPORT DoubleMatrixException : public std::exception
{
public:
    DoubleMatrixException(unsigned int msgCode = 0) NOEXCEPT;
    DoubleMatrixException(const DoubleMatrixException&) NOEXCEPT;
    DoubleMatrixException& operator= (const DoubleMatrixException&) NOEXCEPT;
    //virtual ~DoubleMatrixException()  {}
    virtual const char* what() const NOEXCEPT;

private:
    unsigned int msgCode;
};

class MINIMUMSHARED_EXPORT DoubleMatrix
{
public:
    explicit DoubleMatrix(unsigned int rows=0, unsigned int cols=0, double value=0.0);
    DoubleMatrix(const DoubleMatrix &matrix);
    DoubleMatrix(const DoubleVector &vector);
    virtual ~DoubleMatrix();

    unsigned int rows() const;
    unsigned int cols() const;

    bool empty() const;
    void clear();
    void resize(unsigned int rows, unsigned int cols, double value=0.0);

    double& operator()(unsigned int row, unsigned int col);
    const double& operator()(unsigned int row, unsigned int col) const;

    double& at(unsigned int row, unsigned int col);
    const double& at(unsigned int row, unsigned int col) const;

    double* operator [](unsigned int row) const;
    double* operator [](unsigned int row);

    DoubleVector row(unsigned int r) const;
    DoubleVector col(unsigned int c) const;

    void setColumn(unsigned int c, const DoubleVector& col);
    void setRow(unsigned int r, const DoubleVector& row);

    DoubleMatrix& operator =(const DoubleMatrix &m);
    DoubleMatrix& operator =(const DoubleVector &v);
    //DoubleMatrix& operator +(const DoubleMatrix &m);
    DoubleMatrix& operator +=(const DoubleMatrix &m);
    DoubleMatrix& operator -=(const DoubleMatrix &m);
    DoubleMatrix& operator *=(const DoubleMatrix &m);

    //DoubleMatrix& operator-(const DoubleMatrix &matrix);
    //DoubleMatrix& operator/(const DoubleMatrix &matrix);
    //DoubleMatrix& operator*(const double scalar);
    //DoubleMatrix& operator*(const DoubleVector &vector);
    //DoubleMatrix& operator*(const DoubleMatrix &matrix);

    void print();
    void randomData();

    bool equals(const DoubleMatrix& matrix) const;
    bool dimEquals(const DoubleMatrix& matrix) const;

    double min() const;
    double max() const;

    double determinant2() const;
    double determinant() const;
    void transpose();
    void inverse();
    DoubleMatrix minor(unsigned int row, unsigned int col) const;

    double** data() const;
    double** data();

    void switchRows(unsigned int row1, unsigned int row2);
    void switchCols(unsigned int col1, unsigned int col2);

    bool isIdentityMatrix() const;
    bool isDiagonalMatrix() const;
    bool isSquareMatrix() const;

    friend MINIMUMSHARED_EXPORT DoubleMatrix operator +(DoubleMatrix m1, const DoubleMatrix& m2);
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator -(DoubleMatrix m1, const DoubleMatrix& m2);
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator *(const DoubleMatrix&, const DoubleMatrix&);
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator *(double, DoubleMatrix);
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator *(DoubleMatrix, double);
//    friend MINIMUMSHARED_EXPORT DoubleMatrix operator *(const DoubleMatrix&, const DoubleVector&);
//    friend MINIMUMSHARED_EXPORT DoubleMatrix operator ~(const DoubleMatrix&);

private:
    unsigned int mRows;
    unsigned int mCols;
    double **mData;
};

void MINIMUMSHARED_EXPORT GaussianElimination(DoubleMatrix m, DoubleVector b, DoubleVector &x);

#endif // MATRIX2D_H
