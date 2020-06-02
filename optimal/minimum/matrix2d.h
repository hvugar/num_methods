#ifndef MATRIX2D_H
#define MATRIX2D_H

#include "global.h"
#include "vector2d.h"

class MINIMUMSHARED_EXPORT DoubleMatrix
{
public:
    static DoubleMatrix DiagonalMatrix(const DoubleVector& vector);
    static DoubleMatrix IdentityMatrix(size_t n);
    static DoubleMatrix ZeroMatrix(size_t rows, size_t cols);
    static DoubleMatrix HilbertMatrix(size_t rows, size_t cols);

    explicit DoubleMatrix(size_t rows=0, size_t cols=0, double value=0.0);
    DoubleMatrix(const DoubleMatrix &matrix);
    DoubleMatrix(const DoubleVector &vector);
    DoubleMatrix(const double* const data, size_t length, size_t rows, size_t cols);
    DoubleMatrix(double** data, size_t rows, size_t cols);
    virtual ~DoubleMatrix();

    size_t rows() const;
    size_t cols() const;

    bool empty() const;
    void clear();
    void resize(size_t rows, size_t cols, double value=0.0);
    void reset(double value=0.0);

    double& at(size_t row, size_t col);
    const double& at(size_t row, size_t col) const;

    DoubleVector row(size_t r) const;
    DoubleVector col(size_t c) const;

    void setColumn(size_t c, const DoubleVector& col);
    void setRow(size_t r, const DoubleVector& row);

    bool dioqonalMatrix() const;
    bool identityMatrix() const;
    bool zeroMatrix() const;
    bool squareMatrix() const;

    /*************************************************************************************
     *                        Basic matrix operations
     * **********************************************************************************/

    void Transpose();
    void ConjugateTranspose();
    void Inverse();
    void Det();
    void Minors();
    void Tr();
    void MatrixRank();

    /*************************************************************************************
     *                        Matrix operations
     * **********************************************************************************/

    DoubleVector Diagonal();
    void LowerTriangularize();
    void UpperTriangularize();
    void LUDecomposition();
    void Band();

    double trace() const;
    double det() const;
    DoubleVector eigenValues() const;

    bool equals(const DoubleMatrix& matrix) const;
    bool dimEquals(const DoubleMatrix& matrix) const;

    double min() const;
    double max() const;

    double determinant2() const;
    double determinant() const;
    void transpose();
    void inverse();
    DoubleMatrix minor(size_t row, size_t col) const;

    double** data() const;
    double** data();

    void switchRows(size_t row1, size_t row2);
    void switchCols(size_t col1, size_t col2);

    double* operator [](size_t row) const;
    double* operator [](size_t row);

    double& operator ()(size_t row, size_t col);
    const double& operator ()(size_t row, size_t col) const;

    DoubleMatrix& operator =(const DoubleMatrix& m);
    DoubleMatrix& operator =(const DoubleVector& m);
    DoubleMatrix& operator +=(const DoubleMatrix& m);
    DoubleMatrix& operator -=(const DoubleMatrix& m);
    DoubleMatrix& operator *=(const DoubleMatrix& m);
    DoubleMatrix& operator *=(double scalar);
    DoubleMatrix& operator *=(const DoubleVector& v);

    friend MINIMUMSHARED_EXPORT DoubleMatrix operator +(DoubleMatrix m1, const DoubleMatrix& m2);
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator -(DoubleMatrix m1, const DoubleMatrix& m2);
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator *(const DoubleMatrix& matrix1, const DoubleMatrix& matrix2);
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator *(double scalar, DoubleMatrix m);
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator *(DoubleMatrix m, double scalar);
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator *(DoubleMatrix m, const DoubleVector& v);

    friend MINIMUMSHARED_EXPORT bool operator ==(const DoubleMatrix& matrix1, const DoubleMatrix& matrix2);
    friend MINIMUMSHARED_EXPORT bool operator !=(const DoubleMatrix& matrix1, const DoubleMatrix& matrix2);

    friend MINIMUMSHARED_EXPORT DoubleMatrix operator ~(const DoubleMatrix&); // transose matrix
    friend MINIMUMSHARED_EXPORT DoubleMatrix operator !(const DoubleMatrix&); // inverse

    friend class DoubleVector;
private:
    size_t mRows;
    size_t mCols;
    double **mData;
};

#endif // MATRIX2D_H
