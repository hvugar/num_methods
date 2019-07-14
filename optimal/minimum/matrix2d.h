#ifndef MATRIX2D_H
#define MATRIX2D_H

#include "global.h"
#include "vector2d.h"

class MINIMUMSHARED_EXPORT DoubleMatrix
{
public:
    static DoubleMatrix DiagonalMatrix(const DoubleVector& vector);
    static DoubleMatrix IdentityMatrix(unsigned int n);
    static DoubleMatrix ZeroMatrix(unsigned int rows, unsigned int cols);
    static DoubleMatrix HilbertMatrix(unsigned int rows, unsigned int cols);

    explicit DoubleMatrix(unsigned int rows=0, unsigned int cols=0, double value=0.0);
    DoubleMatrix(const DoubleMatrix &matrix);
    DoubleMatrix(const DoubleVector &vector);
    virtual ~DoubleMatrix();

    unsigned int rows() const;
    unsigned int cols() const;

    bool empty() const;
    void clear();
    void resize(unsigned int rows, unsigned int cols, double value=0.0);
    void reset(double value=0.0);

    double& at(unsigned int row, unsigned int col);
    const double& at(unsigned int row, unsigned int col) const;

    DoubleVector row(unsigned int r) const;
    DoubleVector col(unsigned int c) const;

    void setColumn(unsigned int c, const DoubleVector& col);
    void setRow(unsigned int r, const DoubleVector& row);

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
    DoubleMatrix minor(unsigned int row, unsigned int col) const;

    double** data() const;
    double** data();

    void switchRows(unsigned int row1, unsigned int row2);
    void switchCols(unsigned int col1, unsigned int col2);

    double* operator [](unsigned int row) const;
    double* operator [](unsigned int row);

    double& operator ()(unsigned int row, unsigned int col);
    const double& operator ()(unsigned int row, unsigned int col) const;

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
    unsigned int mRows;
    unsigned int mCols;
    double **mData;
};

#endif // MATRIX2D_H
