#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <global.h>
#include <vector2d.h>
#include <vector>

class MINIMUMSHARED_EXPORT DoubleMatrix
{
public:
    explicit DoubleMatrix(unsigned int rows=0, unsigned int cols=0, double value=0.0);
    virtual ~DoubleMatrix();

    unsigned int rows() const;
    unsigned int cols() const;
    unsigned int size() const;

    void clear();
    void resize(unsigned int rows, unsigned int cols, double value=0.0);

    double& operator()(unsigned int row, unsigned int col);
    const double& operator()(unsigned int row, unsigned int col) const;

    double& at(unsigned int row, unsigned int col);
    const double& at(unsigned int row, unsigned int col) const;

    DoubleVector operator [](unsigned int row) const;

    DoubleVector row(unsigned int r) const;
    DoubleVector col(unsigned int c) const;

    DoubleMatrix& operator=(const DoubleMatrix &matrix);
    DoubleMatrix& operator+(const DoubleMatrix &matrix);
    DoubleMatrix& operator-(const DoubleMatrix &matrix);
    DoubleMatrix& operator*(const DoubleMatrix &matrix);
    DoubleMatrix& operator/(const DoubleMatrix &matrix);
    DoubleMatrix& operator*(const double scalar);
    DoubleMatrix& operator*(const DoubleVector &vector);

    void print();
    void randomData();

    bool equals(const DoubleMatrix& matrix) const;
    bool dimEquals(const DoubleMatrix& matrix) const;

    double min() const { return 0.0; }
    double max() const { return 0.0; }

    double determinant();
    void transpose();
    void inverse();
    DoubleMatrix minor(unsigned int row, unsigned int col);

    double** data() const;
    double** data();


private:
    unsigned int mRows;
    unsigned int mCols;
    double **mData;
};

#endif // MATRIX2D_H
