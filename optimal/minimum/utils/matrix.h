#ifndef MATRIX_H
#define MATRIX_H

#include <global.h>

template <typename T>
class MINIMUMSHARED_EXPORT Matrix
{
public:
    Matrix(unsigned int rows=0, unsigned int cols=0);
    Matrix(const Matrix& matrix);
    //Matrix(const Vector<T>& vector);
    virtual ~Matrix();

    unsigned int rows() const;
    unsigned int cols() const;

    bool empty() const;
    void clear();
    void resize(unsigned int rows, unsigned int cols, double value=0.0);

    T& operator()(unsigned int row, unsigned int col);
    const T& operator()(unsigned int row, unsigned int col) const;

    T& at(unsigned int row, unsigned int col);
    const T& at(unsigned int row, unsigned int col) const;

    T* operator [](unsigned int row) const;
    T* operator [](unsigned int row);

    //Vector<T> row(unsigned int r) const;
    //Vector<T> col(unsigned int c) const;

    //void setColumn(unsigned int c, const Vector<T>& col);
    //void setRow(unsigned int r, const Vector<T>& row);

    Matrix<T>& operator=(const Matrix<T> &matrix);
    //Matrix<T>& operator=(const Vector<T> &vector);

    bool equals(const Matrix<T>& matrix) const;
    bool dimEquals(const Matrix<T>& matrix) const;

    double min() const;
    double max() const;

    T** data() const;
    T** data();

private:
    unsigned int mRows;
    unsigned int mCols;
    double **mData;
};

#endif // MATRIX_H
