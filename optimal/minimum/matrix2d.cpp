#include "matrix2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include "printer.h"

DoubleMatrixException::DoubleMatrixException(unsigned int msgCode) NOEXCEPT : msgCode(msgCode) {}

DoubleMatrixException::DoubleMatrixException(const DoubleMatrixException & e) NOEXCEPT : msgCode(e.msgCode) {}

DoubleMatrixException& DoubleMatrixException::operator =(const DoubleMatrixException& e) NOEXCEPT
{
    msgCode = e.msgCode;
    return *this;
}

const char* DoubleMatrixException::what() const NOEXCEPT
{
    if (msgCode == 1) return "Dimension of matrixs do not matches!";
    if (msgCode == 2) return "Matrix is not square matrix!";
    if (msgCode == 3) return "Matrix columns and row are not equals!";
    return "";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix::DoubleMatrix(unsigned int rows, unsigned int cols, double value) : mRows(0), mCols(0), mData(NULL)
{
    if (rows > 0 && cols > 0)
    {
        mRows = rows;
        mCols = cols;
        mData = (double**)(malloc(sizeof(double*)*rows));
        for (unsigned int i=0; i<rows; i++)
        {
            mData[i] = (double*)malloc(sizeof(double)*cols);
            for (unsigned int j=0; j<cols; j++) mData[i][j] = value;
        }
    }
}

DoubleMatrix::DoubleMatrix(const DoubleMatrix &matrix) : mRows(0), mCols(0), mData(NULL)
{
    //puts("DoubleMatrix::DoubleMatrix(const DoubleMatrix &matrix)");
    if (matrix.mRows > 0 && matrix.mCols > 0)
    {
        mRows = matrix.mRows;
        mCols = matrix.mCols;
        mData = (double**) (malloc(sizeof(double*)*mRows));
        for (unsigned int i=0; i<mRows; i++)
        {
            mData[i] = (double*)malloc(sizeof(double)*mCols);
            memcpy(mData[i], matrix.mData[i], sizeof(double)*mCols);
        }
    }
}

DoubleMatrix::DoubleMatrix(const DoubleVector &vector) : mRows(0), mCols(0), mData(NULL)
{
    //puts("DoubleMatrix::DoubleMatrix(const DoubleVector &vector) : mRows(0), mCols(0), mData(NULL)");

    if (vector.length() > 0)
    {
        mRows = vector.length();
        mCols = 1;
        mData = (double**) malloc(sizeof(double*)*mRows);
        for (unsigned int i=0; i<mRows; i++)
        {
            mData[i] = (double*) malloc(sizeof(double)*mCols);
            mData[i][0] = vector.at(i);
        }
    }
}

DoubleMatrix::~DoubleMatrix()
{
    //puts("DoubleMatrix::~DoubleMatrix()1");
    clear();
    //puts("DoubleMatrix::~DoubleMatrix()2");
}

unsigned int DoubleMatrix::rows() const
{
    return mRows;
}

unsigned int DoubleMatrix::cols() const
{
    return mCols;
}

bool DoubleMatrix::empty() const
{
    return (mRows == 0 || mCols == 0 || mData == NULL);
}

void DoubleMatrix::clear()
{
    if (mData != NULL)
    {
        for (unsigned int row=0; row<mRows; row++)
        {
            free(mData[row]);
            mData[row] = NULL;
        }

        free(mData);
        mData = NULL;
        mRows = 0;
        mCols = 0;
    }
}

void DoubleMatrix::resize(unsigned int rows, unsigned int cols, double value)
{
    if (rows == 0 && cols == 0) clear();

    if (rows == 0 || cols == 0) return;

    if (rows > 0 && cols > 0)
    {
        if (mData == NULL)
        {
            mData = (double**) malloc(sizeof(double*)*rows);
            for (unsigned int row=0; row<rows; row++)
            {
                mData[row] = (double*) malloc(sizeof(double)*cols);
                for (unsigned int col=0; col<cols; col++) mData[row][col] = value;
            }
            mRows = rows;
            mCols = cols;
        }
        else
        {
            double **ptr = (double **) malloc(sizeof(double*) * rows);
            for (unsigned int row=0; row<rows; row++)
            {
                ptr[row] = (double*) malloc(sizeof(double)*cols);

                if (row < mRows)
                {
                    for (unsigned int col=0; col<cols; col++)
                        if (col < mCols) ptr[row][col] = mData[row][col]; else ptr[row][col] = value;
                }
                else
                {
                   for (unsigned int col=0; col<cols; col++)
                       ptr[row][col] = value;
                }
            }

            clear();

            mData = ptr;
            mRows = rows;
            mCols = cols;
        }
    }
}

double& DoubleMatrix::at(unsigned int row, unsigned int col)
{
    if (row>=mRows)
    {
        throw std::out_of_range("row index out of range");
    }
    if (col>=mCols)
    {
        throw std::out_of_range("column index out of range");
    }
    return mData[row][col];
}

const double& DoubleMatrix::at(unsigned int row, unsigned int col) const
{
    if (row>=mRows)
    {
        throw std::out_of_range("row index out of range");
    }
    if (col>=mCols)
    {
        throw std::out_of_range("column index out of range");
    }
    return mData[row][col];
}

DoubleVector DoubleMatrix::row(unsigned int r) const
{
    return DoubleVector(mData[r], mCols);
}

DoubleVector DoubleMatrix::col(unsigned int c) const
{
    DoubleVector v(mRows);
    for (unsigned int r=0; r<mRows; r++) v[r] = mData[r][c];
    return v;
}

void DoubleMatrix::setColumn(unsigned int c, const DoubleVector &col)
{
    if (mRows == col.length() && c < mCols)
    {
        for (unsigned int r=0; r<mRows; r++) mData[r][c] = col.at(r);
    }
    else
    {
        throw DoubleMatrixException(0);
    }
}

void DoubleMatrix::setRow(unsigned int r, const DoubleVector &row)
{
    if (mCols == row.length() && r < mRows)
    {
        for (unsigned int c=0; c<mCols; c++) mData[r][c] = row.at(c);
    }
    else
    {
        throw DoubleMatrixException(0);
    }
}

double** DoubleMatrix::data() const
{
    return mData;
}

double** DoubleMatrix::data()
{
    return mData;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix DoubleMatrix::DiagonalMatrix(const DoubleVector& v)
{
    unsigned int n = v.length();
    DoubleMatrix m(n, n, 0.0);
    for (unsigned int i=0; i<n; i++) m.mData[i][i] = v[i];
    return m;
}

DoubleMatrix DoubleMatrix::IdentityMatrix(unsigned int n)
{
    DoubleMatrix m(n, n, 0.0);
    for (unsigned int i=0; i<n; i++) m.mData[i][i] = 1.0;
    return m;
}

DoubleMatrix DoubleMatrix::ZeroMatrix(unsigned int rows, unsigned int cols)
{
    DoubleMatrix matrix(rows, cols, 0.0);
    return matrix;
}

DoubleMatrix DoubleMatrix::HilbertMatrix(unsigned int rows, unsigned int cols)
{
    DoubleMatrix matrix(rows, cols);
    for (unsigned int row=0; row<rows; row++)
        for (unsigned int col=0; col<cols; col++)
            matrix.mData[row][col] = 1.0/(double)(row+col+1);
    return matrix;
}

///////////////////////////////////////////////////////////////////////////////////////////////

bool DoubleMatrix::dioqonalMatrix() const
{
    for (unsigned int row=0; row < mRows; row++)
    {
        for (unsigned int col=0; col < mCols; col++)
        {
            if (row != col && fabs(mData[row][col]) > DBL_EPSILON) return false;
        }
    }
    return true;
}

bool DoubleMatrix::identityMatrix() const
{
    for (unsigned int row=0; row < mRows; row++)
    {
        for (unsigned int col=0; col < mCols; col++)
        {
            if (row != col && fabs(mData[row][col]) > DBL_EPSILON) return false;
            if (row == col && fabs(mData[row][col]) != 1.0) return false;
        }
    }
    return true;
}

bool DoubleMatrix::zeroMatrix() const
{
    for (unsigned int row=0; row < mRows; row++)
        for (unsigned int col=0; col < mCols; col++)
            if (fabs(mData[row][col]) > DBL_EPSILON) return false;
    return true;
}

bool DoubleMatrix::squareMatrix() const
{
    return mRows == mCols;
}

///////////////////////////////////////////////////////////////////////////////////////////////

bool DoubleMatrix::dimEquals(const DoubleMatrix &matrix) const
{
    return (mRows == matrix.mRows && mCols == matrix.mCols);
}

bool DoubleMatrix::equals(const DoubleMatrix &matrix) const
{
    if (dimEquals(matrix))
    {
        bool equals = true;
        for (unsigned int i=0; i<mRows; i++)
        {
            for (unsigned int j=0; j<mCols; j++)
            {
                if (mData[i][j] != matrix.mData[i][j])
                {
                    equals = false;
                    break;
                }
            }
            if (!equals) break;
        }
        return equals;
    }
    else
        return false;
}

double DoubleMatrix::min() const
{
    double minimum = DBL_MAX;
    for (unsigned int row=0; row<mRows; row++)
    {
        for (unsigned int col=0; col<mCols; col++)
        {
            if (mData[row][col] < minimum) minimum = mData[row][col];
        }
    }
    return minimum;
}

double DoubleMatrix::max() const
{
    double maximum = DBL_MIN;
    for (unsigned int row=0; row<mRows; row++)
    {
        for (unsigned int col=0; col<mCols; col++)
        {
            if (mData[row][col] > maximum) maximum = mData[row][col];
        }
    }
    return maximum;
}

double DoubleMatrix::determinant2() const
{
    double det = 0.0;

    if (mRows != mCols)
    {
        throw DoubleMatrixException(2);
    }

    if (mRows == 1 && mCols == 1)
        return mData[0][0];
    else if (mRows == 2 && mCols == 2)
        return mData[0][0]*mData[1][1] - mData[0][1]*mData[1][0];
    else
    {
        unsigned int i;
        for (i=0; i<mCols; i++)
        {
            DoubleMatrix mnr = minor(0, i);
            det += (i%2==0 ? +1 : -1) * mData[0][i] * mnr.determinant2();
        }
    }

    return det;
}

double DoubleMatrix::determinant() const
{
    double det = 0.0;

    if (mRows != mCols)
    {
        throw DoubleMatrixException(2);
    }

    DoubleMatrix T = *this;

    for (unsigned k=0; k < mRows; k++)
    {
        if (fabs(T.at(k,k)) <= DBL_EPSILON)
        {
            bool swiched = false;
            for (unsigned int p=k+1; p<mRows; p++)
            {
                if (fabs(T.at(p,k)) > DBL_EPSILON)
                {
                    T.switchRows(k, p);
                    swiched = true;
                    //break;
                }
                if (swiched==true) break;
            }
        }

        for (unsigned int j=k+1; j<mRows; j++)
        {
            double c = -(T.at(j,k)/T.at(k,k));
            for (unsigned int i=k; i<mCols; i++)
            {
                T.at(j,i) = T.at(j,i) + T.at(k,i) * c;
            }
        }
    }
    det = 1.0;
    for (unsigned int i=0; i<mRows; i++) det *= T.at(i,i);
    T.clear();

    return det;
}

void DoubleMatrix::transpose()
{
    if (empty()) return;

    unsigned int rows = mCols;
    unsigned int cols = mRows;
    double **data = (double**)(malloc(sizeof(double*)*rows));
    for (unsigned int i=0; i<rows; i++)
    {
        data[i] = (double*)malloc(sizeof(double)*cols);
        for (unsigned int j=0; j<cols; j++) data[i][j] = mData[j][i];
    }
    clear();
    mRows = rows;
    mCols = cols;
    mData = data;
}

void DoubleMatrix::inverse()
{
    double idet = 1.0/determinant();

    DoubleMatrix m = *this;

    for (unsigned int r=0; r<rows(); r++)
    {
        for (unsigned int c=0; c<cols(); c++)
        {
            DoubleMatrix minor = m.minor(r,c);
            mData[r][c] = minor.determinant();
            if ((r+c)%2==1) mData[r][c] *= -1.0;
            minor.clear();
        }
    }

    transpose();

    for (unsigned int r=0; r<rows(); r++)
    {
        for (unsigned int c=0; c<cols(); c++)
        {
            mData[r][c] *= idet;
        }
    }

    m.clear();
}

DoubleMatrix DoubleMatrix::minor(unsigned int row, unsigned int col) const
{
    DoubleMatrix m(mRows-1, mCols-1);
    for (unsigned int i=0; i<mRows; i++)
    {
        for (unsigned int j=0; j<mCols; j++)
        {
            if (i == row && j == col) continue;
            if (i < row  && j <  col) { m.mData[i][j]     = mData[i][j]; continue; }
            if (i < row  && j >  col) { m.mData[i][j-1]   = mData[i][j]; continue; }
            if (i > row  && j <  col) { m.mData[i-1][j]   = mData[i][j]; continue; }
            if (i > row  && j >  col) { m.mData[i-1][j-1] = mData[i][j]; continue; }
        }
    }
    return m;
}

///////////////////////////////////////////////////////////////////////////////////////////////

double* DoubleMatrix::operator [](unsigned int row) const
{
    return mData[row];
}

double* DoubleMatrix::operator [](unsigned int row)
{
    return mData[row];
}

double& DoubleMatrix::operator()(unsigned int row, unsigned int col)
{
    return mData[row][col];
}

const double& DoubleMatrix::operator()(unsigned int row, unsigned int col) const
{
    return mData[row][col];
}

///////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix& DoubleMatrix::operator =(const DoubleMatrix &m)
{
    // the matrix object holds reusable storage, such as a heap-allocated buffer mData

    if (this != &m) // self-assignment check expected
    {
        clear();

        if (m.mRows > 0 && m.mCols > 0)
        {
            mRows = m.mRows;
            mCols = m.mCols;
            mData = (double**) malloc(sizeof(double*)*mRows);
            for (unsigned int i=0; i<mRows; i++)
            {
                mData[i] = (double *) malloc(sizeof(double)*mCols);
                memcpy(mData[i], m.mData[i], sizeof(double)*mCols);
            }
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator =(const DoubleVector &v)
{
    if (v.length() > 0)
    {
        clear();
        mRows = v.length();
        mCols = 1;
        mData = (double**) malloc(sizeof(double*)*mRows);
        for (unsigned int i=0; i<mRows; i++)
        {
            mData[i] = (double*) malloc(sizeof(double)*mCols);
            mData[i][0] = v.at(i);
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator +=(const DoubleMatrix &m)
{
    if (!dimEquals(m))
    {
        throw std::out_of_range("dimension dont match.");
    }

    for (unsigned int rw=0; rw < mRows; rw++)
    {
        for (unsigned int cl=0; cl < mCols; cl++)
        {
            mData[rw][cl] += m.mData[rw][cl];
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator -=(const DoubleMatrix &m)
{
    if (!dimEquals(m))
    {
        throw std::out_of_range("dimension dont match.");
    }

    for (unsigned int rw=0; rw < mRows; rw++)
    {
        for (unsigned int cl=0; cl < mCols; cl++)
        {
            mData[rw][cl] -= m.mData[rw][cl];
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator *=(const DoubleMatrix &m)
{
    if (mCols != m.mRows)
    {
        printf("DoubleMatrix& DoubleMatrix::operator *=(const DoubleMatrix &m) this:%d %d m:%d %d\n",
               mRows, mCols, m.rows(), m.cols());
        throw DoubleMatrixException(3);
    }

    //    DoubleMatrix m;
    //    m.resize(m1.rows(), m2.cols());

    //    for (unsigned int i=0; i<m.rows(); i++)
    //    {
    //        for (unsigned int j=0; j<m.cols(); j++)
    //        {
    //            double sum = 0.0;
    //            for (unsigned int k=0; k<m1.cols(); k++) sum += m1.at(i,k)*m2.at(k,j);
    //            m.at(i,j) = sum;
    //        }
    //    }
    //    return m;

    //    for (unsigned int rw=0; rw < mRows; rw++)
    //    {
    //        for (unsigned int cl=0; cl < mCols; cl++)
    //        {
    //            mData[rw][cl] -= m.mData[rw][cl];
    //        }
    //    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator *=(double scalar)
{
    for (unsigned int row=0; row<mRows; row++)
    {
        for (unsigned int col=0; col<mCols; col++)
        {
            mData[row][col] *= scalar;
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator *=(const DoubleVector& v)
{
    if (mCols != v.mLength)
    {
        printf("DoubleMatrix& DoubleMatrix::operator *=(const DoubleVector& v) this:%d %d v:%d\n", mRows, mCols, v.mLength);
        throw DoubleMatrixException(3);
    }

    DoubleMatrix m = *this;

    resize(mRows, 1);
    for (unsigned int row=0; row<mRows; row++)
    {
        mData[row][0] = 0.0;
        for (unsigned int col=0; col<m.mCols; col++) mData[row][0] += m.mData[row][col]*v[col];
    }
    m.clear();

    return *this;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix operator +(DoubleMatrix m1, const DoubleMatrix& m2)
{
    if (!m1.dimEquals(m2))
    {
        printf("DoubleMatrix operator +(const DoubleMatrix& m1, const DoubleMatrix& m2) %d %d %d %d\n",
               m1.rows(), m1.cols(), m2.rows(), m2.cols());
        throw DoubleMatrixException(1);
    }

    for (unsigned int row=0; row<m1.mRows; row++)
    {
        for (unsigned int col=0; col<m1.mCols; col++)
        {
            m1.mData[row][col] += m2.mData[row][col];
        }
    }
    return m1;
}

DoubleMatrix operator -(DoubleMatrix m1, const DoubleMatrix& m2)
{
    if (!m1.dimEquals(m2))
    {
        printf("DoubleMatrix operator -(DoubleMatrix m1, const DoubleMatrix& m2) %d %d %d %d\n",
               m1.rows(), m1.cols(), m2.rows(), m2.cols());
        throw DoubleMatrixException(1);
    }

    for (unsigned int row=0; row<m1.mRows; row++)
    {
        for (unsigned int col=0; col<m1.mCols; col++)
        {
            m1.mData[row][col] -= m2.mData[row][col];
        }
    }
    return m1;
}

DoubleMatrix operator *(const DoubleMatrix &m1, const DoubleMatrix &m2)
{
    if (m1.cols() != m2.rows())
    {
        printf("DoubleMatrix& DoubleMatrix::operator *(const DoubleMatrix &matrix) %d %d %d %d\n",
               m1.rows(), m1.cols(), m2.rows(), m2.cols());
        throw DoubleMatrixException(3);
    }

    DoubleMatrix m;
    m.resize(m1.rows(), m2.cols());

    for (unsigned int i=0; i<m.rows(); i++)
    {
        for (unsigned int j=0; j<m.cols(); j++)
        {
            double sum = 0.0;
            for (unsigned int k=0; k<m1.cols(); k++) sum += m1.at(i,k)*m2.at(k,j);
            m.at(i,j) = sum;
        }
    }
    return m;
}

DoubleMatrix operator *(double scalar, DoubleMatrix m)
{
    unsigned int m_row = m.rows();
    unsigned int m_col = m.cols();

    for (unsigned int row=0; row < m_row; row++)
    {
        for (unsigned int col=0; col < m_col; col++)
        {
            m.mData[row][col] *= scalar;
        }
    }

    return m;
}

DoubleMatrix operator *(DoubleMatrix m, double scalar)
{
    return scalar*m;
}

DoubleMatrix operator *(DoubleMatrix m, const DoubleVector& v)
{
    if (m.mCols != v.length())
    {
        printf("DoubleMatrix DoubleMatrix::operator *=(const DoubleVector& v) this:%d %d v:%d\n", m.mRows, m.mCols, v.length());
        throw DoubleMatrixException(3);
    }

    DoubleMatrix m1(m.mRows, 1);

    for (unsigned int row=0; row<m.mRows; row++)
    {
        m1.mData[row][0] = 0.0;
        for (unsigned int col=0; col<m.mCols; col++) m1[row][0] += m.mData[row][col]*v[col];
    }

    return m1;
}

///////////////////////////////////////////////////////////////////////////////////////////////

bool operator ==(const DoubleMatrix &matrix1, const DoubleMatrix& matrix2)
{
    if (!matrix1.dimEquals(matrix2)) return false;

    for (unsigned int row = 0; row < matrix1.mRows; row++)
    {
        for (unsigned int col = 0; col < matrix1.mCols; col++)
        {
            if (matrix1.mData[row][col] != matrix2.mData[row][col]) return false;
        }
    }
    return true;
}

bool operator !=(const DoubleMatrix &matrix1, const DoubleMatrix& matrix2)
{
    return !(matrix1==matrix2);
}

///////////////////////////////////////////////////////////////////////////////////////////////

void DoubleMatrix::switchRows(unsigned int row1, unsigned int row2)
{
    double *row = (double*)malloc(sizeof(double) * mCols);
    memcpy(row, mData[row1], sizeof(double)*mCols);
    memcpy(mData[row1], mData[row2], sizeof(double)*mCols);
    memcpy(mData[row2], row, sizeof(double)*mCols);
    free(row);
}

void DoubleMatrix::switchCols(unsigned int col1 UNUSED_PARAM, unsigned int col2 UNUSED_PARAM)
{
}
