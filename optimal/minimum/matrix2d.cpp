#include "matrix2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

DoubleMatrixException::DoubleMatrixException(unsigned int msgCode) noexcept : msgCode(msgCode) {}

DoubleMatrixException::DoubleMatrixException(const DoubleMatrixException & e) noexcept : msgCode(e.msgCode) {}

DoubleMatrixException& DoubleMatrixException::operator =(const DoubleMatrixException& e) noexcept
{
    msgCode = e.msgCode;
    return *this;
}
DoubleMatrixException::~DoubleMatrixException() {}

const char* DoubleMatrixException::what() const noexcept
{
    if (msgCode == 1) return "Dimension of matrixs do not matches!";
    if (msgCode == 2) return "Matrix is not square matrix!";
    if (msgCode == 3) return "Matrix columns and row are not equals!";
    return "";
}


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

    if (vector.size() > 0)
    {
        mRows = vector.size();
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
        for (unsigned int i=0; i<mRows; i++)
        {
            free(mData[i]);
            mData[i] = NULL;
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

    if (rows > 0 && cols > 0)
    {
        if (mData == NULL)
        {
            mData = (double**) malloc(sizeof(double*)*rows);
            for (unsigned int i=0; i<rows; i++)
            {
                mData[i] = (double*) malloc(sizeof(double)*cols);
                for (unsigned int j=0; j<cols; j++) mData[i][j] = value;
            }
            mRows = rows;
            mCols = cols;
        }
        else
        {
            if (rows != mRows)
            {
                double **ptr = (double **) realloc(mData, sizeof(double*) * rows);
                if (cols != mCols)
                {
                    for (unsigned int i=0; i<mRows; i++)
                    {
                        double *pRow = (double *) realloc(ptr[i], sizeof(double) * cols);
                        for (unsigned int j=mCols; j<cols; j++) pRow[j] = value;
                        ptr[i] = pRow;
                    }

                    for (unsigned int i=mRows; i<rows; i++)
                    {
                        double *pRow = (double *) realloc(ptr[i], sizeof(double) * cols);
                        for (unsigned int j=0; j<cols; j++) pRow[j] = value;
                        ptr[i] = pRow;
                    }

                    mCols = cols;
                }
                mRows = rows;
                mData = ptr;
            }
        }
    }
}

double& DoubleMatrix::operator()(unsigned int row, unsigned int col)
{
    return mData[row][col];
}

const double& DoubleMatrix::operator()(unsigned int row, unsigned int col) const
{
    return mData[row][col];
}

double& DoubleMatrix::at(unsigned int row, unsigned int col)
{
    return mData[row][col];
}

const double& DoubleMatrix::at(unsigned int row, unsigned int col) const
{
    return mData[row][col];
}

double* DoubleMatrix::operator [](unsigned int row) const
{
    return mData[row];
}

double* DoubleMatrix::operator [](unsigned int row)
{
    return mData[row];
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

DoubleMatrix& DoubleMatrix::operator= (const DoubleMatrix &m)
{
    if (this != &m)
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

DoubleMatrix& DoubleMatrix::operator= (const DoubleVector &v)
{
    puts("DoubleMatrix& DoubleMatrix::operator= (const DoubleVector &v)");
    if (v.size() > 0)
    {
        clear();
        mRows = v.size();
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
    double _min = DBL_MAX;
    for (unsigned int i=0; i<mRows; i++)
    {
        for (unsigned int j=0; j<mCols; j++)
        {
            if (mData[i][j] < _min) _min = mData[i][j];
        }
    }
    return _min;
}

double DoubleMatrix::max() const
{
    double _max = DBL_MIN;
    for (unsigned int i=0; i<mRows; i++)
    {
        for (unsigned int j=0; j<mCols; j++)
        {
            if (mData[i][j] > _max) _max = mData[i][j];
        }
    }
    return _max;
}

double DoubleMatrix::determinant() const
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
            det += (i%2==0 ? +1 : -1) * mData[0][i] * mnr.determinant();
        }
    }

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

//DoubleMatrix& DoubleMatrix::operator +(const DoubleMatrix &m)
//{
//    if (!dimEquals(m))
//    {
//        printf("DoubleMatrix& DoubleMatrix::operator +(const DoubleMatrix &matrix) %d %d %d %d\n", rows(), cols(), m.rows(), m.cols());
//        throw DoubleMatrixException(1);
//    }

//    for (unsigned int i=0; i<mRows; i++)
//    {
//        for (unsigned int j=0; j<mCols; j++)
//        {
//            mData[i][j] += m.mData[i][j];
//        }
//    }
//    return *this;
//}

DoubleMatrix operator+(const DoubleMatrix& m1, const DoubleMatrix& m2)
{
    if (!m1.dimEquals(m2))
    {
        printf("DoubleMatrix operator+(const DoubleMatrix& m1, const DoubleMatrix& m2) %d %d %d %d\n", m1.rows(), m1.cols(), m2.rows(), m2.cols());
        throw DoubleMatrixException(1);
    }

    DoubleMatrix m;
    m.resize(m1.rows(), m2.cols());
    for (unsigned int i=0; i<m.mRows; i++)
    {
        for (unsigned int j=0; j<m.mCols; j++)
        {
            m.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
        }
    }
    return m;
}

DoubleMatrix operator-(const DoubleMatrix& m1, const DoubleMatrix& m2)
{
    if (!m1.dimEquals(m2))
    {
        printf("DoubleMatrix operator+(const DoubleMatrix& m1, const DoubleMatrix& m2) %d %d %d %d\n", m1.rows(), m1.cols(), m2.rows(), m2.cols());
        throw DoubleMatrixException(1);
    }

    DoubleMatrix m;
    m.resize(m1.rows(), m2.cols());
    for (unsigned int i=0; i<m.mRows; i++)
    {
        for (unsigned int j=0; j<m.mCols; j++)
        {
            m.mData[i][j] = m1.mData[i][j]-m2.mData[i][j];
        }
    }
    return m;
}

//DoubleMatrix& DoubleMatrix::operator *(const DoubleMatrix &m)
//{
//    if (mCols != m.mRows)
//    {
//        printf("DoubleMatrix& DoubleMatrix::operator *(const DoubleMatrix &matrix) %d %d %d %d\n", rows(), cols(), m.rows(), m.cols());
//        throw DoubleMatrixException(3);
//    }


//    unsigned int rows1 = mRows;
//    unsigned int cols1 = m.mCols;
//    double** mdata = (double**)(malloc(sizeof(double*)*rows1));

//    for (unsigned int i=0; i<rows1; i++)
//    {
//        mdata[i] = (double*)malloc(sizeof(double)*cols1);
//    }

//    for (unsigned int i=0; i<mRows; i++)
//    {
//        for (unsigned int j=0; j<m.mCols; j++)
//        {
//            double sum = 0.0;
//            for (unsigned int k=0; k<mRows; k++)
//            {
//                sum += mData[i][k] * m.mData[k][j];
//            }
//            mdata[i][j] = sum;
//        }
//    }

//    clear();

//    this->mData = mdata;
//    this->mRows = rows1;
//    this->mCols = cols1;

//    return *this;
//}

DoubleMatrix operator*(const DoubleMatrix &m1, const DoubleMatrix &m2)
{
    if (m1.cols() != m2.rows())
    {
        printf("DoubleMatrix& DoubleMatrix::operator *(const DoubleMatrix &matrix) %d %d %d %d\n", m1.rows(), m1.cols(), m2.rows(), m2.cols());
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

void DoubleMatrix::print()
{
    for (unsigned int i=0; i<mRows; i++)
    {
        for (unsigned int j=0; j<mCols; j++)
        {
            printf("%10.6f ", mData[i][j]);
        }
        puts("");
    }
}

void DoubleMatrix::randomData()
{
    for (unsigned int i=0; i<mRows; i++)
    {
        for (unsigned int j=0; j<mCols; j++)
        {
            mData[i][j] = (rand() % 100000) * 0.00001;
        }
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

bool DoubleMatrix::isIdentityMatrix() const
{
    return false;
}

#include "printer.h"
void GaussianElimination(DoubleMatrix A, DoubleVector b, DoubleVector &x)
{
    const unsigned int ui = (unsigned)0-1;

    puts("....");
    IPrinter::print(A,12,12,18,14);

    unsigned int n = x.size();
    for (unsigned k=0; k < n-1; k++)
    {
        if (fabs(A.at(k,k)) <= DBL_EPSILON)
        {
            for (unsigned int p = k+1; p < n; p++)
            {
                if (fabs(A[k][p]) <= DBL_EPSILON)
                {
                    A.switchRows(k, p);
                    break;
                }
            }
        }

        for (unsigned int j=(k+1); j<n; j++)
        {
            double c = A.at(j,k)/A.at(k,k);
            for (unsigned int i=k; i<n; i++)
            {
                A.at(j,i) = A.at(j,i) - A.at(k,i) * c;
            }
            b[j] = b[j] - b[k] *c;
            puts("....");
            IPrinter::print(A,12,12,18,14);
            puts("....");
            int aaa;
            scanf("%d", &aaa);
        }
    }
    //puts("....");
    //IPrinter::print(A,12,12,18,14);
    puts("....");

    for (unsigned int i=(n-1); i!=ui; i--)
    {
        for (unsigned int j=(n-1); j>i; j--) b[i] -= (A.at(i,j) * x[j]);
        x[i] = b[i] / A.at(i,i);
    }
}
