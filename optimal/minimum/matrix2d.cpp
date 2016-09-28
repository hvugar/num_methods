#include "matrix2d.h"
#include "exceptions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <exception>

MatrixException::MatrixException(int messageType) : messageType(messageType)
{}

const char* MatrixException::what() const noexcept
{
    return "Message!";
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

DoubleMatrix::~DoubleMatrix()
{
    clear();
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

DoubleMatrix& DoubleMatrix::operator= (const DoubleMatrix &matrix)
{
    if (this != &matrix)
    {
        if (matrix.mRows > 0 && matrix.mCols > 0)
        {
            mRows = matrix.mRows;
            mCols = matrix.mCols;
            mData = (double**)malloc(sizeof(double*)*mRows);
            for (unsigned int i=0; i<mRows; i++)
            {
                mData[i] = (double *)malloc(sizeof(double)*mCols);
                memcpy(mData[i], matrix.mData[i], sizeof(double)*mCols);
            }
        }
        else
        {
            mRows = mCols = 0;
            mData = NULL;
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
    return 0.0;
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

DoubleMatrix DoubleMatrix::minor(unsigned int row, unsigned int col)
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

DoubleMatrix& DoubleMatrix::operator +(const DoubleMatrix &matrix)
{
    if (!dimEquals(matrix))
    {
        throw MatrixException(0);
    }
    else
    {
        for (unsigned int i=0; i<mRows; i++)
        {
            for (unsigned int j=0; j<mCols; j++)
            {
                mData[i][j] += matrix.mData[i][j];
            }
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator *(const DoubleMatrix &matrix)
{
    if (mCols == 0 && mCols != matrix.mRows) return *this;

    return *this;
}

DoubleMatrix mult(const DoubleMatrix &m1, const DoubleMatrix &m2)
{
    DoubleMatrix m;
    m.resize(m1.rows(), m2.cols());

    for (unsigned int i=0; i<m.rows(); i++)
    {
        for (unsigned int j=0; j<m.cols(); j++)
        {
            double sum = 0.0;
            for (unsigned int p=0; p<m1.cols(); p++) sum += m1.at(i,p)*m2.at(p,j);
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
            mData[i][j] = (rand() % 100) * 0.01;
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

void GaussianElimination(DoubleMatrix A, DoubleVector b, DoubleVector &x)
{
    const unsigned int ui = (unsigned)0-1;

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
            for (unsigned int i=k; i<n; i++) A.at(j,i) = A.at(j,i) - A.at(k,i) * c;
            b[j] = b[j] - b[k] *c;
        }
    }

    for (unsigned int i=(n-1); i!=ui; i--)
    {
        for (unsigned int j=(n-1); j>i; j--) b[i] -= (A.at(i,j) * x[j]);
        x[i] = b[i] / A.at(i,i);
    }
}
