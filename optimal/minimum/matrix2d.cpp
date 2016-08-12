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
        for (unsigned int j=0; j<rows; j++)
        {
            mData[j] = (double*)malloc(sizeof(double)*cols);
            for (unsigned int i=0; i<cols; i++) mData[j][i] = value;
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
        for (unsigned int j=0; j<mRows; j++)
        {
            mData[j] = (double*)malloc(sizeof(double)*mCols);
            memcpy(mData[j], matrix.mData[j], sizeof(double)*mCols);
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

//unsigned int DoubleMatrix::size() const
//{
//    return rows();
//}

void DoubleMatrix::clear()
{
    if (mData != NULL)
    {
        for (unsigned int i=0; i < mRows; i++)
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
            for (unsigned int j=0; j<rows; j++)
            {
                mData[j] = (double*) malloc(sizeof(double)*cols);
                for (unsigned int i=0; i<cols; i++) mData[j][i] = value;
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
                    for (unsigned int j=0; j<mRows; j++)
                    {
                        double *pRow = (double *) realloc(ptr[j], sizeof(double) * cols);
                        for (unsigned int i=mCols; i<cols; i++) pRow[i] = value;
                        ptr[j] = pRow;
                    }

                    for (unsigned int j=mRows; j<rows; j++)
                    {
                        double *pRow = (double *) realloc(ptr[j], sizeof(double) * cols);
                        for (unsigned int i=0; i<cols; i++) pRow[i] = value;
                        ptr[j] = pRow;
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
            for (unsigned int j=0; j<mRows; j++)
            {
                mData[j] = (double *)malloc(sizeof(double)*mCols);
                memcpy(mData[j], matrix.mData[j], sizeof(double)*mCols);
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
        for (unsigned int j=0; j<mRows; j++)
        {
            for (unsigned int i=0; i<mCols; i++)
            {
                if (mData[j][i] != matrix.mData[j][i])
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
    for (unsigned int j=0; j<mRows; j++)
    {
        for (unsigned int i=0; i<mCols; i++)
        {
            if (mData[j][i] < _min) _min = mData[j][i];
        }
    }
    return _min;
}

double DoubleMatrix::max() const
{
    double _max = DBL_MIN;
    for (unsigned int j=0; j<mRows; j++)
    {
        for (unsigned int i=0; i<mCols; i++)
        {
            if (mData[j][i] > _max) _max = mData[j][i];
        }
    }
    return _max;
}


DoubleMatrix& DoubleMatrix::operator +(const DoubleMatrix &matrix)
{
    if (!dimEquals(matrix))
    {
        throw MatrixException(0);
    }
    else
    {
        for (unsigned int j=0; j<mRows; j++)
        {
            for (unsigned int i=0; i<mCols; i++)
            {
                mData[j][i] += matrix.mData[j][i];
            }
        }
    }
    return *this;
}

void DoubleMatrix::print()
{
    for (unsigned int j=0; j<mRows; j++)
    {
        for (unsigned int i=0; i<mCols; i++)
        {
            printf("%10.6f ", mData[j][i]);
        }
        puts("");
    }
}

void DoubleMatrix::randomData()
{
    for (unsigned int j=0; j<mRows; j++)
    {
        for (unsigned int i=0; i<mCols; i++)
        {
            mData[j][i] = (rand() % 100) * 0.01;
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

void DoubleMatrix::changeRows(unsigned int i, unsigned int j)
{
    double *row = (double*)malloc(sizeof(double) * mCols);
    memcpy(row, mData[i], sizeof(double)*mCols);
    memcpy(mData[i], mData[j], sizeof(double)*mCols);
    memcpy(mData[j], row, sizeof(double)*mCols);
    free(row);
}

void GaussianElimination(DoubleMatrix m, DoubleVector b, DoubleVector &x)
{
    const unsigned int ui = (unsigned)0-1;

    unsigned int n = x.size();
    for (unsigned k=0; k < n-1; k++)
    {
        if (fabs(m.at(k,k)) <= DBL_EPSILON)
        {
            for (unsigned int p = k+1; p < n; p++)
            {
                if (m[k][p] != 0.0)
                {
                    m.changeRows(k, p);
                    break;
                }
            }
        }

        for (unsigned int j=(k+1); j<n; j++)
        {
            double c = m.at(j,k)/m.at(k,k);
            for (unsigned int i=k; i<n; i++) m.at(j,i) = m.at(j,i) - m.at(k,i) * c;
            b[j] = b[j] - b[k] *c;
        }
    }

    for (unsigned int i=(n-1); i!=ui; i--)
    {
        for (unsigned int j=(n-1); j>i; j--) b[i] -= (m.at(i,j) * x[j]);
        x[i] = b[i] / m.at(i,i);
    }
}
