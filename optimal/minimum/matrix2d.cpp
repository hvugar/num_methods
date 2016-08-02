#include "matrix2d.h"
#include "exceptions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <exception>

MatrixException::MatrixException(int messageType) : messageType(messageType)
{}

const char* MatrixException::what() const noexcept
{
    return "Message!";
}

DoubleMatrix::DoubleMatrix(unsigned int rows, unsigned int cols, double value) : mRows(rows), mCols(cols), mData(NULL)
{
    if (rows > 0 && cols > 0)
    {
        mData = (double**)(malloc(sizeof(double*)*rows));
        for (unsigned int i=0; i<rows; i++)
        {
            mData[i] = (double*)malloc(sizeof(double)*cols);
            for (unsigned int j=0; j<cols; j++) mData[i][j] = value;
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

unsigned int DoubleMatrix::size() const
{
    return rows();
}

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

DoubleVector DoubleMatrix::operator [](unsigned int row) const
{
    DoubleVector vector(mCols);
    for (unsigned int i=0; i<mCols; i++) vector[i] = mData[row][i];
    return vector;
}

DoubleMatrix& DoubleMatrix::operator =(const DoubleMatrix &matrix)
{
    if (this != &matrix)
    {
        this->resize(matrix.rows(), matrix.cols());
        for (unsigned int j=0; j<mRows; j++)
        {
            for (unsigned int i=0; i<mCols; i++)
            {
                mData[j][i] = matrix.mData[j][i];
            }
        }
    }
    return *this;
}

bool DoubleMatrix::dimEquals(const DoubleMatrix &matrix) const
{
    return (rows() == matrix.rows() && cols() == matrix.cols());
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
            mData[j][i] = rand() % 100;
        }
    }
}
