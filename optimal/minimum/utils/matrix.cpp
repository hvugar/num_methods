#include "matrix.h"
#include <cstring>
#include <cstdlib>
#include <exception>

template <typename T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols) : mRows(rows), mCols(cols), mData(NULL)
{
    if (rows > 0 && cols > 0)
    {
        mRows = rows;
        mCols = cols;
        mData = (T**) malloc(sizeof(T*) * rows);
        for (unsigned int r=0; r<rows; r++) mData[r] = (T*) malloc(sizeof(T)*cols);
    }
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& matrix) : mRows(0), mCols(0), mData(NULL)
{
    if (matrix.mRows > 0 && matrix.mCols > 0)
    {
        mRows = matrix.mRows;
        mCols = matrix.mCols;
        mData = (T**) malloc(sizeof(T*) * mRows);
        for (unsigned int r=0; r<mRows; r++)
        {
            mData[r] = (T*) malloc(sizeof(T) * mCols);
            memcpy(mData[r], matrix.mData[r], sizeof(T)*mCols);
        }
    }
}

template <typename T>
unsigned int Matrix<T>::rows() const
{
    return mRows;
}

template <typename T>
unsigned int Matrix<T>::cols() const
{
    return mCols;
}

template <typename T>
bool Matrix<T>::empty() const
{
    return (mRows == 0 || mCols == 0 || mData == NULL);
}

template <typename T>
void Matrix<T>::clear()
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

template <typename T>
void Matrix<T>::resize(unsigned int rows, unsigned int cols, double value)
{
    if (rows == 0 && cols == 0) clear();

    if (rows > 0 && cols > 0)
    {
        if (mData == NULL)
        {
            mData = (T**) malloc(sizeof(T*)*rows);
            for (unsigned int i=0; i<rows; i++)
            {
                mData[i] = (T*) malloc(sizeof(T)*cols);
                for (unsigned int j=0; j<cols; j++) mData[i][j] = value;
            }
            mRows = rows;
            mCols = cols;
        }
        else
        {
            if (rows != mRows)
            {
                T **ptr = (T **) realloc(mData, sizeof(T*) * rows);
                if (cols != mCols)
                {
                    for (unsigned int i=0; i<mRows; i++)
                    {
                        T *pRow = (T *) realloc(ptr[i], sizeof(T) * cols);
                        for (unsigned int j=mCols; j<cols; j++) pRow[j] = value;
                        ptr[i] = pRow;
                    }

                    for (unsigned int i=mRows; i<rows; i++)
                    {
                        T *pRow = (T *) realloc(ptr[i], sizeof(T) * cols);
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

template <typename T>
T& Matrix<T>::operator()(unsigned int row, unsigned int col)
{
    return mData[row][col];
}

template <typename T>
const T& Matrix<T>::operator()(unsigned int row, unsigned int col) const
{
    return mData[row][col];
}

template <typename T>
T& Matrix<T>::at(unsigned int row, unsigned int col)
{
    if (row>=mRows)
    {
        //throw std::out_of_range("row index out of range");
    }
    if (col>=mCols)
    {
        //throw std::out_of_range("column index out of range");
    }
    return mData[row][col];
}

template <typename T>
const T& Matrix<T>::at(unsigned int row, unsigned int col) const
{
    if (row>=mRows)
    {
        //throw std::out_of_range("row index out of range");
    }
    if (col>=mCols)
    {
        //throw std::out_of_range("column index out of range");
    }
    return mData[row][col];
}

template <typename T>
T* Matrix<T>::operator [](unsigned int row) const
{
    return mData[row];
}

template <typename T>
T* Matrix<T>::operator [](unsigned int row)
{
    return mData[row];
}
