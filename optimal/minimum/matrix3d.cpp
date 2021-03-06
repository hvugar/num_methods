#include "matrix3d.h"
#include <stdio.h>
#include <stdlib.h>

DoubleCube::DoubleCube(unsigned int depth, unsigned int rows, unsigned int cols, double value)
    : mDepth(depth), mRows(rows), mCols(cols), pData(nullptr)
{
    if (mDepth != 0 && mRows != 0 && mCols != 0)
    {
        pData = static_cast<double***>(malloc(sizeof(double**)*depth));
        for (unsigned int k=0; k<depth; k++)
        {
            pData[k] = static_cast<double**>(malloc(sizeof(double*)*rows));
            for (unsigned int j=0; j<rows; j++)
            {
                pData[k][j] = static_cast<double*>(malloc(sizeof(double)*cols));
                for (unsigned int i=0; i<cols; i++) pData[k][j][i] = value;
            }
        }
    }
}

DoubleCube::~DoubleCube()
{
    clear();
}

void DoubleCube::clear()
{
    if (pData != nullptr)
    {
        for (unsigned int k=0; k<mDepth; k++)
        {
            for (unsigned int j=0; j<mRows; j++)
            {
                free(pData[k][j]);
                pData[k][j] = nullptr;
            }
            free(pData[k]);
            pData[k] = nullptr;
        }
        free(pData);
        pData = nullptr;

        mCols = 0;
        mRows = 0;
        mDepth = 0;
    }
}

void DoubleCube::resize(unsigned int depth, unsigned int rows, unsigned int cols, double value)
{
    C_UNUSED(depth);
    C_UNUSED(rows);
    C_UNUSED(cols);
    C_UNUSED(value);

    if (depth == 0 && rows == 0 && cols == 0) clear();

    if (depth > 0 && rows > 0 && cols > 0)
    {
        if (pData == nullptr)
        {
            pData = static_cast<double***>(malloc(sizeof(double**)*depth));
            for (unsigned int k=0; k<depth; k++)
            {
                pData[k] = static_cast<double**>(malloc(sizeof(double*)*rows));
                for (unsigned int j=0; j<rows; j++)
                {
                    pData[k][j] = static_cast<double*>(malloc(sizeof(double)*cols));
                    for (unsigned int i=0; i<cols; i++) pData[k][j][i] = value;
                }

            }
            mDepth = depth;
            mRows = rows;
            mCols = cols;
        }
        else
        {
//            if (depth != mDepth)
//            {
//                double **ptr = (double **)realloc(mData, sizeof(double*) * rows);
//                if (cols != mCols)
//                {
//                    for (unsigned int j=0; j<mRows; j++)
//                    {
//                        double *pRow = (double *) realloc(ptr[j], sizeof(double) * cols);
//                        for (unsigned int i=mCols; i<cols; i++) pRow[i] = value;
//                        ptr[j] = pRow;
//                    }

//                    for (unsigned int j=mRows; j<rows; j++)
//                    {
//                        double *pRow = (double *) realloc(ptr[j], sizeof(double) * cols);
//                        for (unsigned int i=0; i<cols; i++) pRow[i] = value;
//                        ptr[j] = pRow;
//                    }

//                    mCols = cols;
//                }
//                mRows = rows;
//                mData = ptr;
//            }
        }
    }
}

//DoubleMatrix DoubleCube::operator [](unsigned int depth) const
//{
//    DoubleMatrix matrix(mRows, mCols);
//    for (unsigned int j=0; j<mRows; j++)
//        for (unsigned int i=0; i<mCols; i++)
//            matrix.at(j, i) = pData[depth][j][i];
//    return matrix;
//}

DoubleMatrix DoubleCube::matrix(unsigned int z) const
{
    DoubleMatrix mx(mRows, mCols);
    for (unsigned int j=0; j<mRows; j++)
        for (unsigned int i=0; i<mCols; i++)
            mx.at(j, i) = pData[z][j][i];
    return mx;
}

unsigned int DoubleCube::depth() const
{
    return mDepth;
}

unsigned int DoubleCube::rows() const
{
    return mRows;
}

unsigned int DoubleCube::cols() const
{
    return mCols;
}

double& DoubleCube::operator ()(unsigned int depth, unsigned int row, unsigned int col)
{
    return pData[depth][row][col];
}

const double& DoubleCube::operator ()(unsigned int depth, unsigned int row, unsigned int col) const
{
    return pData[depth][row][col];
}

double& DoubleCube::at(unsigned int k, unsigned int j, unsigned int i)
{
    return pData[k][j][i];
}

const double& DoubleCube::at(unsigned int k, unsigned int j, unsigned int i) const
{
    return pData[k][j][i];
}
