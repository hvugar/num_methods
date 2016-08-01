#include "matrix3d.h"
#include <stdio.h>
#include <stdlib.h>

//DoubleCube::DoubleCube() : std::vector<DoubleMatrix>()
//{}

//DoubleCube::~DoubleCube()
//{}

//void DoubleCube::Resize(unsigned int Nz, unsigned int Ny, unsigned Nx)
//{
//    Clear();

//    resize(Nz);
//    for (unsigned int k=0; k<size(); k++)
//    {
//        this[k].resize(Ny);
//        for (unsigned int m=0; m<this[k].size(); m++)
//        {
//            //this[k][m].resize(Nx);
//        }
//    }
//}

//void DoubleCube::Clear()
//{
//    for (unsigned int k=0; k<size(); k++)
//    {
//        for (unsigned int m=0; m<this[k].size(); m++)
//        {
//            this[k][m].clear();
//        }
//        this[k].clear();
//    }
//    this->clear();
//}

DoubleCube::DoubleCube(unsigned int z, unsigned int rows, unsigned int cols, double value) : mZ(z), mRows(rows), mCols(cols), pData(NULL)
{
    if (mZ != 0 && mRows != 0 && mCols != 0)
    {
        pData = (double***) malloc(sizeof(double**) * z);
        for (unsigned int k=0; k<z; k++)
        {
            pData[k] = (double**) malloc(sizeof(double*) * rows);
            for (unsigned int j=0; j<rows; j++)
            {
                pData[k][j] = (double*) malloc(sizeof(double) * cols);
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
    if (pData == NULL)
    {
        for (unsigned int k=0; k<mZ; k++)
        {
            for (unsigned int j=0; j<mRows; j++)
            {
                free(pData[k][j]);
                pData[k][j] = NULL;
            }
            free(pData[k]);
            pData[k] = NULL;
        }
        free(pData);
        pData = NULL;

        mCols = 0;
        mRows = 0;
        mZ = 0;
    }
}

void DoubleCube::resize(unsigned int z, unsigned int rows, unsigned int cols, double value)
{

}

DoubleMatrix DoubleCube::operator [](unsigned int z) const
{
}
