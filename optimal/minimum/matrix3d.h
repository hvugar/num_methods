#ifndef MATRIX3D_H
#define MATRIX3D_H

#include "global.h"
#include "matrix2d.h"

//class MINIMUMSHARED_EXPORT DoubleCube : public std::vector<DoubleMatrix>
//{
//public:
//    DoubleCube();
//    virtual ~DoubleCube();

//    void Resize(unsigned int Nz, unsigned int Ny, unsigned Nx);
//    void Clear();
//};

class MINIMUMSHARED_EXPORT DoubleCube
{
public:
    explicit DoubleCube(unsigned int z=0, unsigned int rows=0, unsigned int cols=0, double value=0.0);
    virtual ~DoubleCube();

    void Resize(unsigned int Nz, unsigned int Ny, unsigned Nx);
    void Clear();

    void clear();
    void resize(unsigned int z, unsigned int rows, unsigned int cols, double value=0.0);

    DoubleMatrix operator[] (unsigned int z) const;

    unsigned int depth() const;
    unsigned int rows() const;
    unsigned int cols() const;

private:
    unsigned int mDepth;
    unsigned int mRows;
    unsigned int mCols;
    double ***pData;
};

#endif // MATRIX3D_H
