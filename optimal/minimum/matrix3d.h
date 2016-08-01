#ifndef MATRIX3D_H
#define MATRIX3D_H

#include "global.h"
#include "matrix2d.h"

class MINIMUMSHARED_EXPORT DoubleCube : public std::vector<DoubleMatrix>
{
public:
    DoubleCube();
    virtual ~DoubleCube();

    void Resize(unsigned int Nz, unsigned int Ny, unsigned Nx);
    void Clear();
};

#endif // MATRIX3D_H
