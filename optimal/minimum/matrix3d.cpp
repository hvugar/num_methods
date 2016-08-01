#include "matrix3d.h"

DoubleCube::DoubleCube() : std::vector<DoubleMatrix>()
{}

DoubleCube::~DoubleCube()
{}

void DoubleCube::Resize(unsigned int Nz, unsigned int Ny, unsigned Nx)
{
    Clear();

    resize(Nz);
    for (unsigned int k=0; k<size(); k++)
    {
        this[k].resize(Ny);
        for (unsigned int m=0; m<this[k].size(); m++)
        {
            //this[k][m].resize(Nx);
        }
    }
}

void DoubleCube::Clear()
{
    for (unsigned int k=0; k<size(); k++)
    {
        for (unsigned int m=0; m<this[k].size(); m++)
        {
            this[k][m].clear();
        }
        this[k].clear();
    }
    this->clear();
}
