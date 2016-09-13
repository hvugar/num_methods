#include "r1functionvector.h"

R1FunctionVector::R1FunctionVector(unsigned int size)
{
    sz = size;
    pfx = new R1Function* [sz];
}

R1FunctionVector::~R1FunctionVector()
{
    if (pfx != NULL) delete [] pfx;
}

R1Function* R1FunctionVector::at(unsigned int i) const
{
    return pfx[i];
}
