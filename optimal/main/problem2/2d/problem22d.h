#ifndef PROBLEM22D_H
#define PROBLEM22D_H

#include "iproblem2forward2d.h"
#include "iproblem2backward2d.h"

#include <function.h>
#include <gradient.h>

class Problem22D
{
public:
    static void Main(int argc, char* argv[]);

private:
    Dimension mTimeDimension;
    Dimension mSpaceDimensionX;
    Dimension mSpaceDimensionY;
};

#endif // PROBLEM22D_H
