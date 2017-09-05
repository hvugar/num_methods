#ifndef RANDOM_H
#define RANDOM_H

#include "../global.h"
#include "../matrix2d.h"

class MINIMUMSHARED_EXPORT Random
{
public:
    static double value(int min, int max, unsigned int precition);
    static void fillMatrix(DoubleMatrix& m, int min, int max, unsigned int precition);
    static void fillVector(DoubleVector& v, int min, int max, unsigned int precition);
};

#endif // RANDOM_H
