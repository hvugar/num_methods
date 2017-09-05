#include "random.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double Random::value(int min, int max, unsigned int precition)
{
    int a = rand() % (max-min);
    double b = pow(10.0, precition);
    double f = (rand() % (int)b) / b;
    return (min+a)+f;
}

void Random::fillMatrix(DoubleMatrix &m, int min, int max, unsigned int precition)
{
    unsigned int rows = m.rows();
    unsigned int cols = m.cols();
    for (unsigned int row = 0; row < rows; row++)
    {
        for (unsigned int col = 0; col < cols; col++)
        {
            m[row][col] = value(min, max, precition);
        }
    }
}

void Random::fillVector(DoubleVector &v, int min, int max, unsigned int precition)
{
    unsigned int length = v.length();
    for (unsigned int i = 0; i < length; i++)
    {
        v[i] = value(min, max, precition);
    }
}
