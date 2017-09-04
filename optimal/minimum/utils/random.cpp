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
