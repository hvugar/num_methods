#ifndef RANDOM_H
#define RANDOM_H

#include "../global.h"

class MINIMUMSHARED_EXPORT Random
{
public:
    static double value(int min, int max, unsigned int precition);
};

#endif // RANDOM_H
