#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>

#ifdef OS_UNIX
#define MINIMUMSHARED_EXPORT
#else
#if defined(MINIMUM_LIBRARY)
#  define MINIMUMSHARED_EXPORT __declspec(dllexport)
#else
#  define MINIMUMSHARED_EXPORT __declspec(dllimport)
#endif
#endif

#define C_UNUSED(x) (void)x;

enum Boundary
{
    Left = 0,
    Right = 0
};

enum BoundaryCondition
{
    Dirichlet = 0,
    Neumann = 1,
    Robin = 2,
    Mixed = 3
};

#endif // GLOBAL

