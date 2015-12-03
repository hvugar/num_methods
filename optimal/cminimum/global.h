#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>

#if defined(C_MINIMUM_LIBRARY)
#  define C_MINIMUMSHARED_EXPORT __declspec(dllexport)
#else
#  define C_MINIMUMSHARED_EXPORT __declspec(dllimport)
#endif

#endif // GLOBAL
