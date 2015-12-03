#ifndef CGLOBAL_H
#define CGLOBAL_H

#include <stdio.h>

#if defined(C_MINIMUM_LIBRARY)
#  define C_MINIMUMSHARED_EXPORT __declspec(dllexport)
#else
#  define C_MINIMUMSHARED_EXPORT __declspec(dllimport)
#endif

#endif // CGLOBAL_H
