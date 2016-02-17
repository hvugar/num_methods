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

#endif // GLOBAL

