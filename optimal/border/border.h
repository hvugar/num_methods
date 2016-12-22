#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>

#ifdef _WIN32
#if defined(BORDER_LIBRARY)
#  define BORDERSHARED_EXPORT __declspec(dllexport)
#else
#  define BORDERSHARED_EXPORT __declspec(dllimport)
#endif
#else
#define BORDERSHARED_EXPORT
#endif

#define C_UNUSED(x) (void)x;

#ifdef __GNUC__
#define UNUSED_PARAM __attribute__ ((unused))
#else
#define UNUSED_PARAM
#endif

#endif // GLOBAL

