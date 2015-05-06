#ifndef GLOBAL_H
#define GLOBAL_H

#if defined(MINIMUM_LIBRARY)
#  define MINIMUMSHARED_EXPORT __declspec(dllexport)
#else
#  define MINIMUMSHARED_EXPORT __declspec(dllimport)
#endif

#endif // GLOBAL

