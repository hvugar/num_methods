#ifndef PROBLEM2H_GLOBAL_H
#define PROBLEM2H_GLOBAL_H

#if defined(PROBLEM2H_LIBRARY)
#  define PROBLEM2HSHARED_EXPORT __declspec(dllexport)
#else
#  define PROBLEM2HSHARED_EXPORT __declspec(dllimport)
#endif

#endif // PROBLEM2H_GLOBAL_H
