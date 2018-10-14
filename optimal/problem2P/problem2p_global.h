#ifndef PROBLEM2P_GLOBAL_H
#define PROBLEM2P_GLOBAL_H

#ifdef _WIN32
#if defined(PROBLEM2P_LIBRARY)
#  define PROBLEM2PSHARED_EXPORT __declspec(dllexport)
#else
#define PROBLEM2PSHARED_EXPORT __declspec(dllimport)
#endif
#else
#define PROBLEM2PSHARED_EXPORT
#endif

#endif // PROBLEM2P_GLOBAL_H
