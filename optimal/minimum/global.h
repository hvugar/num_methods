#ifndef GLOBAL_H
#define GLOBAL_H

#ifdef _WIN32
#if defined(MINIMUM_LIBRARY)
#  define MINIMUMSHARED_EXPORT __declspec(dllexport)
#else
#  define MINIMUMSHARED_EXPORT __declspec(dllimport)
#endif
#else
#define MINIMUMSHARED_EXPORT
#endif

#define C_UNUSED(x) (void)x;

#ifdef __GNUC__
#define UNUSED_PARAM __attribute__ ((unused))
#else
#define UNUSED_PARAM
#define M_PI 3.14
#endif

#if (defined(__cplusplus))
#if (__cplusplus == 201103)
#define NOEXCEPT noexcept
#else
#define NOEXCEPT throw()
#endif
#endif

#endif // GLOBAL
