#ifndef PROBLEM2H_GLOBAL_H
#define PROBLEM2H_GLOBAL_H

#ifdef _WIN32
#if defined(PROBLEM2H_LIBRARY)
#  define PROBLEM2HSHARED_EXPORT __declspec(dllexport)
#else
#  define PROBLEM2HSHARED_EXPORT __declspec(dllimport)
#endif
#else
#define PROBLEM2HSHARED_EXPORT
#endif

#define SAVE_TO_IMG_

#endif // PROBLEM2H_GLOBAL_H
