#ifndef PROBLEM1P_GLOBAL_H
#define PROBLEM1P_GLOBAL_H

#ifdef _WIN32
#if defined(PROBLEM1P_LIBRARY)
#   define PROBLEM1PSHARED_EXPORT __declspec(dllexport)
#else
#   define PROBLEM1PSHARED_EXPORT __declspec(dllimport)
#endif
#else
#   define PROBLEM1PSHARED_EXPORT
#endif

#ifdef USE_LIB_IMAGING
#include <QPixmap>
#include <QGuiApplication>
#include <imaging.h>
#endif

#endif // PROBLEM1P_GLOBAL_H
