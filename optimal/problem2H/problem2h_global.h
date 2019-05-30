#ifndef PROBLEM2H_GLOBAL_H
#define PROBLEM2H_GLOBAL_H

#ifdef _WIN32
#if defined(PROBLEM2H_LIBRARY)
#   define PROBLEM2HSHARED_EXPORT __declspec(dllexport)
#else
#   define PROBLEM2HSHARED_EXPORT __declspec(dllimport)
#endif
#else
#   define PROBLEM2HSHARED_EXPORT
#endif

#define SAVE_TO_IMG_
#define OLD_VERSION_
#define NEW_VERSION
//#define TIME_DISCRETE_H

#ifdef USE_IMAGING
#include <QPixmap>
#include <QGuiApplication>
#include <imaging.h>
#endif

#define DISCRETE_DELTA_TIME_2

#endif // PROBLEM2H_GLOBAL_H
