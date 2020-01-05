#ifndef PROBLEM1H_GLOBAL_H
#define PROBLEM1H_GLOBAL_H

#ifdef _WIN32
#if defined(PROBLEM1H_LIBRARY)
#   define PROBLEM1HSHARED_EXPORT __declspec(dllexport)
#else
#   define PROBLEM1HSHARED_EXPORT __declspec(dllimport)
#endif
#else
#   define PROBLEM1HSHARED_EXPORT
#endif

#ifdef USE_IMAGING
#include <QPixmap>
#include <QGuiApplication>
#include <imaging.h>
#endif

#define SAVE_TO_IMG_
#define OLD_VERSION_
#define NEW_VERSION
//#define TIME_DISCRETE_H
#define DISCRETE_DELTA_TIME

#endif // PROBLEM2H_GLOBAL_H
