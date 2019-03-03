#ifndef PROBLEM0H_GLOBAL_H
#define PROBLEM0H_GLOBAL_H

#ifdef _WIN32
#if defined(PROBLEM2H_LIBRARY)
#   define PROBLEM2HSHARED_EXPORT __declspec(dllexport)
#else
#   define PROBLEM0HSHARED_EXPORT __declspec(dllimport)
#endif
#else
#   define PROBLEM0HSHARED_EXPORT
#endif

#define SAVE_TO_IMG_

#ifdef USE_IMAGING
#include <QPixmap>
#include <QGuiApplication>
#include <imaging.h>
#endif

#endif // PROBLEM0H_GLOBAL_H
