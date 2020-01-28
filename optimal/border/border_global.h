#ifndef BORDER_GLOBAL_H
#define BORDER_GLOBAL_H

#ifdef _WIN32
#if defined(BORDER_LIBRARY)
#  define BORDERSHARED_EXPORT __declspec(dllexport)
#else
#  define BORDERSHARED_EXPORT __declspec(dllimport)
#endif
#else
#define BORDERSHARED_EXPORT
#endif

#ifdef USE_LIB_IMAGING
#include <QPixmap>
#include <QGuiApplication>
#include <imaging.h>
#endif

#include <grid/hibvp.h>
#include <grid/pibvp.h>
#include <grid/pibvpX.h>
#include <matrix2d.h>
#include <benchmark.h>
#include <time.h>

#endif // BORDER_GLOBAL_H
