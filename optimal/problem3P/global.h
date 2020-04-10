#ifndef PROBLEM3P_GLOBAL_H
#define PROBLEM3P_GLOBAL_H

#ifdef _WIN32
#if defined(PROBLEM3P_LIBRARY)
#   define PROBLEM3P_SHARED_EXPORT __declspec(dllexport)
#else
#   define PROBLEM3P_SHARED_EXPORT __declspec(dllimport)
#endif
#else
#   define PROBLEM3P_SHARED_EXPORT
#endif

//#define NDEBUG
#include <assert.h>

#include <cfloat>
#include <climits>
#include <functional>

#include <gradient.h>
#include <projection.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>
#include <grid/pibvp.h>
#include <ode/lode1o.h>
#include <benchmark.h>
#include <deltagrid.h>
#include <r1minimize.h>
#include <function.h>

#ifdef USE_LIB_IMAGING
#include <QPixmap>
#include <QGuiApplication>
#include <QFile>
#include <QDir>
#include <imaging.h>
#endif

#ifdef USE_LIB_XLSX_WRITER
#include <xlsxwriter.h>
#endif

#endif // PROBLEM3P_GLOBAL_H
