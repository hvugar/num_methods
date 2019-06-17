#ifndef PROBLEM0H_GLOBAL_H
#define PROBLEM0H_GLOBAL_H

#ifdef _WIN32
#if defined(PROBLEM0H_LIBRARY)
#   define PROBLEM0HSHARED_EXPORT __declspec(dllexport)
#else
#   define PROBLEM0HSHARED_EXPORT __declspec(dllimport)
#endif
#else
#   define PROBLEM0HSHARED_EXPORT
#endif

#define SAVE_TO_IMG_

#include <function.h>
#include <gradient.h>
#include <projection.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>
#include <grid/hibvp.h>

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

#endif // PROBLEM0H_GLOBAL_H
