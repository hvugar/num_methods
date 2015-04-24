#ifndef QMINIMUM_GLOBAL_H
#define QMINIMUM_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(QMINIMUM_LIBRARY)
#  define QMINIMUMSHARED_EXPORT Q_DECL_EXPORT
#else
#  define QMINIMUMSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // QMINIMUM_GLOBAL_H
