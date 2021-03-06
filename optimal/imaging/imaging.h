#ifndef IMAGING_H
#define IMAGING_H

#include <global.h>

#include <math.h>
#include <float.h>

#include <QtGui/QPixmap>
#include <QtGui/QPainter>
#include <QtGui/QLinearGradient>

#include <vector2d.h>
#include <matrix2d.h>

void MINIMUMSHARED_EXPORT visualizeVectorHeat(const DoubleVector &v, double min, double max, QPixmap &img, unsigned int w=0, unsigned int h=0);
void MINIMUMSHARED_EXPORT visualizeMatrixHeat(const DoubleMatrix &m, double min, double max, QPixmap &img, unsigned int w=0, unsigned int h=0);
void MINIMUMSHARED_EXPORT visualHeatColorGradinet1(QPixmap& img, int w=256, int h=10);
void MINIMUMSHARED_EXPORT visualHeatColorGradinet2(QPixmap& img, int w=256, int h=10);
void MINIMUMSHARED_EXPORT visualGrayScale(const DoubleMatrix &m, double min, double max, QPixmap &pxm, size_t w=0, size_t h=0);

void MINIMUMSHARED_EXPORT visualString(const DoubleVector& v, double min, double max, int w, int h, QPixmap &pxm,
                                       QColor bg = Qt::transparent, QColor fg = Qt::black, const QString &filename = QString());
void MINIMUMSHARED_EXPORT visualString1(const DoubleVector& v, double min, double max, int w, int h, QColor bg = Qt::transparent,
                                       QColor fg = Qt::black, const QString &filename = QString());

#endif // IMAGING_H
