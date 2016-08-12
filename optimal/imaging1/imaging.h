#ifndef IMAGING_H
#define IMAGING_H

#include <QtGui/QPixmap>
#include <QtGui/QPainter>

#include <vector2d.h>
#include <matrix2d.h>
#include <global.h>

#ifdef __cplusplus
extern "C" {
#endif

void visualizeVectorHeat(const DoubleVector &v, double min, double max, QPixmap &img, unsigned int w=0, unsigned int h=10);
void visualizeMatrixHeat(const DoubleMatrix& m, double min, double max, QPixmap &img, unsigned int w=0, unsigned int h=0);

#ifdef __cplusplus
}
#endif



#endif // IMAGING_H
