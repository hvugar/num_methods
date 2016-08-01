#ifndef HEATIMAGER_H
#define HEATIMAGER_H

#include <QtGui/QPixmap>
#include <QPainter>
#include <QTextStream>
#include <QFile>
#include <QStringList>
#include <QDebug>
#include <QPen>
#include <stdlib.h>
#include <QVector>
#include <math.h>
//#include <doublevector.h>
#include <vector2d.h>
#include <matrix2d.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef QVector<double> QDoubleVector;
typedef QVector<QDoubleVector> QDoubleMatrix;

void createHeatImage1(int argc, char *argv[]);
void createHeatImage2(int argc, char *argv[]);

#endif // HEATIMAGER_H
