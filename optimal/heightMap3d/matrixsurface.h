#ifndef MATRIXSURFACE_H
#define MATRIXSURFACE_H

#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QValue3DAxis>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtDataVisualization/Q3DCamera>
#include <QtDataVisualization/QHeightMapSurfaceDataProxy>

#include <QDebug>

using namespace QtDataVisualization;

class MatrixSurface : public Q3DSurface
{
    Q_OBJECT
public:
    explicit MatrixSurface(QObject *parent = nullptr);

    QSurfaceDataProxy *m_Proxy;
    QSurface3DSeries  *m_Series;

    float minX;
    float maxX;
    float minZ;
    float maxZ;
    float minY;
    float maxY;

    float rotationX;
    float rotationZ;
    float rotationY;

    int countX;
    int countZ;
    int countY;

    void fillMatrix(const QString &filename, int w, int h);
    void fillMatrix1(const QString &filename, int w, int h);

    int count = 0;

signals:

public slots:
    void timeout();
};

#endif // MATRIXSURFACE_H
