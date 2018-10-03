#ifndef HEIGHTMAPSURFACE_H
#define HEIGHTMAPSURFACE_H

#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QValue3DAxis>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtDataVisualization/Q3DCamera>
#include <QtDataVisualization/QHeightMapSurfaceDataProxy>

using namespace QtDataVisualization;

class HeightMapSurface : public Q3DSurface
{
    Q_OBJECT
public:
    explicit HeightMapSurface(QObject *parent = nullptr);

private:
    QHeightMapSurfaceDataProxy *m_Proxy;
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

    //int countX;
    //int countZ;

signals:

public slots:
};

#endif // HEIGHTMAPSURFACE_H
