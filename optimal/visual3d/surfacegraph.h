#ifndef SURFACEGRAPH_H
#define SURFACEGRAPH_H

#include <QtCore/QObject>
#include <QtGui/QLinearGradient>
#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QValue3DAxis>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>

using namespace QtDataVisualization;

class SurfaceGraph : public QObject
{
    Q_OBJECT
public:
    explicit SurfaceGraph(Q3DSurface *surface);
    virtual ~SurfaceGraph();

    void fillSqrtSinProxy();

signals:

public slots:

private:
    Q3DSurface *m_graph;
    QSurfaceDataProxy *m_sqrtSinProxy;
    QSurface3DSeries *m_sqrtSinSeries;

    double sampleMax;
    double sampleMin;

    int sampleCountX;
    int sampleCountZ;
};

#endif // SURFACEGRAPH_H
