#include "heightmapsurface.h"

HeightMapSurface::HeightMapSurface(QObject *parent) : Q3DSurface()
{
    minX = +0.0f;
    maxX = +1000.0f;
    minZ = +0.0f;
    maxZ = +1000.0f;
    minY = 0.0f;
    maxY = 1.0f;

    rotationX = 30.0f;
    rotationY = 90.0f;
    rotationZ = 30.0f;

    //countX = 1000;
    //countZ = 1000;

    setReflection(true);
    setSelectionMode(QAbstract3DGraph::SelectionNone);
    setShadowQuality(QAbstract3DGraph::ShadowQualityNone);

    setAxisX(new QValue3DAxis);
    setAxisY(new QValue3DAxis);
    setAxisZ(new QValue3DAxis);

    axisX()->setLabelFormat("x: %.2f");
    axisZ()->setLabelFormat("z: %.2f");
    axisY()->setLabelFormat("y: %.2f");
    axisX()->setRange(minX, maxX);
    axisY()->setRange(minY, maxY);
    axisZ()->setRange(minZ, maxZ);

    axisX()->setLabelAutoRotation(rotationX);
    axisY()->setLabelAutoRotation(rotationY);
    axisZ()->setLabelAutoRotation(rotationZ);


    QImage heightMapImage("E:/image1000.png");
    //QImage heightMapImage("E:/ASTGDEMV2_0N38E045.jpg");
    m_Proxy = new QHeightMapSurfaceDataProxy(heightMapImage);
    m_Proxy->setValueRanges(0.0f, 1000.0f, 0.0f, 1000.0f);

    m_Series = new QSurface3DSeries(m_Proxy);
    m_Series->setItemLabelFormat(QStringLiteral("(@xLabel, @zLabel): @yLabel"));
    m_Series->setDrawMode(QSurface3DSeries::DrawSurface);
    m_Series->setFlatShadingEnabled(true);

    QLinearGradient gr;
    gr.setColorAt(0.0, Qt::darkBlue);
    gr.setColorAt(0.1, Qt::blue);
    gr.setColorAt(0.2, Qt::darkGreen);
    gr.setColorAt(0.3, Qt::green);
    gr.setColorAt(0.4, Qt::yellow);
    gr.setColorAt(0.5, Qt::darkYellow);
    gr.setColorAt(0.6, 0xFFC700);
    gr.setColorAt(0.7, 0xFF7700);
    gr.setColorAt(1.0, Qt::darkRed);

    m_Series->setBaseGradient(gr);
    m_Series->setColorStyle(Q3DTheme::ColorStyleRangeGradient);

    addSeries(m_Series);
}
