#include "matrixsurface.h"
#include <matrix2d.h>
#include <printer.h>

MatrixSurface::MatrixSurface(QObject *parent) : Q3DSurface()
{
    minX = +0.0f;
    maxX = +100.0f;
    minZ = +0.0f;
    maxZ = +100.0f;
    minY = -0.456183f*100.0f;
    maxY = +1.255503f*100.0f;

    rotationX = 30.0f;
    rotationY = 90.0f;
    rotationZ = 30.0f;

    countX = 100;
    countZ = 100;

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
    axisY()->setSegmentCount(1);

    //setAspectRatio(10.0);
    setFlipHorizontalGrid(false);

    axisX()->setLabelAutoRotation(rotationX);
    axisY()->setLabelAutoRotation(rotationY);
    axisZ()->setLabelAutoRotation(rotationZ);

    m_Proxy = new QSurfaceDataProxy();

    fillMatrix("e:/data/txt/image200.txt", 101, 101);

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
    gr.setColorAt(0.10, Qt::darkRed);

    m_Series->setBaseGradient(gr);
    m_Series->setColorStyle(Q3DTheme::ColorStyleRangeGradient);

    addSeries(m_Series);
}

void MatrixSurface::fillMatrix(const QString &filename, int w, int h)
{
    DoubleMatrix m(w, h);

    QFile file(filename);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);

    unsigned int j=0;
    QString line = in.readLine();
    while (!line.isNull() && !line.isEmpty())
    {
        unsigned int i=0;
        QStringList list = line.split(" ");
        for (int k=0; k<list.size(); k++)
        {
            QString str = list[k];
            if (!str.isNull() && !str.isEmpty())
            {
                double u = str.toDouble();
                m[j][i] = u;
                i++;
            }
        }
        line = in.readLine();
        j++;
    }
    file.close();

    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve(countZ+1);

    for (int i = 0 ; i <= countZ ; i++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(countX+1);
        float z = i;
        int index = 0;
        for (int j = 0; j <= countX; j++) {
            float x = j;
            float y = (float) (m[i][j]*100.0f);
            (*newRow)[index++].setPosition(QVector3D(x, y, z));
        }
        *dataArray << newRow;
    }
    m_Proxy->resetArray(dataArray);
}
