#include "matrixsurface.h"
#include <matrix2d.h>
#include <printer.h>
#include "player.h"
#include <QTimer>

MatrixSurface::MatrixSurface(QObject *parent) : Q3DSurface()
{
    minX = +0.0f;
    maxX = +2880.0f;

    minZ = +0.0f;
    maxZ = +1920.0f;

    minY = -36.0f;
    maxY = +4780.0f;

//    maxY = +1.255582;
//    minY = -0.452912;

//    maxY = +0.115524;
//    minY = -0.052888;

    rotationX = 30.0f;
    rotationY = 90.0f;
    rotationZ = 30.0f;

    countX = 2880;
    countZ = 1920;

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

//    setAspectRatio(2.5);
//    setFlipHorizontalGrid(false);
    setAspectRatio(50.0f);

    axisX()->setLabelAutoRotation(rotationX);
    axisY()->setLabelAutoRotation(rotationY);
    axisZ()->setLabelAutoRotation(rotationZ);

    m_Proxy = new QSurfaceDataProxy();
    m_Series = new QSurface3DSeries(m_Proxy);
    m_Series->setItemLabelFormat(QStringLiteral("(@xLabel, @yLabel): @xLabel"));
    m_Series->setDrawMode(QSurface3DSeries::DrawSurface);
    m_Series->setFlatShadingEnabled(true);

//    QLinearGradient gr;
//    gr.setColorAt(0.00, Qt::darkBlue);
//    gr.setColorAt(0.05, Qt::blue);
//    //gr.setColorAt(0.05, Qt::white);
//    gr.setColorAt(0.10, Qt::darkGreen);
//    gr.setColorAt(0.20, Qt::green);
//    gr.setColorAt(0.30, Qt::yellow);
//    gr.setColorAt(0.40, Qt::darkYellow);
//    gr.setColorAt(0.50, 0xFFC700);
//    gr.setColorAt(0.60, 0xFF7700);
//    gr.setColorAt(0.80, Qt::red);
//    gr.setColorAt(1.00, Qt::darkRed);


    QLinearGradient gr;
    gr.setColorAt(0.00, 0x4cafea);
    gr.setColorAt(0.0001, 0x4cafea);
    gr.setColorAt(0.0002, 0x58c2e9);
    gr.setColorAt(0.0005, 0x9fd5f7);
    gr.setColorAt(0.0020, 0x8bc7e9);
    gr.setColorAt(0.0040, 0xc4e0e4);
    gr.setColorAt(0.0050, 0x8fc963);
    gr.setColorAt(0.0060, 0x98cf71);
    gr.setColorAt(0.0080, 0xc5e274);
    gr.setColorAt(0.0500, 0xffdaa3);
    gr.setColorAt(0.0600, 0xf2ca85);
    gr.setColorAt(0.0700, 0xf0ba64);
    gr.setColorAt(0.0900, 0xe9a844);
    gr.setColorAt(0.4900, 0xd38c22);
    gr.setColorAt(0.5950, 0xb1772d);
    gr.setColorAt(1.0000, 0xa35e23);

    m_Series->setBaseGradient(gr);
    m_Series->setColorStyle(Q3DTheme::ColorStyleRangeGradient);

    removeSeries(m_Series);
    addSeries(m_Series);

    fillMatrix1(":/maps/azerb_dem", 2880, 1920);

//    fillMatrix("e:/data/txt/image400.txt", 101, 101);
//    QTimer *timer = new QTimer;
//    QTimer::connect(timer, SIGNAL(timeout()), this, SLOT(timeout()));
//    timer->setInterval(100);
//    timer->start();
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
            float y = (float) (m[i][j]);
            (*newRow)[index++].setPosition(QVector3D(x, y, z));
        }
        *dataArray << newRow;
    }
    m_Proxy->resetArray(dataArray);

//    QWidget::
}

void MatrixSurface::fillMatrix1(const QString &filename, int w, int h)
{
    DoubleMatrix m(h, w);

    QFile file(filename);
    file.open(QIODevice::ReadOnly);
qDebug() << "OK1";
    unsigned int row = 0;
    unsigned int col = 0;
    while (!file.atEnd())
    {

        QDataStream in;
        in.setDevice(&file);
        QByteArray data = file.read(24);
        double *xyz = (double*)data.data();
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];

        m[row][col] = z;

        row++;

//        qDebug() << x << y << z;

        if (row % 1920 == 0)
        {
//            qDebug() << row << col;
            row = 0;
            col++;

//            if (col == 2)
//            break;
        }
    }

//    qDebug() << m.rows() << m.cols();

//    IPrinter::printMatrix(m);

    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve(countZ);

    for (int i = 0 ; i < countZ ; i++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(countX);
        float z = i;
        int index = 0;
        for (int j = 0; j < countX; j++) {
            float x = j;
            float y = (float) (m[i][j]);
            (*newRow)[index++].setPosition(QVector3D(x, y, z));
        }
        *dataArray << newRow;
    }
    m_Proxy->resetArray(dataArray);

    qDebug() << "OK3";
}

void MatrixSurface::timeout()
{
    count += 2;
    //qDebug() << QString("e:/data/txt/image%1.txt").arg(count);
    fillMatrix(QString("e:/data/txt/image%1.txt").arg(count), 1001, 1001);
}
