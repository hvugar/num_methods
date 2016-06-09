#include "sample1.h"

QSample1::QSample1(QWidget *parent) : QWidget(parent)
{
    resize(1280, 720);
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);

    for (unsigned int i=0; i<=2000; i++)
    {
        QPixmap pixmap1 = QPixmap(QSize(1280, 720));
        pixmap1.fill(Qt::white);

        QPixmap pixmap2;
        pixmap2.load(QString("data3/image%1.png").arg(i));

        QPixmap pixmap3;
        pixmap3.load(QString("data4/image%1.png").arg(i));

        QPixmap pixmap4;
        pixmap4.load(QString("data5/image%1.png").arg(i));

        QPainter painter(&pixmap1);
        painter.drawPixmap(0,0,pixmap2);
        painter.drawPixmap(0,0,pixmap3);
        painter.drawPixmap(0,0,pixmap4);

        pixmap1.save(QString("data6/image%1.png").arg(i), "PNG");
    }

    return;

    QFile file("matrix.txt");
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);

    double maximum = -9999.9;
    double minimum = +9999.9;

    DoubleMatrix m;
    m.resize(2001);

    unsigned int j=0;
    QString line = in.readLine();
    while (!line.isNull() && !line.isEmpty())
    {
        unsigned int i=0;
        QStringList list = line.split(" ");
        m[j].resize(1001);
        for (int k=0; k<list.size(); k++)
        {
            QString str = list[k];
            if (!str.isNull() && !str.isEmpty())
            {
                double u = str.toDouble();
                m[j][i] = u;
                if (u<=minimum) minimum = u;
                if (u>=maximum) maximum = u;
                i++;
            }
        }
        line = in.readLine();
        j++;
    }
    file.close();

    minimum = -1.2627;
    maximum = +1.2627;

    qDebug() << minimum << maximum;

    for (unsigned int j=0; j<m.size(); j++)
    {
        const DoubleVector &v = m[j];
        int v_size = v.size();
        QPixmap pixmap = QPixmap(QSize(1280, 720));
        pixmap.fill(Qt::transparent);
        QPainter painter(&pixmap);
        painter.setPen(Qt::green);
        painter.setRenderHint(QPainter::Antialiasing, true);
        for (int i=0; i<v_size-1; i++)
        {
            QPoint point1(140+i+0, 360-600/(maximum-minimum)*m[j][i+0]);
            QPoint point2(140+i+1, 360-600/(maximum-minimum)*m[j][i+1]);
            painter.drawLine(point1, point2);
        }
        painter.setPen(Qt::black);
        painter.setBrush(Qt::black);
        QPoint center1(140, 360);
        painter.drawEllipse(center1, 6, 6);
        QPoint center2(1140, 360);
        painter.drawEllipse(center2, 6, 6);
        pixmap.save(QString("data5/image%1.png").arg(j), "PNG");
    }
}

void QSample1::paintEvent(QPaintEvent *e)
{
}
