#include <QApplication>
#include <QtGui/QPixmap>
#include <QPainter>
#include <QTextStream>
#include <QFile>
#include <QStringList>
#include <QDebug>
#include <QPen>
#include <stdlib.h>
#include <QVector>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef QVector<double> QDoubleVector;
typedef QVector<QDoubleVector> QDoubleMatrix;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv, false);

    int width = 101;
    int height = 101;
    QPixmap pixmap(QSize(width, height));
    pixmap.fill(Qt::white);
    QPainter painter(&pixmap);
    painter.drawPoint(0,0);

    QFile file("optimal.txt");
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);

    double maximum = -999999999999.9;
    double minimum = +999999999999.9;

    QDoubleMatrix m;
    m.resize(height);

    unsigned int j=0;
    QString line = in.readLine();
    while (!line.isNull())
    {
        QStringList list = line.split(" ");
        m[j].resize(width);
        for (int i=0; i<list.size(); i++)
        {
            double u = list[i].toDouble();
            m[j][i] = u;
            if (u<minimum) minimum = u;
            if (u>maximum) maximum = u;
        }
        line = in.readLine();
        j++;
    }
    file.close();

    qDebug() << minimum << maximum << m.size() << m[0].size() << width << height;

//    minimum = 1.0;
//    maximum = 3.0;
    for (int j=0; j<m.size(); j++)
    {
        for (int i=0; i<m[j].size(); i++)
        {
            double u = m[j][i];
            double ratio = 2.0 * (u-minimum) / (maximum - minimum);
            int b = int(MAX(0, 255*(1 - ratio)));
            int r = int(MAX(0, 255*(ratio - 1)));
            int g = 255 - b - r;
            QColor c(r, g, b);
            painter.setPen(c);
            painter.drawPoint(i,height-j-1);
        }
    }


    pixmap.save(file.fileName()+".png", "PNG");

    qDebug() << "end";
    return a.exec();
}
