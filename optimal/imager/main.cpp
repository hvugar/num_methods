#include <QApplication>
#include <QPixmap>
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

    int width = 1001;
    int height = 1001;
    QPixmap pixmap(QSize(width, height));
    pixmap.fill(Qt::white);

    QFile file("data.txt");
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);
    QPainter painter(&pixmap);

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
            painter.drawPoint(i,height-j);
        }
    }


    pixmap.save("image.png");

    return 0;
}
