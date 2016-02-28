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
#include "widget2.h"
#include <math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef QVector<double> QDoubleVector;
typedef QVector<QDoubleVector> QDoubleMatrix;

void createHeatImage(int width, int height, const QString &inFile, const QString &outFile);

int main(int argc, char *argv[])
{
    QApplication a(argc, argv, false);

    if (argc < 9)
    {
        printf("Usage: images.exe -w 100 -h 100 -i filename.txt -o filename.png\n");
    }
    else
    {
        createHeatImage(QString(argv[2]).toInt(), QString(argv[4]).toInt(), QString(argv[6]), QString(argv[8]));
    }

    return 0;
}

void createHeatImage(int width, int height, const QString &inFile, const QString &outFile)
{
    QPixmap pixmap(QSize(width, height));
    pixmap.fill(Qt::white);
    QPainter painter(&pixmap);
    painter.drawPoint(0,0);

    QFile file(inFile);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);

    double maximum = -999999999999.9;
    double minimum = +999999999999.9;

    QDoubleMatrix m;
    m.resize(height);

    unsigned int j=0;
    QString line = in.readLine();
    while (!line.isNull() && !line.isEmpty())
    {
        unsigned int i=0;
        QStringList list = line.split(" ");
        m[j].resize(width);
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

    printf("Minimum: %.10f Maximum: %.10f width: %d height %d\n", minimum, maximum, width, height);
    minimum = -0.09900679249999990000;
    maximum = +0.12911580289999900000;

//    FILE *file1 = fopen("minmax.txt1", "a");
//    fprintf(file1, "%.20f %.20f\n", minimum, maximum);
//    fclose(file1);

    double sum = 0.0;
    for (int j=0; j<m.size(); j++)
    {
        for (int i=0; i<m[j].size(); i++)
        {
            double u = m[j][i];

            //unsigned int c = (u-minimum)*(maximum-minimum)*0xffffff;
            double ratio = 0.0;
            if (minimum!=maximum)
                ratio = 2.0 * (u-minimum) / (maximum - minimum);
            int b = int(MAX(0, 255*(1 - ratio)));
            int r = int(MAX(0, 255*(ratio - 1)));
            int g = 255 - b - r;
            QColor c(r, g, b);
            painter.setPen(c);
            painter.drawPoint(i,height-j-1);

            sum += u*u;
        }
    }
    sum = sqrt(sum);
    painter.setPen(Qt::black);
    painter.drawText(10, 10, QString::number(sum, 'f', 6));


    pixmap.save(outFile, "PNG");

    //qDebug() << "end";
}
