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

void createHeatImage1(int argc, char *argv[]);
void createHeatImage2(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    QApplication a(argc, argv, true);
    createHeatImage1(argc, argv);
    //createHeatImage2(argc, argv);

    return 0;
}

void createHeatImage1(int argc, char *argv[])
{
    if (argc < 9)
    {
        printf("Usage: images.exe -w 100 -h 100 -i filename.txt -o filename.png\n");
        return;
    }

    unsigned int width = QString(argv[2]).toInt();
    unsigned int height = QString(argv[4]).toInt();
    QString inFile = QString(argv[6]);
    QString outFile = QString(argv[8]);


    QPixmap pixmap(QSize(width, height));
    pixmap.fill(Qt::white);
    QPainter painter(&pixmap);
    //painter.drawPoint(0,0);

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

    //minimum = 4.0;
    printf("Minimum: %.10f Maximum: %.10f width: %d height %d\n", minimum, maximum, width, height);
    //minimum = -0.15603702370000000000;
    //maximum = +0.19309054760000000000;

    //FILE *file1 = fopen("00minmax.txt1", "a");
    //fprintf(file1, "%.20f %.20f\n", minimum, maximum);
    //fclose(file1);

    //minimum = -0.19900679249999980000;
    //maximum = +0.22911580289999700000;

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
        }
    }

    //    double sum = 0.0;
    //    for (unsigned int j=0; j<=100; j++)
    //    {
    //        for (unsigned int i=0; i<=100; i++)
    //        {
    //            double k = 1.0;
    //            if (i==0 || i==100) k *= 0.5;
    //            if (j==0 || j==100) k *= 0.5;
    //            sum = sum + k * (m[j][i]) * (m[j][i]);
    //        }
    //    }
    //    double h1 = 0.010;
    //    double h2 = 0.010;
    //    sum = h1*h2*sum;
    //    painter.setPen(Qt::black);
    //    painter.drawText(10, 10, QString::number(sum, 'f', 6));

    pixmap.save(outFile, "PNG");
}

void createHeatImage2(int argc, char *argv[])
{
    if (argc < 7)
    {
        printf("Usage: images.exe -w 100 -i filename.txt -o filename.png\n");
        return;
    }

    unsigned int width = QString(argv[2]).toInt();
    QString inFile = QString(argv[4]);
    QString outFile = QString(argv[6]);

    QPixmap pixmap(QSize(width, 10));
    pixmap.fill(Qt::white);
    QPainter painter(&pixmap);

    QFile file(inFile);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);

    double maximum = -999999999999.9;
    double minimum = +999999999999.9;

    QDoubleVector m;
    m.resize(width);

    QString line = in.readLine();
    QStringList list = line.split(" ");
    unsigned int i=0;
    for (int j=0; j<list.size(); j++)
    {
        QString str = list[j];
        if (!str.isNull() && !str.isEmpty())
        {
            double u = str.toDouble();
            m[i] = u;
            if (u<=minimum) minimum = u;
            if (u>=maximum) maximum = u;
            i++;
        }
    }
    file.close();

    minimum = 0.0;
    printf("Minimum: %.10f Maximum: %.10f width: %d\n", minimum, maximum, width);

    for (int i=0; i<m.size(); i++)
    {
        double u = m[i];

        //unsigned int c = (u-minimum)*(maximum-minimum)*0xffffff;
        double ratio = 0.0;
        if (minimum!=maximum)
            ratio = 2.0 * (u-minimum) / (maximum - minimum);
        int b = int(MAX(0, 255*(1 - ratio)));
        int r = int(MAX(0, 255*(ratio - 1)));
        int g = 255 - b - r;
        QColor c(r, g, b);
        painter.setPen(c);
        painter.drawLine(i,0,i,10);
    }

    pixmap.save(outFile, "PNG");
}
