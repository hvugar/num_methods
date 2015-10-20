#include <QApplication>
#include <QPixmap>
#include <QPainter>
#include <QTextStream>
#include <QFile>
#include <QStringList>
#include <QDebug>
#include <QPen>
#include <stdlib.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int main(int argc, char *argv[])
{
    QApplication a(argc, argv, false);

    QPixmap pixmap(QSize(2000, 2000));
    pixmap.fill(Qt::white);

    QFile file("data2000.txt");
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);
    QPainter painter(&pixmap);

    QString line = in.readLine();
    unsigned int j=0;

    double maximum = 3.0;
    double minimum = 1.0;

    while (!line.isNull())
    {
        QStringList list = line.split(" ");

        for (int i=0; i<list.size(); i++)
        {
            double u = list[i].toDouble();

            double ratio = 2.0 * (u-minimum) / (maximum - minimum);
            int b = int(MAX(0, 255*(1 - ratio)));
            int r = int(MAX(0, 255*(ratio - 1)));
            int g = 255 - b - r;

            QColor c(r, g, b);// = (unsigned int)((u-2.0)*(0xFFFFFF));
            painter.setPen(c);

//            if (list[i].toDouble() >= 2.00 && list[i].toDouble() < 2.05) painter.setPen(QPen(0xFFBE7A));
//            if (list[i].toDouble() >= 2.05 && list[i].toDouble() < 2.10) painter.setPen(QPen(0xFFB473));
//            if (list[i].toDouble() >= 2.10 && list[i].toDouble() < 2.15) painter.setPen(QPen(0xFFAA6D));
//            if (list[i].toDouble() >= 2.15 && list[i].toDouble() < 2.20) painter.setPen(QPen(0xFFA066));
//            if (list[i].toDouble() >= 2.20 && list[i].toDouble() < 2.25) painter.setPen(QPen(0xF09660));
//            if (list[i].toDouble() >= 2.25 && list[i].toDouble() < 2.30) painter.setPen(QPen(0xE08C5A));
//            if (list[i].toDouble() >= 2.30 && list[i].toDouble() < 2.35) painter.setPen(QPen(0xD08253));
//            if (list[i].toDouble() >= 2.35 && list[i].toDouble() < 2.40) painter.setPen(QPen(0xC0784D));
//            if (list[i].toDouble() >= 2.40 && list[i].toDouble() < 2.45) painter.setPen(QPen(0xB06E46));
//            if (list[i].toDouble() >= 2.45 && list[i].toDouble() < 2.50) painter.setPen(QPen(0xB06E46));
//            if (list[i].toDouble() >= 2.50 && list[i].toDouble() < 2.55) painter.setPen(QPen(0x905A3A));
//            if (list[i].toDouble() >= 2.55 && list[i].toDouble() < 2.60) painter.setPen(QPen(0x805033));
//            if (list[i].toDouble() >= 2.60 && list[i].toDouble() < 2.65) painter.setPen(QPen(0x805033));
//            if (list[i].toDouble() >= 2.65 && list[i].toDouble() < 2.70) painter.setPen(QPen(0x603C26));
//            if (list[i].toDouble() >= 2.70 && list[i].toDouble() < 2.75) painter.setPen(QPen(0x50321F));
//            if (list[i].toDouble() >= 2.75 && list[i].toDouble() < 2.80) painter.setPen(QPen(0x40281A));
//            if (list[i].toDouble() >= 2.80 && list[i].toDouble() < 2.85) painter.setPen(QPen(0x301E13));
//            if (list[i].toDouble() >= 2.85 && list[i].toDouble() < 2.90) painter.setPen(QPen(0x1F130D));
//            if (list[i].toDouble() >= 2.90 && list[i].toDouble() < 2.95) painter.setPen(QPen(0x100906));
//            if (list[i].toDouble() >= 2.95 && list[i].toDouble() < 3.00) painter.setPen(QPen(0x000001));

//            painter.setPen((int)(list[i].toDouble()*0x02ffff));

            painter.drawPoint(i,j);
        }

        line = in.readLine();
        j++;

    }
    file.close();

    pixmap.save("image.png");

    return 0;
}
