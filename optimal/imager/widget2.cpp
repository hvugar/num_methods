#include "widget2.h"

Widget2::Widget2(QWidget *parent) :
    QWidget(parent)
{
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);
    resize(1000, 600);
    painted = false;

    counter = 0;

    p = new DataPainter;
    p->w = this;
    p->start();
}

void Widget2::paintEvent(QPaintEvent *e)
{
    Q_UNUSED(e);

    QPainter painter(this);

    painter.setPen(QPen(Qt::red, 1.0));
    painter.setRenderHint(QPainter::Antialiasing, true);
    for (int i=1; i<points.size()-1; i++) {
        double y1 = height()/2 - (points.at(i-1)) * 20;
        double y2 = height()/2 - (points.at(i-0)) * 20;
        QLine line((i-1), y1, i, y2);

        painter.drawLine(line);
    }
    painter.drawText(20, 20, QString::number(time, 'f', 4));
    painter.drawText(20, 40, QString::number(integral, 'f', 8));

    painter.end();
    //painted = true;

    {
        QPixmap p(1000, 600);
        p.fill(Qt::white);
        QPainter painter1(&p);
        painter1.setPen(QPen(Qt::red, 1.0));
        painter1.setRenderHint(QPainter::Antialiasing, true);
        for (int i=1; i<points.size()-1; i++) {
            double y1 = height()/2 - (points.at(i-1)) * 20;
            double y2 = height()/2 - (points.at(i-0)) * 20;
            QLine line((i-1), y1, i, y2);
            painter1.drawLine(line);
        }
        painter1.drawText(20, 20, QString::number(time, 'f', 4));
        painter1.drawText(20, 40, QString::number(integral, 'f', 8));
        p.save("pages/page"+QString::number(++counter)+".png", "png");
    }
}

void DataPainter::run()
{
    QFile f("data1.txt");
    f.open(QFile::Text | QFile::ReadOnly);
    QTextStream in(&f);
    QString line = in.readLine();
    while (!line.isNull())
    {
        QStringList list = line.split(" ");
        w->points.clear();
        for (int i=2; i<list.size()-1; i++) {
            w->points.append(list.at(i).toDouble());
        }
        w->time = list.at(0).toDouble();
        w->integral = list.at(1).toDouble();
        //        qDebug() << list.size();
        line = in.readLine();
        w->update();

        msleep(50);
    }
    f.close();
}
