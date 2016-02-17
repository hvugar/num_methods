#include "widget.h"

Widget::Widget(QWidget *parent) : QWidget(parent)
{
    setAutoFillBackground(true);
    setBackgroundRole(QPalette::Base);
    resize(600, 400);
}

Widget::~Widget()
{

}

void Widget::paintEvent(QPaintEvent *evt)
{
    QPainter painter(this);
    //painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen(QPen(Qt::black, 1.0, Qt::SolidLine));
    painter.drawLine(100, 30, 100, 370);
    painter.drawLine(50, 350, 550, 350);

    painter.setPen(QPen(Qt::black, 2.0, Qt::SolidLine));
    QVector<QPoint> points;
    points.append(QPoint(50, 350));
    points.append(QPoint(250, 350));
    points.append(QPoint(250, 150));
    points.append(QPoint(400, 150));
    points.append(QPoint(400, 350));
    points.append(QPoint(550, 350));
    painter.drawLines(points);
    painter.drawLine(250, 150, 250, 350);
    painter.drawLine(400, 150, 400, 350);

    QImage img("E:/project/hvugar/num_methods/trunk/optimal/grapher/images/arrow.png", "PNG");
    painter.save();
    painter.rotate(-30);
    painter.drawImage(10, 10, img);
}
