#ifndef WIDGET2_H
#define WIDGET2_H

#include <QWidget>
#include <QPainter>
#include <QPaintEvent>
#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QPainterPath>
#include <QDebug>
#include <QThread>
#include <QRunnable>
#include <QVector>

class DataPainter;

class Widget2 : public QWidget
{
    Q_OBJECT
public:
    explicit Widget2(QWidget *parent = 0);

    virtual void paintEvent(QPaintEvent* e);

    QVector<double> points;
    double time;
    double integral;
    DataPainter *p;
    int counter;
private:
    bool painted;

signals:

public slots:

};

class DataPainter : public QThread
{
public:
    virtual void run();
    Widget2* w;
    int counter;
    double t;
};

#endif // WIDGET2_H
