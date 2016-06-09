#ifndef SAMPLE1_H
#define SAMPLE1_H

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
#include <doublevector.h>

class QSample1 : public QWidget
{
    Q_OBJECT
public:
    explicit QSample1(QWidget *parent = 0);
    virtual ~QSample1() {}

    virtual void paintEvent(QPaintEvent* e);
};

#endif // SAMPLE1_H
