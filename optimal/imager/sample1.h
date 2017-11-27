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
#include <QTimer>
#include <QRunnable>
#include <QVector>
#include <vector2d.h>
#include <matrix2d.h>
#include <pde_old/hyperbolicequation.h>

class QSample1 : public QWidget, public IHyperbolicEquation
{
    Q_OBJECT
public:
    explicit QSample1(QWidget *parent = 0);
    virtual ~QSample1() {}

    virtual void paintEvent(QPaintEvent* e);

public slots:
    void onTimeout();

protected:
    virtual double f(unsigned int i, unsigned int j) const;
    virtual double initial1(unsigned int i) const;
    virtual double initial2(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;

private:
    double x0;
    double x1;
    double t0;
    double t1;
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double a;
    DoubleMatrix m;

    QTimer timer;
    unsigned int j;
};

#endif // SAMPLE1_H
