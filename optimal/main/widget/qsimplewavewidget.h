#ifndef QSIMPLEWAVEWIDGET_H
#define QSIMPLEWAVEWIDGET_H

#include <QObject>
#include <QWidget>
#include <QTimer>
#include <QPainter>
#include <QPixmap>
#include <QApplication>
#include <hyperbolicequation.h>

class QSimpleWaveWidget : public QWidget, public IHyperbolicEquation
{
    Q_OBJECT
public:
    static int main(int argc, char ** argv);

    explicit QSimpleWaveWidget(QWidget *parent = 0);
    virtual ~QSimpleWaveWidget() {}

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

signals:

public slots:
};

#endif // QSIMPLEWAVEWIDGET_H
