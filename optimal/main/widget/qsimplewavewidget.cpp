#include "qsimplewavewidget.h"

int QSimpleWaveWidget::main(int argc, char **argv)
{
    QApplication app(argc, argv);
    QSimpleWaveWidget w;
    w.show();
    return app.exec();
}

QSimpleWaveWidget::QSimpleWaveWidget(QWidget *parent) : QWidget(parent)
{
    resize(1280, 720);
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);

    QPixmap background;
    background.load("simple-yellow.png", "PNG");

    QPixmap equation;
    equation.load("simple_wave/eqaution1.png", "PNG");

    for (unsigned int i=0; i<=1000; i++)
    {
        QString filename;
        if (i<10000) filename = QString("simple_wave/layer1/0000%1.png").arg(i);
        if (i<1000) filename = QString("simple_wave/layer1/00000%1.png").arg(i);
        if (i<100) filename = QString("simple_wave/layer1/000000%1.png").arg(i);
        if (i<10) filename = QString("simple_wave/layer1/0000000%1.png").arg(i);

        QPixmap pixmap;
        pixmap.load(filename, "PNG");

        QPixmap pixmap1(1280, 720);
        QPainter painter(&pixmap1);
        painter.drawPixmap(0, 0, background);
        painter.drawPixmap(1280-350, 20, equation);
        painter.drawPixmap(0, 0, pixmap);

        if (i<10000) filename = QString("simple_wave/full/0000%1.png").arg(i);
        if (i<1000) filename = QString("simple_wave/full/00000%1.png").arg(i);
        if (i<100) filename = QString("simple_wave/full/000000%1.png").arg(i);
        if (i<10) filename = QString("simple_wave/full/0000000%1.png").arg(i);

        pixmap1.save(filename, "PNG");
     }

    return;

    a = 1.0;

    x0 = 0.0;
    t0 = 0.0;
    x1 = 1.0;
    t1 = 1.0;

    hx = 0.001;
    ht = 0.001;
    N  = 1000;
    M  = 1000;

    calculateU(m, hx, ht, M, N);

    j = 0;
    timer.setInterval(40);
    connect(&timer, SIGNAL(timeout()), this, SLOT(onTimeout()));
    timer.start();
}

void QSimpleWaveWidget::paintEvent(QPaintEvent *e)
{
    //    double minimum = m.min();
    //    double maximum = m.max();

    //    const DoubleVector &v = m[j];
    //    int v_size = v.size();
    //    QPixmap pixmap = QPixmap(QSize(1280, 720));
    //    pixmap.fill(Qt::transparent);
    //    QPainter painter(&pixmap);
    //    painter.setPen(Qt::blue);
    //    painter.setRenderHint(QPainter::Antialiasing, true);
    //    for (int i=0; i<v_size-1; i++)
    //    {
    //        QPoint point1(140+i+0, 360-500/(maximum-minimum)*m[j][i+0]);
    //        QPoint point2(140+i+1, 360-500/(maximum-minimum)*m[j][i+1]);
    //        painter.drawLine(point1, point2);
    //    }
    //    painter.setPen(Qt::black);
    //    painter.setBrush(Qt::black);
    //    QPoint center1(140, 360);
    //    painter.drawEllipse(center1, 6, 6);
    //    QPoint center2(1140, 360);
    //    painter.drawEllipse(center2, 6, 6);

    //    QString filename;
    //    if (j<10000) filename = QString("simple_wave/layer1/0000%1.png").arg(j);
    //    if (j<1000) filename = QString("simple_wave/layer1/00000%1.png").arg(j);
    //    if (j<100) filename = QString("simple_wave/layer1/000000%1.png").arg(j);
    //    if (j<10) filename = QString("simple_wave/layer1/0000000%1.png").arg(j);

    //    pixmap.save(filename, "PNG");

    //    QPainter p(this);
    //    p.drawPixmap(0,0,pixmap);


    /////////////////////////////////////////////////////////////////////////
}

void QSimpleWaveWidget::onTimeout()
{
    update();
    j++;
    //if (j > M) timerNumber = 0;
    if (j > M)
        exit(0);
}

double QSimpleWaveWidget::f(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double QSimpleWaveWidget::initial1(unsigned int i) const
{
    double x = i*hx;
    //    return -4.0*(x-0.5)*(x-0.5)+1.0;
    return sin(2.0*M_PI*x);
}

double QSimpleWaveWidget::initial2(unsigned int i) const
{
    return 0.0;
}

double QSimpleWaveWidget::boundary(Boundary type, unsigned int j) const
{
    return 0.0;
}


