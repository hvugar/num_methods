#include "sample1.h"

QSample1::QSample1(QWidget *parent) : QWidget(parent)
{
    resize(1280, 720);
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);

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

    //    double minimum = m.min();
    //    double maximum = m.max();

    //    for (unsigned int j=0; j<m.size(); j++)
    //    {
    //        const DoubleVector &v = m[j];
    //        int v_size = v.size();
    //        QPixmap pixmap = QPixmap(QSize(1280, 720));
    //        pixmap.fill(Qt::transparent);
    //        QPainter painter(&pixmap);
    //        painter.setPen(Qt::green);
    //        painter.setRenderHint(QPainter::Antialiasing, true);
    //        for (int i=0; i<v_size-1; i++)
    //        {
    //            QPoint point1(140+i+0, 360-600/(maximum-minimum)*m[j][i+0]);
    //            QPoint point2(140+i+1, 360-600/(maximum-minimum)*m[j][i+1]);
    //            painter.drawLine(point1, point2);
    //        }
    //        painter.setPen(Qt::black);
    //        painter.setBrush(Qt::black);
    //        QPoint center1(140, 360);
    //        painter.drawEllipse(center1, 6, 6);
    //        QPoint center2(1140, 360);
    //        painter.drawEllipse(center2, 6, 6);
    //        pixmap.save(QString("data1/image%1.png").arg(j), "PNG");
    //    }

    //    connect(this, SIGNAL(timeout()), this, SLOT(onTimeout()));

    //    for (unsigned int i=0; i<=M; i++)
    //    {
    //        QPixmap pixmap1 = QPixmap(QSize(1280, 720));
    //        pixmap1.fill(Qt::white);

    //        QPixmap pixmap2;
    //        pixmap2.load(QString("data3/image%1.png").arg(i));

    //        QPixmap pixmap3;
    //        pixmap3.load(QString("data4/image%1.png").arg(i));

    //        QPixmap pixmap4;
    //        pixmap4.load(QString("data5/image%1.png").arg(i));

    //        QPainter painter(&pixmap1);
    //        painter.drawPixmap(0,0,pixmap2);
    //        painter.drawPixmap(0,0,pixmap3);
    //        painter.drawPixmap(0,0,pixmap4);

    //        pixmap1.save(QString("data6/image%1.png").arg(i), "PNG");
    //    }

    //    return;

    //    QFile file("matrix.txt");
    //    file.open(QIODevice::ReadOnly | QIODevice::Text);
    //    QTextStream in(&file);

    //    double maximum = -9999.9;
    //    double minimum = +9999.9;

    //    DoubleMatrix m;
    //    m.resize(2001);

    //    unsigned int j=0;
    //    QString line = in.readLine();
    //    while (!line.isNull() && !line.isEmpty())
    //    {
    //        unsigned int i=0;
    //        QStringList list = line.split(" ");
    //        m[j].resize(1001);
    //        for (int k=0; k<list.size(); k++)
    //        {
    //            QString str = list[k];
    //            if (!str.isNull() && !str.isEmpty())
    //            {
    //                double u = str.toDouble();
    //                m[j][i] = u;
    //                if (u<=minimum) minimum = u;
    //                if (u>=maximum) maximum = u;
    //                i++;
    //            }
    //        }
    //        line = in.readLine();
    //        j++;
    //    }
    //    file.close();

    //    minimum = -1.2627;
    //    maximum = +1.2627;

    //    qDebug() << minimum << maximum;


}

void QSample1::paintEvent(QPaintEvent *e)
{
    Q_UNUSED(e);

    double minimum = m.min();
    double maximum = m.max();

    const DoubleVector &v = m.row(j);
    int v_size = v.length();
    QPixmap pixmap = QPixmap(QSize(1280, 720));
    pixmap.fill(Qt::white);
    QPainter painter(&pixmap);
    painter.setPen(Qt::blue);
    painter.setRenderHint(QPainter::Antialiasing, true);
    for (int i=0; i<v_size-1; i++)
    {
        QPoint point1(140+i+0, 360-600/(maximum-minimum)*m[j][i+0]);
        QPoint point2(140+i+1, 360-600/(maximum-minimum)*m[j][i+1]);
        painter.drawLine(point1, point2);
    }
    painter.setPen(Qt::black);
    painter.setBrush(Qt::black);
    QPoint center1(140, 360);
    painter.drawEllipse(center1, 6, 6);
    QPoint center2(1140, 360);
    painter.drawEllipse(center2, 6, 6);

    QString filename;
    if (j<10000) filename = QString("data1/0000%1.png").arg(j);
    if (j<1000) filename = QString("data1/00000%1.png").arg(j);
    if (j<100) filename = QString("data1/000000%1.png").arg(j);
    if (j<10) filename = QString("data1/0000000%1.png").arg(j);

    pixmap.save(filename, "PNG");

    QPainter p(this);
    p.drawPixmap(0,0,pixmap);
}

double QSample1::initial1(unsigned int i) const
{
    double x = i*hx;
    //    return -4.0*(x-0.5)*(x-0.5)+1.0;
    return sin(2.0*M_PI*x);
}

double QSample1::initial2(unsigned int i) const
{
    C_UNUSED(i);
    //double x = i*hx;
    return 0.0;
}

double QSample1::boundary(Boundary type, unsigned int j) const
{
    C_UNUSED(type);
    C_UNUSED(j);
    return 0.0;
}

double QSample1::f(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

void QSample1::onTimeout()
{
    update();
    j++;
    //if (j > M) timerNumber = 0;
    if (j > M)
        exit(0);
}

