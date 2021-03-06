#include "imaging.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void visualizeVectorHeat(const DoubleVector &v, double min, double max, QPixmap &img, unsigned int w, unsigned int h)
{
    C_UNUSED(w);
    int size = static_cast<int>(v.length());
    img = QPixmap(size, h);
    img.fill(Qt::transparent);
    QPainter painter(&img);
    for (unsigned int i=0; i<v.length(); i++)
    {
        double u = v[i];
        double ratio = 0.0;
        if (min != max)
        {
            ratio = 2.0 * (u-min) / (max - min);
        }
        int b = int(MAX(0, 255*(1 - ratio)));
        int r = int(MAX(0, 255*(ratio - 1)));
        int g = 255 - b - r;
        QColor c(r, g, b);
        painter.setPen(c);
        painter.drawLine(i,0,i,h);
    }
}

void visualizeMatrixHeat(const DoubleMatrix& m, double min, double max, QPixmap &img, unsigned int w, unsigned int h)
{
    C_UNUSED(w);
    C_UNUSED(h);
    unsigned int rows = m.rows();
    unsigned int cols = m.cols();

    img = QPixmap(static_cast<int>(cols), static_cast<int>(rows));
    img.fill(Qt::transparent);
    QPainter painter(&img);
    painter.setRenderHint(QPainter::Antialiasing, true);
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            double u = m.at(j,i);
            double ratio = 0.0;
            if (fabs(max-min) >= 0.0000000001) ratio = 2.0 * (u-min) / (max - min);
            int b = int(MAX(0, 255*(1 - ratio)));
            int r = int(MAX(0, 255*(ratio - 1)));
            int g = 255 - b - r;
            QColor c(r, g, b);
            painter.setPen(c);
            QPoint pnt(static_cast<int>(i), static_cast<int>(rows-j-1));
            painter.drawPoint(pnt);
        }
    }
}

void visualHeatColorGradinet1(QPixmap &img, int w, int h)
{
    QLinearGradient gradient(0,0,w,h);
    gradient.setColorAt(0.00, Qt::blue);
    gradient.setColorAt(0.25, Qt::cyan);
    gradient.setColorAt(0.50, Qt::green);
    gradient.setColorAt(0.75, Qt::yellow);
    gradient.setColorAt(1.00, Qt::red);

    img = QPixmap(w,h);
    img.fill(Qt::transparent);
    QPainter painter(&img);
    painter.fillRect(0,0,w,h,QBrush(gradient));
    painter.end();
}

void visualHeatColorGradinet2(QPixmap &img, int w, int h)
{
    QLinearGradient gradient(0,0,w,h);
    gradient.setColorAt(0.000, Qt::black);
    gradient.setColorAt(0.167, Qt::blue);
    gradient.setColorAt(0.333, Qt::cyan);
    gradient.setColorAt(0.500, Qt::green);
    gradient.setColorAt(0.667, Qt::yellow);
    gradient.setColorAt(0.833, Qt::red);
    gradient.setColorAt(1.000, Qt::white);

    img = QPixmap(w,h);
    img.fill(Qt::transparent);
    QPainter painter(&img);
    painter.fillRect(0,0,w,h,QBrush(gradient));
    painter.end();
}

void visualGrayScale(const DoubleMatrix &m, double min, double max, QPixmap &pxm, size_t, size_t)
{
    unsigned int rows = m.rows();
    unsigned int cols = m.cols();

    pxm = QPixmap(static_cast<int>(cols), static_cast<int>(rows));
    pxm.fill(Qt::transparent);
    QPainter painter(&pxm);
    painter.setRenderHint(QPainter::Antialiasing, true);
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            int gray = static_cast<int>((1.0 - (m[j][i]-min)/(max-min))*255);
            QColor color(gray, gray, gray);
            painter.setPen(color);
            painter.drawPoint(static_cast<int>(i), static_cast<int>(rows-(j+1)));
        }
    }
}

void visualString(const DoubleVector& v, double min, double max, int w, int h, QPixmap &pxm, QColor bg, QColor fg, const QString &filename)
{
    pxm = QPixmap(w, h);
    pxm.fill(bg);

    double hx = static_cast<double>(w) / static_cast<double>(v.length());
    double hy = static_cast<double>(h) / static_cast<double>(max-min);

    QPainter painter(&pxm);

    painter.drawLine(0, h/2, w, h/2);

    painter.setPen(fg);
    painter.setRenderHint(QPainter::Antialiasing, true);
    for (unsigned int i=0; i<=v.length()-2;i++)
    {
        QPointF p1;
        p1.rx() = i*hx;
        p1.ry() = (max-v[i])*hy;
        QPointF p2;
        p2.rx() = (i+1)*hx;
        p2.ry() = (max-v[i+1])*hy;
        painter.drawLine(p1, p2);
    }

    if (!filename.isNull()) pxm.save(filename, "PNG");
}

void visualString1(const DoubleVector& v, double min, double max, int w, int h, QColor bg, QColor fg, const QString &filename)
{
    QPixmap pxm(w, h);
    pxm.fill(bg);

    double hx = static_cast<double>(w) / static_cast<double>(v.length());
    double hy = static_cast<double>(h) / static_cast<double>(max-min);

    QPainter painter(&pxm);
    painter.setPen(fg);
    painter.setRenderHint(QPainter::Antialiasing, true);
    for (unsigned int i=0; i<=v.length()-2;i++)
    {
        QPointF p1;
        p1.rx() = i*hx;
        p1.ry() = (max-v[i])*hy;
        QPointF p2;
        p2.rx() = (i+1)*hx;
        p2.ry() = (max-v[i+1])*hy;
        painter.drawLine(p1, p2);
    }
    if (!filename.isNull()) pxm.save(filename, "PNG");
}
