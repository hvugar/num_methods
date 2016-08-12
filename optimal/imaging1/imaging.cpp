#include "imaging.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void visualizeVectorHeat(const DoubleVector &v, double min, double max, QPixmap &img, unsigned int w, unsigned int h)
{
    C_UNUSED(w);
    unsigned int size = v.size();
    img = QPixmap(size, h);
    img.fill(Qt::transparent);
    QPainter painter(&img);
    for (unsigned int i=0; i<v.size(); i++)
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
    img = QPixmap(cols, rows);
    img.fill(Qt::transparent);
    QPainter painter(&img);
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            double u = m.at(j,i);
            double ratio = 0.0;
            if (min!=max) ratio = 2.0 * (u-min) / (max - min);
            int b = int(MAX(0, 255*(1 - ratio)));
            int r = int(MAX(0, 255*(ratio - 1)));
            int g = 255 - b - r;
            QColor c(r, g, b);
            painter.setPen(c);
            painter.drawPoint(i, rows-j-1);
        }
    }
}
