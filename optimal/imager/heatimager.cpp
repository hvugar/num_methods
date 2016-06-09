#include "heatimager.h"

void MatrixHeatImaging(const DoubleMatrix &m, double minimum, double maximum, QPixmap &pixmap, int width, int height)
{
    pixmap = QPixmap(QSize(width, height));
    pixmap.fill(Qt::white);
    QPainter painter(&pixmap);

    int m_size = m.size();

    for (int j=0; j<m_size; j++)
    {
        const DoubleVector &v = m[j];
        int v_size = v.size();
        for (int i=0; i<v_size; i++)
        {
            double u = v[i];
            double ratio = 0.0;
            if (minimum!=maximum)
                ratio = 2.0 * (u-minimum) / (maximum - minimum);
            int b = int(MAX(0, 255*(1 - ratio)));
            int r = int(MAX(0, 255*(ratio - 1)));
            int g = 255 - b - r;
            QColor c(r, g, b);
            painter.setPen(c);
            painter.drawPoint(i,height-j-1);
        }
    }
}

void createHeatImage1(int argc, char *argv[])
{
    if (argc < 9)
    {
        printf("Usage: imager.exe -w 100 -h 100 -i filename.txt -o filename.png -min 0.0 -max 5.0\n");
        return;
    }

    unsigned int width = QString(argv[2]).toInt();
    unsigned int height = QString(argv[4]).toInt();
    QString inFile = QString(argv[6]);
    QString outFile = QString(argv[8]);
    double min = QString(argv[10]).toDouble();
    double max = QString(argv[12]).toDouble();

    //    QPixmap pixmap(QSize(width, height));
    //    pixmap.fill(Qt::white);
    //    QPainter painter(&pixmap);
    //    painter.drawPoint(0,0);

    QFile file(inFile);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);

    double maximum = -9999.9;
    double minimum = +9999.9;

    DoubleMatrix m;
    m.resize(height);

    unsigned int j=0;
    QString line = in.readLine();
    while (!line.isNull() && !line.isEmpty())
    {
        unsigned int i=0;
        QStringList list = line.split(" ");
        m[j].resize(width);
        for (int k=0; k<list.size(); k++)
        {
            QString str = list[k];
            if (!str.isNull() && !str.isEmpty())
            {
                double u = str.toDouble();
                m[j][i] = u;
                if (u<=minimum) minimum = u;
                if (u>=maximum) maximum = u;
                i++;
            }
        }
        line = in.readLine();
        j++;
    }
    file.close();

    if (min != max)
    {
        minimum = min;
        maximum = max;
    }
    printf("File: %s Minimum: %.10f Maximum: %.10f width: %d height %d\n", inFile.toAscii().data(), minimum, maximum, width, height);

    QPixmap pixmap;
    MatrixHeatImaging(m, minimum, maximum, pixmap, width, height);

//    for (int j=0; j<m.size(); j++)
//    {
//        for (int i=0; i<m[j].size(); i++)
//        {
//            double u = m[j][i];
//            double ratio = 0.0;
//            if (minimum!=maximum)
//                ratio = 2.0 * (u-minimum) / (maximum - minimum);
//            int b = int(MAX(0, 255*(1 - ratio)));
//            int r = int(MAX(0, 255*(ratio - 1)));
//            int g = 255 - b - r;
//            QColor c(r, g, b);
//            painter.setPen(c);
//            painter.drawPoint(i,height-j-1);
//        }
//    }

//    double sum = 0.0;
//    for (unsigned int j=0; j<height; j++)
//    {
//        for (unsigned int i=0; i<width; i++)
//        {
//            //            double k = 1.0;
//            //            if (i==0 || i==width-1)  k *= 0.5;
//            //            if (j==0 || j==height-1) k *= 0.5;
//            //            sum = sum + k * (m[j][i]) * (m[j][i]);
//            sum = sum + m[j][i];
//        }
//    }
//    double h1 = 1.0 / (height-1);
//    double h2 = 1.0 / (width-1);
//    //    sum = h1*h2*sum;
//    painter.setPen(Qt::black);
//    painter.drawText(10, 10, QString::number(sum, 'f', 6));

    pixmap.save(outFile, "PNG");
}

void createHeatImage2(int argc, char *argv[])
{
    if (argc < 7)
    {
        printf("Usage: images.exe -w 100 -i filename.txt -o filename.png\n");
        return;
    }

    unsigned int width = QString(argv[2]).toInt();
    QString inFile = QString(argv[4]);
    QString outFile = QString(argv[6]);

    QPixmap pixmap(QSize(width, 10));
    pixmap.fill(Qt::white);
    QPainter painter(&pixmap);

    QFile file(inFile);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);

    double maximum = 6.6724519545;
    double minimum = 0.0;

    QDoubleVector m;
    m.resize(width);

    QString line = in.readLine();
    QStringList list = line.split(" ");
    unsigned int i=0;
    for (int j=0; j<list.size(); j++)
    {
        QString str = list[j];
        if (!str.isNull() && !str.isEmpty())
        {
            double u = str.toDouble();
            m[i] = u;
            if (u<=minimum) minimum = u;
            if (u>=maximum) maximum = u;
            i++;
        }
    }
    file.close();

    printf("Minimum: %.10f Maximum: %.10f width: %d output: %s\n", minimum, maximum, width, outFile.toAscii().data());
    //minimum = 1.0;

    for (int i=0; i<m.size(); i++)
    {
        double u = m[i];

        //unsigned int c = (u-minimum)*(maximum-minimum)*0xffffff;
        double ratio = 0.0;
        if (minimum!=maximum)
            ratio = 2.0 * (u-minimum) / (maximum - minimum);
        int b = int(MAX(0, 255*(1 - ratio)));
        int r = int(MAX(0, 255*(ratio - 1)));
        int g = 255 - b - r;
        QColor c(r, g, b);
        painter.setPen(c);
        painter.drawLine(i,0,i,10);
    }

    pixmap.save(outFile, "PNG");
}

