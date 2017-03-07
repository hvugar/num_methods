#include "heatimager.h"

void MatrixHeatImaging(const DoubleMatrix &m, double minimum, double maximum, QPixmap &pixmap, int width, int height)
{
    pixmap = QPixmap(QSize(width, height));
    pixmap.fill(Qt::white);
    QPainter painter(&pixmap);

    unsigned int rows = m.rows();
    unsigned int cols = m.cols();

    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            double u = m.at(j,i);
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

    qDebug() << width << height << inFile << outFile << min << max;

    QFile file(inFile);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);

    double maximum = -9999.9;
    double minimum = +9999.9;

    DoubleMatrix m(height, width);

    unsigned int j=0;
    QString line = in.readLine();
    while (!line.isNull() && !line.isEmpty())
    {
        unsigned int i=0;
        QStringList list = line.split(" ");
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

    minimum = -0.3892016224094016;
    maximum = +0.3951541064046917;

    printf("File: %s Minimum: %.10f Maximum: %.10f width: %d height %d\n", inFile.toLatin1().data(), minimum, maximum, width, height);


    QPixmap pixmap;
    MatrixHeatImaging(m, minimum, maximum, pixmap, width, height);
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

    printf("Minimum: %.10f Maximum: %.10f width: %d output: %s\n", minimum, maximum, width, outFile.toLatin1().data());
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

