#include <QApplication>
#include <QPixmap>
#include <QPainter>
#include <QTextStream>
#include <QFile>
#include <QStringList>
#include <QDebug>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv, false);

    QPixmap pixmap(QSize(2000, 2000));
    pixmap.fill(Qt::white);

    QFile file("data.txt");
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);
    QPainter painter(&pixmap);

    QString line = in.readLine();
    unsigned int j=0;
    while (!line.isNull())
    {
        QStringList list = line.split(" ");

        for (unsigned int i=0; i<list.size(); i++)
        {
            painter.setPen((int)(list[i].toDouble()*0x02ffff));
            painter.drawPoint(i,j);

        }

        line = in.readLine();
        j++;

    }
    file.close();

    pixmap.save("image.png");

    return 0;
}
