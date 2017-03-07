#include <QApplication>

#include "widget2.h"
#include "heatimager.h"
#include "sample1.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    for (int i=0; i<argc; i++)
    {
        if (QString(argv[i]).compare("dim2") == 0)
        {
            createHeatImage1(argc, argv);
        }
        if (QString(argv[i]).compare("dim1") == 0)
        {
            createHeatImage2(argc, argv);
        }
    }
//    QSample1 smp;
//    smp.show();

    return 0;
}


