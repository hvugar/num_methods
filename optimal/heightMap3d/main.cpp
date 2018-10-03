#include "widget.h"
#include <QApplication>
#include "heightmapsurface.h"
#include "matrixsurface.h"
#include <QMessageBox>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    //Q3DSurface *surface = new HeightMapSurface();
    Q3DSurface *surface = new MatrixSurface();

    if (!surface->hasContext()) {
        QMessageBox msgBox;
        msgBox.setText("Couldn't initialize the OpenGL context.");
        msgBox.exec();
        return -1;
    }

    QWidget *w = QWidget::createWindowContainer(surface);
    w->show();

    return a.exec();
}
