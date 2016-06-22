#include "headers.h"

#include "widget/qsimplewavewidget.h"
#include <matrix.h>
#include <time.h>

int main(int argc, char ** argv)
{
    srand(time(NULL));

    SampleBorderHyperBolic::main(argc, argv);
    QSimpleWaveWidget::main(argc, argv);

    return 0;
}
