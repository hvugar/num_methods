#include "headers.h"

#include "widget/qsimplewavewidget.h"
#include <matrix.h>

int main(int argc, char ** argv)
{
//    SampleBorderHyperBolic::main(argc, argv);
//    QSimpleWaveWidget::main(argc, argv);

    srand(3000);

    struct Matrix *m1 = matrix_new(2,5);
    struct Matrix *m2 = matrix_new(5,4);
    matrix_rand(m1);
    matrix_rand(m2);
    matrix_print(m1);
    puts("-----");
    matrix_print(m2);
    puts("-----");

    struct Matrix *mult = matrix_mult(m1, m2);
    matrix_print(mult);
    matrix_free(mult);

    matrix_free(m2);
    matrix_free(m1);

    return 0;
}
