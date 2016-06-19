#include "headers.h"

#include "widget/qsimplewavewidget.h"
#include <matrix.h>
#include <time.h>

int main(int argc, char ** argv)
{
    srand(time(NULL));

    typedef union {
      unsigned long long i64;
      double d64;
    } dbl_64;

    double value = 0.0;
    dbl_64 s;
    s.i64 = 0x00000001;
    s.d64 = 0.0000000000000000001;
    //printf("%llu\n", s.i64);
    printf("%.55f\n", s.d64);

//    SampleBorderHyperBolic::main(argc, argv);
//    QSimpleWaveWidget::main(argc, argv);


//    CVector v(100);
//    v[0] = 1.0;
//    printf("%f\n", v[0]);

//    struct Matrix *m1 = matrix_new(2,5);
//    struct Matrix *m2 = matrix_new(5,4);
//    matrix_rand(m1);
//    matrix_rand(m2);
//    matrix_print(m1);
//    puts("-----");
//    matrix_print(m2);
//    puts("-----");

//    struct Matrix *mult = matrix_mult(m1, m2);
//    matrix_print(mult);
//    matrix_free(mult);

//    matrix_free(m2);
//    matrix_free(m1);

    return 0;
}
