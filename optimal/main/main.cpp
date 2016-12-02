#include "problem1.h"
#include "loadedsystems.h"
#include "example1.h"
#include "example2.h"
#include "example3.h"
#include "example4.h"
#include "problem1.h"
#include <cmethods.h>
#include <bordertest1.h>

#include <QtGui/QGuiApplication>
//#include <imaging.h>

#include <float.h>
#include <time.h>

#include "bordertest.h"

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    //srand(time(NULL));

    //BorderTest::Main(argc, argv);
    //BorderTest1::Main(argc, argv);
    //Example2::Main(argc, argv);
    //Problem1::Main(argc, argv);
    //Example3 e;
    Example4::Main(argc, argv);

//    DoubleVector v;
//    v << 1.1 << 2.2 << 3.3;
//    DoubleMatrix m(v);
//    printf("%f %f %f\n", v.at(0), v.at(1), v.at(2));
//    printf("%d %d\n", m.rows(), m.cols());
//    printf("%f %f %f\n", m.at(1,0), m.at(2,0), m.at(2,0));

    return 0;
}
