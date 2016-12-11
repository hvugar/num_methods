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
    srand(time(NULL));

    //BorderTest::Main(argc, argv);
    //BorderTest1::Main(argc, argv);
    //Example2::Main(argc, argv);
    //Problem1::Main(argc, argv);
    //Example3 e;
    Example4::Main(argc, argv);

//    DoubleMatrix m(5,5);
//    m.at(0,0) = 0.0; m.at(0,1) = 0.0; m.at(0,2) = 0.0; m.at(0,3) = 1.0; m.at(0,4) = 4.0;
//    m.at(1,0) = 1.0; m.at(1,1) = 3.0; m.at(1,2) = 2.0; m.at(1,3) = 0.0; m.at(1,4) = 1.0;
//    m.at(2,0) = 5.0; m.at(2,1) = 0.0; m.at(2,2) = 1.0; m.at(2,3) = 2.0; m.at(2,4) = 0.0;
//    m.at(3,0) = 0.0; m.at(3,1) = 5.0; m.at(3,2) = 1.0; m.at(3,3) = 4.0; m.at(3,4) = 0.0;
//    m.at(4,0) = 2.0; m.at(4,1) = 0.0; m.at(4,2) = 0.0; m.at(4,3) = 1.0; m.at(4,4) = 2.0;
//    double d1 = m.determinant1();
//    printf("%14.10f\n", d1);
    return 0;
}
