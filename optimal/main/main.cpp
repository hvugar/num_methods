#include "problem1.h"
#include "loadedsystems.h"
#include "example1.h"
#include "example2.h"
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
    Problem1::Main(argc, argv);

    return 0;
}
