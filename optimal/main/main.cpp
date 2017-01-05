#include "problem1.h"
#include "loadedsystems.h"
#include "example1.h"
#include "example2.h"
#include "example3.h"
#include "example4.h"
#include "example5.h"
#include "problem1.h"
#include <cmethods.h>

#include <../border/borderparabolic.h>
#include <../border/borderparabolic2d.h>
#include <../border/borderhyperbolic2d.h>

#include <bordertest1.h>
#include <sampleboundaryproblem1.h>

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
    //Example4::Main(argc, argv);
    //BorderParabolic::Main(argc, argv);
    //BorderParabolic2D::Main(argc, argv);
    //BorderHyperbolic2D::Main(argc, argv);
    //Example5 e5;
    SampleBoundaryProblem1::Main(argc, argv);
    return 0;
}
