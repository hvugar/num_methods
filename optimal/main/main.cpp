#include "problem1.h"
#include "problem1k.h"
#include "problem1z.h"
#include "problem1x.h"
#include "problem1x1.h"
#include "problem1kz.h"
#include "problem1kzx.h"

#include "problem3.h"
#include "loadedsystems.h"
#include "example1.h"
#include "example2.h"
#include "example3.h"
#include <cmethods.h>

#include <QtGui/QGuiApplication>
//#include <imaging.h>

#include <float.h>
#include <time.h>

int main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    srand(time(NULL));

//    Problem1KZX::Main(argc, argv);
//    Example3::Main(argc, argv);
    Example2::Main(argc, argv);
//    Problem1K::Main(argc, argv);

    return 0;
}
