#include "problem1.h"
#include "problem1k.h"
#include "problem1z.h"
#include "problem1kz.h"
#include "problem3.h"

#include <QtGui/QGuiApplication>
#include <imaging.h>

#include <float.h>

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    //srand(time(NULL));
    Problem1K::Main(argc, argv);
    return 0;
}
