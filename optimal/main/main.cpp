#include "headers.h"
#include <time.h>
#include "gaussianelimination.h"
#include "problem1.h"
#include "problem2.h"

#include <QtGui/QGuiApplication>
#include <imaging.h>

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    srand(time(NULL));

    Problem2::Main(argc, argv);

    return 0;
}
