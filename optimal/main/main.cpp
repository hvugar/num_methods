#include "headers.h"
#include <time.h>
#include "gaussianelimination.h"
#include "problem1.h"
#include "problem2.h"

#include <QtGui/QGuiApplication>
#include <imaging.h>

union Comp
{
    float d;
    unsigned int u;
};

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    srand(time(NULL));

//    Problem2::Main(argc, argv);

    Comp c;
    for (unsigned int i=0x3F0FFFFF; i>0x3F0FFFF0; i--)
    {
        c.u = i;
        printf("%.8f %X %u\n", c.d, c.u, c.u);
    }

    return 0;
}
