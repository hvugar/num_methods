#include "example31.h"
#include "example32.h"
#include "example33.h"
#include "example34.h"
#include "example35.h"
#include "example36.h"
#include "example37.h"
#include "example38.h"
#include "example335.h"
#include "example336.h"

#define SAMPLE_31
//#define SAMPLE_32
//#define SAMPLE_33
//#define SAMPLE_34
//#define SAMPLE_35
//#define SAMPLE_36
//#define SAMPLE_37
//#define SAMPLE_371
//#define SAMPLE_335
//#define SAMPLE_336


int main(int argc, char *argv[])
{
#ifdef SAMPLE_31
    puts("Sample 3.1");
    BorderParabolic1D31::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
#ifdef SAMPLE_32
    puts("Sample 3.2");
    BorderParabolic1D32::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
#ifdef SAMPLE_33
    puts("Sample 3.3");
    BorderParabolic2D33::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
#ifdef SAMPLE_34
    puts("Sample 3.4");
    BorderParabolic2D34::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
#ifdef SAMPLE_35
    puts("Sample 3.5");
    Parabolic1DControl35::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
#ifdef SAMPLE_36
    puts("Sample 3.6");
    Parabolic1DControl36::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
#ifdef SAMPLE_37
    puts("Sample 3.7");
    Parabolic2DControl37::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
#ifdef SAMPLE_371
    puts("Sample 3.71");
    Parabolic1DControl38::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
#ifdef SAMPLE_335
    puts("Sample 335");
    Parabolic1DControl335::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
#ifdef SAMPLE_336
    puts("Sample 336");
    Parabolic1DControl336::main(argc, argv);
    puts("----------------------------------------------------------------------------");
#endif
}
