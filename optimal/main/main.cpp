#include "headers.h"

#include "control/example1.h"
#include "control/example2.h"

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    Example2::main(argc, argv);

//    ControlFunction4::main(argc, argv);
//    SampleMain();
//    SampleLoaderBorder::main();
//    HyperbolicControl2D21::main(argc, argv);
//    Rosenbrock::main(argc, argv);
    return 0;
}
