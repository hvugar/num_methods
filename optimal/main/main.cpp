#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gradient_prj.h>
#include <methods.h>
#include <gridmethod.h>

#include "cfunction1.h"
#include "cfunction2.h"
#include "cfunction.h"
#include "rosenbrock.h"

int main()
{
//    Rosenbrock::main();
//    CFunction1::main();
//    CFunction2::main();
    ControlFunction::main();

    return 0;
}

