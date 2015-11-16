#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gridmethod.h>
#include <rungekutta.h>
#include <doublevector.h>
#include <parabolicequation.h>
#include <hyperbolicequation.h>
#include <r1minimize.h>

#include <iostream>
#include <stdexcept>

#include "cfunction1.h"
#include "cfunction2.h"
#include "cfunction3.h"
#include "cfunction.h"
#include "rosenbrock.h"
#include "bealesfunction.h"
#include "boothfunction.h"
#include "pointcontrol.h"
#include "pointcontrol1.h"
#include "pointcontrol2.h"
#include "utils.h"

#include "heat/1d/heatcontrol.h"
#include "heat/1d/heatcontrol1.h"
#include "heat/1d/heatcontroldeltaf.h"
#include "heat/1d/heatcontroldeltax.h"
#include "heat/2d/heatcontrol2d.h"
#include "heat/2d/heatcontrol2delta.h"
#include "heat/2d/heatcontrol2deltaf.h"
#include "heat/2d/heatcontrol2deltax.h"

#include "hyperbolic/hyperbolic1dx.h"
#include "hyperbolic/hyperboliccontrol1d2.h"
#include "hyperbolic/hyperboliccontrol1d3.h"
#include "hyperbolic/hyperboliccontrol1dt.h"

struct MyFunc : public R1Function {
    double fx(double x) { return x*x; }
};

int main()
{
    MyFunc fx;

    double a,b,x;
    try {
        stranghLineSearch(4.0, 0.001, a, b, &fx);
        printf("a: %f b: %f\n", a, b);
        goldenSectionSearch(a, b, x, &fx, 0.00001);
        printf("a: %f b: %f x: %f\n", a, b, x);
    } catch (std::invalid_argument &ex) {
        std::cout << "Invalid argument: " << ex.what() << std::endl;
    } catch (std::runtime_error &ex) {
        std::cout << "Runtime error: " << ex.what() << std::endl;
    }

    std::cout << "(after exception)\n";


//    std::vector<unsigned char> v;
//    v.resize(10);
//    printf("0x%X %d 0x%X\n", v.data(), v.size(), &v);
//    v.resize(20);
//    printf("0x%X %d 0x%X\n", v.data(), v.size(), &v);
//    v.resize(30);
//    printf("0x%X %d 0x%X\n", v.data(), v.size(), &v);
//    printf("%u\n", v.max_size());
    //DblVector v(10);
//    printf("Size: %d\n", v.size());
//    for (unsigned int i=0; i<v.size(); i++)
//    {
//        v.data()[i] = i+1.0;
//        printf("%2.1f\n", v.at(i));
//    }

//    puts("add");
//    v.add(11.0);
//    printf("Size: %d\n", v.size());
//    v.data()[10]=11.0;
//    for (unsigned int i=0; i<v.size(); i++)
//    {
//        printf("%2.1f\n", v.at(i));
//    }

//    puts("insert");
//    v.insert(11, 12.0);
//    printf("Size: %d\n", v.size());
//    for (unsigned int i=0; i<v.size(); i++)
//    {
//        printf("%2.1f\n", v.at(i));
//    }

//    puts("delete");
//    v.remove(11);
//    printf("Size: %d\n", v.size());
//    for (unsigned int i=0; i<v.size(); i++)
//    {
//        printf("%2.1f\n", v.at(i));
//    }

    HyperbolicControl1DT::main();
    return 0;
}



