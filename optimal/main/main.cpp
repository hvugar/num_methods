#include "headers.h"


class A : public R1Function
{
    virtual double fx(double x);
};

double A::fx(double x) {
    return -200*(x-0.3)*(x-0.3)*(x-0.3)*(x-0.3)+0.2;
}

int main(int argc, char ** argv)
{
//    HeatControl2DeltaF::main(argc, argv);
//    HeatControlDeltaF::main();
//    HeatControlDeltaX::main(argc, argv);
//    DiscreteHyperbolic1::main();
//    HyperbolicControlH::main();
//    HyperbolicControl2D1::main();
//    puts("--------------------------------------");
//    HyperbolicControl2DMX::main();
//    puts("--------------------------------------");
//    HyperbolicControl2DMV::main();
//    puts("--------------------------------------");
//    BorderHyperbolic::main();
//    BorderParabolic2D::main();
    HyperbolicControl2D24::main(argc, argv);

//    double a = 0.14;
//    double b = 0.41;
//    double x = 0.0;
//    A f;
//    goldenSectionSearch(a, b, x, &f, 0.0001);
//    printf("%.8f %.8f %.8f\n", a, b, x);

    return 0;
}
