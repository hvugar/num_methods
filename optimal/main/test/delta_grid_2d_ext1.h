#ifndef DELTA_GRID_2D_EXT1_H
#define DELTA_GRID_2D_EXT1_H

#include <deltagrid.h>
#include <integral.h>

class DeltaGrid1DExt1
{
public:
    static void Main(int argc, char ** argv);
    static void example1();
    static void example2();
    static void example3();

    double fx(double x) const;
    double dx(double x) const;
};

class DeltaGrid2DExt1
{
public:
    static void Main(int argc, char ** argv);
    static void example1();
    static void example2();

    double fx(double x, double y) const;
    double dx(double x, double y) const;
    double dy(double x, double y) const;
};

#endif // DELTA_GRID_2D_EXT1_H
