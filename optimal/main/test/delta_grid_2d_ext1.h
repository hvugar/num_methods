#ifndef DELTA_GRID_2D_EXT1_H
#define DELTA_GRID_2D_EXT1_H

#include <deltagrid.h>

class DeltaGrid1DExt1
{
public:
    static void Main(int argc, char ** argv);

    double fx(double x) const;
    double dx(double x) const;
};

class DeltaGrid2DExt1
{
public:
    static void Main(int argc, char ** argv);

    double fx(double x, double y) const;
    double dx(double x, double y) const;
    double dy(double x, double y) const;
};

#endif // DELTA_GRID_2D_EXT1_H
