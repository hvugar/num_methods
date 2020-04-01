#ifndef GRID_H
#define GRID_H

#include <global.h>
#include <vector>
#include <math.h>
#include <cstdio>
#include <float.h>
#include <math.h>
#include <limits>
#include <exception>
#include <stdexcept>
#include "../cmethods.h"
#include "../linearequation.h"
#include "../printer.h"
#include "../matrix2d.h"
#include "../matrix3d.h"


enum class FiniteDifference
{
    Explicit = 0,
    Implicit = 1
};

struct MINIMUMSHARED_EXPORT SpacePoint
{
    SpacePoint(double x = 0.0, double y = 0.0, double z = 0.0);

    double x;
    double y;
    double z;
};

struct MINIMUMSHARED_EXPORT TimeMoment
{
    TimeMoment(double t = 0.0);
    double t;
};

struct MINIMUMSHARED_EXPORT SpaceNodePDE : public SpacePoint
{
    SpaceNodePDE(int i = 0, double x = 0.0, int j = 0, double y = 0.0, int k = 0, double z = 0.0);
    SpaceNodePDE(const SpaceNodePDE &sn);

    int i;
    int j;
    int k;
};

struct MINIMUMSHARED_EXPORT TimeNodePDE : public TimeMoment
{
    explicit TimeNodePDE(unsigned int i=0, double t = 0.0);
    TimeNodePDE(const TimeNodePDE &tn);
    unsigned int i;
};

struct MINIMUMSHARED_EXPORT PointODE
{
    PointODE(double x = 0.0);

    double x;
};

struct MINIMUMSHARED_EXPORT PointNodeODE : public PointODE
{
    PointNodeODE(double x = 0.0, int i = 0) : PointODE(x), i(i) {}

    int i;
};

class MINIMUMSHARED_EXPORT Dimension
{
public:
//    enum SpaceDimension
//    {
//        DimensionX = 1,
//        DimensionY = 2,
//        DimensionZ = 3
//    };

    Dimension(double step = 0.01, int min = 0, int max = 100);
    Dimension(const Dimension &);
    Dimension & operator =(const Dimension &);

    double step() const;
    void setStep(double step);
    int min() const;
    void setMin(int min);
    int max() const;
    void setMax(int max);
    unsigned int size() const;

    Dimension& incMax();
    Dimension& decMax();
    Dimension& addMinNode();
    Dimension& delMinNode();

protected:
    double _step;
    int _min;
    int _max;
};

#endif // GRID_H
