#ifndef GRID_H
#define GRID_H

#include <global.h>
#include <vector>

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
    SpaceNodePDE(int i=0, double x=0.0, int j=0, double y=0.0, int k=0, double z=0.0);
    SpaceNodePDE(const SpaceNodePDE &sn);

    int i;
    int j;
    int k;
};

struct MINIMUMSHARED_EXPORT TimeNodePDE : public TimeMoment
{
    TimeNodePDE(unsigned int i=0, double t = 0.0);
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
    inline PointNodeODE(double x = 0.0, int i = 0) : PointODE(x), i(i) {}

    int i;
};

class MINIMUMSHARED_EXPORT Dimension
{
public:
    enum SpaceDimension
    {
        DimensionX = 1,
        DimensionY = 2,
        DimensionZ = 3
    };

    Dimension(double step = 0.01, int min = 0, int max = 100);

    auto step() const -> double;
    auto setStep(double step) -> void;
    auto min() const -> int;
    auto setMin(int min) -> void;
    auto max() const -> int;
    auto setMax(int max) -> void;
    auto size() const -> int;

protected:
    double _step;
    int _min;
    int _max;
};

#endif // GRID_H
