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

struct MINIMUMSHARED_EXPORT SpaceNodePDE : public SpacePoint
{
    SpaceNodePDE(int i=0, double x=0.0, int j=0, double y=0.0, int k=0, double z=0.0);
    SpaceNodePDE(const SpaceNodePDE &sn);

    int i;
    int j;
    int k;
};

struct MINIMUMSHARED_EXPORT TimeMoment
{
    TimeMoment(double t = 0.0);
    double t;
};

struct MINIMUMSHARED_EXPORT TimeNodePDE : public TimeMoment
{
    TimeNodePDE(unsigned int i=0, double t = 0.0);
    TimeNodePDE(const TimeNodePDE &tn);
    unsigned int i;
};

struct MINIMUMSHARED_EXPORT GridNodeODE
{
    inline GridNodeODE(double x, int i) : x(x), i(i) {}
    double x;
    int i;
};


class MINIMUMSHARED_EXPORT Dimension
{
public:
    enum SpaceDimension
    {
        DimensionX = 0,
        DimensionY = 1,
        DimensionZ = 2
    };

    Dimension(double step = 0.01, int min = 0, int max = 100);

    auto step() const -> double;
    auto min() const -> int;
    auto max() const -> int;
    auto size() const -> int;

protected:
    double m_step;
    int m_min;
    int m_max;
};

class MINIMUMSHARED_EXPORT UniformODEGrid
{
public:
    UniformODEGrid();
    explicit UniformODEGrid(const Dimension &dimension);
    virtual ~UniformODEGrid();

    const Dimension &dimension() const;

private:
    Dimension _dimension;
};

class MINIMUMSHARED_EXPORT UniformPDEGrid
{
public:
    UniformPDEGrid();
    UniformPDEGrid(const Dimension &timeDimension, std::vector<Dimension> &spaceDimensions);

    void setTimeDimension(const Dimension &timeDimension);
    void addSpaceDimension(const Dimension &spaceDimension);

    const Dimension& timeDimension() const;
    const Dimension& spaceDimension(Dimension::SpaceDimension dimension) const;
private:
    Dimension mtimeDimension;
    std::vector<Dimension> mspaceDimensions;
};

#endif // GRID_H
