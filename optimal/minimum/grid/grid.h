#ifndef GRID_H
#define GRID_H

#include <global.h>
#include <vector>

struct MINIMUMSHARED_EXPORT GridNodeODE
{
    inline GridNodeODE(double x, int i) : x(x), i(i) {}
    double x;
    int i;
};

struct MINIMUMSHARED_EXPORT SpaceNodePDE
{
    double x;
    double y;
    double z;

    unsigned int i;
    unsigned int j;
    unsigned int k;
};

struct MINIMUMSHARED_EXPORT TimeNodePDE
{
    double t;
    unsigned int i;
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

    Dimension(double step=0.01, int minN = 0, int maxN=100);

    virtual double step() const;
    virtual int minN() const;
    virtual int maxN() const;
    virtual int sizeN() const;

protected:
    double mstep;
    int mminN;
    int mmaxN;
};

class MINIMUMSHARED_EXPORT UniformODEGrid
{
public:
    UniformODEGrid(double step = 0.0, int min = 0, int max = 0);

    double step() const;
    int minN() const;
    int maxN() const;
    int sizeN() const;
private:
    double mstep;
    int mminN;
    int mmaxN;
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
