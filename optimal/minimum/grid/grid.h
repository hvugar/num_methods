#ifndef GRID_H
#define GRID_H

#include <global.h>
#include <vector>

struct MINIMUMSHARED_EXPORT SpaceNode
{
    unsigned int i;
    unsigned int j;
    unsigned int k;

    double x;
    double y;
    double z;
};

struct MINIMUMSHARED_EXPORT TimeNode
{
    unsigned int i;
    double t;
};

class MINIMUMSHARED_EXPORT Dimension
{
public:
    enum SpaceDimension
    {
        Dim1 = 0,
        Dim2 = 1,
        Dim3 = 2
    };

    Dimension(double step=0.01, unsigned int maxN=100, unsigned int minN = 0);

    virtual double step() const;
    virtual unsigned int minN() const;
    virtual unsigned int maxN() const;
    virtual unsigned int sizeN() const;

protected:
    double mstep;
    unsigned int mmaxN;
    unsigned int mminN;
};

class MINIMUMSHARED_EXPORT ODEGrid
{
public:
    ODEGrid();
    explicit ODEGrid(const Dimension &dimension);

    const Dimension &dimension() const;
    void setDimension(const Dimension &dimension);
private:
    Dimension mdimension;
};

class MINIMUMSHARED_EXPORT PDEGrid
{
public:
    PDEGrid(const Dimension &timeDimension, std::vector<Dimension> &spaceDimension);

    void setTimeDimension(const Dimension &timeDimension);
    void addSpaceDimension(const Dimension &spaceDimension);

    const Dimension &timeDimension() const;
    const Dimension &spaceDimension() const;
private:
    Dimension mtimeDimension;
    std::vector<Dimension> mspaceDimensions;
};

#endif // GRID_H
