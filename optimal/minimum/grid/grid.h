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

struct MINIMUMSHARED_EXPORT SpacePoint
{
    double x;
    double y;
    double z;
};

//struct MINIMUMSHARED_EXPORT SpaceGridNode
//{
//    unsigned int i;
//    unsigned int j;
//    unsigned int k;
//};

struct MINIMUMSHARED_EXPORT SpaceNodePDE
{
    SpaceNodePDE();
    SpaceNodePDE(unsigned int i, double x);
    SpaceNodePDE(unsigned int i, unsigned int j, double x, double y);
    SpaceNodePDE(unsigned int i, unsigned int j, unsigned int k, double x, double y, double z);
    SpaceNodePDE(const SpaceNodePDE& node);

    void setPoint(double x, double y);

    unsigned int i;
    double x;
    unsigned int j;
    double y;
    unsigned int k;
    double z;
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
    //Dimension(const Dimension &dim);
    //virtual ~Dimension() {}

    //Dimension& operator=(const Dimension &dim) const;

    double step() const;
    int minN() const;
    int maxN() const;
    int sizeN() const;

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

class UniformPDEGrid
{
public:
    MINIMUMSHARED_EXPORT UniformPDEGrid();
    MINIMUMSHARED_EXPORT UniformPDEGrid(const Dimension &timeDimension, std::vector<Dimension> &spaceDimensions);

    MINIMUMSHARED_EXPORT void setTimeDimension(const Dimension &timeDimension);
    MINIMUMSHARED_EXPORT void addSpaceDimension(const Dimension &spaceDimension);

    MINIMUMSHARED_EXPORT const Dimension& timeDimension() const;
    MINIMUMSHARED_EXPORT const Dimension& spaceDimension(Dimension::SpaceDimension dimension) const;
private:
    Dimension mtimeDimension;
    std::vector<Dimension> mspaceDimensions;
};

#endif // GRID_H
