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
    Dimension(double step=0.01, unsigned int maxN=100, unsigned int minN = 0);

    virtual double step() const;
    virtual double min() const;
    virtual double max() const;
    virtual unsigned int minN() const;
    virtual unsigned int maxN() const;
    virtual unsigned int sizeN() const;

protected:
    double mstep;
    unsigned int mmaxN;
    unsigned int mminN;
    double mmin;
    double mmax;
};

class MINIMUMSHARED_EXPORT SpaceDimension : public Dimension
{
public:
    enum DimensionSize
    {
        Dim1 = 0,
        Dim2 = 1,
        Dim3 = 2
    };

    SpaceDimension(double step=0.01, unsigned int maxN=100, unsigned int minN = 0);
};

class MINIMUMSHARED_EXPORT TimeDimension : public Dimension
{
public:
    TimeDimension(double step=0.01, unsigned int maxN=100, unsigned int minN = 0);
};

class MINIMUMSHARED_EXPORT GridPDE
{
public:
    const TimeDimension& timeDimension() const;
    void setTimeDimension(const TimeDimension &timeDimension);

    void addSpaceDimension(const SpaceDimension &spaceDimSize);
    unsigned int spaceDimSize() const;
    const SpaceDimension& spaceDimensions(SpaceDimension::DimensionSize dimension) const;

private:
    std::vector<SpaceDimension> mSpaceDimensions;
    TimeDimension mTimeDimension;
};

class MINIMUMSHARED_EXPORT GridODE
{
public:
    const Dimension& dimension() const;
    void setDimension(const Dimension &dimension);
private:
    Dimension mdimension;
};


#endif // GRID_H
