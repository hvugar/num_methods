#ifndef GRID_H
#define GRID_H

#include <global.h>
#include <vector>

class MINIMUMSHARED_EXPORT SpaceDimension
{
public:
    enum Dimension
    {
        Dim1 = 1,
        Dim2 = 2,
        Dim3 = 3
    };

    SpaceDimension(double hx=0.01, double x1=0.0, double x2=1.0, unsigned int N1 = 0, unsigned int N2=100);

    double hx() const;
    double x1() const;
    double x2() const;
    unsigned int N1() const;
    unsigned int N2() const;

private:
    double _hx;
    double _x1;
    double _x2;
    unsigned int _N1;
    unsigned int _N2;
};

class MINIMUMSHARED_EXPORT TimeDimension
{
public:
    TimeDimension(double ht=0.01, double t1=0.0, double t2=1.0, unsigned int M1 = 0, unsigned int M2=100);

    double ht() const;
    double t1() const;
    double t2() const;
    unsigned int M1() const;
    unsigned int M2() const;

private:
    double _ht;
    double _t1;
    double _t2;
    unsigned int _M1;
    unsigned int _M2;
};

class MINIMUMSHARED_EXPORT GridPDE
{
public:
    const TimeDimension& timeDimension() const;
    void setTimeDimension(const TimeDimension &timeDimension);

    void addSpaceDimension(const SpaceDimension &spaceDimSize);
    unsigned int spaceDimSize() const;
    const SpaceDimension& spaceDimensions(SpaceDimension::Dimension dimension) const;

private:
    std::vector<SpaceDimension> mSpaceDimensions;
    TimeDimension mTimeDimension;
};

class MINIMUMSHARED_EXPORT GridODE
{

};


#endif // GRID_H
