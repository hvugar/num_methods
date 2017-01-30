#include "grid.h"

SpaceDimension::SpaceDimension(double h, double x1, double x2, unsigned int N1, unsigned int N2)
    : _hx(h), _x1(x1), _x2(x2), _N1(N1), _N2(N2)
{}

double SpaceDimension::hx() const
{
    return _hx;
}

double SpaceDimension::x1() const
{
    return _x1;
}

double SpaceDimension::x2() const
{
    return _x2;
}

unsigned int SpaceDimension::N1() const
{
    return _N1;
}

unsigned int SpaceDimension::N2() const
{
    return _N2;
}

TimeDimension::TimeDimension(double ht, double t1, double t2, unsigned int M1, unsigned int M2)
    : _ht(ht), _t1(t1), _t2(t2), _M1(M1), _M2(M2)
{}

double TimeDimension::ht() const
{
    return _ht;
}

double TimeDimension::t1() const
{
    return _t1;
}

double TimeDimension::t2() const
{
    return _t2;
}

unsigned int TimeDimension::M1() const
{
    return _M1;
}

unsigned int TimeDimension::M2() const
{
    return _M2;
}

const TimeDimension& GridPDE::timeDimension() const
{
    return mTimeDimension;
}

void GridPDE::setTimeDimension(const TimeDimension &timeDimension)
{
    mTimeDimension = timeDimension;
}

const SpaceDimension& GridPDE::spaceDimensions(SpaceDimension::Dimension dim) const
{
    return mSpaceDimensions.at(dim-1);
}

void GridPDE::addSpaceDimension(const SpaceDimension &spaceDimension)
{
    mSpaceDimensions.push_back(spaceDimension);
}

unsigned int GridPDE::spaceDimSize() const
{
    return mSpaceDimensions.size();
}
