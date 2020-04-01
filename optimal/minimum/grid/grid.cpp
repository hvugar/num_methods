#include "grid.h"

SpacePoint::SpacePoint(double x, double y, double z) : x(x), y(y), z(z) {}

SpaceNodePDE::SpaceNodePDE(int i, double x, int j, double y, int k, double z) : i(i), j(j), k(k)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

SpaceNodePDE::SpaceNodePDE(const SpaceNodePDE &sn) : SpacePoint (sn.x, sn.y, sn.z), i(sn.i), j(sn.j), k(sn.k) {}

/////////////////////////////////////////////////////////////////////////////////////////////////////

TimeMoment::TimeMoment(double t) : t(t) {}

TimeNodePDE::TimeNodePDE(unsigned int i, double t) : TimeMoment (t), i(i) {}

TimeNodePDE::TimeNodePDE(const TimeNodePDE &tn) : TimeMoment (tn.t), i(tn.i) {}

/////////////////////////////////////////////////////////////////////////////////////////////////////

PointODE::PointODE(double x) : x(x) {}

/////////////////////////////////////////////////////////////////////////////////////////////////////

Dimension::Dimension(double step, int min, int max) : _step(step), _min(min), _max(max) {}

Dimension::Dimension(const Dimension &dimension) : _step(dimension._step), _min(dimension._min), _max(dimension._max) {}

Dimension & Dimension::operator =(const Dimension &other)
{
    this->_step = other._step;
    this->_min = other._min;
    this->_max = other._max;
    return *this;
}

double Dimension::step() const { return _step; }

void Dimension::setStep(double step) { _step = step; }

int Dimension::min() const { return _min; }

void Dimension::setMin(int min) { _min = min; }

int Dimension::max() const { return _max; }

void Dimension::setMax(int max) { _max = max; }

unsigned int Dimension::size() const { return static_cast<unsigned int>(_max-_min+1); }

Dimension& Dimension::incMax()
{
    _max++;
    return *this;
}

Dimension& Dimension::decMax()
{
    _max--;
    return *this;
}

Dimension& Dimension::addMinNode()
{
    _min++;
    return *this;
}

Dimension& Dimension::delMinNode()
{
    _min--;
    return *this;
}
