#include "grid.h"
#include <cstdio>

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

auto Dimension::step() const -> double { return _step; }

auto Dimension::setStep(double step) -> void { _step = step; }

auto Dimension::min() const -> int { return _min; }

auto Dimension::setMin(int min) -> void { _min = min; }

auto Dimension::max() const -> int { return _max; }

auto Dimension::setMax(int max) -> void { _max = max; }

auto Dimension::size() const -> int { return _max-_min; }
