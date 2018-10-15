#include "problem2h_common.h"

SpacePointInfo::SpacePointInfo()
{
    init(0);
}

SpacePointInfo::SpacePointInfo(unsigned int length)
{
    init(length);
}

SpacePointInfo::~SpacePointInfo()
{
    clear();
}

void SpacePointInfo::init(unsigned int length)
{
    this->length = length;
    vl.resize(length);
    dx.resize(length);
    dy.resize(length);
}
void SpacePointInfo::clear()
{
    dy.clear();
    dx.clear();
    vl.clear();
    length = 0;
}

auto ExtendedSpacePoint::contains(int nx, int ny) const -> bool
{
    return (minX <= nx && nx <= maxX && minY <= ny && ny <= maxY);
}
