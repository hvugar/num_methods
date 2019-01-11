#include "problem2h_common.h"

auto ExtendedSpacePointH::contains(int nx, int ny) const -> bool
{
    return (minX <= nx && nx <= maxX && minY <= ny && ny <= maxY);
}

SpacePointInfoH::SpacePointInfoH()
{
    init(0);
}

SpacePointInfoH::SpacePointInfoH(unsigned int length)
{
    init(length);
}

SpacePointInfoH::~SpacePointInfoH()
{
    clear();
}

void SpacePointInfoH::init(unsigned int length)
{
    this->length = length;
    vl.resize(length);
    dx.resize(length);
    dy.resize(length);
}
void SpacePointInfoH::clear()
{
    dy.clear();
    dx.clear();
    vl.clear();
    length = 0;
}
