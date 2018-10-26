#include "problem2p_common.h"

SpacePointInfoP::SpacePointInfoP()
{
    init(0);
}

SpacePointInfoP::SpacePointInfoP(unsigned int length)
{
    init(length);
}

SpacePointInfoP::~SpacePointInfoP()
{
    clear();
}

void SpacePointInfoP::init(unsigned int length)
{
    this->length = length;
    vl.resize(length);
    dx.resize(length);
    dy.resize(length);
}
void SpacePointInfoP::clear()
{
    dy.clear();
    dx.clear();
    vl.clear();
    length = 0;
}

auto ExtendedSpacePointP::contains(int nx, int ny) const -> bool
{
    return (minX <= nx && nx <= maxX && minY <= ny && ny <= maxY);
}
