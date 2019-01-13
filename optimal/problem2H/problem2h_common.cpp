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

EquationParameterH1::EquationParameterH1(unsigned int Nc, unsigned int No, unsigned int Ns, double a, double alpha)
{
    this->Nc = Nc;
    this->No = No;
    this->Ns = Ns;
    this->a = a;
    this->alpha = alpha;

    op.k.resize(Nc, No); op.z.resize(Nc, No); op.xi.resize(No); op.eta.resize(Nc);
    rp.k.resize(Nc, No); rp.z.resize(Nc, No); rp.xi.resize(No); rp.eta.resize(Nc);
    pulseVector.resize(Ns);
}
