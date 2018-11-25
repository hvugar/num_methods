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

//EquationParameterHE::EquationParameterHE()
//{
//    EquationParameterHE(0,0,0,0,0);
//}

//EquationParameterHE::EquationParameterHE(unsigned int Nc, unsigned int No, unsigned int Nd, unsigned int length, unsigned int gw)
//{
//    this->Nc = Nc;
//    this->No = No;
//    this->Nd = Nd;
//    k.resize(Nc, No, 0.0);
//    z.resize(Nc, No, 0.0);
//    xi.resize(No);
//    eta.resize(Nc);
//    q.resize(Nd);
//    theta.resize(Nd);

//    xi_ext.resize(No);  for (unsigned int j=0; j<No; j++) xi_ext[j].length = length;
//    eta_ext.resize(Nc); for (unsigned int i=0; i<No; i++) eta_ext[i].length = length;
//    theta.resize(Nd);
//}


