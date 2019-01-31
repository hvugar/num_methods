#include "problem2h_common.h"

InitialPulse2D::InitialPulse2D() : theta(SpacePoint()), q(0.0) {}

InitialPulse2D::InitialPulse2D(const SpacePoint &sp, double q) : theta(sp), q(q) {}

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

auto EquaParameter2H::initPulseParemeters(unsigned int Np) -> void
{
    this->Np = Np;
    pulses.resize(Np);
}

auto EquaParameter2H::initParemeters(unsigned int Nt, unsigned int Nc, unsigned int No) -> void
{
    this->Nt = Nt;
    this->Nc = Nc;
    this->No = No;

    opt.k.resize(Nt);
    opt.z.resize(Nt);
    reg.k.resize(Nt);
    reg.z.resize(Nt);

    for (unsigned int s=0; s<Nt; s++)
    {
        opt.k[s].resize(Nc, No, 0.0);
        opt.z[s].resize(Nc, No, 0.0);
        reg.k[s].resize(Nc, No, 0.0);
        reg.z[s].resize(Nc, No, 0.0);
    }

    opt.ksi.resize(No);
    reg.ksi.resize(No);

    opt.eta.resize(Nc);
    reg.eta.resize(Nc);
}

auto EquaParameter2H::OptimalParameterFromVector(const DoubleVector &x) -> void
{
    unsigned int index = 0;

    for (unsigned int s=0; s<opt.k.size(); s++) opt.k[s].clear();
    for (unsigned int s=0; s<opt.z.size(); s++) opt.z[s].clear();

    opt.k.resize(Nt);
    opt.z.resize(Nt);
    for (unsigned int s=0; s<Nt; s++)
    {
        opt.k[s].clear(); opt.k[s].resize(Nc, No);
        opt.z[s].clear(); opt.z[s].resize(Nc, No);
    }
    opt.ksi.clear(); opt.ksi.resize(No);
    opt.eta.clear(); opt.eta.resize(Nc);

    for (unsigned int s=0; s<Nt; s++)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                opt.k[s][i][j] = x[index++];
            }
        }
    }

    for (unsigned int s=0; s<Nt; s++)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                opt.z[s][i][j] = x[index++];
            }
        }
    }

    for (unsigned int j=0; j<No; j++)
    {
        opt.ksi[j].x = x[index++];
        opt.ksi[j].y = x[index++];
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        opt.eta[i].x = x[index++];
        opt.eta[i].y = x[index++];
    }
}

auto EquaParameter2H::OptimalParameterToVector(DoubleVector &x) const -> void
{
    x.clear();
    x.resize(2*Nc*No*Nt+2*No+2*Nc);

    unsigned int index = 0;
    for (unsigned int s=0; s<Nt; s++)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                x[index++] = opt.k[s][i][j];
            }
        }
    }
    for (unsigned int s=0; s<Nt; s++)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                x[index++] = opt.z[s][i][j];
            }
        }
    }

    for (unsigned int j=0; j<No; j++)
    {
        x[index++] = opt.ksi[j].x;
        x[index++] = opt.ksi[j].y;
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        x[index++] = opt.eta[i].x;
        x[index++] = opt.eta[i].y;
    }
}
