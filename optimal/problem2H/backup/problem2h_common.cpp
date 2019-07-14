#include "problem2h_common.h"

//InitialPulse2D::InitialPulse2D() : theta(SpacePoint()), q(0.0) {}

//InitialPulse2D::InitialPulse2D(const SpacePoint &sp, double q) : theta(sp), q(q) {}

//InitialPulse2D::InitialPulse2D(const InitialPulse2D &pulse) : theta(pulse.theta), q(pulse.q) {}

//auto ExtendedSpacePointH::contains(int nx, int ny) const -> bool
//{
//    return (minX <= nx && nx <= maxX && minY <= ny && ny <= maxY);
//}

//SpacePointInfoH::SpacePointInfoH()
//{
//    init(0);
//}

//SpacePointInfoH::SpacePointInfoH(unsigned int length)
//{
//    init(length);
//}

//SpacePointInfoH::~SpacePointInfoH()
//{
//    clear();
//}

//void SpacePointInfoH::init(unsigned int length)
//{
//    this->length = length;
//    vl.resize(length);
//    dx.resize(length);
//    dy.resize(length);
//    vx.resize(length);
//    vy.resize(length);
//}

//void SpacePointInfoH::clear()
//{
//    vy.clear();
//    vx.clear();
//    dy.clear();
//    dx.clear();
//    vl.clear();
//    length = 0;
//}

//auto EquaParameter2H::initPulseParemeters(unsigned int Np) -> void
//{
//    this->Np = Np;
//    pulses.resize(Np);
//}

//auto EquaParameter2H::initParemeters(unsigned int Nt, unsigned int Nc, unsigned int No) -> void
//{
//    this->Nt = Nt;
//    this->Nc = Nc;
//    this->No = No;

//#if defined (DISCRETE_DELTA_TIME_1)
//    opt.k.resize(Nt);
//    opt.z.resize(Nt);
//    reg.k.resize(Nt);
//    reg.z.resize(Nt);

//    for (unsigned int s=0; s<Nt; s++)
//    {
//        opt.k[s].resize(Nc, No, 0.0);
//        opt.z[s].resize(Nc, No, 0.0);
//        reg.k[s].resize(Nc, No, 0.0);
//        reg.z[s].resize(Nc, No, 0.0);
//    }
//#endif

//#if defined (DISCRETE_DELTA_TIME_2)
//    opt.k.resize(Nc, No, 0.0);
//    opt.z.resize(Nc, No, 0.0);
//    reg.k.resize(Nc, No, 0.0);
//    reg.z.resize(Nc, No, 0.0);
//#endif

//    opt.ksi.resize(No);
//    reg.ksi.resize(No);

//    opt.eta.resize(Nc);
//    reg.eta.resize(Nc);
//}

//auto EquaParameter2H::OptimalParameterFromVector(const DoubleVector &x) -> void
//{
//    unsigned int index = 0;

//#if defined (DISCRETE_DELTA_TIME_1)
//    for (unsigned int s=0; s<opt.k.size(); s++) opt.k[s].clear();
//    for (unsigned int s=0; s<opt.z.size(); s++) opt.z[s].clear();

//    opt.k.resize(Nt);
//    opt.z.resize(Nt);
//    for (unsigned int s=0; s<Nt; s++)
//    {
//        opt.k[s].clear(); opt.k[s].resize(Nc, No);
//        opt.z[s].clear(); opt.z[s].resize(Nc, No);
//    }
//    opt.ksi.clear(); opt.ksi.resize(No);
//    opt.eta.clear(); opt.eta.resize(Nc);

//    for (unsigned int s=0; s<Nt; s++)
//    {
//        for (unsigned int i=0; i<Nc; i++)
//        {
//            for (unsigned int j=0; j<No; j++)
//            {
//                opt.k[s][i][j] = x[index++];
//            }
//        }
//    }

//    for (unsigned int s=0; s<Nt; s++)
//    {
//        for (unsigned int i=0; i<Nc; i++)
//        {
//            for (unsigned int j=0; j<No; j++)
//            {
//                opt.z[s][i][j] = x[index++];
//            }
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        opt.ksi[j].x = x[index++];
//        opt.ksi[j].y = x[index++];
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        opt.eta[i].x = x[index++];
//        opt.eta[i].y = x[index++];
//    }
//#endif

//#if defined (DISCRETE_DELTA_TIME_2)
//    opt.k.clear(); opt.k.resize(Nc, No);
//    opt.z.clear(); opt.z.resize(Nc, No);
//    opt.ksi.clear(); opt.ksi.resize(No);
//    opt.eta.clear(); opt.eta.resize(Nc);

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            opt.k[i][j] = x[index++];
//        }
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            opt.z[i][j] = x[index++];
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        opt.ksi[j].x = x[index++];
//        opt.ksi[j].y = x[index++];
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        opt.eta[i].x = x[index++];
//        opt.eta[i].y = x[index++];
//    }
//#endif
//}

//auto EquaParameter2H::OptimalParameterToVector(DoubleVector &x) const -> void
//{
//#if defined (DISCRETE_DELTA_TIME_1)
//    x.clear();
//    x.resize(2*Nc*No*Nt+2*No+2*Nc);

//    unsigned int index = 0;
//    for (unsigned int s=0; s<Nt; s++)
//    {
//        for (unsigned int i=0; i<Nc; i++)
//        {
//            for (unsigned int j=0; j<No; j++)
//            {
//                x[index++] = opt.k[s][i][j];
//            }
//        }
//    }
//    for (unsigned int s=0; s<Nt; s++)
//    {
//        for (unsigned int i=0; i<Nc; i++)
//        {
//            for (unsigned int j=0; j<No; j++)
//            {
//                x[index++] = opt.z[s][i][j];
//            }
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        x[index++] = opt.ksi[j].x;
//        x[index++] = opt.ksi[j].y;
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        x[index++] = opt.eta[i].x;
//        x[index++] = opt.eta[i].y;
//    }
//#endif

//#if defined (DISCRETE_DELTA_TIME_2)
//    x.clear();
//    x.resize(2*Nc*No+2*No+2*Nc);

//    unsigned int index = 0;
//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            x[index++] = opt.k[i][j];
//        }
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            x[index++] = opt.z[i][j];
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        x[index++] = opt.ksi[j].x;
//        x[index++] = opt.ksi[j].y;
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        x[index++] = opt.eta[i].x;
//        x[index++] = opt.eta[i].y;
//    }
//#endif
//}

//auto EquaParameter2H::RegularParameterToVector(DoubleVector &x) const -> void
//{
//#if defined (DISCRETE_DELTA_TIME_1)
//    x.clear();
//    x.resize(2*Nc*No*Nt+2*No+2*Nc);

//    unsigned int index = 0;
//    for (unsigned int s=0; s<Nt; s++)
//    {
//        for (unsigned int i=0; i<Nc; i++)
//        {
//            for (unsigned int j=0; j<No; j++)
//            {
//                x[index++] = reg.k[s][i][j];
//            }
//        }
//    }
//    for (unsigned int s=0; s<Nt; s++)
//    {
//        for (unsigned int i=0; i<Nc; i++)
//        {
//            for (unsigned int j=0; j<No; j++)
//            {
//                x[index++] = reg.z[s][i][j];
//            }
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        x[index++] = reg.ksi[j].x;
//        x[index++] = reg.ksi[j].y;
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        x[index++] = reg.eta[i].x;
//        x[index++] = reg.eta[i].y;
//    }
//#endif

//#if defined (DISCRETE_DELTA_TIME_2)
//    x.clear();
//    x.resize(2*Nc*No+2*No+2*Nc);

//    unsigned int index = 0;
//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            x[index++] = reg.k[i][j];
//        }
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            x[index++] = reg.z[i][j];
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        x[index++] = reg.ksi[j].x;
//        x[index++] = reg.ksi[j].y;
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        x[index++] = reg.eta[i].x;
//        x[index++] = reg.eta[i].y;
//    }
//#endif
//}

//auto EquaParameter2H::RegularParameterFromVector(const DoubleVector &rx) -> void
//{
//#if defined (DISCRETE_DELTA_TIME_1)
//    unsigned int index = 0;

//    for (unsigned int s=0; s<reg.k.size(); s++) reg.k[s].clear();
//    for (unsigned int s=0; s<reg.z.size(); s++) reg.z[s].clear();

//    reg.k.resize(Nt);
//    reg.z.resize(Nt);
//    for (unsigned int s=0; s<Nt; s++)
//    {
//        reg.k[s].clear(); reg.k[s].resize(Nc, No);
//        reg.z[s].clear(); reg.z[s].resize(Nc, No);
//    }
//    reg.ksi.clear(); reg.ksi.resize(No);
//    reg.eta.clear(); reg.eta.resize(Nc);

//    for (unsigned int s=0; s<Nt; s++)
//    {
//        for (unsigned int i=0; i<Nc; i++)
//        {
//            for (unsigned int j=0; j<No; j++)
//            {
//                reg.k[s][i][j] = rx[index++];
//            }
//        }
//    }

//    for (unsigned int s=0; s<Nt; s++)
//    {
//        for (unsigned int i=0; i<Nc; i++)
//        {
//            for (unsigned int j=0; j<No; j++)
//            {
//                reg.z[s][i][j] = rx[index++];
//            }
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        reg.ksi[j].x = rx[index++];
//        reg.ksi[j].y = rx[index++];
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        reg.eta[i].x = rx[index++];
//        reg.eta[i].y = rx[index++];
//    }
//#endif

//#if defined (DISCRETE_DELTA_TIME_2)
//    unsigned int index = 0;

//    reg.k.clear(); reg.k.resize(Nc, No);
//    reg.z.clear(); reg.z.resize(Nc, No);
//    reg.ksi.clear(); reg.ksi.resize(No);
//    reg.eta.clear(); reg.eta.resize(Nc);

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            reg.k[i][j] = rx[index++];
//        }
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            reg.z[i][j] = rx[index++];
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        reg.ksi[j].x = rx[index++];
//        reg.ksi[j].y = rx[index++];
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        reg.eta[i].x = rx[index++];
//        reg.eta[i].y = rx[index++];
//    }
//#endif
//}

//auto EquaParameter2H::printOptimalParemeters() const -> void
//{
//    //    printf("k : "); IPrinter::print(x.mid(00, 19), x.mid(00, 19).length(), 9, 6);
//    //    printf("k : "); IPrinter::print(x.mid(20, 39), x.mid(20, 39).length(), 9, 6);
//    //    printf("z : "); IPrinter::print(x.mid(40, 59), x.mid(40, 59).length(), 9, 6);
//    //    printf("z : "); IPrinter::print(x.mid(60, 79), x.mid(60, 79).length(), 9, 6);
//    //    printf("xy: "); IPrinter::print(x.mid(80, 87), x.mid(80, 87).length(), 9, 6);
//}

//auto EquaParameter2H::printRegularParemeters() const -> void
//{

//}
