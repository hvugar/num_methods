#include "problem2h_ibvp.h"

Problem2HNDirichletForward1::Problem2HNDirichletForward1() : CdIHyperbolicIBVP ()
{}

double Problem2HNDirichletForward1::initial1(const SpaceNodePDE &) const { return 0.0; }

double Problem2HNDirichletForward1::initial2(const SpaceNodePDE &sn) const
{
    return ixv[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

double Problem2HNDirichletForward1::boundary(const SpaceNodePDE&, const TimeNodePDE&) const { return 0.0; }

double Problem2HNDirichletForward1::f(const SpaceNodePDE &sn, const TimeNodePDE &) const
{
    return fxv[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HNDirichletForward1::setParameters(const EquationParameterH &equationParameter, const OptimizeParameterH &optimizeParameter, unsigned int LD)
{
    clear();
    mEquationParameter = equationParameter;
    mOptimizeParameter = optimizeParameter;

    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionY);
    const Dimension &time = timeDimension();
    const unsigned int N = static_cast<unsigned int>(dimX.size());
    const unsigned int M = static_cast<unsigned int>(dimY.size());
    const unsigned int L = static_cast<unsigned int>(time.size());
    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = dimY.step();

    const unsigned int No = equationParameter.No;
    const unsigned int Nc = equationParameter.Nc;
    const unsigned int Ns = equationParameter.Ns;

    for (unsigned int s=0; s<tetaGridList.size(); s++) tetaGridList[s].cleanGrid(); tetaGridList.clear();
    for (unsigned int j=0; j<msrmGridList.size(); j++) msrmGridList[j].cleanGrid(); msrmGridList.clear();
    for (unsigned int i=0; i<cntrGridList.size(); i++) cntrGridList[i].cleanGrid(); cntrGridList.clear();

    msrmGridList.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        const SpacePoint &sp = optimizeParameter.xi[j];
        msrmGridList[j].initGrid(N, hx, M, hy);
        msrmGridList[j].distributeGauss(sp);
    }

    cntrGridList.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        const SpacePoint &sp = optimizeParameter.eta[i];
        cntrGridList[i].initGrid(N, hx, M, hy);
        cntrGridList[i].distributeGauss(sp);
    }

    tetaGridList.resize(Ns);
    for (unsigned int s=0; s<Ns; s++)
    {
        const SpacePoint &sp = equationParameter.pulses[s].theta;
        tetaGridList[s].initGrid(N, hx, M, hy);
        tetaGridList[s].distributeGauss(sp, 8, 8);
    }

    ixv.resize(M+1, N+1, 0.0);
    fxv.resize(M+1, N+1, 0.0);

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            ixv[m][n] = 0.0;
            for (unsigned int s=0; s<Ns; s++) ixv[m][n] += equationParameter.pulses[s].q*tetaGridList[s].weight(n,m);
        }
    }

    vu.resize(LD+1);
    for (unsigned int i=0; i<=LD; i++)
    {
        vu[i].resize(M+1, N+1);
    }

    u_info.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        SpacePointInfoH &inf = u_info[j];
        const SpacePoint &sp = optimizeParameter.xi[j];
        u_info[j].init(L+1);
        inf.x = sp.x;
        inf.y = sp.y;
        inf.init(L+1);
    }
}

void Problem2HNDirichletForward1::clear()
{
    for (unsigned int j=0; j<msrmGridList.size(); j++) msrmGridList[j].cleanGrid(); msrmGridList.clear();
    for (unsigned int i=0; i<cntrGridList.size(); i++) cntrGridList[i].cleanGrid(); cntrGridList.clear();
    for (unsigned int s=0; s<tetaGridList.size(); s++) tetaGridList[s].cleanGrid(); tetaGridList.clear();
    ixv.clear();
    fxv.clear();
    for (unsigned int ln=0; ln<vu.size(); ln++) vu[ln].clear(); vu.clear();
    for (unsigned int j=0; j<u_info.size(); j++) u_info[j].clear(); u_info.clear();
}

void Problem2HNDirichletForward1::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const static Dimension &dimY = spaceDimension(Dimension::DimensionY);
    const static Dimension &time = timeDimension();
    const static unsigned int N = static_cast<unsigned int>(dimX.size());
    const static unsigned int M = static_cast<unsigned int>(dimY.size());
    const static unsigned int L = static_cast<unsigned int>(time.size());
    const static double hx = dimX.step();
    const static double hy = dimY.step();
    //const static double ht = dimY.step();

    const static unsigned int No = mEquationParameter.No;
    const static unsigned int Nc = mEquationParameter.Nc;
    const static unsigned int Ns = mEquationParameter.Ns;

    Problem2HNDirichletForward1 *pf = const_cast<Problem2HNDirichletForward1*>(this);

    if (ln==500)
    {
        IPrinter::printSeperatorLine("*layerInfo: ln == 2");
        IPrinter::printMatrix(u);
    }

    if (ln >= L-LD)
    {
        pf->vu[ln-(L-LD)] = u;
    }

    if (ln > 0)
    {
        double *_u = new double[No];
        double *_v = new double[Nc];

        for (unsigned int j=0; j<No; j++)
        {
            _u[j] = 0.0;
            const  DeltaGrid2D &mdg = msrmGridList[j];
            for (unsigned int m=mdg.minY(); m<=mdg.maxY(); m++)
            {
                for (unsigned int n=mdg.minX(); n<=mdg.maxX(); n++)
                {
                    _u[j] += u[m][n] * mdg.weight(n,m) * (hx*hy);
                }
            }
            //_u[j] *= (1.0 + noise * (rand()%2==0 ? +1.0 : -1.0));
        }

        for (unsigned int i=0; i<Nc; i++)
        {
            _v[i] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                _v[i] += mOptimizeParameter.k[i][j] * (_u[j] - mOptimizeParameter.z[i][j]);
            }
        }

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                fxv[m][n] = 0.0;
                for (unsigned int i=0; i<Nc; i++)
                {
                    fxv[m][n] += _v[i] * cntrGridList[i].weight(n,m);
                }
            }
        }

        delete [] _u;
        delete [] _v;

        for (unsigned int j=0; j<No; j++)
        {
            SpacePointInfoH &ui = pf->u_info[j];
            const DeltaGrid2D &mdg = msrmGridList[j];
            for (unsigned int m=mdg.minY(); m<=mdg.maxY(); m++)
            {
                for (unsigned int n=mdg.minX(); n<=mdg.maxX(); n++)
                {
                    ui.vl[ln] += u[m][n] * mdg.weight(n,m) * (hx*hy);
                }
            }

            double px = mdg.p().x;
            double py = mdg.p().y;
            unsigned int rx = mdg.rx();
            unsigned int ry = mdg.rx();

            ui.dx[ln] = (u[ry][rx+1] - u[ry][rx-1])/(2.0*hx);
            ui.dy[ln] = (u[ry+1][rx] - u[ry-1][rx])/(2.0*hy);

            ui.dx[ln] += ((px-rx*hx)/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
            ui.dy[ln] += ((py-ry*hy)/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);

            //ui.dxx[ln] = (1.0/(hx*hx))*(u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1]);
            //ui.dyy[ln] = (1.0/(hy*hy))*(u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx]);
        }

    }
}

//**********************************************************************************************************//

Problem2HNDirichletBackward1::Problem2HNDirichletBackward1(const Problem2HNDirichletForward1 &fw) : _fw(fw)
{
    setTimeDimension(fw.timeDimension());
    addSpaceDimension(fw.spaceDimension(Dimension::DimensionX));
    addSpaceDimension(fw.spaceDimension(Dimension::DimensionY));
}

double Problem2HNDirichletBackward1::initial1(const SpaceNodePDE &) const { return 0.0; }

double Problem2HNDirichletBackward1::initial2(const SpaceNodePDE &) const { return 0.0; }

double Problem2HNDirichletBackward1::boundary(const SpaceNodePDE&, const TimeNodePDE&) const { return 0.0; }

double Problem2HNDirichletBackward1::f(const SpaceNodePDE &sn, const TimeNodePDE &) const
{
    return fxv[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HNDirichletBackward1::layerInfo(const DoubleMatrix &p, unsigned int ln) const
{
//    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
//    const Dimension &dimY = spaceDimension(Dimension::DimensionY);
//    const Dimension &time = timeDimension();
//    const unsigned int N = static_cast<unsigned int>(dimX.size());
//    const unsigned int M = static_cast<unsigned int>(dimY.size());
//    const unsigned int L = static_cast<unsigned int>(time.size());
//    const double hx = dimX.step();
//    const double hy = dimY.step();
//    const double ht = dimY.step();

//    const EquationParameterH &mep = _fw.mEquationParameter;
//    const OptimizeParameterH &mop = _fw.mOptimizeParameter;

//    const unsigned int No = mep.No;
//    const unsigned int Nc = mep.Nc;

//    if (ln==500)
//    {
//        IPrinter::printSeperatorLine("*layerInfo: ln == 2");
//        IPrinter::printMatrix(p);
//    }

//    if (ln == L)
//    {
//        for (unsigned int m=0; m<=M; m++)
//        {
//            for (unsigned int n=0; n<=N; n++)
//            {
//                ixv[m][n] = 0.0;
//            }
//        }
//    }

//    if (ln < L)
//    {
//        double *_p = new double[Nc];
//        for (unsigned int i=0; i<Nc; i++)
//        {
//            _p[i] = 0.0;
//            const DeltaGrid &mdg = _fw.cntrGridList[i];
//            for (unsigned int m=mdg.minY(); m<=mdg.maxY(); m++)
//            {
//                for (unsigned int n=mdg.minX(); n<=mdg.maxX(); n++)
//                {
//                    _p[i] += p[m][n] * mdg.weight(n,m) * (hx*hy);
//                }
//            }
//            //_u[j] *= (1.0 + noise * (rand()%2==0 ? +1.0 : -1.0));
//        }

//        double *_w = new double[No];
//        for (unsigned int j=0; j<No; j++)
//        {
//            _w[j] = 0.0;
//            for (unsigned int i=0; i<Nc; i++)
//            {
//                _w[j] += mop.k[i][j] * (_p[i]);
//            }
//        }
//        delete [] _p;

//        for (unsigned int m=0; m<=M; m++)
//        {
//            for (unsigned int n=0; n<=N; n++)
//            {
//                fxv[m][n] = 0.0;
//                for (unsigned int j=0; j<No; j++)
//                {
//                    fxv[m][n] += _w[j] * _fw.msrmGridList[j].weight(n,m);
//                }

//                if (ln < L-LD)
//                {
//                } else {
//                    fxv[m][n] -= 2.0*vu.at(ln)[m][n];
//                }
//            }
//        }
//    }
}

//**********************************************************************************************************//

//**********************************************************************************************************//
