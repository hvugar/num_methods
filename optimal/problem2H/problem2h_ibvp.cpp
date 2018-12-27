#include "problem2h_ibvp.h"

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

void Problem2HNDirichletForward1::setEquationParameters(const EquationParameterH &e_prm, const OptimizeParameterH &o_prm, unsigned int N, double hx, unsigned int M, double hy)
{
    IPrinter::printSeperatorLine();
    this->e_prm = &e_prm;
    this->o_prm = &o_prm;

    const unsigned int No = e_prm.No;
    const unsigned int Nc = e_prm.Nc;
    const unsigned int Ns = e_prm.Ns;

    for (unsigned int s=0; s<thetaGridList.size(); s++) thetaGridList[s].cleanGrid();
    thetaGridList.clear();
    for (unsigned int j=0; j<msrmGridList.size(); j++) msrmGridList[j].cleanGrid();
    msrmGridList.clear();
    for (unsigned int i=0; i<cntrGridList.size(); i++) cntrGridList[i].cleanGrid();
    cntrGridList.clear();

    msrmGridList.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        msrmGridList[j].initGrid(N, hx, M, hy);
        msrmGridList[j].setPoint(o_prm.xi[j], 1, 1);
    }

    cntrGridList.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        cntrGridList[i].initGrid(N, hx, M, hy);
        cntrGridList[i].setPoint(o_prm.eta[i], 1, 1);
    }

    thetaGridList.resize(Ns);
    for (unsigned int s=0; s<Ns; s++)
    {
        thetaGridList[s].initGrid(N, hx, M, hy);
        thetaGridList[s].setPoint(e_prm.theta[s], 8, 8);
    }

    ixv.resize(M+1, N+1, 0.0);
    fxv.resize(M+1, N+1, 0.0);

    IPrinter::printSeperatorLine();
}

void Problem2HNDirichletForward1::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const static Dimension &dimY = spaceDimension(Dimension::DimensionY);
    const static unsigned int N = static_cast<unsigned int>(dimX.size());
    const static unsigned int M = static_cast<unsigned int>(dimY.size());
    const static double hx = dimX.step();
    const static double hy = dimY.step();

    const static unsigned int No = e_prm->No;
    const static unsigned int Nc = e_prm->Nc;
    const static unsigned int Ns = e_prm->Ns;

    if (ln==500)
    {
        IPrinter::printSeperatorLine("*layerInfo: ln == 2");
        IPrinter::printMatrix(u);
    }

    if (ln == 0)
    {
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                ixv[m][n] = 0.0;
                for (unsigned int s=0; s<Ns; s++) ixv[m][n] += e_prm->q[s]*thetaGridList[s].weight(n,m);
            }
        }
    }

    if (ln > 0)
    {
        double *_u = new double[No];
        for (unsigned int j=0; j<No; j++)
        {
            _u[j] = 0.0;
            const  DeltaGrid &mdg = msrmGridList[j];
            for (unsigned int m=mdg.minY(); m<=mdg.maxY(); m++)
            {
                for (unsigned int n=mdg.minX(); n<=mdg.maxX(); n++)
                {
                    _u[j] += u[m][n] * mdg.weight(n,m) * (hx*hy);
                }
            }
            //_u[j] *= (1.0 + noise * (rand()%2==0 ? +1.0 : -1.0));
        }

        double *_v = new double[Nc];
        for (unsigned int i=0; i<Nc; i++)
        {
            _v[i] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                _v[i] += o_prm->k[i][j] * (_u[j] - o_prm->z[i][j]);
            }
        }
        delete [] _u;

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
    }
}
