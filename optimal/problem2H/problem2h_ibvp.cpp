#include "problem2h_ibvp.h"

double Problem2HNDirichletForward1::initial1(const SpaceNodePDE &) const { return 0.0; }

double Problem2HNDirichletForward1::initial2(const SpaceNodePDE &sn) const
{
    static unsigned int Ns = e_prm->Ns;
    static std::vector<DeltaGrid> thetaGridList;

    const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const static Dimension &dimY = spaceDimension(Dimension::DimensionY);
    const static unsigned int N = static_cast<unsigned int>(dimX.size());
    const static unsigned int M = static_cast<unsigned int>(dimY.size());
    const static double hx = dimX.step();
    const static double hy = dimY.step();

    if (!_initialCalculation)
    {
        thetaGridList.resize(Ns);
        for (unsigned int s=0; s<Ns; s++)
        {
            thetaGridList[s].initGrid(N, hx, M, hy);
            thetaGridList[s].setPoint(e_prm->theta[s], 5, 5);
        }
        const_cast<Problem2HNDirichletForward1*>(this)->_initialCalculation = true;
    }

    double sum = 0.0;
    for (unsigned int s=0; s<Ns; s++) sum += e_prm->q[s]*thetaGridList[s].weight(sn);
    return sum;
}

double Problem2HNDirichletForward1::boundary(const SpaceNodePDE&, const TimeNodePDE&) const { return 0.0; }

double Problem2HNDirichletForward1::f(const SpaceNodePDE &sn, const TimeNodePDE &) const
{
    return fx_value[sn.j][sn.i];
}

void Problem2HNDirichletForward1::setEquationParameters(const EquationParameterH &e_prm,
                                                        const OptimizeParameterH &o_prm,
                                                        unsigned int N, double hx,
                                                        unsigned int M, double hy)
{
    this->e_prm = &e_prm;
    this->o_prm = &o_prm;

    const unsigned int No = e_prm.No;
    const unsigned int Nc = e_prm.Nc;

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

    msrmGridList.resize(No);
    for (unsigned int i=0; i<Nc; i++)
    {
        cntrGridList[i].initGrid(N, hx, M, hy);
        cntrGridList[i].setPoint(o_prm.eta[i], 1, 1);
    }
}

void Problem2HNDirichletForward1::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const static Dimension &dimY = spaceDimension(Dimension::DimensionY);
    const static unsigned int N = static_cast<unsigned int>(dimX.size());
    const static unsigned int M = static_cast<unsigned int>(dimY.size());
    const static double hx = dimX.step();
    const static double hy = dimY.step();

    if (ln > 0)
    {
        const unsigned int No = e_prm->No;
        const unsigned int Nc = e_prm->Nc;

        double *_u = new double[No];
        for (size_t j=0; j<No; j++)
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
            //u_xi[j] *= (1.0 + noise * (rand()%2==0 ? +1.0 : -1.0));
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

        for (unsigned int m=1; m<=M-1; m++)
        {
            for (unsigned int n=1; n<=N-1; n++)
            {
                double _fx = 0.0;
                for (unsigned int i=0; i<Nc; i++)
                {
                    _fx += _v[i] * cntrGridList[i].weight(n,m);
                }
                fx_value[m][n] = _fx;
            }
        }
    }
}
