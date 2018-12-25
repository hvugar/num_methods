#include "problem2h_ibvp.h"

double Problem2HNDirichletForward1::initial1(const SpaceNodePDE &) const { return 0.0; }

double Problem2HNDirichletForward1::initial2(const SpaceNodePDE &sn) const
{
    static unsigned int Ns = e_prm->Ns;
    static std::vector<DeltaGrid> thetaGridList;

    if (!_initialCalculation)
    {
        thetaGridList.resize(Ns);
        for (unsigned int s=0; s<Ns; s++)
        {
            thetaGridList[s].initGrid(grid.N, grid.hx, grid.M, grid.hy);
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
    static unsigned int No = e_prm->No;
    static unsigned int Nc = e_prm->Nc;
    static unsigned int Ns = e_prm->Ns;

    for (unsigned int i=0; i<Nc; i++)
    {
        double vi = 0.0;
        for (unsigned int j=0; j<No; j++)
        {
        }
    }

    return 0.0;
}

void Problem2HNDirichletForward1::setU10(const DoubleMatrix &u10)
{
    const size_t No = msrmGridList.size();
    u_xi.clear();
    u_xi.resize(No);

    for (size_t j=0; j<No; j++)
    {
        u_xi[j] = 0.0;
        const  DeltaGrid &mdg = msrmGridList[j];
        for (unsigned int m=mdg.minY(); m<=mdg.maxY(); m++)
        {
            for (unsigned int n=mdg.minX(); n<=mdg.maxX(); n++)
            {
                u_xi[j] += u10[m][n] * mdg.weight(n,m) * (grid.hx*grid.hy);
            }
        }
        //u_xi[j] *= (1.0 + noise * (rand()%2==0 ? +1.0 : -1.0));
    }
}

void Problem2HNDirichletForward1::setEquationParameters(const EquationParameterH &e_prm, const OptimizeParameterH& o_prm)
{
    this->e_prm = &e_prm;
    this->o_prm = &o_prm;

    unsigned int No = e_prm.No;
    for (unsigned int j=0; j<msrmGridList.size(); j++) msrmGridList[j].cleanGrid();
    for (unsigned int j=0; j<No; j++)
    {
        msrmGridList[j].initGrid(grid.N, grid.hx, grid.M, grid.hy);
        msrmGridList[j].setPoint(o_prm.xi[j], 1, 1);
    }
}
