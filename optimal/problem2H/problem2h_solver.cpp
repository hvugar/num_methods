#include "problem2h_solver.h"

void Problem2HSolver::Main(int argc, char **argv)
{
    Problem2HSolver ps;
    Dimension dimensionX(0.01, 0, 100);
    Dimension dimensionY(0.01, 0, 100);
    Dimension dimensionT(0.01, 0, 1000);

    ps.setTimeDimension(dimensionT);
    ps.addSpaceDimension(dimensionX);
    ps.addSpaceDimension(dimensionY);
    ps.initMatrixes(dimensionX, dimensionY, dimensionT, 2, 2, 2);

    ps.setWaveSpeed(1.0);
    ps.setWaveDissipation(0.01);

    unsigned int Nq = 2;
    SpacePoint* zta = new SpacePoint[2];
    zta[0] = SpacePoint(0.25, 0.25);
    zta[1] = SpacePoint(0.75, 0.75);
    double* q = new double[Nq];
    q[0] = 0.05;
    q[1] = 0.05;
    ps.setInitialConditionMatrix(zta, q, Nq);

    std::vector<SpacePoint> eta;
    eta.push_back(SpacePoint(0.65, 0.34));
    eta.push_back(SpacePoint(0.25, 0.75));

    std::vector<SpacePoint> ksi;
    ksi.push_back(SpacePoint(0.22, 0.54));
    ksi.push_back(SpacePoint(0.82, 0.27));

    ps.fw().initControlMeasurementDeltaGrid(eta, ksi);
    ps.fw().implicit_calculate_D2V1();

    delete [] q;
    delete [] zta;
}

//Problem2HForward& Problem2HSolver::forward()
//{
//    return *this;
//}

//Problem2HBackward& Problem2HSolver::backward()
//{
//    return *this;
//}

void Problem2HCommon::initMatrixes(const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension UNUSED_PARAM,
                                   unsigned int Nq, unsigned int Nc, unsigned No)
{
    uint32_t N = static_cast<uint32_t>(dimensionX.size());
    uint32_t M = static_cast<uint32_t>(dimensionY.size());
    uint32_t L = static_cast<uint32_t>(timeDimension.size());
    double hx = dimensionX.step();
    double hy = dimensionY.step();
    uint32_t length = 2*L + 1;

    f_initialMatrix.clear();
    f_initialMatrix.resize(M+1, N+1, 0.0);

    f_layerMatrix.clear();
    f_layerMatrix.resize(M+1, N+1, 0.0);

    b_initialMatrix.clear();
    b_initialMatrix.resize(M+1, N+1, 0.0);

    b_layerMatrix.clear();
    b_layerMatrix.resize(M+1, N+1, 0.0);

//    this->Nq = Nq;
//    this->zta.resize(Nq);
//    for (unsigned int i=0; i<Nq; i++)
//    {
//        this->zta[i].cleanGrid();
//        this->zta[i].initGrid(N, hx, M, hy);
//    }

    this->Nc = Nc;
    this->eta.resize(Nc);
    this->p_info.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        this->eta[i].cleanGrid();
        this->eta[i].initGrid(N, hx, M, hy);
        this->p_info[i].vl.resize(length);
        this->p_info[i].dx.resize(length);
        this->p_info[i].dy.resize(length);
    }

    this->No = No;
    this->ksi.resize(No);
    this->u_info.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        this->ksi[j].cleanGrid();
        this->ksi[j].initGrid(N, hx, M, hy);
        this->u_info[j].vl.resize(length);
        this->u_info[j].dx.resize(length);
        this->u_info[j].dy.resize(length);
    }

    k.resize(Nc, No, -1.0);
    z.resize(Nc, No, +0.0);
}

void Problem2HCommon::clearMatrixes()
{
    f_initialMatrix.clear();
    f_layerMatrix.clear();
    b_initialMatrix.clear();
    b_layerMatrix.clear();

//    for (unsigned int i=0; i<Nq; i++) zta[i].cleanGrid(); Nq = 0;

    for (unsigned int i=0; i<Nc; i++)
    {
        eta[i].cleanGrid();
        p_info[i].vl.clear();
        p_info[i].dx.clear();
        p_info[i].dy.clear();
        Nc = 0;
    }
    eta.clear();
    p_info.clear();

    for (unsigned int j=0; j<No; j++)
    {
        ksi[j].cleanGrid();
        u_info[j].vl.clear();
        u_info[j].dx.clear();
        u_info[j].dy.clear();
        No = 0;
    }
    ksi.clear();
    u_info.clear();
}

///////////////////////////////////

double Problem2HForward::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    if (condition == InitialCondition::FirstDerivative)
        return f_initialMatrix[static_cast<uint32_t>(sn.j)][static_cast<uint32_t>(sn.i)];
    else
        return 0.0;
}

double Problem2HForward::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    return 0.0;
}

double Problem2HForward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return f_layerMatrix[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HForward::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    Problem2HForward* pf = const_cast<Problem2HForward*>(this);
    for (unsigned int j=0; j<No; j++)
    {
        double u_vl, u_dx, u_dy;
        u_vl = ksi[j].lumpPointGauss(u, u_dx, u_dy);
        pf->u_info[j].vl[tn.i] = u_vl;
        pf->u_info[j].dx[tn.i] = u_dx;
        pf->u_info[j].dy[tn.i] = u_dy;
    }
    const_cast<Problem2HForward*>(this)->prepareLayerMatrix(u, tn);

    saveToTextF(u, tn);
}

void Problem2HForward::saveToTextF(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn) const
{
    static double MIN = +100000.0;
    static double MAX = -100000.0;

    std::string txt_number = std::to_string(tn.i);
    std::string filename = std::string("data/problem2H/f/txt/f_") +
            std::string(4 - txt_number.length(), '0') + txt_number + std::string(".txt");
    IPrinter::print(u, filename.c_str());
    if (MIN > u.min()) MIN = u.min();
    if (MAX < u.max()) MAX = u.max();
    printf("Forward: %4d %0.3f %10.8f %10.8f %10.8f %10.8f %4d %4d\n", tn.i, tn.t, u.min(), u.max(), MIN, MAX, 0, 0);
}

void Problem2HForward::setInitialConditionMatrix(const SpacePoint* zta, const double* q, unsigned int Nq)
{
    const Dimension &dimensionX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimensionY = spaceDimension(Dimension::DimensionY);
    const unsigned int N = static_cast<unsigned int>(dimensionX.size());
    const unsigned int M = static_cast<unsigned int>(dimensionY.size());
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();

    for (unsigned int i=0; i<Nq; i++)
    {
        DeltaGrid2D deltaGrid;
        deltaGrid.initGrid(N, hx, M, hy);
        deltaGrid.distributeGauss(zta[i], 5, 5);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                f_initialMatrix[m][n] += q[i] * deltaGrid.weight(n, m);
            }
        }

        deltaGrid.cleanGrid();
    }
}

void Problem2HForward::clrInitialConditionMatrix()
{
    f_initialMatrix.clear();
}

void Problem2HForward::initControlMeasurementDeltaGrid(std::vector<SpacePoint> &eta, std::vector<SpacePoint> &ksi)
{
    this->Nc = static_cast<unsigned int>(eta.size());
    this->eta.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        this->eta[i].resetGrid();
        this->eta[i].distributeGauss(eta[i], 1, 1);
    }

    this->No = static_cast<unsigned int>(ksi.size());
    this->ksi.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        this->ksi[j].resetGrid();
        this->ksi[j].distributeGauss(ksi[j], 1, 1);
    }
}

void Problem2HForward::prepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn)
{
    const Dimension &dimensionX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimensionY = spaceDimension(Dimension::DimensionY);
    const Dimension &time = timeDimension();

    unsigned int N = static_cast<unsigned int>(dimensionX.size());
    unsigned int M = static_cast<unsigned int>(dimensionY.size());
    unsigned int L = static_cast<unsigned int>(time.size());
    double ht = time.step();

    for (unsigned int m=0; m<=M; m++) for (unsigned int n=0; n<=N; n++) f_layerMatrix[m][n] = 0.0;

    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        if (tn.i == 2*times[s].i) { wt = 2.0/ht; } else { continue; }
        //    }

        double* v = new double[Nc];
        //for (unsigned int i=0; i<Nc; i++) v[i] = 0.1;

        for (unsigned int i=0; i<Nc; i++)
        {
            v[i] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                v[i] += k.at(i,j) * (u_info[j].vl[tn.i]-z.at(i,j));
            }
        }

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                for (unsigned int i=0; i<Nc; i++)
                {
                    f_layerMatrix[m][n] += v[i] * eta[i].weight(n, m);
                }
            }
        }

        delete [] v;

    }
}

///////////////////////////////////

double Problem2HBackward::initial(const SpaceNodePDE &, InitialCondition) const
{
    return 0.0;
}

double Problem2HBackward::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    return 0.0;
}

double Problem2HBackward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return b_layerMatrix[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HBackward::layerInfo(const DoubleMatrix &, const TimeNodePDE &) const
{

}
