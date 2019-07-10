#include "problem2h_solver.h"

void Problem2HSolver::Main(int argc UNUSED_PARAM, char* argv[] UNUSED_PARAM)
{
    Problem2HForward ps;
    Dimension dimensionX(0.01, 0, 100);
    Dimension dimensionY(0.01, 0, 100);
    Dimension dimensionT(0.01, 0, 1000);

    //ps.setTimeDimension(dimensionT);
    //ps.setSpaceDimensionX(dimensionX);
    //ps.setSpaceDimensionY(dimensionY);
    ps.initMatrixes(dimensionX, dimensionY, dimensionT, 2, 2, 2);

    ps.setDimension(dimensionX, dimensionY, dimensionT);

    ps.setWaveSpeed(1.0);
    ps.setWaveDissipation(0.01);

    ps.k.resize(2, 2, -1.0);
    ps.z.resize(2, 2, +0.0);

    unsigned int Nq = 2;
    SpacePoint* zta = new SpacePoint[2];
    zta[0] = SpacePoint(0.25, 0.25);
    zta[1] = SpacePoint(0.75, 0.75);
    double* q = new double[Nq];
    q[0] = 0.05;
    q[1] = 0.05;
    ps.setInitialConditionMatrix(zta, q, Nq);

    SpacePoint* eta = new SpacePoint[2];
    eta[0] = SpacePoint(0.65, 0.34);
    eta[1] = SpacePoint(0.25, 0.75);

    SpacePoint* ksi = new SpacePoint[2];
    ksi[0] = SpacePoint(0.22, 0.54);
    ksi[1] = SpacePoint(0.82, 0.27);

    ps.initControlMeasurementDeltaGrid(2, 2);
    ps.distributeControlMeasurementDeltaGrid(eta, 2, ksi, 2);
    ps.implicit_calculate_D2V1();
    ps.clearControlMeasurementDeltaGrid();
    puts("OK");

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

//    f_initialMatrix.clear();
//    f_initialMatrix.resize(M+1, N+1, 0.0);

//    f_layerMatrix.clear();
//    f_layerMatrix.resize(M+1, N+1, 0.0);

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
//    f_initialMatrix.clear();
//    f_layerMatrix.clear();
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

double Problem2HForward::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    return 0.0;
}

double Problem2HForward::f(const SpaceNodePDE &sn, const TimeNodePDE &) const
{
    return f_crLayerMatrix[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HForward::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    if (tn.i == 0 || tn.i == 1 || tn.i == 2 || tn.i == 3)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
    }

    Problem2HForward* pf = const_cast<Problem2HForward*>(this);
    for (unsigned int j=0; j<No; j++)
    {
        double u_vl, u_dx, u_dy;
        u_vl = _deltaGridMeasurement[j].lumpPointGauss(u, u_dx, u_dy);
        pf->u__info[j].vl[tn.i] = u_vl;
        pf->u__info[j].dx[tn.i] = u_dx;
        pf->u__info[j].dy[tn.i] = u_dy;
    }
    const_cast<Problem2HForward*>(this)->prepareLayerMatrix(u, tn);

    //save2TextFile(u, tn);
}

void Problem2HForward::save2TextFile(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn) const
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

void Problem2HForward::setDimension(const Dimension& dimensionX, const Dimension& dimensionY, const Dimension &timeDimension)
{
    setSpaceDimensionX(dimensionX);
    setSpaceDimensionY(dimensionY);
    setTimeDimension(timeDimension);

    f_initialMatrix.resize(static_cast<unsigned int>(dimensionX.size())+1,
                           static_cast<unsigned int>(dimensionY.size())+1);
    f_crLayerMatrix.resize(static_cast<unsigned int>(dimensionX.size())+1,
                           static_cast<unsigned int>(dimensionY.size())+1);
}

void Problem2HForward::setInitialConditionMatrix(const SpacePoint* zta, const double* q, unsigned int Nq)
{
    const unsigned int N = static_cast<unsigned int>(_spaceDimensionX.size());
    const unsigned int M = static_cast<unsigned int>(_spaceDimensionY.size());
    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();

    for (unsigned int s=0; s<Nq; s++)
    {
        DeltaGrid2D deltaGrid;
        deltaGrid.initGrid(N, hx, M, hy);
        deltaGrid.distributeGauss(zta[s], 5, 5);

        const unsigned minX = deltaGrid.minX();
        const unsigned maxX = deltaGrid.maxX();
        const unsigned minY = deltaGrid.minY();
        const unsigned maxY = deltaGrid.maxY();

        for (unsigned int m=minY; m<=maxY; m++)
        {
            for (unsigned int n=minX; n<=maxX; n++)
            {
                f_initialMatrix[m][n] += q[s] * deltaGrid.weight(n, m);
            }
        }

        deltaGrid.cleanGrid();
    }
}

void Problem2HForward::clrInitialConditionMatrix()
{
    f_initialMatrix.clear();
}

void Problem2HForward::initControlMeasurementDeltaGrid(unsigned int Nc, unsigned int No)
{
    const unsigned int N = static_cast<unsigned int>(_spaceDimensionX.size());
    const unsigned int M = static_cast<unsigned int>(_spaceDimensionY.size());
    const unsigned int L = static_cast<unsigned int>(_timeDimension.size());
    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();
    const unsigned int length = 2*L+1;

    _deltaGridControl = new DeltaGrid2D[Nc];
    for (unsigned int i=0; i<Nc; i++) _deltaGridControl[i].initGrid(N, hx, M, hy);

    _deltaGridMeasurement = new DeltaGrid2D[No];
    u__info = new SpacePointInfo[No];
    for (unsigned int j=0; j<No; j++)
    {
        _deltaGridMeasurement[j].initGrid(N, hx, M, hy);

        u__info[j].vl.resize(length);
        u__info[j].dx.resize(length);
        u__info[j].dy.resize(length);
    }
}

void Problem2HForward::distributeControlMeasurementDeltaGrid(const SpacePoint *eta, unsigned int Nc, const SpacePoint *ksi, unsigned int No)
{
    for (unsigned int i=0; i<Nc; i++)
    {
        _deltaGridControl[i].resetGrid();
        _deltaGridControl[i].distributeGauss(eta[i], 1, 1);
    }
    for (unsigned int j=0; j<No; j++)
    {
        _deltaGridMeasurement[j].resetGrid();
        _deltaGridMeasurement[j].distributeGauss(ksi[j], 1, 1);
    }
}

void Problem2HForward::prepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn)
{
    //const Dimension &dimensionX = spaceDimension(Dimension::DimensionX);
    //const Dimension &dimensionY = spaceDimension(Dimension::DimensionY);
    const Dimension &time = timeDimension();

    //unsigned int N = static_cast<unsigned int>(dimensionX.size());
    //unsigned int M = static_cast<unsigned int>(dimensionY.size());
    //unsigned int L = static_cast<unsigned int>(time.size());
    double ht = time.step();

    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        if (tn.i == 2*times[s].i) { wt = 2.0/ht; } else { continue; }

        double* v = new double[Nc];

        for (unsigned int i=0; i<Nc; i++)
        {
            v[i] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                v[i] += k.at(i,j) * (u__info[j].vl[tn.i]-z.at(i,j));
            }
        }

        for (unsigned int i=0; i<Nc; i++)
        {
            const unsigned minX = _deltaGridControl[i].minX();
            const unsigned maxX = _deltaGridControl[i].maxX();
            const unsigned minY = _deltaGridControl[i].minY();
            const unsigned maxY = _deltaGridControl[i].maxY();

            for (unsigned int m=minY; m<=maxY; m++)
            {
                for (unsigned int n=minX; n<=maxX; n++)
                {
                    if (i==0) f_crLayerMatrix[m][n] = 0.0;
                    f_crLayerMatrix[m][n] += v[i] * eta[i].weight(n, m);
                }
            }
        }

        delete [] v;
    }
}

void Problem2HForward::clearControlMeasurementDeltaGrid()
{
//    for (unsigned int i=0; i<Nc; i++) _deltaGridControl[i].cleanGrid();s
//    for (unsigned int j=0; j<No; j++) _deltaGridMeasurement[j].cleanGrid();

//    delete [] _deltaGridControl;
//    delete [] _deltaGridMeasurement;
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
