#include "problem2h_solver.h"

void Problem2HSolver::Main(int argc UNUSED_PARAM, char* argv[] UNUSED_PARAM)
{
    Problem2HSolver ps;
    ps.setDimensions(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, 300));
    ps.setEquationParameters(1.0, 0.0);

    ps.Nt = 10;
    ps.times.resize(ps.Nt);
    ps.times[0] = TimeNodePDE(30, 0.3);
    ps.times[1] = TimeNodePDE(60, 0.6);
    ps.times[2] = TimeNodePDE(90, 0.9);
    ps.times[3] = TimeNodePDE(120, 1.2);
    ps.times[4] = TimeNodePDE(150, 1.5);
    ps.times[5] = TimeNodePDE(180, 1.8);
    ps.times[6] = TimeNodePDE(210, 2.1);
    ps.times[7] = TimeNodePDE(240, 2.4);
    ps.times[8] = TimeNodePDE(270, 2.7);
    ps.times[9] = TimeNodePDE(300, 3.0);

    unsigned int Nq = 2;
    SpacePoint* zta = new SpacePoint[2];
    zta[0] = SpacePoint(0.25, 0.25);
    zta[1] = SpacePoint(0.75, 0.75);
    double* q = new double[Nq];
    q[0] = 0.05;
    q[1] = 0.05;
    ps.setInitialConditionMatrix(zta, q, Nq);

    std::vector<SpacePoint> eta(2);
    eta[0] = SpacePoint(0.65, 0.34);
    eta[1] = SpacePoint(0.25, 0.75);

    std::vector<SpacePoint> ksi(2);
    ksi[0] = SpacePoint(0.22, 0.54);
    ksi[1] = SpacePoint(0.82, 0.27);

    DoubleMatrix k(2, 2, -1.0);
    DoubleMatrix z(2, 2, +0.0);

    ps.setParameterCounts(2, 2, ps.Problem2HWaveEquationIBVP::spaceDimensionX(),
                                ps.Problem2HWaveEquationIBVP::spaceDimensionY(),
                                ps.Problem2HWaveEquationIBVP::timeDimension());
    ps.setOptimizedParameters(k, z, eta, ksi);
    ps.distributeControlDeltaGrid();
    ps.distributeMeasurementDeltaGrid();

    ps.Problem2HWaveEquationIBVP::implicit_calculate_D2V1();

    delete [] q;
    delete [] zta;

    puts("Finished");
}

void Problem2HSolver::setDimensions(const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension)
{
    Problem2HWaveEquationIBVP::setTimeDimension(timeDimension);
    Problem2HWaveEquationIBVP::setSpaceDimensions(dimensionX, dimensionY);

    Problem2HConjugateWaveEquationIBVP::setTimeDimension(timeDimension);
    Problem2HConjugateWaveEquationIBVP::setSpaceDimensions(dimensionX, dimensionY);
}

void Problem2HSolver::setEquationParameters(double waveSpeed, double waveDissipation)
{
    Problem2HWaveEquationIBVP::setWaveSpeed(waveSpeed);
    Problem2HWaveEquationIBVP::setWaveDissipation(waveDissipation);

    Problem2HConjugateWaveEquationIBVP::setWaveSpeed(waveSpeed);
    Problem2HConjugateWaveEquationIBVP::setWaveDissipation(waveDissipation);
}

void Problem2HWaveEquationIBVP::setSpaceDimensions(const Dimension &dimensionX, const Dimension &dimensionY)
{
    IWaveEquationIBVP::setSpaceDimensions(dimensionX, dimensionY);
    f_initialMatrix.resize(static_cast<unsigned int>(dimensionY.size())+1, static_cast<unsigned int>(dimensionX.size())+1, 0.0);
    f_crLayerMatrix.resize(static_cast<unsigned int>(dimensionY.size())+1, static_cast<unsigned int>(dimensionX.size())+1, 0.0);
}

void Problem2HConjugateWaveEquationIBVP::setSpaceDimensions(const Dimension &dimensionX, const Dimension &dimensionY)
{
    InitialBoundaryValueProblemPDE::setSpaceDimensions(dimensionX, dimensionY);
    b_initialMatrix.resize(static_cast<unsigned int>(dimensionY.size())+1, static_cast<unsigned int>(dimensionX.size())+1, 0.0);
    b_crLayerMatrix.resize(static_cast<unsigned int>(dimensionY.size())+1, static_cast<unsigned int>(dimensionX.size())+1, 0.0);
}

Problem2HCommon::~Problem2HCommon() {}

void Problem2HCommon::setParameterCounts(unsigned int Nc, unsigned int No, const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension)
{
    this->Nc = Nc;
    this->No = No;

    const unsigned int N = static_cast<unsigned int>(dimensionX.size());
    const unsigned int M = static_cast<unsigned int>(dimensionY.size());
    const unsigned int L = static_cast<unsigned int>(timeDimension.size());
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();
    const unsigned int length = 2*L+1;

    _deltaGridControl = new DeltaGrid2D[Nc];
    p_info = new SpacePointInfo[Nc];
    for (unsigned int i=0; i<Nc; i++)
    {
        _deltaGridControl[i].initGrid(N, hx, M, hy);
        p_info[i].vl.resize(length);
        p_info[i].dx.resize(length);
        p_info[i].dy.resize(length);
    }

    _deltaGridMeasurement = new DeltaGrid2D[No];
    u_info = new SpacePointInfo[No];
    for (unsigned int j=0; j<No; j++)
    {
        _deltaGridMeasurement[j].initGrid(N, hx, M, hy);
        u_info[j].vl.resize(length);
        u_info[j].dx.resize(length);
        u_info[j].dy.resize(length);
    }
}

void Problem2HCommon::setOptimizedParameters(const DoubleMatrix &k, const DoubleMatrix &z, const std::vector<SpacePoint> &eta, const std::vector<SpacePoint> &ksi)
{
    this->k = k;
    this->z = z;
    this->eta = eta;
    this->ksi = ksi;
}

void Problem2HCommon::distributeControlDeltaGrid()
{
    for (unsigned int i=0; i<Nc; i++)
    {
        _deltaGridControl[i].resetGrid();
        _deltaGridControl[i].distributeGauss(eta[i], 1, 1);
    }
}

void Problem2HCommon::distributeMeasurementDeltaGrid()
{
    for (unsigned int j=0; j<No; j++)
    {
        _deltaGridMeasurement[j].resetGrid();
        _deltaGridMeasurement[j].distributeGauss(ksi[j], 1, 1);
    }
}

///////////////////////////////////

double Problem2HWaveEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    if (condition == InitialCondition::FirstDerivative)
        return f_initialMatrix[static_cast<uint32_t>(sn.j)][static_cast<uint32_t>(sn.i)];
    else
        return 0.0;
}

double Problem2HWaveEquationIBVP::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    return 0.0;
}

double Problem2HWaveEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &) const
{
    return f_crLayerMatrix[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HWaveEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    if (tn.i == 0 || tn.i == 1 || tn.i == 2 || tn.i == 3)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
    }

    const_cast<Problem2HWaveEquationIBVP*>(this)->layerInfoPrepareLayerMatrix(u, tn);
    layerInfoSave2TextFile(u, tn);
}

void Problem2HWaveEquationIBVP::layerInfoSave2TextFile(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn) const
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

void Problem2HWaveEquationIBVP::setInitialConditionMatrix(const SpacePoint* zta, const double* q, unsigned int Nq)
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

void Problem2HWaveEquationIBVP::clrInitialConditionMatrix()
{
    f_initialMatrix.clear();
}

void Problem2HWaveEquationIBVP::layerInfoPrepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn)
{
    //const Dimension &dimensionX = spaceDimension(Dimension::DimensionX);
    //const Dimension &dimensionY = spaceDimension(Dimension::DimensionY);
    const Dimension &time = timeDimension();

    //unsigned int N = static_cast<unsigned int>(dimensionX.size());
    //unsigned int M = static_cast<unsigned int>(dimensionY.size());
    //unsigned int L = static_cast<unsigned int>(time.size());
    double ht = time.step();

    for (unsigned int j=0; j<No; j++)
    {
        double u_vl, u_dx, u_dy;
        u_vl = _deltaGridMeasurement[j].lumpPointGauss(u, u_dx, u_dy);
        u_info[j].vl[tn.i] = u_vl;
        u_info[j].dx[tn.i] = u_dx;
        u_info[j].dy[tn.i] = u_dy;
    }

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
                v[i] += k.at(i,j) * (u_info[j].vl[tn.i]-z.at(i,j));
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
                    f_crLayerMatrix[m][n] += v[i] * _deltaGridControl[i].weight(n, m);
                }
            }
        }

        delete [] v;
    }
}

double Problem2HConjugateWaveEquationIBVP::initial(const SpaceNodePDE &, InitialCondition) const
{
    return 0.0;
}

double Problem2HConjugateWaveEquationIBVP::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    return 0.0;
}

double Problem2HConjugateWaveEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return b_crLayerMatrix[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HConjugateWaveEquationIBVP::layerInfo(const DoubleMatrix &, const TimeNodePDE &) const
{

}
