#include "problem2h_solver.h"

void Problem2HSolver::Main(int argc, char **argv)
{
    Problem2HForward ps;
    ps.setTimeDimension(Dimension(0.01, 0, 1000));
    ps.addSpaceDimension(Dimension(0.01, 0, 100));
    ps.addSpaceDimension(Dimension(0.01, 0, 100));
    ps.setWaveSpeed(1.0);
    ps.setWaveDissipation(0.01);

    std::vector<double> q;
    q.push_back(0.05);
    q.push_back(0.05);

    std::vector<SpacePoint> theta;
    theta.push_back(SpacePoint(0.25, 0.25));
    theta.push_back(SpacePoint(0.75, 0.75));

    ps.initInitialConditionMatrix(2, q, theta);

    std::vector<SpacePoint> eta;
    eta.push_back(SpacePoint(0.65, 0.34));
    eta.push_back(SpacePoint(0.25, 0.75));

    std::vector<SpacePoint> ksi;
    ksi.push_back(SpacePoint(0.22, 0.54));
    ksi.push_back(SpacePoint(0.82, 0.27));

    ps.initControlMeasurementDeltaGrid(eta, ksi);
    ps.implicit_calculate_D2V1();
}

//Problem2HForward& Problem2HSolver::forward()
//{
//    return *this;
//}

//Problem2HBackward& Problem2HSolver::backward()
//{
//    return *this;
//}

void Problem2HCommon::initMatrixes(const Dimension &dimensionX, const Dimension &dimensionY)
{
    f_initialMatrix.clear();
    f_initialMatrix.resize(static_cast<uint32_t>(dimensionY.size())+1,
                           static_cast<uint32_t>(dimensionX.size())+1, 0.0);
    f_layerMatrix.clear();
    f_layerMatrix.resize(static_cast<uint32_t>(dimensionY.size())+1,
                         static_cast<uint32_t>(dimensionX.size())+1, 0.0);
    b_initialMatrix.clear();
    b_initialMatrix.resize(static_cast<uint32_t>(dimensionY.size())+1,
                           static_cast<uint32_t>(dimensionX.size())+1, 0.0);
    b_layerMatrix.clear();
    b_layerMatrix.resize(static_cast<uint32_t>(dimensionY.size())+1,
                         static_cast<uint32_t>(dimensionX.size())+1, 0.0);
}

double Problem2HForward::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    if (condition == InitialCondition::FirstDerivative)
        return f_initialMatrix[static_cast<uint32_t>(sn.j)][static_cast<uint32_t>(sn.i)];
    else
        return 0.0;
}

double Problem2HForward::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 0.0;
}

double Problem2HForward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return f_layerMatrix[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HForward::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
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

void Problem2HForward::initInitialConditionMatrix(uint32_t Nt, const std::vector<double> &q, const std::vector<SpacePoint> &theta)
{
    const Dimension &dimensionX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimensionY = spaceDimension(Dimension::DimensionY);
    unsigned int N = static_cast<unsigned int>(dimensionX.size());
    unsigned int M = static_cast<unsigned int>(dimensionY.size());

    this->q = q;
    this->theta = theta;

    thetaDeltaGrid.resize(Nt);
    for (unsigned int i=0; i<Nt; i++)
    {
        thetaDeltaGrid[i].cleanGrid();
        thetaDeltaGrid[i].initGrid(static_cast<uint32_t>(dimensionX.size()), dimensionX.step(),
                                   static_cast<uint32_t>(dimensionY.size()), dimensionX.step());
        thetaDeltaGrid[i].distributeGauss(theta[i], 5, 5);
    }

    initialMatrix.resize(static_cast<uint32_t>(dimensionY.size())+1, static_cast<uint32_t>(dimensionX.size())+1, 0.0);
    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            for (unsigned int i=0; i<Nt; i++)
            {
                initialMatrix[m][n] += q[i] * thetaDeltaGrid[i].weight(n, m);
            }
        }
    }
}

void Problem2HForward::clearInitialConditionMatrix()
{
    this->q.clear();
    this->theta.clear();

    for (unsigned int i=0; i<Nt; i++) thetaDeltaGrid[i].cleanGrid();
    initialMatrix.clear();
}

void Problem2HForward::initControlMeasurementDeltaGrid(std::vector<SpacePoint> &eta, std::vector<SpacePoint> &ksi)
{
    const Dimension &time = timeDimension();
    const Dimension &dimensionX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimensionY = spaceDimension(Dimension::DimensionY);

    const unsigned int length = 2*static_cast<unsigned int>(time.size())+1;

    this->Nc = static_cast<unsigned int>(eta.size());
    this->eta.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        this->eta[i].vl.resize(length);
        this->eta[i].dx.resize(length);
        this->eta[i].dy.resize(length);
        this->eta[i].deltaGrid.cleanGrid();
        this->eta[i].deltaGrid.initGrid(static_cast<uint32_t>(dimensionX.size()), dimensionX.step(),
                                        static_cast<uint32_t>(dimensionY.size()), dimensionX.step());
        this->eta[i].deltaGrid.distributeGauss(eta[i]);
    }

    this->No = static_cast<unsigned int>(ksi.size());
    this->ksi.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        this->ksi[j].vl.resize(length);
        this->ksi[j].dx.resize(length);
        this->ksi[j].dy.resize(length);
        this->ksi[j].deltaGrid.cleanGrid();
        this->ksi[j].deltaGrid.initGrid(static_cast<uint32_t>(dimensionX.size()), dimensionX.step(),
                                        static_cast<uint32_t>(dimensionY.size()), dimensionX.step());
        this->ksi[j].deltaGrid.distributeGauss(ksi[j]);
    }

    k.resize(Nc, No, -0.1);
    z.resize(Nc, No, 0.001);

    layerMatrix.resize(static_cast<uint32_t>(dimensionY.size())+1,
                       static_cast<uint32_t>(dimensionX.size())+1, 0.0);
}

void Problem2HForward::prepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn)
{
    const Dimension &dimensionX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimensionY = spaceDimension(Dimension::DimensionY);
    unsigned int N = static_cast<unsigned int>(dimensionX.size());
    unsigned int M = static_cast<unsigned int>(dimensionY.size());

    double u_vl, u_dx, u_dy, a1, b1;
    for (unsigned int j=0; j<No; j++)
    {
        u_vl = ksi[j].deltaGrid.lumpPointGauss(u, u_dx, u_dy, a1, b1);
        ksi[j].vl[tn.i] = u_vl;
        ksi[j].dx[tn.i] = u_dx;
        ksi[j].dy[tn.i] = u_dy;
    }

    double* v = new double[Nc];
    for (unsigned int i=0; i<Nc; i++)
    {
        v[i] = 0.0;
        for (unsigned int j=0; j<No; j++)
        {
            v[i] += k.at(i,j) * (ksi[j].vl[tn.i]-z.at(i,j));
        }
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            layerMatrix[m][n] = 0.0;
            for (unsigned int i=0; i<Nc; i++)
            {
                layerMatrix[m][n] += v[i] * eta[i].deltaGrid.weight(n, m);
            }
        }
    }

    delete [] v;
}

double Problem2HBackward::initial(const SpaceNodePDE &, InitialCondition) const
{
    return 0.0;
}

double Problem2HBackward::boundary(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2HBackward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return b_layerMatrix[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HBackward::layerInfo(const DoubleMatrix &, const TimeNodePDE &) const
{

}
