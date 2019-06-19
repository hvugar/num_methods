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
    q.push_back(0.03);
    q.push_back(0.03);

    std::vector<SpacePoint> theta;
    theta.push_back(SpacePoint(0.25, 0.25));
    theta.push_back(SpacePoint(0.75, 0.75));

    ps.initInitialConditionMatrix(2, q, theta);
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

double Problem2HForward::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    if (condition == InitialCondition::FirstDerivative)
        return initialMatrix[static_cast<uint32_t>(sn.j)][static_cast<uint32_t>(sn.i)];
    else
        return 0.0;
}

double Problem2HForward::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 0.0;
}

double Problem2HForward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 0.0;
}

void Problem2HForward::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
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
    unsigned int M = static_cast<unsigned int>(dimensionX.size());

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

double Problem2HBackward::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    return 0.0;
}

double Problem2HBackward::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 0.0;
}

double Problem2HBackward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 0.0;
}

void Problem2HBackward::layerInfo(const DoubleMatrix &, const TimeNodePDE &) const
{

}
