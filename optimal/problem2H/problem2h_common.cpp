#include "problem2h_common.h"

auto ExtendedSpacePointH::contains(int nx, int ny) const -> bool
{
    return (minX <= nx && nx <= maxX && minY <= ny && ny <= maxY);
}

SpacePointInfoH::SpacePointInfoH()
{
    init(0);
}

SpacePointInfoH::SpacePointInfoH(unsigned int length)
{
    init(length);
}

SpacePointInfoH::~SpacePointInfoH()
{
    clear();
}

void SpacePointInfoH::init(unsigned int length)
{
    this->length = length;
    vl.resize(length);
    dx.resize(length);
    dy.resize(length);
}
void SpacePointInfoH::clear()
{
    dy.clear();
    dx.clear();
    vl.clear();
    length = 0;
}

//EquationParameterHE::EquationParameterHE()
//{
//    EquationParameterHE(0,0,0,0,0);
//}

//EquationParameterHE::EquationParameterHE(unsigned int Nc, unsigned int No, unsigned int Nd, unsigned int length, unsigned int gw)
//{
//    this->Nc = Nc;
//    this->No = No;
//    this->Nd = Nd;
//    k.resize(Nc, No, 0.0);
//    z.resize(Nc, No, 0.0);
//    xi.resize(No);
//    eta.resize(Nc);
//    q.resize(Nd);
//    theta.resize(Nd);

//    xi_ext.resize(No);  for (unsigned int j=0; j<No; j++) xi_ext[j].length = length;
//    eta_ext.resize(Nc); for (unsigned int i=0; i<No; i++) eta_ext[i].length = length;
//    theta.resize(Nd);
//}

DeltaGrid::DeltaGrid() : _N(0), _M(0), _hx(0.0), _hy(0.0) {}

DeltaGrid::~DeltaGrid() { cleanGrid(); }

auto DeltaGrid::initGrid(unsigned int N, double hx, unsigned int M, double hy, const SpacePoint& p, unsigned int sigmaXN, unsigned int sigmaYN) -> void
{
    this->_N = N;
    this->_M = M;
    this->_hx = hx;
    this->_hy = hy;
    this->_p = p;


    unsigned int kx = 4*sigmaXN;
    unsigned int ky = 4*sigmaYN;
    double sigmaX = hx*sigmaXN;
    double sigmaY = hy*sigmaYN;

    nwmx = new double*[M+1];
    for (unsigned int m=0; m<=M; m++) nwmx[m] = new double[N+1];

    _rx = static_cast<unsigned int>(round(p.x*N));
    _ry = static_cast<unsigned int>(round(p.y*M));

    _minX = _rx - kx; _maxX = _rx + kx;
    _minY = _ry - ky; _maxY = _ry + ky;

    double sumX = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++) sumX += exp(-((n*hx-p.x)*(n*hx-p.x))/(2.0*sigmaX*sigmaX));
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=_minY; m<=_maxY; m++) sumY += exp(-((m*hy-p.y)*(m*hy-p.y))/(2.0*sigmaY*sigmaY));
    sumY *= hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/((2.0*M_PI)*sigma);

    SpaceNodePDE sn;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=_minX; n<=_maxX; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            nwmx[m][n] = factor*exp(-0.5*(((sn.x-p.x)*(sn.x-p.x))/(sigmaX*sigmaX)+((sn.y-p.y)*(sn.y-p.y))/(sigmaY*sigmaY)));
        }
    }
}

auto DeltaGrid::cleanGrid() -> void
{
    if (_M==0 || _N == 0) return;
    for (unsigned int m=0; m<=_M; m++) delete [] nwmx[m];
    delete [] nwmx;
    _N = 0;
    _M = 0;
}

auto DeltaGrid::isCenter(const SpaceNodePDE &sn) const -> bool
{
    return (sn.i == static_cast<int>(_rx)) && (sn.j == static_cast<int>(_ry));
}

auto DeltaGrid::isCenter(unsigned int n, unsigned int m) const -> bool
{
    return (n == _rx) && (m == _ry);
}

auto DeltaGrid::isContains(const SpaceNodePDE &sn) const -> bool
{
    return (sn.i >= static_cast<int>(_minX)) && (sn.i <= static_cast<int>(_maxX)) &&
           (sn.j >= static_cast<int>(_minY)) && (sn.j <= static_cast<int>(_maxY));
}

auto DeltaGrid::isContains(unsigned int n, unsigned int m) const -> bool
{
    return (n >= _minX) && (n <= _maxX) && (m >= _minY) && (m <= _maxY);
}

auto DeltaGrid::weight(const SpaceNodePDE &sn) const -> double
{
    return nwmx[sn.j][sn.i];
}

auto DeltaGrid::weight(unsigned int n, unsigned int m) const -> double
{
    return nwmx[m][n];
}


