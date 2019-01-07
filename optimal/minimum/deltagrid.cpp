#include "deltagrid.h"

const char* DeltaGridException::what() const noexcept
{
    return "delta grid exception:";
}

DeltaGrid2D::DeltaGrid2D() : mN(0), mM(0), mhx(0.0), mhy(0.0) {}

DeltaGrid2D::~DeltaGrid2D() { cleanGrid(); }

auto DeltaGrid2D::initGrid(unsigned int N, double hx, unsigned int M, double hy) -> void
{
    this->mN = N;
    this->mM = M;
    this->mhx = hx;
    this->mhy = hy;

    m_nodes = new double*[M+1];
    for (unsigned int m=0; m<=M; m++)
    {
        m_nodes[m] = new double[N+1];
        for (unsigned int n=0; n<=N; n++) m_nodes[m][n] = 0.0;
    }

    _rows = new bool[M+1]; for (unsigned int m=0; m<=M; m++) _rows[m] = false;
    _cols = new bool[N+1]; for (unsigned int n=0; n<=N; n++) _cols[n] = false;
}

auto DeltaGrid2D::distributeGauss(const SpacePoint& sp, unsigned sigmaXNum, unsigned int sigmaYNum) -> void
{
    unsigned int kx = 3*sigmaXNum;
    unsigned int ky = 3*sigmaYNum;
    double sigmaX = mhx*sigmaXNum;
    double sigmaY = mhy*sigmaYNum;

    _rx = static_cast<unsigned int>( round(sp.x*mN) );
    _ry = static_cast<unsigned int>( round(sp.y*mM) );

    if (_rx < kx or _ry < ky or _rx > mN-kx or _ry > mM-ky) throw DeltaGridException();

    mp = sp;

    _minX = _rx - kx; _maxX = _rx + kx;
    _minY = _ry - ky; _maxY = _ry + ky;

    double sumX = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++) sumX += exp(-((n*mhx-sp.x)*(n*mhx-sp.x))/(2.0*sigmaX*sigmaX));
    sumX *= mhx;

    double sumY = 0.0;
    for (unsigned int m=_minY; m<=_maxY; m++) sumY += exp(-((m*mhy-sp.y)*(m*mhy-sp.y))/(2.0*sigmaY*sigmaY));
    sumY *= mhy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/(2.0*M_PI*sigma);

    SpaceNodePDE sn;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*mhy;
        for (unsigned int n=_minX; n<=_maxX; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*mhx;
            m_nodes[m][n] = factor*exp(-(((sn.x-sp.x)*(sn.x-sp.x))/(2.0*sigmaX*sigmaX) + ((sn.y-sp.y)*(sn.y-sp.y))/(2.0*sigmaY*sigmaY)));
        }
    }

    for (unsigned int m=_minY; m<=_maxY; m++) _rows[m] = true;
    for (unsigned int n=_minX; n<=_maxX; n++) _cols[n] = true;
}

auto DeltaGrid2D::distributeSigle(const SpacePoint& sp) -> void
{
    _rx = static_cast<unsigned int>( round(sp.x*mN) );
    _ry = static_cast<unsigned int>( round(sp.y*mM) );

    mp = sp;

    _minX = _maxX = _rx;
    _minY = _maxY = _ry;

    m_nodes[_ry][_rx] = 1.0/(mhx*mhy);
}

auto DeltaGrid2D::distributeRect4(const SpacePoint& sp) -> void
{}

auto DeltaGrid2D::cleanGrid() -> void
{
    if (mM==0 || mN == 0) return;
    for (unsigned int m=0; m<=mM; m++) delete [] m_nodes[m];
    delete [] m_nodes;
    mN = 0;
    mM = 0;
}

auto DeltaGrid2D::isCenter(const SpaceNodePDE &sn) const -> bool
{
    return (sn.i == static_cast<int>(_rx)) && (sn.j == static_cast<int>(_ry));
}

auto DeltaGrid2D::isCenter(unsigned int n, unsigned int m) const -> bool
{
    return (n == _rx) && (m == _ry);
}

auto DeltaGrid2D::isContains(const SpaceNodePDE &sn) const -> bool
{
    return (sn.i >= static_cast<int>(_minX)) && (sn.i <= static_cast<int>(_maxX)) &&
           (sn.j >= static_cast<int>(_minY)) && (sn.j <= static_cast<int>(_maxY));
}

auto DeltaGrid2D::isContains(unsigned int n, unsigned int m) const -> bool
{
    return (n >= _minX) && (n <= _maxX) &&
           (m >= _minY) && (m <= _maxY);
}

auto DeltaGrid2D::weight(const SpaceNodePDE &sn) const -> double
{
    return m_nodes[sn.j][sn.i];
}

auto DeltaGrid2D::weight(unsigned int n, unsigned int m) const -> double
{
    return m_nodes[m][n];
}

auto DeltaGrid2D::rx() const -> unsigned int { return _rx; }
auto DeltaGrid2D::ry() const -> unsigned int { return _ry; }

auto DeltaGrid2D::p() const -> const SpacePoint& { return mp; }
auto DeltaGrid2D::p() -> SpacePoint& { return mp; }

auto DeltaGrid2D::minX() const -> unsigned int { return _minX; }
auto DeltaGrid2D::maxX() const -> unsigned int { return _maxX; }
auto DeltaGrid2D::minY() const -> unsigned int { return _minY; }
auto DeltaGrid2D::maxY() const -> unsigned int { return _maxY; }
