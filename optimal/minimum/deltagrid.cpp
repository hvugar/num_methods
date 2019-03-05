#include "deltagrid.h"

DeltaGridException::DeltaGridException(const std::string &msg) : message(msg) {}

const char* DeltaGridException::what() const noexcept
{
    std::string msg = "DeltaGridException: " + message;
    return msg.data();
}

DeltaGrid2D::DeltaGrid2D() : _N(0), _M(0), _hx(0.0), _hy(0.0) {}

DeltaGrid2D::~DeltaGrid2D() { cleanGrid(); }

auto DeltaGrid2D::initGrid(unsigned int N, double hx, unsigned int M, double hy) -> void
{
    this->_N = N;
    this->_M = M;
    this->_hx = hx;
    this->_hy = hy;

    m_nodes = new double*[M+1];
    for (unsigned int m=0; m<=M; m++)
    {
        m_nodes[m] = new double[N+1];
        for (unsigned int n=0; n<=N; n++) m_nodes[m][n] = 0.0;
    }

    _rows = new bool[M+1]; for (unsigned int m=0; m<=M; m++) _rows[m] = false;
    _cols = new bool[N+1]; for (unsigned int n=0; n<=N; n++) _cols[n] = false;
}

auto DeltaGrid2D::distributeGauss(const SpacePoint& sp, unsigned int sigmaXNum, unsigned int sigmaYNum) -> void
{
    unsigned int kx = 3*sigmaXNum;
    unsigned int ky = 3*sigmaYNum;
    double sigmaX = _hx*sigmaXNum;
    double sigmaY = _hy*sigmaYNum;

    _rx = static_cast<unsigned int>( round(sp.x*_N) );
    _ry = static_cast<unsigned int>( round(sp.y*_M) );

    if (_rx < kx or _ry < ky or _rx > _N-kx or _ry > _M-ky) throw DeltaGridException("Point:["+std::to_string(sp.x)+","+std::to_string(sp.y)+"]");

    _p = sp;

    _minX = _rx - kx; _maxX = _rx + kx;
    _minY = _ry - ky; _maxY = _ry + ky;

    double sumX = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++) sumX += exp(-((n*_hx-sp.x)*(n*_hx-sp.x))/(2.0*sigmaX*sigmaX));
    sumX *= _hx;

    double sumY = 0.0;
    for (unsigned int m=_minY; m<=_maxY; m++) sumY += exp(-((m*_hy-sp.y)*(m*_hy-sp.y))/(2.0*sigmaY*sigmaY));
    sumY *= _hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/(2.0*M_PI*sigma);

    SpaceNodePDE sn;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*_hy;
        for (unsigned int n=_minX; n<=_maxX; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*_hx;
            m_nodes[m][n] = factor*exp(-(((sn.x-sp.x)*(sn.x-sp.x))/(2.0*sigmaX*sigmaX) +
                                         ((sn.y-sp.y)*(sn.y-sp.y))/(2.0*sigmaY*sigmaY)));
        }
    }

    for (unsigned int m=_minY; m<=_maxY; m++) _rows[m] = true;
    for (unsigned int n=_minX; n<=_maxX; n++) _cols[n] = true;
}

auto DeltaGrid2D::distributeSigle(const SpacePoint& sp) -> void
{
    _rx = static_cast<unsigned int>( round(sp.x*_N) );
    _ry = static_cast<unsigned int>( round(sp.y*_M) );

    _p = sp;

    _minX = _maxX = _rx;
    _minY = _maxY = _ry;

    m_nodes[_ry][_rx] = 1.0/(_hx*_hy);
}

auto DeltaGrid2D::distributeRect4(const SpacePoint&) -> void
{}

auto DeltaGrid2D::consentrateInPoint(const DoubleMatrix &u) const -> double
{
//    double pu = 0.0;
//    for (unsigned int m=minY(); m<=maxY(); m++)
//    {
//        for (unsigned int n=minX(); n<=maxX(); n++)
//        {
//            pu += u[m][n] * weight(n,m) * _hx * _hy;
//        }
//    }
//    return pu;

//    double pu = 0.0;
//    const unsigned int rx = static_cast<unsigned int>(_rx);
//    const unsigned int ry = static_cast<unsigned int>(_ry);
//    const double px = p().x;
//    const double py = p().y;
//    double x0, x1, x2; x0 = (rx-1)*_hx; x1 = rx*_hx; x2 = (rx+1)*_hx;
//    double y0, y1, y2; y0 = (ry-1)*_hy; y1 = ry*_hy; y2 = (ry+1)*_hy;
//    pu = (((px-x1)*(px-x2))/((x0-x1)*(x0-x2))) * ( (((py-y1)*(py-y2))/((y0-y1)*(y0-y2)))*u[ry-1][rx-1]
//                                                 + (((py-y2)*(py-y0))/((y1-y2)*(y1-y0)))*u[ry+0][rx-1] +
//                                                   (((py-y0)*(py-y1))/((y2-y0)*(y2-y1)))*u[ry+1][rx-1] )
//       + (((px-x2)*(px-x0))/((x1-x2)*(x1-x0))) * ( (((py-y1)*(py-y2))/((y0-y1)*(y0-y2)))*u[ry-1][rx+0]
//                                                 + (((py-y2)*(py-y0))/((y1-y2)*(y1-y0)))*u[ry+0][rx+0] +
//                                                   (((py-y0)*(py-y1))/((y2-y0)*(y2-y1)))*u[ry+1][rx+0] )
//       + (((px-x0)*(px-x1))/((x2-x0)*(x2-x1))) * ( (((py-y1)*(py-y2))/((y0-y1)*(y0-y2)))*u[ry-1][rx+1]
//                                                 + (((py-y2)*(py-y0))/((y1-y2)*(y1-y0)))*u[ry+0][rx+1] +
//                                                   (((py-y0)*(py-y1))/((y2-y0)*(y2-y1)))*u[ry+1][rx+1] );
//    return pu;

    double pu = 0.0;
    const unsigned int rx = static_cast<unsigned int>(_rx);
    const unsigned int ry = static_cast<unsigned int>(_ry);
    const double px = p().x;
    const double py = p().y;
    double x0, x1, x2, x3; x0 = (rx-1)*_hx; x1 = rx*_hx; x2 = (rx+1)*_hx; x3 = (rx+2)*_hx;
    double y0, y1, y2, y3; y0 = (ry-1)*_hy; y1 = ry*_hy; y2 = (ry+1)*_hy; y3 = (ry+2)*_hy;

    double Lx0 = (px-x1)*(px-x2)*(px-x3); double L00 = (x0-x1)*(x0-x2)*(x0-x3);
    double Lx1 = (px-x0)*(px-x2)*(px-x3); double L11 = (x1-x0)*(x1-x2)*(x1-x3);
    double Lx2 = (px-x0)*(px-x1)*(px-x3); double L22 = (x2-x0)*(x2-x1)*(x2-x3);
    double Lx3 = (px-x0)*(px-x1)*(px-x2); double L33 = (x3-x0)*(x3-x1)*(x3-x2);

    double Ry0 = (py-y1)*(py-y2)*(py-y3); double R00 = (y0-y1)*(y0-y2)*(y0-y3);
    double Ry1 = (py-y0)*(py-y2)*(py-y3); double R11 = (y1-y0)*(y1-y2)*(y1-y3);
    double Ry2 = (py-y0)*(py-y1)*(py-y3); double R22 = (y2-y0)*(y2-y1)*(y2-y3);
    double Ry3 = (py-y0)*(py-y1)*(py-y2); double R33 = (y3-y0)*(y3-y1)*(y3-y2);

    pu = (Lx0/L00) * ( (Ry0/R00)*u[ry-1][rx-1] + (Ry1/R11)*u[ry+0][rx-1] + (Ry2/R22)*u[ry+1][rx-1] + (Ry3/R33)*u[ry+2][rx-1] )
       + (Lx1/L11) * ( (Ry0/R00)*u[ry-1][rx+0] + (Ry1/R11)*u[ry+0][rx+0] + (Ry2/R22)*u[ry+1][rx+0] + (Ry3/R33)*u[ry+2][rx+0] )
       + (Lx2/L22) * ( (Ry0/R00)*u[ry-1][rx+1] + (Ry1/R11)*u[ry+0][rx+1] + (Ry2/R22)*u[ry+1][rx+1] + (Ry3/R33)*u[ry+2][rx+1] )
       + (Lx3/L33) * ( (Ry0/R00)*u[ry-1][rx+2] + (Ry1/R11)*u[ry+0][rx+2] + (Ry2/R22)*u[ry+1][rx+2] + (Ry3/R33)*u[ry+2][rx+2] );
    return pu;
}

auto DeltaGrid2D::consentrateInPoint(const DoubleMatrix &u, double &dx, double &dy) const -> double
{
    const unsigned int rx = static_cast<unsigned int>(_rx);
    const unsigned int ry = static_cast<unsigned int>(_ry);

    const double px = p().x;
    const double py = p().y;

//    dx = (u[ry][rx+1] - u[ry][rx-1])/(2.0*_hx);
//    dy = (u[ry+1][rx] - u[ry-1][rx])/(2.0*_hy);

//    dx += (px-rx*_hx)*((u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1])/(_hx*_hx));
//    dy += (py-ry*_hy)*((u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx])/(_hy*_hy));

//    dx = (u[ry][rx-2] - 8.0*u[ry][rx-1] + 8.0*u[ry][rx+1] - u[ry][rx+2])/(12.0*_hx);
//    dy = (u[ry-2][rx] - 8.0*u[ry-1][rx] + 8.0*u[ry+1][rx] - u[ry+2][rx])/(12.0*_hy);

//    dx += ((px-rx*_hx))*((-2.0*u[ry][rx-2] + 32.0*u[ry][rx-1] - 60.0*u[ry][rx] + 32.0*u[ry][rx+1] - 2.0*u[ry][rx+2])/(24.0*_hx*_hx));
//    dy += ((py-ry*_hy ))*((-2.0*u[ry-2][rx] + 32.0*u[ry-1][rx] - 60.0*u[ry][rx] + 32.0*u[ry+1][rx] - 2.0*u[ry+2][rx])/(24.0*_hy*_hy));

//    double x0, x1, x2; x0 = (rx-1)*_hx; x1 = rx*_hx; x2 = (rx+1)*_hx;
//    double y0, y1, y2; y0 = (ry-1)*_hy; y1 = ry*_hy; y2 = (ry+1)*_hy;
//    dx = (((px-x1)+(px-x2))/((x0-x1)*(x0-x2))) * ( (((py-y1)*(py-y2))/((y0-y1)*(y0-y2)))*u[ry-1][rx-1]
//                                                 + (((py-y2)*(py-y0))/((y1-y2)*(y1-y0)))*u[ry+0][rx-1] +
//                                                   (((py-y0)*(py-y1))/((y2-y0)*(y2-y1)))*u[ry+1][rx-1] )
//       + (((px-x2)+(px-x0))/((x1-x2)*(x1-x0))) * ( (((py-y1)*(py-y2))/((y0-y1)*(y0-y2)))*u[ry-1][rx+0]
//                                                 + (((py-y2)*(py-y0))/((y1-y2)*(y1-y0)))*u[ry+0][rx+0] +
//                                                   (((py-y0)*(py-y1))/((y2-y0)*(y2-y1)))*u[ry+1][rx+0] )
//       + (((px-x0)+(px-x1))/((x2-x0)*(x2-x1))) * ( (((py-y1)*(py-y2))/((y0-y1)*(y0-y2)))*u[ry-1][rx+1]
//                                                 + (((py-y2)*(py-y0))/((y1-y2)*(y1-y0)))*u[ry+0][rx+1] +
//                                                   (((py-y0)*(py-y1))/((y2-y0)*(y2-y1)))*u[ry+1][rx+1] );

//    dy = (((py-y1)+(py-y2))/((y0-y1)*(y0-y2))) * ( (((px-x1)*(px-x2))/((x0-x1)*(x0-x2)))*u[ry-1][rx-1]
//                                                 + (((px-x2)*(px-x0))/((x1-x2)*(x1-x0)))*u[ry-1][rx+0] +
//                                                   (((px-x0)*(px-x1))/((x2-x0)*(x2-x1)))*u[ry-1][rx+1] )
//       + (((py-y2)+(py-y0))/((y1-y2)*(y1-y0))) * ( (((px-x1)*(px-x2))/((x0-x1)*(x0-x2)))*u[ry+0][rx-1]
//                                                 + (((px-x2)*(px-x0))/((x1-x2)*(x1-x0)))*u[ry+0][rx+0] +
//                                                   (((px-x0)*(px-x1))/((x2-x0)*(x2-x1)))*u[ry+0][rx+1] )
//       + (((py-y0)+(py-y1))/((y2-y0)*(y2-y1))) * ( (((px-x1)*(px-x2))/((x0-x1)*(x0-x2)))*u[ry+1][rx-1]
//                                                 + (((px-x2)*(px-x0))/((x1-x2)*(x1-x0)))*u[ry+1][rx+0] +
//                                                   (((px-x0)*(px-x1))/((x2-x0)*(x2-x1)))*u[ry+1][rx+1] );

    double x0, x1, x2, x3; x0 = (rx-1)*_hx; x1 = rx*_hx; x2 = (rx+1)*_hx; x3 = (rx+2)*_hx;
    double y0, y1, y2, y3; y0 = (ry-1)*_hy; y1 = ry*_hy; y2 = (ry+1)*_hy; y3 = (ry+2)*_hy;

    double Lx0 = (px-x1)*(px-x2)*(px-x3); double L00 = (x0-x1)*(x0-x2)*(x0-x3);
    double Lx1 = (px-x0)*(px-x2)*(px-x3); double L11 = (x1-x0)*(x1-x2)*(x1-x3);
    double Lx2 = (px-x0)*(px-x1)*(px-x3); double L22 = (x2-x0)*(x2-x1)*(x2-x3);
    double Lx3 = (px-x0)*(px-x1)*(px-x2); double L33 = (x3-x0)*(x3-x1)*(x3-x2);
    double Lx0x = (px-x1)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x1);
    double Lx1x = (px-x0)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x0);
    double Lx2x = (px-x0)*(px-x1)+(px-x1)*(px-x3)+(px-x3)*(px-x0);
    double Lx3x = (px-x0)*(px-x1)+(px-x1)*(px-x2)+(px-x2)*(px-x0);

    double Ry0 = (py-y1)*(py-y2)*(py-y3); double R00 = (y0-y1)*(y0-y2)*(y0-y3);
    double Ry1 = (py-y0)*(py-y2)*(py-y3); double R11 = (y1-y0)*(y1-y2)*(y1-y3);
    double Ry2 = (py-y0)*(py-y1)*(py-y3); double R22 = (y2-y0)*(y2-y1)*(y2-y3);
    double Ry3 = (py-y0)*(py-y1)*(py-y2); double R33 = (y3-y0)*(y3-y1)*(y3-y2);

    double Ry0y = (py-y1)*(py-y2)+(py-y2)*(py-y3)+(py-y3)*(py-y1);
    double Ry1y = (py-y0)*(py-y2)+(py-y2)*(py-y3)+(py-y3)*(py-y0);
    double Ry2y = (py-y0)*(py-y1)+(py-y1)*(py-y3)+(py-y3)*(py-y0);
    double Ry3y = (py-y0)*(py-y1)+(py-y1)*(py-y2)+(py-y2)*(py-y0);

    dx = (Lx0x/L00) * ( (Ry0/R00)*u[ry-1][rx-1] + (Ry1/R11)*u[ry+0][rx-1] + (Ry2/R22)*u[ry+1][rx-1] + (Ry3/R33)*u[ry+2][rx-1] )
       + (Lx1x/L11) * ( (Ry0/R00)*u[ry-1][rx+0] + (Ry1/R11)*u[ry+0][rx+0] + (Ry2/R22)*u[ry+1][rx+0] + (Ry3/R33)*u[ry+2][rx+0] )
       + (Lx2x/L22) * ( (Ry0/R00)*u[ry-1][rx+1] + (Ry1/R11)*u[ry+0][rx+1] + (Ry2/R22)*u[ry+1][rx+1] + (Ry3/R33)*u[ry+2][rx+1] )
       + (Lx3x/L33) * ( (Ry0/R00)*u[ry-1][rx+2] + (Ry1/R11)*u[ry+0][rx+2] + (Ry2/R22)*u[ry+1][rx+2] + (Ry3/R33)*u[ry+2][rx+2] );

    dy = (Ry0y/R00) * ( (Lx0/L00)*u[ry-1][rx-1] + (Lx1/L11)*u[ry-1][rx+0] + (Lx2/L22)*u[ry-1][rx+1] + (Lx3/L33)*u[ry-1][rx+2] )
       + (Ry1y/R11) * ( (Lx0/L00)*u[ry+0][rx-1] + (Lx1/L11)*u[ry+0][rx+0] + (Lx2/L22)*u[ry+0][rx+1] + (Lx3/L33)*u[ry+0][rx+2] )
       + (Ry2y/R22) * ( (Lx0/L00)*u[ry+1][rx-1] + (Lx1/L11)*u[ry+1][rx+0] + (Lx2/L22)*u[ry+1][rx+1] + (Lx3/L33)*u[ry+1][rx+2] )
       + (Ry3y/R33) * ( (Lx0/L00)*u[ry+2][rx-1] + (Lx1/L11)*u[ry+2][rx+0] + (Lx2/L22)*u[ry+2][rx+1] + (Lx3/L33)*u[ry+2][rx+2] );

    return consentrateInPoint(u);
}

auto DeltaGrid2D::cleanGrid() -> void
{
    if (_M==0 || _N == 0) return;
    for (unsigned int m=0; m<=_M; m++) delete [] m_nodes[m];
    delete [] m_nodes;
    delete [] _rows;
    delete [] _cols;
    _N = 0;
    _M = 0;
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

auto DeltaGrid2D::p() const -> const SpacePoint& { return _p; }
auto DeltaGrid2D::p() -> SpacePoint& { return _p; }

auto DeltaGrid2D::minX() const -> unsigned int { return _minX; }
auto DeltaGrid2D::maxX() const -> unsigned int { return _maxX; }
auto DeltaGrid2D::minY() const -> unsigned int { return _minY; }
auto DeltaGrid2D::maxY() const -> unsigned int { return _maxY; }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DeltaGrid1D::DeltaGrid1D() : _N(0), _hx(0.0) {}

DeltaGrid1D::~DeltaGrid1D() { cleanGrid(); }

auto DeltaGrid1D::initGrid(unsigned int N, double hx) -> void
{
    this->_N = N;
    this->_hx = hx;

    m_nodes = new double[N+1];
    for (unsigned int n=0; n<=N; n++) m_nodes[n] = 0.0;
}

auto DeltaGrid1D::cleanGrid() -> void
{
    if (_N == 0) return;
    delete [] m_nodes;
    _N = 0;
}

auto DeltaGrid1D::distributeGauss(const SpacePoint& sp, unsigned sigmaXNum) -> void
{
    unsigned int kx = 3*sigmaXNum;
    double sigmaX = _hx*sigmaXNum;

    _rx = static_cast<unsigned int>( round(sp.x*_N) );

    if (_rx < kx || _rx > _N-kx) throw DeltaGridException(std::to_string(sp.x));

    _p = sp;
    _minX = _rx - kx; _maxX = _rx + kx;

    double sumX = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++) sumX += exp(-((n*_hx-sp.x)*(n*_hx-sp.x))/(2.0*sigmaX*sigmaX));
    sumX *= _hx;

    double sigma = sumX / sqrt(2.0*M_PI);
    double factor = 1.0 / (sqrt(2.0*M_PI)*sigma);

    SpaceNodePDE sn;
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        sn.i = static_cast<int>(n); sn.x = n*_hx;
        m_nodes[n] = factor*exp(-(((sn.x-sp.x)*(sn.x-sp.x))/(2.0*sigmaX*sigmaX)));
    }
}

auto DeltaGrid1D::distributeSigle(const SpacePoint& sp) -> void
{
    _p = sp;
    _rx = static_cast<unsigned int>( round(sp.x*_N) );
    _minX = _maxX = _rx;
    m_nodes[_rx] = 1.0/_hx;
}

auto DeltaGrid1D::distributeRect4(const SpacePoint&) -> void
{}

auto DeltaGrid1D::consentrateInPoint(const DoubleVector &u) const -> double
{
    double pu = 0.0;
    for (unsigned int n=minX(); n<=maxX(); n++)
    {
        pu += u[n] * weight(n) * _hx;
    }
    return pu;
}

auto DeltaGrid1D::consentrateInPoint(const DoubleVector &u, double &dx) const -> double
{
    const unsigned int rx = static_cast<unsigned int>(_rx);
    const double px = p().x;

    //    dx = (u[ry][rx+1] - u[ry][rx-1])/(2.0*_hx);
    //    dy = (u[ry+1][rx] - u[ry-1][rx])/(2.0*_hy);

    //    dx += (px-rx*_hx)*((u[ry][rx+1] - 2.0*u[ry][rx] + u[ry][rx-1])/(_hx*_hx));
    //    dy += (py-ry*_hy)*((u[ry+1][rx] - 2.0*u[ry][rx] + u[ry-1][rx])/(_hy*_hy));

    dx = (u[rx-2] - 8.0*u[rx-1] + 8.0*u[rx+1] - u[rx+2])/(12.0*_hx);
    dx += ((px-rx*_hx))*((-2.0*u[rx-2] + 32.0*u[rx-1] - 60.0*u[rx] + 32.0*u[rx+1] - 2.0*u[rx+2])/(24.0*_hx*_hx));
    return consentrateInPoint(u);
}

auto DeltaGrid1D::isCenter(const SpaceNodePDE &sn) const -> bool { return (sn.i == static_cast<int>(_rx)); }
auto DeltaGrid1D::isCenter(unsigned int n) const -> bool { return (n == _rx); }
auto DeltaGrid1D::isContains(const SpaceNodePDE &sn) const -> bool { return (sn.i >= static_cast<int>(_minX)) && (sn.i <= static_cast<int>(_maxX)); }
auto DeltaGrid1D::isContains(unsigned int n) const -> bool { return (n >= _minX) && (n <= _maxX); }
auto DeltaGrid1D::weight(const SpaceNodePDE &sn) const -> double { return m_nodes[sn.i]; }
auto DeltaGrid1D::weight(unsigned int n) const -> double { return m_nodes[n]; }
auto DeltaGrid1D::rx() const -> unsigned int { return _rx; }
auto DeltaGrid1D::p() const -> const SpacePoint& { return _p; }
auto DeltaGrid1D::p() -> SpacePoint& { return _p; }
auto DeltaGrid1D::minX() const -> unsigned int { return _minX; }
auto DeltaGrid1D::maxX() const -> unsigned int { return _maxX; }

