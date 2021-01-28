#include "deltagrid.h"

DeltaGrid2D::DeltaGrid2D() : _N(0), _M(0), _hx(0.0), _hy(0.0), _p(SpacePoint()),
    m_nodes(nullptr), m_der_x(nullptr), m_der_y(nullptr), _rx(0), _ry(0),
    _minX(0), _maxX(0), _minY(0), _maxY(0)  {}

DeltaGrid2D::DeltaGrid2D(const DeltaGrid2D &other)
{
    this->_N = other._N;
    this->_M = other._M;
    this->_hx = other._hx;
    this->_hy = other._hy;
    this->_p = other._p;
    this->m_nodes = other.m_nodes;
    this->m_der_x = other.m_der_x;
    this->m_der_y = other.m_der_y;
    this->_minX = other._minX;
    this->_maxX = other._maxX;
    this->_minY = other._minY;
    this->_maxY = other._maxY;
}

DeltaGrid2D& DeltaGrid2D::operator =(const DeltaGrid2D &other)
{
    if (this == &other) { return *this; }

    this->_N = other._N;
    this->_M = other._M;
    this->_hx = other._hx;
    this->_hy = other._hy;
    this->_p = other._p;
    this->m_nodes = other.m_nodes;
    this->m_der_x = other.m_der_x;
    this->m_der_y = other.m_der_y;
    this->_minX = other._minX;
    this->_maxX = other._maxX;
    this->_minY = other._minY;
    this->_maxY = other._maxY;
    return *this;
}

DeltaGrid2D::~DeltaGrid2D() { cleanGrid(); }

auto DeltaGrid2D::initGrid(unsigned int N, double hx, unsigned int M, double hy) -> void
{
    this->_N = N;
    this->_M = M;
    this->_hx = hx;
    this->_hy = hy;

    m_nodes = new double*[M+1];
    m_der_x = new double*[M+1];
    m_der_y = new double*[M+1];
    for (unsigned int m=0; m<=M; m++)
    {
        m_nodes[m] = new double[N+1];
        m_der_x[m] = new double[N+1];
        m_der_y[m] = new double[N+1];
        for (unsigned int n=0; n<=N; n++)
        {
            m_nodes[m][n] = m_der_x[m][n] = m_der_y[m][n] = 0.0;
        }
    }
}

auto DeltaGrid2D::initGrid(const Dimension &dimensionX, const Dimension &dimensionY) -> void
{
    initGrid(static_cast<unsigned int>(dimensionX.size()), dimensionX.step(),
             static_cast<unsigned int>(dimensionY.size()), dimensionY.step());
}

auto DeltaGrid2D::cleanGrid() -> void
{
    if (_M==0 || _N == 0) return;
    for (unsigned int m=0; m<=_M; m++)
    {
        delete [] m_nodes[m];
        delete [] m_der_x[m];
        delete [] m_der_y[m];
    }
    delete [] m_nodes;
    delete [] m_der_x;
    delete [] m_der_y;
    _N = 0;
    _hx = 0.0;
    _M = 0;
    _hy = 0.0;
    _p.x = _p.y = 0.0;
}

auto DeltaGrid2D::reset() -> void
{
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        for (unsigned int n=_minX; n<=_maxX; n++)
        {
            m_nodes[m][n] = m_der_x[m][n] = m_der_y[m][n] = 0.0;
        }
    }
}

auto DeltaGrid2D::resetAll() -> void
{
    for (unsigned int m=0; m<=_M; m++)
    {
        for (unsigned int n=0; n<=_N; n++)
        {
            m_nodes[m][n] = m_der_x[m][n] = m_der_y[m][n] = 0.0;
        }
    }
}

auto DeltaGrid2D::gaussWeight(const SpacePoint &sp, const SpacePoint &mu, double sigmaX, double sigmaY) const -> double
{
    const double factor1 = 1.0/(2.0*M_PI*sigmaX*sigmaY);
    const double sigmax2 = 1.0/(2.0*sigmaX*sigmaX);
    const double sigmay2 = 1.0/(2.0*sigmaY*sigmaY);
    return factor1 * exp(-(sigmax2*(sp.x-mu.x)*(sp.x-mu.x)+sigmay2*(sp.y-mu.y)*(sp.y-mu.y)));
}

auto DeltaGrid2D::distributeGauss(const SpacePoint& sp, unsigned int nodeX_per_sigmaX, unsigned int nodeY_per_sigmaY) -> void
{
    reset();

    const unsigned int sigma_count = 4;
    const unsigned int kx = sigma_count * nodeX_per_sigmaX;
    const unsigned int ky = sigma_count * nodeY_per_sigmaY;
    const double sigmaX = _hx * nodeX_per_sigmaX;
    const double sigmaY = _hy * nodeY_per_sigmaY;

    _rx = static_cast<unsigned int>( round(sp.x*_N) );
    _ry = static_cast<unsigned int>( round(sp.y*_M) );

    if (_rx < kx || _ry < ky || _rx > _N-kx || _ry > _M-ky)
        throw delta_grid_exception(std::string("Point:["+std::to_string(sp.x)+","+std::to_string(sp.y)+"]"));

    _p = sp;

    _minX = _rx - kx; _maxX = _rx + kx;
    _minY = _ry - ky; _maxY = _ry + ky;

    double sumX = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        double k = 1.0;
        //if (n==_minX || n==_maxX) k = 0.5;
        sumX += k * exp(-((n*_hx-sp.x)*(n*_hx-sp.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= _hx;

    double sumY = 0.0;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        double k = 1.0;
        //if (m==_minY || m==_maxY) k = 0.5;
        sumY += k * exp(-((m*_hy-sp.y)*(m*_hy-sp.y))/(2.0*sigmaY*sigmaY));
    }
    sumY *= _hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/(2.0*M_PI*sigma);
    //factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);

    SpaceNodePDE sn;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*_hy;
        for (unsigned int n=_minX; n<=_maxX; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*_hx;
            m_nodes[m][n] = factor*exp(-(((sn.x-sp.x)*(sn.x-sp.x))/(2.0*sigmaX*sigmaX) +
                                         ((sn.y-sp.y)*(sn.y-sp.y))/(2.0*sigmaY*sigmaY)));
            m_der_x[m][n] =  m_nodes[m][n] * (-(sn.x-sp.x)/(sigmaX*sigmaX));
            m_der_y[m][n] =  m_nodes[m][n] * (-(sn.y-sp.y)/(sigmaY*sigmaY));
        }
    }
}

auto DeltaGrid2D::distributeGauss1(const SpacePoint& sp, unsigned int nodeX_per_sigmaX, unsigned int nodeY_per_sigmaY) -> void
{
    reset();

    const unsigned int sigma_count = 4;
    const unsigned int kx = sigma_count * nodeX_per_sigmaX;
    const unsigned int ky = sigma_count * nodeY_per_sigmaY;
    const double sigmaX = _hx * nodeX_per_sigmaX;
    const double sigmaY = _hy * nodeY_per_sigmaY;

    _rx = static_cast<unsigned int>( round(sp.x*_N) );
    _ry = static_cast<unsigned int>( round(sp.y*_M) );

    if (_rx < kx || _ry < ky || _rx > _N-kx || _ry > _M-ky)
        throw delta_grid_exception(std::string("Point:["+std::to_string(sp.x)+","+std::to_string(sp.y)+"]"));

    _p = sp;

    _minX = _rx - kx; _maxX = _rx + kx;
    _minY = _ry - ky; _maxY = _ry + ky;

    double sumX = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        double k = 1.0;
        //if (n==_minX || n==_maxX) k = 0.5;
        sumX += k * exp(-((n*_hx-sp.x)*(n*_hx-sp.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= _hx;

    double sumY = 0.0;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        double k = 1.0;
        //if (m==_minY || m==_maxY) k = 0.5;
        sumY += k * exp(-((m*_hy-sp.y)*(m*_hy-sp.y))/(2.0*sigmaY*sigmaY));
    }
    sumY *= _hy;

    //double sigma = (sumX*sumY) / (2.0*M_PI);
    //double factor = 1.0/(2.0*M_PI*sigma);
    double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);

    SpaceNodePDE sn;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*_hy;
        for (unsigned int n=_minX; n<=_maxX; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*_hx;
            m_nodes[m][n] = factor*exp(-(((sn.x-sp.x)*(sn.x-sp.x))/(2.0*sigmaX*sigmaX) +
                                         ((sn.y-sp.y)*(sn.y-sp.y))/(2.0*sigmaY*sigmaY)));
            m_der_x[m][n] =  m_nodes[m][n] * (-(sn.x-sp.x)/(sigmaX*sigmaX));
            m_der_y[m][n] =  m_nodes[m][n] * (-(sn.y-sp.y)/(sigmaY*sigmaY));
        }
    }
}

auto DeltaGrid2D::lumpPointGauss(const DoubleMatrix &u) const -> double
{
    double pu = 0.0;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        for (unsigned int n=_minX; n<=_maxX; n++)
        {
            pu += m_nodes[m][n] * _hx * _hy * u[m][n];
        }
    }
    return pu;
}

auto DeltaGrid2D::lumpPointGauss(const DoubleMatrix &u, double &ux, double &uy) const -> double
{
    double pu = 0.0;
    ux = 0.0;
    uy = 0.0;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        for (unsigned int n=_minX; n<=_maxX; n++)
        {
            //double k = 1.0;
            //if (m == _minY || m == _maxY) k *= 0.5;
            //if (n == _minX || n == _maxX) k *= 0.5;

            double _w = m_nodes[m][n];
            pu += _w * _hx * _hy * u[m][n];
            ux += _w * _hx * _hy * (u[m][n+1]-u[m][n-1])/(2.0*_hx);
            uy += _w * _hx * _hy * (u[m+1][n]-u[m-1][n])/(2.0*_hy);
        }
    }
    return pu;
}

auto DeltaGrid2D::lumpPointGauss1(const DoubleMatrix &u, double &ux, double &uy) const -> double
{
    double pu = 0.0;
    ux = 0.0;
    uy = 0.0;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        for (unsigned int n=_minX; n<=_maxX; n++)
        {
            pu += _hx * _hy * u[m][n] * m_nodes[m][n];
            ux -= _hx * _hy * u[m][n] * m_der_x[m][n];
            uy -= _hx * _hy * u[m][n] * m_der_y[m][n];
        }
    }
    return pu;
}

auto DeltaGrid2D::distributeSigle(const SpacePoint &sp) -> void
{
    _rx = static_cast<unsigned int>( round(sp.x*_N) );
    _ry = static_cast<unsigned int>( round(sp.y*_M) );

    _p = sp;

    _minX = _maxX = _rx;
    _minY = _maxY = _ry;

    m_nodes[_ry][_rx] = 1.0/(_hx*_hy);
}

auto DeltaGrid2D::distributeRect4(const SpacePoint &) -> void
{}

auto DeltaGrid2D::consentrateInPoint(const DoubleMatrix &u, unsigned int v) const -> double
{
    double pu = 0.0;

    if (v==3)
    {
        const double px = _p.x;
        const double py = _p.y;
        const unsigned int rx = static_cast<unsigned int>(floor(px/_hx));
        const unsigned int ry = static_cast<unsigned int>(floor(py/_hy));
        const double x0 = static_cast<double>(rx+0)*_hx;
        const double x1 = static_cast<double>(rx+1)*_hx;
        const double y0 = static_cast<double>(ry+0)*_hy;
        const double y1 = static_cast<double>(ry+1)*_hy;

        const double Lx0 = (px-x1); const double Lx00 = (x0-x1);
        const double Lx1 = (px-x0); const double Lx11 = (x1-x0);

        const double Ly0 = (py-y1); const double Ly00 = (y0-y1);
        const double Ly1 = (py-y0); const double Ly11 = (y1-y0);

        pu = 0.0;
        pu += (Ly0/Ly00) * ( (Lx0/Lx00)*u[ry+0][rx+0] + (Lx1/Lx11)*u[ry+0][rx+1] );
        pu += (Ly1/Ly11) * ( (Lx0/Lx00)*u[ry+1][rx+0] + (Lx1/Lx11)*u[ry+1][rx+1] );
        return pu;
    }

    if (v==4)
    {
        const double px = _p.x;
        const double py = _p.y;
        const unsigned int rx = static_cast<unsigned int>(round(px/_hx));
        const unsigned int ry = static_cast<unsigned int>(round(py/_hy));
        const double x0 = static_cast<double>(rx-1)*_hx;
        const double x1 = static_cast<double>(rx+0)*_hx;
        const double x2 = static_cast<double>(rx+1)*_hx;
        const double y0 = static_cast<double>(ry-1)*_hy;
        const double y1 = static_cast<double>(ry+0)*_hy;
        const double y2 = static_cast<double>(ry+1)*_hy;

        const double Lx0 = (px-x1)*(px-x2); const double Lx00 = (x0-x1)*(x0-x2);
        const double Lx1 = (px-x0)*(px-x2); const double Lx11 = (x1-x0)*(x1-x2);
        const double Lx2 = (px-x0)*(px-x1); const double Lx22 = (x2-x0)*(x2-x1);

        const double Ly0 = (py-y1)*(py-y2); const double Ly00 = (y0-y1)*(y0-y2);
        const double Ly1 = (py-y0)*(py-y2); const double Ly11 = (y1-y0)*(y1-y2);
        const double Ly2 = (py-y0)*(py-y1); const double Ly22 = (y2-y0)*(y2-y1);

        pu = 0.0;
        pu += (Ly0/Ly00) * ( (Lx0/Lx00)*u[ry-1][rx-1] + (Lx1/Lx11)*u[ry-1][rx+0] + (Lx2/Lx22)*u[ry-1][rx+1] );
        pu += (Ly1/Ly11) * ( (Lx0/Lx00)*u[ry+0][rx-1] + (Lx1/Lx11)*u[ry+0][rx+0] + (Lx2/Lx22)*u[ry+0][rx+1] );
        pu += (Ly2/Ly22) * ( (Lx0/Lx00)*u[ry+1][rx-1] + (Lx1/Lx11)*u[ry+1][rx+0] + (Lx2/Lx22)*u[ry+1][rx+1] );
        return pu;
    }

    if (v==5)
    {
        const double px = _p.x;
        const double py = _p.y;
        const unsigned int rx = static_cast<unsigned int>(floor(px/_hx));
        const unsigned int ry = static_cast<unsigned int>(floor(py/_hy));
        const double x0 = static_cast<double>(rx-1)*_hx;
        const double x1 = static_cast<double>(rx+0)*_hx;
        const double x2 = static_cast<double>(rx+1)*_hx;
        const double x3 = static_cast<double>(rx+2)*_hx;
        const double y0 = static_cast<double>(ry-1)*_hy;
        const double y1 = static_cast<double>(ry+0)*_hy;
        const double y2 = static_cast<double>(ry+1)*_hy;
        const double y3 = static_cast<double>(ry+2)*_hy;

        const double Lx0 = (px-x1)*(px-x2)*(px-x3); const double Lx00 = (x0-x1)*(x0-x2)*(x0-x3);
        const double Lx1 = (px-x0)*(px-x2)*(px-x3); const double Lx11 = (x1-x0)*(x1-x2)*(x1-x3);
        const double Lx2 = (px-x0)*(px-x1)*(px-x3); const double Lx22 = (x2-x0)*(x2-x1)*(x2-x3);
        const double Lx3 = (px-x0)*(px-x1)*(px-x2); const double Lx33 = (x3-x0)*(x3-x1)*(x3-x2);

        const double Ly0 = (py-y1)*(py-y2)*(py-y3); const double Ly00 = (y0-y1)*(y0-y2)*(y0-y3);
        const double Ly1 = (py-y0)*(py-y2)*(py-y3); const double Ly11 = (y1-y0)*(y1-y2)*(y1-y3);
        const double Ly2 = (py-y0)*(py-y1)*(py-y3); const double Ly22 = (y2-y0)*(y2-y1)*(y2-y3);
        const double Ly3 = (py-y0)*(py-y1)*(py-y2); const double Ly33 = (y3-y0)*(y3-y1)*(y3-y2);

        pu = 0.0;
        pu += (Ly0/Ly00) * ( (Lx0/Lx00)*u[ry-1][rx-1] + (Lx1/Lx11)*u[ry-1][rx+0] + (Lx2/Lx22)*u[ry-1][rx+1] + (Lx3/Lx33)*u[ry-1][rx+2] );
        pu += (Ly1/Ly11) * ( (Lx0/Lx00)*u[ry+0][rx-1] + (Lx1/Lx11)*u[ry+0][rx+0] + (Lx2/Lx22)*u[ry+0][rx+1] + (Lx3/Lx33)*u[ry+0][rx+2] );
        pu += (Ly2/Ly22) * ( (Lx0/Lx00)*u[ry+1][rx-1] + (Lx1/Lx11)*u[ry+1][rx+0] + (Lx2/Lx22)*u[ry+1][rx+1] + (Lx3/Lx33)*u[ry+1][rx+2] );
        pu += (Ly3/Ly33) * ( (Lx0/Lx00)*u[ry+2][rx-1] + (Lx1/Lx11)*u[ry+2][rx+0] + (Lx2/Lx22)*u[ry+2][rx+1] + (Lx3/Lx33)*u[ry+2][rx+2] );
        return pu;
    }

    return pu;
}

auto DeltaGrid2D::derivativesInPoint(const DoubleMatrix &u, double &dx, double &dy, unsigned int v) const -> void
{
    if (v == 0)
    {
        const double px = _p.x;
        const double py = _p.y;
        const unsigned int rx = static_cast<unsigned int>(round(px/_hx));
        const unsigned int ry = static_cast<unsigned int>(round(py/_hy));
        dx = (u[ry][rx+1]-u[ry][rx-1])/(2.0*_hx);
        dy = (u[ry+1][rx]-u[ry-1][rx])/(2.0*_hy);
    }

    if (v == 3)
    {
        const double px = _p.x;
        const double py = _p.y;
        const unsigned int rx = static_cast<unsigned int>(floor(px/_hx));
        const unsigned int ry = static_cast<unsigned int>(floor(py/_hy));

        const double x0 = static_cast<double>(rx+0)*_hx;
        const double x1 = static_cast<double>(rx+1)*_hx;
        const double y0 = static_cast<double>(ry+0)*_hy;
        const double y1 = static_cast<double>(ry+1)*_hy;

        const double Lx0 = (px-x1); const double Lx00 = (x0-x1); const double Lx0x = 1.0;
        const double Lx1 = (px-x0); const double Lx11 = (x1-x0); const double Lx1x = 1.0;

        const double Ly0 = (py-y1); const double Ly00 = (y0-y1); const double Ly0y = 1.0;
        const double Ly1 = (py-y0); const double Ly11 = (y1-y0); const double Ly1y = 1.0;

        dx = 0.0;
        dx += (Lx0x/Lx00) * ( (Ly0/Ly00)*u[ry+0][rx+0] + (Ly1/Ly11)*u[ry+1][rx+0] );
        dx += (Lx1x/Lx11) * ( (Ly0/Ly00)*u[ry+0][rx+1] + (Ly1/Ly11)*u[ry+1][rx+1] );

        dy = 0.0;
        dy += (Ly0y/Ly00) * ( (Lx0/Lx00)*u[ry+0][rx+0] + (Lx1/Lx11)*u[ry+0][rx+1] );
        dy += (Ly1y/Ly11) * ( (Lx0/Lx00)*u[ry+1][rx+0] + (Lx1/Lx11)*u[ry+1][rx+1] );
    }

    if (v == 4)
    {
        const double px = _p.x;
        const double py = _p.y;
        const unsigned int rx = static_cast<unsigned int>(round(px/_hx));
        const unsigned int ry = static_cast<unsigned int>(round(py/_hy));

        const double x0 = static_cast<double>(rx-1)*_hx;
        const double x1 = static_cast<double>(rx+0)*_hx;
        const double x2 = static_cast<double>(rx+1)*_hx;
        const double y0 = static_cast<double>(ry-1)*_hy;
        const double y1 = static_cast<double>(ry+0)*_hy;
        const double y2 = static_cast<double>(ry+1)*_hy;

        const double Lx0 = (px-x1)*(px-x2); const double Lx00 = (x0-x1)*(x0-x2); const double Lx0x = (px-x1)+(px-x2);
        const double Lx1 = (px-x0)*(px-x2); const double Lx11 = (x1-x0)*(x1-x2); const double Lx1x = (px-x0)+(px-x2);
        const double Lx2 = (px-x0)*(px-x1); const double Lx22 = (x2-x0)*(x2-x1); const double Lx2x = (px-x0)+(px-x1);

        const double Ly0 = (py-y1)*(py-y2); const double Ly00 = (y0-y1)*(y0-y2); const double Ly0y = (py-y1)+(py-y2);
        const double Ly1 = (py-y0)*(py-y2); const double Ly11 = (y1-y0)*(y1-y2); const double Ly1y = (py-y0)+(py-y2);
        const double Ly2 = (py-y0)*(py-y1); const double Ly22 = (y2-y0)*(y2-y1); const double Ly2y = (py-y0)+(py-y1);

        dx = 0.0;
        dx += (Lx0x/Lx00) * ( (Ly0/Ly00)*u[ry-1][rx-1] + (Ly1/Ly11)*u[ry+0][rx-1] + (Ly2/Ly22)*u[ry+1][rx-1] );
        dx += (Lx1x/Lx11) * ( (Ly0/Ly00)*u[ry-1][rx+0] + (Ly1/Ly11)*u[ry+0][rx+0] + (Ly2/Ly22)*u[ry+1][rx+0] );
        dx += (Lx2x/Lx22) * ( (Ly0/Ly00)*u[ry-1][rx+1] + (Ly1/Ly11)*u[ry+0][rx+1] + (Ly2/Ly22)*u[ry+1][rx+1] );

        dy = 0.0;
        dy += (Ly0y/Ly00) * ( (Lx0/Lx00)*u[ry-1][rx-1] + (Lx1/Lx11)*u[ry-1][rx+0] + (Lx2/Lx22)*u[ry-1][rx+1] );
        dy += (Ly1y/Ly11) * ( (Lx0/Lx00)*u[ry+0][rx-1] + (Lx1/Lx11)*u[ry+0][rx+0] + (Lx2/Lx22)*u[ry+0][rx+1] );
        dy += (Ly2y/Ly22) * ( (Lx0/Lx00)*u[ry+1][rx-1] + (Lx1/Lx11)*u[ry+1][rx+0] + (Lx2/Lx22)*u[ry+1][rx+1] );
    }

    if (v == 5)
    {
        const double px = _p.x;
        const double py = _p.y;
        const unsigned int rx = static_cast<unsigned int>(floor(px/_hx));
        const unsigned int ry = static_cast<unsigned int>(floor(py/_hy));

        const double x0 = static_cast<double>(rx-1)*_hx;
        const double x1 = static_cast<double>(rx+0)*_hx;
        const double x2 = static_cast<double>(rx+1)*_hx;
        const double x3 = static_cast<double>(rx+2)*_hx;
        const double y0 = static_cast<double>(ry-1)*_hy;
        const double y1 = static_cast<double>(ry+0)*_hy;
        const double y2 = static_cast<double>(ry+1)*_hy;
        const double y3 = static_cast<double>(ry+2)*_hy;

        const double Lx0 = (px-x1)*(px-x2)*(px-x3); const double Lx00 = (x0-x1)*(x0-x2)*(x0-x3); const double Lx0x = (px-x1)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x1);
        const double Lx1 = (px-x0)*(px-x2)*(px-x3); const double Lx11 = (x1-x0)*(x1-x2)*(x1-x3); const double Lx1x = (px-x0)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x0);
        const double Lx2 = (px-x0)*(px-x1)*(px-x3); const double Lx22 = (x2-x0)*(x2-x1)*(x2-x3); const double Lx2x = (px-x0)*(px-x1)+(px-x1)*(px-x3)+(px-x3)*(px-x0);
        const double Lx3 = (px-x0)*(px-x1)*(px-x2); const double Lx33 = (x3-x0)*(x3-x1)*(x3-x2); const double Lx3x = (px-x0)*(px-x1)+(px-x1)*(px-x2)+(px-x2)*(px-x0);

        const double Ly0 = (py-y1)*(py-y2)*(py-y3); const double Ly00 = (y0-y1)*(y0-y2)*(y0-y3); const double Ly0y = (py-y1)*(py-y2)+(py-y2)*(py-y3)+(py-y3)*(py-y1);
        const double Ly1 = (py-y0)*(py-y2)*(py-y3); const double Ly11 = (y1-y0)*(y1-y2)*(y1-y3); const double Ly1y = (py-y0)*(py-y2)+(py-y2)*(py-y3)+(py-y3)*(py-y0);
        const double Ly2 = (py-y0)*(py-y1)*(py-y3); const double Ly22 = (y2-y0)*(y2-y1)*(y2-y3); const double Ly2y = (py-y0)*(py-y1)+(py-y1)*(py-y3)+(py-y3)*(py-y0);
        const double Ly3 = (py-y0)*(py-y1)*(py-y2); const double Ly33 = (y3-y0)*(y3-y1)*(y3-y2); const double Ly3y = (py-y0)*(py-y1)+(py-y1)*(py-y2)+(py-y2)*(py-y0);

        dx = 0.0;
        dx += (Lx0x/Lx00) * ( (Ly0/Ly00)*u[ry-1][rx-1] + (Ly1/Ly11)*u[ry+0][rx-1] + (Ly2/Ly22)*u[ry+1][rx-1] + (Ly3/Ly33)*u[ry+2][rx-1] );
        dx += (Lx1x/Lx11) * ( (Ly0/Ly00)*u[ry-1][rx+0] + (Ly1/Ly11)*u[ry+0][rx+0] + (Ly2/Ly22)*u[ry+1][rx+0] + (Ly3/Ly33)*u[ry+2][rx+0] );
        dx += (Lx2x/Lx22) * ( (Ly0/Ly00)*u[ry-1][rx+1] + (Ly1/Ly11)*u[ry+0][rx+1] + (Ly2/Ly22)*u[ry+1][rx+1] + (Ly3/Ly33)*u[ry+2][rx+1] );
        dx += (Lx3x/Lx33) * ( (Ly0/Ly00)*u[ry-1][rx+2] + (Ly1/Ly11)*u[ry+0][rx+2] + (Ly2/Ly22)*u[ry+1][rx+2] + (Ly3/Ly33)*u[ry+2][rx+2] );

        dy = 0.0;
        dy += (Ly0y/Ly00) * ( (Lx0/Lx00)*u[ry-1][rx-1] + (Lx1/Lx11)*u[ry-1][rx+0] + (Lx2/Lx22)*u[ry-1][rx+1] + (Lx3/Lx33)*u[ry-1][rx+2] );
        dy += (Ly1y/Ly11) * ( (Lx0/Lx00)*u[ry+0][rx-1] + (Lx1/Lx11)*u[ry+0][rx+0] + (Lx2/Lx22)*u[ry+0][rx+1] + (Lx3/Lx33)*u[ry+0][rx+2] );
        dy += (Ly2y/Ly22) * ( (Lx0/Lx00)*u[ry+1][rx-1] + (Lx1/Lx11)*u[ry+1][rx+0] + (Lx2/Lx22)*u[ry+1][rx+1] + (Lx3/Lx33)*u[ry+1][rx+2] );
        dy += (Ly3y/Ly33) * ( (Lx0/Lx00)*u[ry+2][rx-1] + (Lx1/Lx11)*u[ry+2][rx+0] + (Lx2/Lx22)*u[ry+2][rx+1] + (Lx3/Lx33)*u[ry+2][rx+2] );
    }
}

auto DeltaGrid2D::consentrateInPoint(const DoubleMatrix &u, double &dx, double &dy, unsigned int v) const -> double
{
    derivativesInPoint(u, dx, dy, v);
    return consentrateInPoint(u, v);
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
    return (n >= _minX) && (n <= _maxX) && (m >= _minY) && (m <= _maxY);
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

auto DeltaGrid2D::hx() const -> double { return _hx; }
auto DeltaGrid2D::hy() const -> double { return _hy; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DeltaGrid1D::DeltaGrid1D() : _N(0), _hx(0.0), _p(SpacePoint()),
    m_nodes(nullptr), m_der_x(nullptr), _rx(0), _minX(0), _maxX(0) {}

DeltaGrid1D::DeltaGrid1D(const DeltaGrid1D &other)
{
    this->_N = other._N;
    this->_hx = other._hx;
    this->_p = other._p;
    this->m_nodes = other.m_nodes;
    this->m_der_x = other.m_der_x;
    this->_minX = other._minX;
    this->_maxX = other._maxX;
}

DeltaGrid1D& DeltaGrid1D::operator =(const DeltaGrid1D &other)
{
    if (this == &other) { return *this; }

    this->_N = other._N;
    this->_hx = other._hx;
    this->_p = other._p;
    this->m_nodes = other.m_nodes;
    this->m_der_x = other.m_der_x;
    this->_minX = other._minX;
    this->_maxX = other._maxX;
    return *this;
}

DeltaGrid1D::~DeltaGrid1D() { cleanGrid(); }

auto DeltaGrid1D::initGrid(unsigned int N, double hx) -> void
{
    this->_N = N;
    this->_hx = hx;

    m_nodes = new double[N+1];
    m_der_x = new double[N+1];
    for (unsigned int n=0; n<=N; n++)
    {
        m_nodes[n] = m_der_x[n] = 0.0;
    }
}

auto DeltaGrid1D::initGrid(const Dimension &dimension) -> void
{
    initGrid(static_cast<unsigned int>(dimension.size()), dimension.step());
}

auto DeltaGrid1D::cleanGrid() -> void
{
    if (_N == 0) return;
    delete [] m_nodes;
    delete [] m_der_x;
    _N = 0;
    _hx = 0.0;
    _p.x = 0.0;
}

auto DeltaGrid1D::reset() -> void
{
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        m_nodes[n] = m_der_x[n] = 0.0;
    }
}

auto DeltaGrid1D::resetAll() -> void
{
    for (unsigned int n=0; n<=_N; n++)
    {
        m_nodes[n] = m_der_x[n] = 0.0;
    }
}

auto DeltaGrid1D::distributeGauss(const SpacePoint &sp, unsigned int nodeX_per_sigmaX) -> void
{
    reset();

    const unsigned int sigma_count = 4;
    const unsigned int kx = sigma_count * nodeX_per_sigmaX;
    const double sigmaX = _hx * nodeX_per_sigmaX;

    _rx = static_cast<unsigned int>( round(sp.x*_N) );

    if (_rx < kx || _rx > _N-kx)
        throw delta_grid_exception(std::string("Point:["+std::to_string(sp.x)+"]"));

    _p = sp;

    _minX = _rx - kx; _maxX = _rx + kx;

    double sumX = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        double k = 1.0;
        //if (n==_minX || n==_maxX) k = 0.5;
        sumX += k * exp(-((n*_hx-sp.x)*(n*_hx-sp.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= _hx;

    double sigma = sumX / sqrt(2.0*M_PI);
    double factor = 1.0 / (sqrt(2.0*M_PI)*sigma);
    //factor = 1.0/(sqrt(2.0*M_PI)*sigmaX);

    SpaceNodePDE sn;
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        sn.i = static_cast<int>(n); sn.x = n*_hx;
        m_nodes[n] = factor*exp(-(((sn.x-sp.x)*(sn.x-sp.x))/(2.0*sigmaX*sigmaX)));
        m_der_x[n] =  m_nodes[n] * (-(sn.x-sp.x)/(sigmaX*sigmaX));
    }
}

auto DeltaGrid1D::lumpPointGauss(const DoubleVector &u) const -> double
{
    double pu = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        pu += m_nodes[n] * _hx * u[n];
    }
    return pu;
}

auto DeltaGrid1D::lumpPointGauss(const DoubleVector &u, double &ux) const -> double
{
    double pu = 0.0;
    ux = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        //double k = 1.0;
        //if (n == _minX || n == _maxX) k *= 0.5;

        double _w = m_nodes[n];
        pu += _w * _hx * u[n];
        ux += _w * _hx * (u[n+1]-u[n-1])/(2.0*_hx);
    }
    return pu;
}

auto DeltaGrid1D::lumpPointGauss1(const DoubleVector &u, double &ux) const -> double
{
    double pu = 0.0;
    ux = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        pu += _hx * u[n] * m_nodes[n];
        ux -= _hx * u[n] * m_der_x[n];
    }
    return pu;
}

auto DeltaGrid1D::distributeSigle(const SpacePoint &sp) -> void
{
    _rx = static_cast<unsigned int>( round(sp.x*_N) );

    _p = sp;

    _minX = _maxX = _rx;

    m_nodes[_rx] = 1.0/_hx;
}

auto DeltaGrid1D::distributeRect4(const SpacePoint &) -> void
{}

auto DeltaGrid1D::consentrateInPoint(const DoubleVector &u, unsigned int v) const -> double
{
    double pu = 0.0;

    if (v==3)
    {
        const double px = _p.x;
        const unsigned int rx = static_cast<unsigned int>(floor(px/_hx));
        const double x0 = (rx+0)*_hx;
        const double x1 = (rx+1)*_hx;

        const double Lx0 = (px-x1); const double Lx00 = (x0-x1);
        const double Lx1 = (px-x0); const double Lx11 = (x1-x0);

        pu = (Lx0/Lx00)*u[rx+0] + (Lx1/Lx11)*u[rx+1];
        return pu;
    }

    if (v==4)
    {
        const double px = _p.x;
        const unsigned int rx = static_cast<unsigned int>(round(px/_hx));
        const double x0 = static_cast<double>(rx-1)*_hx;
        const double x1 = static_cast<double>(rx+0)*_hx;
        const double x2 = static_cast<double>(rx+1)*_hx;

        const double Lx0 = (px-x1)*(px-x2); const double Lx00 = (x0-x1)*(x0-x2);
        const double Lx1 = (px-x0)*(px-x2); const double Lx11 = (x1-x0)*(x1-x2);
        const double Lx2 = (px-x0)*(px-x1); const double Lx22 = (x2-x0)*(x2-x1);

        pu = 0.0;
        pu = (Lx0/Lx00)*u[rx-1] + (Lx1/Lx11)*u[rx+0] + (Lx2/Lx22)*u[rx+1];
        return pu;
    }

    if (v==5)
    {
        const double px = _p.x;
        const unsigned int rx = static_cast<unsigned int>(floor(px/_hx));
        const double x0 = static_cast<double>(rx-1)*_hx;
        const double x1 = static_cast<double>(rx+0)*_hx;
        const double x2 = static_cast<double>(rx+1)*_hx;
        const double x3 = static_cast<double>(rx+2)*_hx;

        const double Lx0 = (px-x1)*(px-x2)*(px-x3); const double Lx00 = (x0-x1)*(x0-x2)*(x0-x3);
        const double Lx1 = (px-x0)*(px-x2)*(px-x3); const double Lx11 = (x1-x0)*(x1-x2)*(x1-x3);
        const double Lx2 = (px-x0)*(px-x1)*(px-x3); const double Lx22 = (x2-x0)*(x2-x1)*(x2-x3);
        const double Lx3 = (px-x0)*(px-x1)*(px-x2); const double Lx33 = (x3-x0)*(x3-x1)*(x3-x2);

        pu = (Lx0/Lx00)*u[rx-1] + (Lx1/Lx11)*u[rx+0] + (Lx2/Lx22)*u[rx+1] + (Lx3/Lx33)*u[rx+2];
        return pu;
    }

    return pu;
}

auto DeltaGrid1D::derivativesInPoint(const DoubleVector &u, double &dx, unsigned int v) const -> void
{
    if (v == 0)
    {
        const double px = _p.x;
        const unsigned int rx = static_cast<unsigned int>(round(px/_hx));
        dx = (u[rx+1]-u[rx-1])/(2.0*_hx);
    }

    if (v == 3)
    {
        const double px = _p.x;
        const unsigned int rx = static_cast<unsigned int>(floor(px/_hx));

        const double x0 = static_cast<double>(rx+0)*_hx;
        const double x1 = static_cast<double>(rx+1)*_hx;

        const double Lx00 = (x0-x1); const double Lx0x = 1.0;
        const double Lx11 = (x1-x0); const double Lx1x = 1.0;

        dx = (Lx0x/Lx00)*u[rx+0] + (Lx1x/Lx11)*u[rx+1];
    }

    if (v == 4)
    {
        const double px = _p.x;
        const unsigned int rx = static_cast<unsigned int>(round(px/_hx));

        const double x0 = static_cast<double>(rx-1)*_hx;
        const double x1 = static_cast<double>(rx+0)*_hx;
        const double x2 = static_cast<double>(rx+1)*_hx;

        const double Lx00 = (x0-x1)*(x0-x2); const double Lx0x = (px-x1)+(px-x2);
        const double Lx11 = (x1-x0)*(x1-x2); const double Lx1x = (px-x0)+(px-x2);
        const double Lx22 = (x2-x0)*(x2-x1); const double Lx2x = (px-x0)+(px-x1);

        dx = (Lx0x/Lx00)*u[rx-1] + (Lx1x/Lx11)*u[rx+0] + (Lx2x/Lx22)*u[rx+1];
    }

    if (v == 5)
    {
        const double px = _p.x;
        const unsigned int rx = static_cast<unsigned int>(floor(px/_hx));

        const double x0 = static_cast<double>(rx-1)*_hx;
        const double x1 = static_cast<double>(rx+0)*_hx;
        const double x2 = static_cast<double>(rx+1)*_hx;
        const double x3 = static_cast<double>(rx+2)*_hx;

        const double Lx00 = (x0-x1)*(x0-x2)*(x0-x3); const double Lx0x = (px-x1)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x1);
        const double Lx11 = (x1-x0)*(x1-x2)*(x1-x3); const double Lx1x = (px-x0)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x0);
        const double Lx22 = (x2-x0)*(x2-x1)*(x2-x3); const double Lx2x = (px-x0)*(px-x1)+(px-x1)*(px-x3)+(px-x3)*(px-x0);
        const double Lx33 = (x3-x0)*(x3-x1)*(x3-x2); const double Lx3x = (px-x0)*(px-x1)+(px-x1)*(px-x2)+(px-x2)*(px-x0);

        dx = (Lx0x/Lx00)*u[rx-1] + (Lx1x/Lx11)*u[rx+0] + (Lx2x/Lx22)*u[rx+1] + (Lx3x/Lx33)*u[rx+2];
    }
}

auto DeltaGrid1D::consentrateInPoint(const DoubleVector &u, double &dx, unsigned int v) const -> double
{
    derivativesInPoint(u, dx, v);
    return consentrateInPoint(u, v);
}

auto DeltaGrid1D::isCenter(const SpaceNodePDE &sn) const -> bool
{
    return (sn.i == static_cast<int>(_rx));
}

auto DeltaGrid1D::isCenter(unsigned int n) const -> bool
{
    return (n == _rx);
}

auto DeltaGrid1D::isContains(const SpaceNodePDE &sn) const -> bool
{
    return (sn.i >= static_cast<int>(_minX)) && (sn.i <= static_cast<int>(_maxX));
}

auto DeltaGrid1D::isContains(unsigned int n) const -> bool
{
    return (n >= _minX) && (n <= _maxX);
}

auto DeltaGrid1D::weight(const SpaceNodePDE &sn) const -> double
{
    return m_nodes[sn.i];
}

auto DeltaGrid1D::weight(unsigned int n) const -> double
{
    return m_nodes[n];
}

auto DeltaGrid1D::rx() const -> unsigned int { return _rx; }

auto DeltaGrid1D::p() const -> const SpacePoint& { return _p; }
auto DeltaGrid1D::p() -> SpacePoint& { return _p; }

auto DeltaGrid1D::minX() const -> unsigned int { return _minX; }
auto DeltaGrid1D::maxX() const -> unsigned int { return _maxX; }

auto DeltaGrid1D::hx() const -> double { return _hx; }

auto DeltaGrid1D::onFlyWeight(const SpaceNodePDE &sn, const SpacePoint &mu, size_t n) const -> double
{
    //const size_t N = _N;
    const double h = _hx;
    const double sigma = h*n;
    const size_t k = 20;
    const size_t kn1 = k*n;
    const size_t kn2 = k*n*2;
    const double lmb = k*n*h;

    const double a = 0.0*h;
    const double b = _N*h;
    const double m = mu.x;

    if (fabs(a-m) >= lmb && fabs(b-m) >= lmb)
    {
        double x = mu.x-lmb, sum = 0.0;
        sum += 0.5*exp(-(((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma)));
        for (size_t i=1; i<kn2; i++)
        {
            x += h;
            sum += exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
        }
        x += h;
        sum += 0.5*exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));

        double factor = 1.0/(sum*h);

        return factor*exp(-((sn.x-mu.x)*(sn.x-mu.x))/(2.0*sigma*sigma));
    }

    if (fabs(a-m) < lmb)
    {
        double x = mu.x, sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;

        size_t s = static_cast<size_t>(floor(fabs(a-m)/h));
        double d = fabs(a-m)-s*h;
        double e = a+d;

        if (d>=0.0)
        {
            sum0 += exp(-((a-mu.x)*(a-mu.x))/(2.0*sigma*sigma));
            sum0 += exp(-((e-mu.x)*(e-mu.x))/(2.0*sigma*sigma));
            sum0 *= 0.5*d;
        }

        if (s>0)
        {
            x = mu.x;

            sum1 += 0.5*exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
            for (size_t i=1; i<s; i++)
            {
                x -= h;
                sum1 += exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
            }
            x -= h;
            sum1 += 0.5*exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
        }

        x = mu.x;
        sum2 += 0.5*exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
        for (size_t i=1; i<kn1; i++)
        {
            x += h;
            sum2 += exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
        }
        x += h;
        sum2 += 0.5*exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
        //printf("x: %f\n", x);

        //printf(">> %f %f %f %d %f %f %14.10f %14.10f %14.10f\n", m, lmb, sn.x, s, d, e, sum0, sum1, sum2);

        double factor = 1.0/((sum1*h)+(sum2*h)+sum0);

        return factor*exp(-((sn.x-mu.x)*(sn.x-mu.x))/(2.0*sigma*sigma));
    }

    if (fabs(b-m) < lmb)
    {
        double x = mu.x, sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;

        size_t s = static_cast<size_t>(fabs(b-m)/h);
        double d = fabs(b-m)-s*h;

        double e = b-d;
        if (d>=0.0)
        {
            sum0 += exp(-((b-mu.x)*(b-mu.x))/(2.0*sigma*sigma));
            sum0 += exp(-((e-mu.x)*(e-mu.x))/(2.0*sigma*sigma));
            sum0 *= 0.5*d;
        }

        if (s>0)
        {
            x = mu.x;

            sum1 += 0.5*exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
            for (size_t i=1; i<s; i++)
            {
                x += h;
                sum1 += exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
            }
            x += h;
            sum1 += 0.5*exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
        }

        x = mu.x;
        sum2 += 0.5*exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
        for (size_t i=1; i<kn1; i++)
        {
            x -= h;
            sum2 += exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));
        }
        x -= h;
        sum2 += 0.5*exp(-((x-mu.x)*(x-mu.x))/(2.0*sigma*sigma));

        //printf(">> %f %f %f %d %f %f %14.10f %14.10f %14.10f\n", m, lmb, sn.x, s, d, e, sum0, sum1, sum2);

        double factor = 1.0/((sum1*h)+(sum2*h)+sum0);

        return factor*exp(-((sn.x-mu.x)*(sn.x-mu.x))/(2.0*sigma*sigma));
    }

    return 0.0;
}

auto DeltaGrid2D::onFlyWeight(const SpacePoint &sp, const SpacePoint &mu, size_t n) const -> double
{
    const double hx = _hx;
    const double hy = _hy;
    const size_t Nx = _N;
    const size_t Ny = _M;
    const double sigmaX = hx*n;
    const double sigmaY = hy*n;
    const size_t k = 2;
    const double lmbX = k*sigmaX;
    const double lmbY = k*sigmaY;

    const double minX = 0.0*hx;
    const double maxX = Nx*hx;
    const double minY = 0.0*hx;
    const double maxY = Ny*hx;
    const double muX = mu.x;
    const double muY = mu.y;
    const size_t kn2 = k*n*2;

    if (fabs(minX-muX) >= lmbX && fabs(maxX-muX) >= lmbX && fabs(minY-muY) >= lmbY && fabs(maxY-muY) >= lmbY)
    {
        double sum = 0.0;
        double x = muX-lmbX;
        double y = muY-lmbY;

        sum += 0.25*exp(-(((lmbX)*(lmbX))/(2.0*sigmaX*sigmaX)+((lmbY)*(lmbY))/(2.0*sigmaY*sigmaY)));
        sum += 0.25*exp(-(((lmbX)*(lmbX))/(2.0*sigmaX*sigmaX)+((lmbY)*(lmbY))/(2.0*sigmaY*sigmaY)));
        sum += 0.25*exp(-(((lmbX)*(lmbX))/(2.0*sigmaX*sigmaX)+((lmbY)*(lmbY))/(2.0*sigmaY*sigmaY)));
        sum += 0.25*exp(-(((lmbX)*(lmbX))/(2.0*sigmaX*sigmaX)+((lmbY)*(lmbY))/(2.0*sigmaY*sigmaY)));

        x = muX-lmbX;
        for (unsigned int i=1; i<kn2; i++)
        {
            x += hx;
            sum += 0.5*exp(-(((x-muX)*(x-muX))/(2.0*sigmaX*sigmaX)+((lmbY)*(lmbY))/(2.0*sigmaY*sigmaY)));
            sum += 0.5*exp(-(((x-muX)*(x-muX))/(2.0*sigmaX*sigmaX)+((lmbY)*(lmbY))/(2.0*sigmaY*sigmaY)));
        }

        y = muY-lmbY;
        for (unsigned int j=1; j<kn2; j++)
        {
            y += hy;
            sum += 0.5*exp(-(((lmbX)*(lmbX))/(2.0*sigmaX*sigmaX)+((y-muY)*(y-muY))/(2.0*sigmaY*sigmaY)));
            sum += 0.5*exp(-(((lmbX)*(lmbX))/(2.0*sigmaX*sigmaX)+((y-muY)*(y-muY))/(2.0*sigmaY*sigmaY)));
        }

        x = muX-lmbX;
        for (unsigned int i=1; i<kn2; i++)
        {
            x += hx;
            y = muY-lmbY;
            for (unsigned int j=1; j<kn2; j++)
            {
                y += hy;
                sum += exp(-(((x-muX)*(x-muX))/(2.0*sigmaX*sigmaX)+((y-muY)*(y-muY))/(2.0*sigmaY*sigmaY)));
            }
        }

        double factor = 1.0/(sum*hx*hy);

        return factor*exp(-(((sp.x-muX)*(sp.x-muX))/(2.0*sigmaX*sigmaX)+((sp.y-muY)*(sp.y-muY))/(2.0*sigmaY*sigmaY)));
    }

    return 0.0;
}

auto DeltaFunction::nearest(const SpaceNodePDE &p, const SpacePoint &m, double hx, double hy, unsigned int Nx, unsigned int Ny) -> double
{
    const double lx = hx*Nx;
    const double ly = hy*Ny;
    const double px = m.x;
    const double py = m.y;
    const unsigned int rx = static_cast<unsigned int>(floor((px/hx)*lx));
    const unsigned int ry = static_cast<unsigned int>(floor((py/hy)*ly));

    double pu = 0.0;
    if (static_cast<unsigned int>(p.i) == rx && static_cast<unsigned int>(p.j) == ry) pu = 1.0/(hx*hy);
    return pu;
}

/**
 * @brief The normal (or Gaussian or Gauss or Laplace–Gauss) distribution
 * @param p Searchin point
 * @param m The mean point paremeter or expectation of the distribution
 * @param sigma The standard deviation
 * @return The value of Gaussian distribution function on point of p
 */
double DeltaFunction::gaussian(double p, double m, double sigma)
{
    const double factor = 1.0/(sqrt(2.0*M_PI)*sigma);
    const double sigma2 = 1.0/(2.0*sigma*sigma);
    return factor * exp(-(sigma2*(p-m)*(p-m)));
}

/**
 * @brief The normal (or Gaussian or Gauss or Laplace–Gauss) distribution
 * @param p Searchin point
 * @param m The mean point paremeter or expectation of the distribution
 * @param sigma The standard deviation
 * @param dx The value of derivative by respect to x paremeter on point of p
 * @return The value of Gaussian distribution function on point of p
 */
double DeltaFunction::gaussian(double p, double m, double sigma, double &dx)
{
    const double factor = 1.0/(sqrt(2.0*M_PI)*sigma);
    const double sigma2 = 1.0/(2.0*sigma*sigma);
    double fx = factor * exp(-(sigma2*(p-m)*(p-m)));
    dx = -fx*(p-m)/(sigma*sigma);
    return fx;
}

/**
 * @brief The multivariate normal (or Gaussian or Gauss or Laplace–Gauss) distribution
 * @param p Searchin point
 * @param m The mean point paremeter or expectation of the distribution : (x,y)
 * @param sigma The standard deviation : (x,y)
 * @return The value of Gaussian distribution function on point of p
 */
double DeltaFunction::gaussian(const SpacePoint &p, const SpacePoint &m, const SpacePoint &sigma)
{
    const double factor1 = 1.0/(2.0*M_PI*sigma.x*sigma.y);
    const double sigmax2 = 1.0/(2.0*sigma.x*sigma.x);
    const double sigmay2 = 1.0/(2.0*sigma.y*sigma.y);
    return factor1 * exp(-(sigmax2*(p.x-m.x)*(p.x-m.x)+sigmay2*(p.y-m.y)*(p.y-m.y)));
}

/**
 * @brief The multivariate normal (or Gaussian or Gauss or Laplace–Gauss) distribution
 * @param p Searchin point
 * @param m The mean point paremeter or expectation of the distribution : (x,y)
 * @param sigma The standard deviation : (x,y)
 * @param dx The value of partial derivative by respect to x paremeter on point of p
 * @param dy The value of partial derivative by respect to y paremeter on point of p
 * @return The value of Gaussian distribution function on point of p
 */
double DeltaFunction::gaussian(const SpacePoint &p, const SpacePoint &m, const SpacePoint &sigma, double &dx, double &dy)
{
    const double factor1 = 1.0/(2.0*M_PI*sigma.x*sigma.y);
    const double sigmax2 = 1.0/(2.0*sigma.x*sigma.x);
    const double sigmay2 = 1.0/(2.0*sigma.y*sigma.y);
    double fx = factor1 * exp(-(sigmax2*(p.x-m.x)*(p.x-m.x)+sigmay2*(p.y-m.y)*(p.y-m.y)));
    dx = -fx*(p.x-m.x)/(sigma.x*sigma.x);
    dy = -fx*(p.y-m.y)/(sigma.y*sigma.y);
    return fx;
}


double DeltaFunction::gaussian(double p, double m, double sigma, size_t k)
{
    const double a = m-k*sigma;
    const double b = m+k*sigma;
    const double sigma2 = 1.0/(2.0*sigma*sigma);
    const size_t N = 1000;

    struct F { static double fx(double x, double m, double b) { return exp(-((x-m)*(x-m))*b); } };

    double A = 0.0;
    {
        double sum = 0.0;
        double h1= fabs(b-a)/N;

        /** Trapezoidal rule **/
        //sum += 0.5*F::fx(a,m,sigma2);
        //for (size_t i=1; i<=999; i++)
        //{
        //   const double x = a+i*h1;
        //   sum += F::fx(x, m, sigma2);
        //}
        //sum += 0.5*F::fx(b,m,sigma2);
        //sum *= h1;

        /** Simpson's rule **/
        for (size_t i=1; i<N; i+=2)
        {
            double x = a+i*h1;
            sum += (F::fx(x-h1, m, sigma2) + 4.0*F::fx(x, m, sigma2) + F::fx(x+h1, m, sigma2));
        }
        sum *= (h1/3.0);

        A = 1.0/sum;
    }

    return A * F::fx(p, m, sigma2);
}

double DeltaFunction::gaussian(const SpacePoint &p, const SpacePoint &m, const SpacePoint &sigma, size_t k)
{
    const double ax = m.x-k*sigma.x;
    const double bx = m.x+k*sigma.x;
    const double ay = m.y-k*sigma.y;
    const double by = m.y+k*sigma.y;
    const double sigmaX2 = 1.0/(2.0*sigma.x*sigma.x);
    const double sigmaY2 = 1.0/(2.0*sigma.y*sigma.y);
    const size_t N = 1000;

    struct F { static inline double fx(double x, double y, double mx, double my, double c2x, double c2y)
        { return exp(-((x-mx)*(x-mx)*c2x + (y-my)*(y-my)*c2y)); } };

    double A = 0.0;
    {
        double sum = 0.0;
        double hx= fabs(bx-ax)/N;
        double hy= fabs(by-ay)/N;

        /** Trapezoidal rule **/
        sum += 0.25*F::fx(ax,ay,m.x,m.y,sigmaX2,sigmaY2);
        sum += 0.25*F::fx(ax,by,m.x,m.y,sigmaX2,sigmaY2);
        sum += 0.25*F::fx(bx,by,m.x,m.y,sigmaX2,sigmaY2);
        sum += 0.25*F::fx(bx,ay,m.x,m.y,sigmaX2,sigmaY2);

        for (size_t i=1; i<N; i++)
        {
            double x = ax+i*hx;
            sum += 0.5*F::fx(x,ay,m.x,m.y,sigmaX2,sigmaY2);
            sum += 0.5*F::fx(x,by,m.x,m.y,sigmaX2,sigmaY2);
        }
        for (size_t j=1; j<N; j++)
        {
            double y = ay+j*hy;
            sum += 0.5*F::fx(ax,y,m.x,m.y,sigmaX2,sigmaY2);
            sum += 0.5*F::fx(bx,y,m.x,m.y,sigmaX2,sigmaY2);
        }
        for (size_t i=1; i<N; i++)
        {
            double x = ax+i*hx;
            for (size_t j=1; j<N; j++)
            {
                double y = ay+j*hy;
                sum += F::fx(x,y,m.x,m.y,sigmaX2,sigmaY2);
            }
        }
        sum *= hx*hy;

        /** Simpson's rule **/
        //for (size_t i=0; i<=N; i++)
        //{
        //    double x = ax+i*hx;
        //    for (size_t j=0; j<=N; j++)
        //    {
        //        double y = ay+j*hy;
        //        double fx = F::fx(x,y,m.x,m.y,sigmaX2,sigmaY2);
        //        if (i==0 || i==N) { fx *= 1.0; } else if (i%2==1) { fx *= 4.0; } else { fx *= 2.0; }
        //        if (j==0 || j==N) { fx *= 1.0; } else if (j%2==1) { fx *= 4.0; } else { fx *= 2.0; }
        //        sum += fx;
        //    }
        //}
        //sum *= ((hx*hy)/9.0);

        A = 1.0/sum;
    }

    return A * F::fx(p.x, p.y, m.x, m.y, sigmaX2, sigmaY2);
}

double DeltaFunction::sinusoid(double /*p*/) { return 0.0; }

double DeltaFunction::sinusoid(const SpacePoint &/*p*/) { return 0.0; }

double DeltaFunction::lumpedPoint2(const DoubleMatrix &m, const SpacePoint &p, const Dimension &dimensionX, const Dimension &dimensionY)
{
    return DeltaFunction::lumpedPoint2(m, p, dimensionX.step(), dimensionY.step(), dimensionX.size()-1, dimensionY.size()-1);
}

double DeltaFunction::lumpedPoint2(const DoubleMatrix &m, const SpacePoint &p, double hx, double hy, unsigned int Nx, unsigned int Ny)
{
    const double lx = hx*Nx;
    const double ly = hy*Ny;
    const double px = p.x;
    const double py = p.y;
    const unsigned int rx = static_cast<unsigned int>(floor((px/hx)*lx));
    const unsigned int ry = static_cast<unsigned int>(floor((py/hy)*ly));
    const double x0 = static_cast<double>(rx+0)*hx;
    const double x1 = static_cast<double>(rx+1)*hx;
    const double y0 = static_cast<double>(ry+0)*hy;
    const double y1 = static_cast<double>(ry+1)*hy;

    const double Lx0 = (px-x1); const double Lx00 = (x0-x1);
    const double Lx1 = (px-x0); const double Lx11 = (x1-x0);

    const double Ly0 = (py-y1); const double Ly00 = (y0-y1);
    const double Ly1 = (py-y0); const double Ly11 = (y1-y0);

    double pu = 0.0;
    pu += (Ly0/Ly00) * ( (Lx0/Lx00)*m[ry+0][rx+0] + (Lx1/Lx11)*m[ry+0][rx+1] );
    pu += (Ly1/Ly11) * ( (Lx0/Lx00)*m[ry+1][rx+0] + (Lx1/Lx11)*m[ry+1][rx+1] );
    return pu;
}

double DeltaFunction::lumpedPoint3(const DoubleMatrix &m, const SpacePoint &p, const Dimension &dimensionX, const Dimension &dimensionY)
{
    return DeltaFunction::lumpedPoint3(m, p, dimensionX.step(), dimensionY.step(), dimensionX.size()-1, dimensionY.size()-1);
}

double DeltaFunction::lumpedPoint3(const DoubleMatrix &m, const SpacePoint &p, double hx, double hy, unsigned int Nx, unsigned int Ny)
{
    const double lx = hx*Nx;
    const double ly = hy*Ny;
    const double px = p.x;
    const double py = p.y;
    const unsigned int rx = static_cast<unsigned int>(round((px/hx)*lx));
    const unsigned int ry = static_cast<unsigned int>(round((py/hy)*ly));
    const double x0 = static_cast<double>(rx-1)*hx;
    const double x1 = static_cast<double>(rx+0)*hx;
    const double x2 = static_cast<double>(rx+1)*hx;
    const double y0 = static_cast<double>(ry-1)*hy;
    const double y1 = static_cast<double>(ry+0)*hy;
    const double y2 = static_cast<double>(ry+1)*hy;

    const double Lx0 = (px-x1)*(px-x2); const double Lx00 = (x0-x1)*(x0-x2);
    const double Lx1 = (px-x0)*(px-x2); const double Lx11 = (x1-x0)*(x1-x2);
    const double Lx2 = (px-x0)*(px-x1); const double Lx22 = (x2-x0)*(x2-x1);

    const double Ly0 = (py-y1)*(py-y2); const double Ly00 = (y0-y1)*(y0-y2);
    const double Ly1 = (py-y0)*(py-y2); const double Ly11 = (y1-y0)*(y1-y2);
    const double Ly2 = (py-y0)*(py-y1); const double Ly22 = (y2-y0)*(y2-y1);

    double pu = 0.0;
    pu += (Ly0/Ly00) * ( (Lx0/Lx00)*m[ry-1][rx-1] + (Lx1/Lx11)*m[ry-1][rx+0] + (Lx2/Lx22)*m[ry-1][rx+1] );
    pu += (Ly1/Ly11) * ( (Lx0/Lx00)*m[ry+0][rx-1] + (Lx1/Lx11)*m[ry+0][rx+0] + (Lx2/Lx22)*m[ry+0][rx+1] );
    pu += (Ly2/Ly22) * ( (Lx0/Lx00)*m[ry+1][rx-1] + (Lx1/Lx11)*m[ry+1][rx+0] + (Lx2/Lx22)*m[ry+1][rx+1] );
    return pu;
}

double DeltaFunction::lumpedPoint4(const DoubleMatrix &m, const SpacePoint &p, const Dimension &dimensionX, const Dimension &dimensionY)
{
    return DeltaFunction::lumpedPoint4(m, p, dimensionX.step(), dimensionY.step(), dimensionX.size()-1, dimensionY.size()-1);
}

double DeltaFunction::lumpedPoint4(const DoubleMatrix &m, const SpacePoint &p, double hx, double hy, unsigned int Nx, unsigned int Ny)
{
    const double lx = hx*Nx;
    const double ly = hy*Ny;
    const double px = p.x;
    const double py = p.y;
    const unsigned int rx = static_cast<unsigned int>(floor((px/hx)*lx));
    const unsigned int ry = static_cast<unsigned int>(floor((py/hy)*ly));
    const double x0 = static_cast<double>(rx-1)*hx;
    const double x1 = static_cast<double>(rx+0)*hx;
    const double x2 = static_cast<double>(rx+1)*hx;
    const double x3 = static_cast<double>(rx+2)*hx;
    const double y0 = static_cast<double>(ry-1)*hy;
    const double y1 = static_cast<double>(ry+0)*hy;
    const double y2 = static_cast<double>(ry+1)*hy;
    const double y3 = static_cast<double>(ry+2)*hy;

    const double Lx0 = (px-x1)*(px-x2)*(px-x3); const double Lx00 = (x0-x1)*(x0-x2)*(x0-x3);
    const double Lx1 = (px-x0)*(px-x2)*(px-x3); const double Lx11 = (x1-x0)*(x1-x2)*(x1-x3);
    const double Lx2 = (px-x0)*(px-x1)*(px-x3); const double Lx22 = (x2-x0)*(x2-x1)*(x2-x3);
    const double Lx3 = (px-x0)*(px-x1)*(px-x2); const double Lx33 = (x3-x0)*(x3-x1)*(x3-x2);

    const double Ly0 = (py-y1)*(py-y2)*(py-y3); const double Ly00 = (y0-y1)*(y0-y2)*(y0-y3);
    const double Ly1 = (py-y0)*(py-y2)*(py-y3); const double Ly11 = (y1-y0)*(y1-y2)*(y1-y3);
    const double Ly2 = (py-y0)*(py-y1)*(py-y3); const double Ly22 = (y2-y0)*(y2-y1)*(y2-y3);
    const double Ly3 = (py-y0)*(py-y1)*(py-y2); const double Ly33 = (y3-y0)*(y3-y1)*(y3-y2);

    double pu = 0.0;
    pu += (Ly0/Ly00) * ( (Lx0/Lx00)*m[ry-1][rx-1] + (Lx1/Lx11)*m[ry-1][rx+0] + (Lx2/Lx22)*m[ry-1][rx+1] + (Lx3/Lx33)*m[ry-1][rx+2] );
    pu += (Ly1/Ly11) * ( (Lx0/Lx00)*m[ry+0][rx-1] + (Lx1/Lx11)*m[ry+0][rx+0] + (Lx2/Lx22)*m[ry+0][rx+1] + (Lx3/Lx33)*m[ry+0][rx+2] );
    pu += (Ly2/Ly22) * ( (Lx0/Lx00)*m[ry+1][rx-1] + (Lx1/Lx11)*m[ry+1][rx+0] + (Lx2/Lx22)*m[ry+1][rx+1] + (Lx3/Lx33)*m[ry+1][rx+2] );
    pu += (Ly3/Ly33) * ( (Lx0/Lx00)*m[ry+2][rx-1] + (Lx1/Lx11)*m[ry+2][rx+0] + (Lx2/Lx22)*m[ry+2][rx+1] + (Lx3/Lx33)*m[ry+2][rx+2] );
    return pu;
}

double DeltaFunction::lumpedPoint4(const DoubleVector &v, double p, double hx, unsigned int Nx)
{
    const double lx = hx*Nx;
    const double px = p;
    const unsigned int rx = static_cast<unsigned int>(floor((px/hx)*lx));
    const double x0 = static_cast<double>(rx-1)*hx;
    const double x1 = static_cast<double>(rx+0)*hx;
    const double x2 = static_cast<double>(rx+1)*hx;
    const double x3 = static_cast<double>(rx+2)*hx;

    const double Lx0 = (px-x1)*(px-x2)*(px-x3); const double Lx00 = (x0-x1)*(x0-x2)*(x0-x3);
    const double Lx1 = (px-x0)*(px-x2)*(px-x3); const double Lx11 = (x1-x0)*(x1-x2)*(x1-x3);
    const double Lx2 = (px-x0)*(px-x1)*(px-x3); const double Lx22 = (x2-x0)*(x2-x1)*(x2-x3);
    const double Lx3 = (px-x0)*(px-x1)*(px-x2); const double Lx33 = (x3-x0)*(x3-x1)*(x3-x2);

    double pu = (Lx0/Lx00)*v[rx-1] + (Lx1/Lx11)*v[rx+0] + (Lx2/Lx22)*v[rx+1] + (Lx3/Lx33)*v[rx+2];
    return pu;
}

double DeltaFunction::lumpedPoint4(const DoubleVector &m, double p, double hx, unsigned int Nx, double &dx)
{
    const double lx = hx*Nx;
    const double px = p;
    const unsigned int rx = static_cast<unsigned int>(floor((px/hx)*lx));
    const double x0 = static_cast<double>(rx-1)*hx;
    const double x1 = static_cast<double>(rx+0)*hx;
    const double x2 = static_cast<double>(rx+1)*hx;
    const double x3 = static_cast<double>(rx+2)*hx;

    const double Lx0 = (px-x1)*(px-x2)*(px-x3); const double Lx00 = (x0-x1)*(x0-x2)*(x0-x3); const double Lx0x = (px-x1)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x1);
    const double Lx1 = (px-x0)*(px-x2)*(px-x3); const double Lx11 = (x1-x0)*(x1-x2)*(x1-x3); const double Lx1x = (px-x0)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x0);
    const double Lx2 = (px-x0)*(px-x1)*(px-x3); const double Lx22 = (x2-x0)*(x2-x1)*(x2-x3); const double Lx2x = (px-x0)*(px-x1)+(px-x1)*(px-x3)+(px-x3)*(px-x0);
    const double Lx3 = (px-x0)*(px-x1)*(px-x2); const double Lx33 = (x3-x0)*(x3-x1)*(x3-x2); const double Lx3x = (px-x0)*(px-x1)+(px-x1)*(px-x2)+(px-x2)*(px-x0);

    double pu = ( (Lx0/Lx00)*m[rx-1] + (Lx1/Lx11)*m[rx+0] + (Lx2/Lx22)*m[rx+1] + (Lx3/Lx33)*m[rx+2] );
    dx = (Lx0x/Lx00) * m[rx-1] + (Lx1x/Lx11) * m[rx+0] + (Lx2x/Lx22) * m[rx+1] + (Lx3x/Lx33) * m[rx+2];

    return pu;
}

double DeltaFunction::lumpedPoint4(const DoubleMatrix &m, const SpacePoint &p, const Dimension &dimensionX, const Dimension &dimensionY, SpacePoint &d)
{
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();
    const double Nx = dimensionX.size()-1;
    const double Ny = dimensionY.size()-1;
    const double lx = hx*Nx;
    const double ly = hy*Ny;
    const double px = p.x;
    const double py = p.y;
    const unsigned int rx = static_cast<unsigned int>(floor((px/hx)*lx));
    const unsigned int ry = static_cast<unsigned int>(floor((py/hy)*ly));
    const double x0 = static_cast<double>(rx-1)*hx;
    const double x1 = static_cast<double>(rx+0)*hx;
    const double x2 = static_cast<double>(rx+1)*hx;
    const double x3 = static_cast<double>(rx+2)*hx;
    const double y0 = static_cast<double>(ry-1)*hy;
    const double y1 = static_cast<double>(ry+0)*hy;
    const double y2 = static_cast<double>(ry+1)*hy;
    const double y3 = static_cast<double>(ry+2)*hy;

    const double Lx0 = (px-x1)*(px-x2)*(px-x3); const double Lx00 = (x0-x1)*(x0-x2)*(x0-x3); const double Lx0x = (px-x1)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x1);
    const double Lx1 = (px-x0)*(px-x2)*(px-x3); const double Lx11 = (x1-x0)*(x1-x2)*(x1-x3); const double Lx1x = (px-x0)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x0);
    const double Lx2 = (px-x0)*(px-x1)*(px-x3); const double Lx22 = (x2-x0)*(x2-x1)*(x2-x3); const double Lx2x = (px-x0)*(px-x1)+(px-x1)*(px-x3)+(px-x3)*(px-x0);
    const double Lx3 = (px-x0)*(px-x1)*(px-x2); const double Lx33 = (x3-x0)*(x3-x1)*(x3-x2); const double Lx3x = (px-x0)*(px-x1)+(px-x1)*(px-x2)+(px-x2)*(px-x0);

    const double Ly0 = (py-y1)*(py-y2)*(py-y3); const double Ly00 = (y0-y1)*(y0-y2)*(y0-y3); const double Ly0y = (py-y1)*(py-y2)+(py-y2)*(py-y3)+(py-y3)*(py-y1);
    const double Ly1 = (py-y0)*(py-y2)*(py-y3); const double Ly11 = (y1-y0)*(y1-y2)*(y1-y3); const double Ly1y = (py-y0)*(py-y2)+(py-y2)*(py-y3)+(py-y3)*(py-y0);
    const double Ly2 = (py-y0)*(py-y1)*(py-y3); const double Ly22 = (y2-y0)*(y2-y1)*(y2-y3); const double Ly2y = (py-y0)*(py-y1)+(py-y1)*(py-y3)+(py-y3)*(py-y0);
    const double Ly3 = (py-y0)*(py-y1)*(py-y2); const double Ly33 = (y3-y0)*(y3-y1)*(y3-y2); const double Ly3y = (py-y0)*(py-y1)+(py-y1)*(py-y2)+(py-y2)*(py-y0);

    double pu = 0.0;
    pu += (Ly0/Ly00) * ( (Lx0/Lx00)*m[ry-1][rx-1] + (Lx1/Lx11)*m[ry-1][rx+0] + (Lx2/Lx22)*m[ry-1][rx+1] + (Lx3/Lx33)*m[ry-1][rx+2] );
    pu += (Ly1/Ly11) * ( (Lx0/Lx00)*m[ry+0][rx-1] + (Lx1/Lx11)*m[ry+0][rx+0] + (Lx2/Lx22)*m[ry+0][rx+1] + (Lx3/Lx33)*m[ry+0][rx+2] );
    pu += (Ly2/Ly22) * ( (Lx0/Lx00)*m[ry+1][rx-1] + (Lx1/Lx11)*m[ry+1][rx+0] + (Lx2/Lx22)*m[ry+1][rx+1] + (Lx3/Lx33)*m[ry+1][rx+2] );
    pu += (Ly3/Ly33) * ( (Lx0/Lx00)*m[ry+2][rx-1] + (Lx1/Lx11)*m[ry+2][rx+0] + (Lx2/Lx22)*m[ry+2][rx+1] + (Lx3/Lx33)*m[ry+2][rx+2] );

    d.x = 0.0;
    d.x += (Lx0x/Lx00) * ( (Ly0/Ly00)*m[ry-1][rx-1] + (Ly1/Ly11)*m[ry+0][rx-1] + (Ly2/Ly22)*m[ry+1][rx-1] + (Ly3/Ly33)*m[ry+2][rx-1] );
    d.x += (Lx1x/Lx11) * ( (Ly0/Ly00)*m[ry-1][rx+0] + (Ly1/Ly11)*m[ry+0][rx+0] + (Ly2/Ly22)*m[ry+1][rx+0] + (Ly3/Ly33)*m[ry+2][rx+0] );
    d.x += (Lx2x/Lx22) * ( (Ly0/Ly00)*m[ry-1][rx+1] + (Ly1/Ly11)*m[ry+0][rx+1] + (Ly2/Ly22)*m[ry+1][rx+1] + (Ly3/Ly33)*m[ry+2][rx+1] );
    d.x += (Lx3x/Lx33) * ( (Ly0/Ly00)*m[ry-1][rx+2] + (Ly1/Ly11)*m[ry+0][rx+2] + (Ly2/Ly22)*m[ry+1][rx+2] + (Ly3/Ly33)*m[ry+2][rx+2] );

    d.y = 0.0;
    d.y += (Ly0y/Ly00) * ( (Lx0/Lx00)*m[ry-1][rx-1] + (Lx1/Lx11)*m[ry-1][rx+0] + (Lx2/Lx22)*m[ry-1][rx+1] + (Lx3/Lx33)*m[ry-1][rx+2] );
    d.y += (Ly1y/Ly11) * ( (Lx0/Lx00)*m[ry+0][rx-1] + (Lx1/Lx11)*m[ry+0][rx+0] + (Lx2/Lx22)*m[ry+0][rx+1] + (Lx3/Lx33)*m[ry+0][rx+2] );
    d.y += (Ly2y/Ly22) * ( (Lx0/Lx00)*m[ry+1][rx-1] + (Lx1/Lx11)*m[ry+1][rx+0] + (Lx2/Lx22)*m[ry+1][rx+1] + (Lx3/Lx33)*m[ry+1][rx+2] );
    d.y += (Ly3y/Ly33) * ( (Lx0/Lx00)*m[ry+2][rx-1] + (Lx1/Lx11)*m[ry+2][rx+0] + (Lx2/Lx22)*m[ry+2][rx+1] + (Lx3/Lx33)*m[ry+2][rx+2] );

    return pu;
}

double DeltaFunction::lumpedPoint4(const DoubleVector &m, double p, const Dimension &dimensionX, double &dx)
{
    const double hx = dimensionX.step();
    const double Nx = dimensionX.size()-1;
    const double lx = hx*Nx;
    const double px = p;
    const unsigned int rx = static_cast<unsigned int>(floor((px/hx)*lx));
    const double x0 = static_cast<double>(rx-1)*hx;
    const double x1 = static_cast<double>(rx+0)*hx;
    const double x2 = static_cast<double>(rx+1)*hx;
    const double x3 = static_cast<double>(rx+2)*hx;

    const double Lx0 = (px-x1)*(px-x2)*(px-x3); const double Lx00 = (x0-x1)*(x0-x2)*(x0-x3); const double Lx0x = (px-x1)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x1);
    const double Lx1 = (px-x0)*(px-x2)*(px-x3); const double Lx11 = (x1-x0)*(x1-x2)*(x1-x3); const double Lx1x = (px-x0)*(px-x2)+(px-x2)*(px-x3)+(px-x3)*(px-x0);
    const double Lx2 = (px-x0)*(px-x1)*(px-x3); const double Lx22 = (x2-x0)*(x2-x1)*(x2-x3); const double Lx2x = (px-x0)*(px-x1)+(px-x1)*(px-x3)+(px-x3)*(px-x0);
    const double Lx3 = (px-x0)*(px-x1)*(px-x2); const double Lx33 = (x3-x0)*(x3-x1)*(x3-x2); const double Lx3x = (px-x0)*(px-x1)+(px-x1)*(px-x2)+(px-x2)*(px-x0);

    double pu = ( (Lx0/Lx00)*m[rx-1] + (Lx1/Lx11)*m[rx+0] + (Lx2/Lx22)*m[rx+1] + (Lx3/Lx33)*m[rx+2] );
    dx = (Lx0x/Lx00) * m[rx-1] + (Lx1x/Lx11) * m[rx+0] + (Lx2x/Lx22) * m[rx+1] + (Lx3x/Lx33) * m[rx+2];

    return pu;
}

/**
 * @brief DeltaFunction::lumpedPointG
 * @param u matrix values of grid
 * @param m searching point
 * @param dimensionX X dimension of the grid
 * @param dimensionY Y dimension of the grid
 * @param nps
 * @param k
 * @return
 */
double DeltaFunction::lumpedPointG(const DoubleMatrix &u, const SpacePoint &m, const Dimension &dimensionX, const Dimension &dimensionY, size_t nps, size_t k)
{
    const double hx = dimensionX.step(), hy = dimensionY.step();
    const unsigned int Nx = dimensionX.size()-1, Ny = dimensionY.size()-1;
    const double lx = hx*Nx;
    const double ly = hy*Ny;
    const double px = m.x;
    const double py = m.y;
    const size_t rx = static_cast<size_t>(round((px/hx)*lx));
    const size_t ry = static_cast<size_t>(round((py/hy)*ly));

    const size_t minX = static_cast<size_t>(rx - nps*k);
    const size_t maxX = static_cast<size_t>(rx + nps*k);
    const size_t minY = static_cast<size_t>(ry - nps*k);
    const size_t maxY = static_cast<size_t>(ry + nps*k);
    const SpacePoint sigma = { nps*hx, nps*hy };

    double uv = 0.0;
    //double ux = 0.0, uy = 0.0;
    SpacePoint p;
    for (size_t j=minY; j<=maxY; j++)
    {
        p.y = j*hy;

        for (size_t i=minX; i<=maxX; i++)
        {
            p.x = i*hx;

            double w = gaussian(p, m, sigma);
            // Trapezoidal rule
            if (i==minX || i==maxX) { w *= 0.5; }
            if (j==minY || j==maxY) { w *= 0.5; }

            // Simpson's rule
            //if (i==minX || i==maxX) { w *= 1.0; } else if (i%2==1) { w *= 4.0; } else { w *= 2.0; }
            //if (j==minY || j==maxY) { w *= 1.0; } else if (j%2==1) { w *= 4.0; } else { w *= 2.0; }
            uv += w * u[j][i];
            //ux += w * (u[j][i+1]-u[j][i-1])/(2.0*hx);
            //uy += w * (u[j+1][i]-u[j-1][i])/(2.0*hy);
        }
    }
    // Trapezoidal rule
    uv *= (hx*hy);

    // Simpson's rule
    //uv *= ((hx*hy)/9.0);

    //ux *= ((hx*hy)/9.0);
    //uy *= ((hx*hy)/9.0);

    return uv;
}

double DeltaFunction::lumpedPointG(const DoubleMatrix &u, const SpacePoint &m, const Dimension &dimensionX, const Dimension &dimensionY, size_t nps, size_t k, double &dx, double &dy)
{
    const double hx = dimensionX.step(), hy = dimensionY.step();
    const unsigned int Nx = dimensionX.size()-1, Ny = dimensionY.size()-1;
    const double lx = hx*Nx;
    const double ly = hy*Ny;
    const double px = m.x;
    const double py = m.y;
    const size_t rx = static_cast<size_t>(round((px/hx)*lx));
    const size_t ry = static_cast<size_t>(round((py/hy)*ly));

    const size_t minX = static_cast<size_t>(rx - nps*k);
    const size_t maxX = static_cast<size_t>(rx + nps*k);
    const size_t minY = static_cast<size_t>(ry - nps*k);
    const size_t maxY = static_cast<size_t>(ry + nps*k);
    const SpacePoint sigma = { nps*hx, nps*hy };

    double uv = 0.0;
    dx = dy = 0.0;
    SpacePoint p;
    for (size_t j=minY; j<=maxY; j++)
    {
        p.y = j*hy;

        for (size_t i=minX; i<=maxX; i++)
        {
            p.x = i*hx;

            double wx, wy;
            double w0 = gaussian(p, m, sigma, wx, wy);
            // Trapezoidal rule
            if (i==minX || i==maxX) { w0 *= 0.5; }
            if (j==minY || j==maxY) { w0 *= 0.5; }

            // Simpson's rule
            //if (i==minX || i==maxX) { w *= 1.0; } else if (i%2==1) { w *= 4.0; } else { w *= 2.0; }
            //if (j==minY || j==maxY) { w *= 1.0; } else if (j%2==1) { w *= 4.0; } else { w *= 2.0; }
            uv += w0 * u[j][i];
            //dx += wx * u[j][i];
            //dy += wy * u[j][i];
            dx += w0 * (u[j][i+1]-u[j][i-1])/(2.0*hx);
            dy += w0 * (u[j+1][i]-u[j-1][i])/(2.0*hy);
        }
    }
    // Trapezoidal rule
    uv *= (hx*hy);
    dx *= (hx*hy);
    dy *= (hx*hy);

    // Simpson's rule
    //uv *= ((hx*hy)/9.0);

    //ux *= ((hx*hy)/9.0);
    //uy *= ((hx*hy)/9.0);

    return uv;
}

