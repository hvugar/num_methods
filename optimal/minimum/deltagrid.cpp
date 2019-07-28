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
            uy += _w * _hx * _hy * (u[m+1][n]-u[m-1][n])/(2.0*_hx);
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
