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
    m_der_x = new double*[M+1];
    m_der_y = new double*[M+1];
    for (unsigned int m=0; m<=M; m++)
    {
        m_nodes[m] = new double[N+1];
        m_der_x[m] = new double[N+1];
        m_der_y[m] = new double[N+1];
        for (unsigned int n=0; n<=N; n++) m_nodes[m][n] = 0.0;
        for (unsigned int n=0; n<=N; n++) m_der_x[m][n] = 0.0;
        for (unsigned int n=0; n<=N; n++) m_der_y[m][n] = 0.0;
    }

    _rows = new bool[M+1]; for (unsigned int m=0; m<=M; m++) _rows[m] = false;
    _cols = new bool[N+1]; for (unsigned int n=0; n<=N; n++) _cols[n] = false;
}

auto DeltaGrid2D::distributeGauss(const SpacePoint& sp, unsigned int nodeX_per_sigmaX, unsigned int nodeY_per_sigmaY) -> void
{
    unsigned int kx = 4 * nodeX_per_sigmaX;
    unsigned int ky = 4 * nodeY_per_sigmaY;

    double sigmaX = _hx * nodeX_per_sigmaX;
    double sigmaY = _hy * nodeY_per_sigmaY;

    _rx = static_cast<unsigned int>( round(sp.x*_N) );
    _ry = static_cast<unsigned int>( round(sp.y*_M) );

    if (_rx < kx || _ry < ky || _rx > _N-kx || _ry > _M-ky)
        throw DeltaGridException("Point:["+std::to_string(sp.x)+","+std::to_string(sp.y)+"]");

    _p = sp;

    _minX = _rx - kx; _maxX = _rx + kx;
    _minY = _ry - ky; _maxY = _ry + ky;

    double sumX = 0.0;
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        double k = 1.0;
        if (n==_minX || n==_maxX) k = 0.5;
        sumX += k * exp(-((n*_hx-sp.x)*(n*_hx-sp.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= _hx;

    double sumY = 0.0;
    for (unsigned int m=_minY; m<=_maxY; m++)
    {
        double k = 1.0;
        if (m==_minY || m==_maxY) k = 0.5;
        sumY += k * exp(-((m*_hy-sp.y)*(m*_hy-sp.y))/(2.0*sigmaY*sigmaY));
    }
    sumY *= _hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/(2.0*M_PI*sigma);
    //    double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);

    for (unsigned int m=0; m<=_M; m++)
    {
        for (unsigned int n=0; n<=_N; n++) m_nodes[m][n] = 0.0;
    }

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

    for (unsigned int m=_minY; m<=_maxY; m++) _rows[m] = true;
    for (unsigned int n=_minX; n<=_maxX; n++) _cols[n] = true;
}

auto DeltaGrid2D::lumpPointGauss(const DoubleMatrix &mx) const -> double
{
    double pu = 0.0;
    for (unsigned int m=minY(); m<=maxY(); m++)
    {
        for (unsigned int n=minX(); n<=maxX(); n++)
        {
            pu += mx[m][n] * weight(n,m) * _hx * _hy;
        }
    }
    return pu;
}

auto DeltaGrid2D::lumpPointGauss(const DoubleMatrix &u, double &ux, double &uy, double &dx, double &dy) const -> double
{
    double pu = 0.0;
    ux = 0.0;
    uy = 0.0;
    dx = 0.0;
    dy = 0.0;
    for (unsigned int m=minY(); m<=maxY(); m++)
    {
        for (unsigned int n=minX(); n<=maxX(); n++)
        {
            //dx += u[m][n] * der_x()[m][n] * _hx * _hy;
            //dy += u[m][n] * der_y()[m][n] * _hx * _hy;
            //pu += mx[m][n] * weight(n,m) * _hx * _hy;

            double _w = weight(n,m);
            pu += _w * _hx * _hy * u[m][n];
            ux += _w * _hx * _hy * (u[m][n+1]-u[m][n-1])/(2.0*_hx);
            uy += _w * _hx * _hy * (u[m+1][n]-u[m-1][n])/(2.0*_hx);
        }
    }
    return pu;
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

auto DeltaGrid2D::consentrateInPoint(const DoubleMatrix &u, unsigned int v) const -> double
{
    double pu = 0.0;

    if (v==3)
    {
        //const unsigned int rx = static_cast<unsigned int>(_rx);
        //const unsigned int ry = static_cast<unsigned int>(_ry);
        const double px = p().x;
        const double py = p().y;
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


    if (v==31)
    {
        const unsigned int rx = static_cast<unsigned int>(_rx);
        const unsigned int ry = static_cast<unsigned int>(_ry);
        const double px = p().x;
        const double py = p().y;
        double x0, x1; x0 = rx*_hx; x1 = (rx+1)*_hx;
        double y0, y1; y0 = ry*_hy; y1 = (ry+1)*_hy;

        double Lx0 = (px-x1); double Lx00 = (x0-x1);
        double Lx1 = (px-x0); double Lx11 = (x1-x0);

        double Ly0 = (py-y1); double Ly00 = (y0-y1);
        double Ly1 = (py-y0); double Ly11 = (y1-y0);

        pu = 0.0;
        pu += (Lx0/Lx00) * ( (Ly0/Ly00)*u[ry+0][rx+0] + (Ly1/Ly11)*u[ry+1][rx+0] );
        pu += (Lx1/Lx11) * ( (Ly0/Ly00)*u[ry+0][rx+1] + (Ly1/Ly11)*u[ry+1][rx+1] );
        return pu;
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    if (v==4)
    {
        //const unsigned int rx = static_cast<unsigned int>(_rx);
        //const unsigned int ry = static_cast<unsigned int>(_ry);
        const double px = p().x;
        const double py = p().y;
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

    if (v==41)
    {
        const unsigned int rx = static_cast<unsigned int>(_rx);
        const unsigned int ry = static_cast<unsigned int>(_ry);
        const double px = p().x;
        const double py = p().y;
        double x0, x1, x2; x0 = (rx-1)*_hx; x1 = rx*_hx; x2 = (rx+1)*_hx;
        double y0, y1, y2; y0 = (ry-1)*_hy; y1 = ry*_hy; y2 = (ry+1)*_hy;

        double Lx0 = (px-x1)*(px-x2); double Lx00 = (x0-x1)*(x0-x2);
        double Lx1 = (px-x0)*(px-x2); double Lx11 = (x1-x0)*(x1-x2);
        double Lx2 = (px-x0)*(px-x1); double Lx22 = (x2-x0)*(x2-x1);

        double Ly0 = (py-y1)*(py-y2); double Ly00 = (y0-y1)*(y0-y2);
        double Ly1 = (py-y0)*(py-y2); double Ly11 = (y1-y0)*(y1-y2);
        double Ly2 = (py-y0)*(py-y1); double Ly22 = (y2-y0)*(y2-y1);

        pu = 0.0;
        pu += (Lx0/Lx00) * ( (Ly0/Ly00)*u[ry-1][rx-1] + (Ly1/Ly11)*u[ry+0][rx-1] + (Ly2/Ly22)*u[ry+1][rx-1] );
        pu += (Lx1/Lx11) * ( (Ly0/Ly00)*u[ry-1][rx+0] + (Ly1/Ly11)*u[ry+0][rx+0] + (Ly2/Ly22)*u[ry+1][rx+0] );
        pu += (Lx2/Lx22) * ( (Ly0/Ly00)*u[ry-1][rx+1] + (Ly1/Ly11)*u[ry+0][rx+1] + (Ly2/Ly22)*u[ry+1][rx+1] );
        return pu;
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    if (v==5)
    {
        //const unsigned int rx = static_cast<unsigned int>(_rx);
        //const unsigned int ry = static_cast<unsigned int>(_ry);
        const double px = p().x;
        const double py = p().y;
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

    if (v==51)
    {
        const unsigned int rx = static_cast<unsigned int>(_rx);
        const unsigned int ry = static_cast<unsigned int>(_ry);
        const double px = p().x;
        const double py = p().y;
        double x0, x1, x2, x3; x0 = (rx-1)*_hx; x1 = rx*_hx; x2 = (rx+1)*_hx; x3 = (rx+2)*_hx;
        double y0, y1, y2, y3; y0 = (ry-1)*_hy; y1 = ry*_hy; y2 = (ry+1)*_hy; y3 = (ry+2)*_hy;

        double Lx0 = (px-x1)*(px-x2)*(px-x3); double Lx00 = (x0-x1)*(x0-x2)*(x0-x3);
        double Lx1 = (px-x0)*(px-x2)*(px-x3); double Lx11 = (x1-x0)*(x1-x2)*(x1-x3);
        double Lx2 = (px-x0)*(px-x1)*(px-x3); double Lx22 = (x2-x0)*(x2-x1)*(x2-x3);
        double Lx3 = (px-x0)*(px-x1)*(px-x2); double Lx33 = (x3-x0)*(x3-x1)*(x3-x2);

        double Ly0 = (py-y1)*(py-y2)*(py-y3); double Ly00 = (y0-y1)*(y0-y2)*(y0-y3);
        double Ly1 = (py-y0)*(py-y2)*(py-y3); double Ly11 = (y1-y0)*(y1-y2)*(y1-y3);
        double Ly2 = (py-y0)*(py-y1)*(py-y3); double Ly22 = (y2-y0)*(y2-y1)*(y2-y3);
        double Ly3 = (py-y0)*(py-y1)*(py-y2); double Ly33 = (y3-y0)*(y3-y1)*(y3-y2);

        pu = 0.0;
        pu += (Lx0/Lx00) * ( (Ly0/Ly00)*u[ry-1][rx-1] + (Ly1/Ly11)*u[ry+0][rx-1] + (Ly2/Ly22)*u[ry+1][rx-1] + (Ly3/Ly33)*u[ry+2][rx-1] );
        pu += (Lx1/Lx11) * ( (Ly0/Ly00)*u[ry-1][rx+0] + (Ly1/Ly11)*u[ry+0][rx+0] + (Ly2/Ly22)*u[ry+1][rx+0] + (Ly3/Ly33)*u[ry+2][rx+0] );
        pu += (Lx2/Lx22) * ( (Ly0/Ly00)*u[ry-1][rx+1] + (Ly1/Ly11)*u[ry+0][rx+1] + (Ly2/Ly22)*u[ry+1][rx+1] + (Ly3/Ly33)*u[ry+2][rx+1] );
        pu += (Lx3/Lx33) * ( (Ly0/Ly00)*u[ry-1][rx+2] + (Ly1/Ly11)*u[ry+0][rx+2] + (Ly2/Ly22)*u[ry+1][rx+2] + (Ly3/Ly33)*u[ry+2][rx+2] );
        return pu;
    }

    return NAN;
}

auto DeltaGrid2D::derivativesInPoint(const DoubleMatrix &u, double &dx, double &dy, unsigned int v) const -> void
{
    if (v == 0)
    {
        const double px = p().x;
        const double py = p().y;
        const unsigned int rx = static_cast<unsigned int>(round(px/_hx));
        const unsigned int ry = static_cast<unsigned int>(round(py/_hy));
        dx = (u[ry][rx+1]-u[ry][rx-1])/(2.0*_hx);
        dy = (u[ry+1][rx]-u[ry-1][rx])/(2.0*_hy);
    }

    if (v == 3)
    {
        const double px = p().x;
        const double py = p().y;
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
        const double px = p().x;
        const double py = p().y;
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
        const double px = p().x;
        const double py = p().y;
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

auto DeltaGrid2D::consentrateInPoint(const DoubleMatrix &m, double &dx, double &dy, unsigned int v) const -> double
{
    derivativesInPoint(m, dx, dy, v);
    return consentrateInPoint(m, v);
}

auto DeltaGrid2D::cleanGrid() -> void
{
    if (_M==0 || _N == 0) return;
    for (unsigned int m=0; m<=_M; m++) delete [] m_nodes[m];
    for (unsigned int m=0; m<=_M; m++) { delete [] m_der_x[m]; delete [] m_der_y[m]; }
    delete [] m_nodes;
    delete [] m_der_x;
    delete [] m_der_y;
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

DeltaGrid1D::DeltaGrid1D() : _N(0), _hx(0.0) {}

DeltaGrid1D::~DeltaGrid1D() { cleanGrid(); }

auto DeltaGrid1D::initGrid(unsigned int N, double hx) -> void
{
    this->_N = N;
    this->_hx = hx;

    m_nodes = new double[N+1];
    for (unsigned int n=0; n<=N; n++) m_nodes[n] = 0.0;
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
    for (unsigned int n=_minX; n<=_maxX; n++)
    {
        double k = 1.0;
        if (n==_minX || n==_maxX) k = 0.5;
        sumX += k * exp(-((n*_hx-sp.x)*(n*_hx-sp.x))/(2.0*sigmaX*sigmaX));
    }
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
    _rx = static_cast<unsigned int>( round(sp.x*_N) );

    _p = sp;

    _minX = _maxX = _rx;

    m_nodes[_rx] = 1.0/_hx;
}

auto DeltaGrid1D::distributeRect4(const SpacePoint&) -> void
{}

auto DeltaGrid1D::consentrateInPoint(const DoubleVector &u, int v) const -> double
{
    double pu = 0.0;

    if (v == 1)
    {
        for (unsigned int n=minX(); n<=maxX(); n++)
        {
            pu += u[n] * weight(n) * _hx;
        }
        return pu;
    }

    if (v == 2)
    {
        double k = 1.0;
        for (unsigned int n=minX(); n<=maxX(); n++)
        {
            if (n == minX() || n == maxX()) k = 0.5; else k = 1.0;
            pu += k * u[n] * weight(n) * _hx;
        }
        return pu;
    }

    if (v == 0)
    {
        const unsigned int rx = static_cast<unsigned int>(_rx);
        const double px = p().x;
        double x0, x1; x0 = rx*_hx; x1 = (rx+1)*_hx;

        double Lx0 = (px-x1); double L00 = (x0-x1);
        double Lx1 = (px-x0); double L11 = (x1-x0);

        pu = (Lx0/L00)*u[rx+0] + (Lx1/L11)*u[rx+1];
        return pu;
    }

    if (v==3)
    {
        const unsigned int rx = static_cast<unsigned int>(_rx);
        const double px = p().x;
        double x0, x1, x2; x0 = (rx-1)*_hx; x1 = rx*_hx; x2 = (rx+1)*_hx;

        double Lx0 = (px-x1)*(px-x2); double L00 = (x0-x1)*(x0-x2);
        double Lx1 = (px-x0)*(px-x2); double L11 = (x1-x0)*(x1-x2);
        double Lx2 = (px-x0)*(px-x1); double L22 = (x2-x0)*(x2-x1);

        pu = (Lx0/L00)*u[rx-1] + (Lx1/L11)*u[rx+0] + (Lx2/L22)*u[rx+1];
        return pu;
    }

    if (v==4)
    {
        const unsigned int rx = static_cast<unsigned int>(_rx);
        const double px = p().x;
        double x0, x1, x2, x3; x0 = (rx-1)*_hx; x1 = rx*_hx; x2 = (rx+1)*_hx; x3 = (rx+2)*_hx;

        double Lx0 = (px-x1)*(px-x2)*(px-x3); double L00 = (x0-x1)*(x0-x2)*(x0-x3);
        double Lx1 = (px-x0)*(px-x2)*(px-x3); double L11 = (x1-x0)*(x1-x2)*(x1-x3);
        double Lx2 = (px-x0)*(px-x1)*(px-x3); double L22 = (x2-x0)*(x2-x1)*(x2-x3);
        double Lx3 = (px-x0)*(px-x1)*(px-x2); double L33 = (x3-x0)*(x3-x1)*(x3-x2);

        pu = (Lx0/L00)*u[rx-1] + (Lx1/L11)*u[rx+0] + (Lx2/L22)*u[rx+1] + (Lx3/L33)*u[rx+2];
        return pu;
    }

    return NAN;
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
    return consentrateInPoint(u, 4);
}

auto DeltaGrid1D::cleanGrid() -> void
{
    if (_N == 0) return;
    delete [] m_nodes;
    _N = 0;
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
