#include "iproblem2backward2d.h"

IProblem2Backward2D::IProblem2Backward2D() {}

void IProblem2Backward2D::setSettings(P2Setting s)
{
    setting = s;
}

void IProblem2Backward2D::calculateMVD(DoubleMatrix &p)
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    Dimension td = timeDimension();

    unsigned int N = xd.sizeN();
    unsigned int M = yd.sizeN();
    unsigned int L = td.sizeN();
    double hx = xd.step();
    double hy = yd.step();
    double ht = td.step();

    double a = setting.a;
    double lambda0 = setting.lambda0;
    double lambda = setting.lambda;
    double theta = setting.theta;
    unsigned int Lc = setting.Lc;
    unsigned int Lo =setting.Lo;

    U.resize(M+1, N+1, 10.0);
    uT.resize(M+1, N+1, 10.0);
    mu.resize(M+1, N+1, 1.0);

    //--------------------------------------------------------------------------------------------//

    std::vector<ControlNode> controlNodes;
    for (unsigned int i=0; i<setting.Lc; i++) extendControlPoint(setting.eta[i], controlNodes, i);

    unsigned int *v1y = new unsigned int[M+1]; for (unsigned int m=0; m<=M; m++) v1y[m] = 0;
    unsigned int *v1x = new unsigned int[N+1]; for (unsigned int n=0; n<=N; n++) v1x[n] = 0;

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;

        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            for (unsigned int j=0; j<setting.Lo; j++)
            {
                double _delta = delta(sn, setting.xi[j], j);

                if (checkDelta(_delta))
                {
                    if ( v1y[m] == 0 ) v1y[m] = 1;
                    if ( v1x[n] == 0 ) v1x[n] = 1;
                }
            }
        }
    }
    //--------------------------------------------------------------------------------------------//

    p.clear();
    p.resize(M+1, N+1);

    DoubleMatrix ph(M+1, N+1);

    //------------------------------------- initial conditions -------------------------------------//
    for (unsigned int m=0; m<=M; m++)
    {
        SpaceNodePDE sn;
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p[m][n] = initial(sn);
        }
    }
    //IPrinter::printMatrix(p);
    //IPrinter::printSeperatorLine();
    //------------------------------------- initial conditions -------------------------------------//

    TimeNodePDE tn;
    for (unsigned int l=L-1; l!=UINT32_MAX; l--)
    {
        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        {
            tn.i = l;
            tn.t = l*ht + 0.5*ht;
            {
                DoubleVector a1(M+1, 0.0);
                DoubleVector b1(M+1, 0.0);
                DoubleVector c1(M+1, 0.0);
                DoubleVector d1(M+1, 0.0);
                DoubleVector x1(M+1, 0.0);

                for (unsigned int n=0; n<=N; n++)
                {
                    SpaceNodePDE sn;
                    sn.i = n; sn.x = n*hx;

                    if (v1x[n] == 0)
                    {
                        for (unsigned int m=0; m<=M; m++)
                        {
                            sn.j = m; sn.y = m*hy;

                            d1[m] = 2.0*p[m][n] - ht*f(sn, tn);

                            if (n==0)       d1[m] += ((a*a*ht)/(hx*hx))*(p[m][0]   - 2.0*p[m][1]   + p[m][2]);
                            if (n>0 && n<N) d1[m] += ((a*a*ht)/(hx*hx))*(p[m][n-1] - 2.0*p[m][n]   + p[m][n+1]);
                            if (n==N)       d1[m] += ((a*a*ht)/(hx*hx))*(p[m][N-2] - 2.0*p[m][N-1] + p[m][N]);

                            if (m == 0)
                            {
                                b1[0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht - (2.0*a*a*lambda*ht)/(hy);
                                c1[0] = -(2.0*a*a*ht)/(hy*hy);

                                d1[0] += ((2.0*a*a*ht)/(hy))*g3(sn, tn);
                            }
                            else if (m == M)
                            {
                                a1[M] = -(2.0*a*a*ht)/(hy*hy);
                                b1[M] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht - (2.0*a*a*lambda*ht)/(hy);

                                d1[M] += ((2.0*a*a*ht)/(hy))*g4(sn, tn);

                            }
                            else
                            {
                                a1[m] = -(a*a*ht)/(hy*hy);
                                b1[m] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht;
                                c1[m] = -(a*a*ht)/(hy*hy);
                            }
                        }

                        tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), x1.data(), M+1);

                        for (unsigned int m=0; m<=M; m++) ph[m][n] = x1[m];
                    }
                }

                a1.clear();
                b1.clear();
                c1.clear();
                d1.clear();
                x1.clear();

            }

            {
                std::vector<unsigned int> cnt;
                for (unsigned int n=0; n<=N; n++) if (v1x[n] != 0) cnt.push_back(n);

                DoubleMatrix w(cnt.size()*(M+1), cnt.size()*(M+1), 0.0);
                DoubleVector d(cnt.size()*(M+1), 0.0);
                DoubleVector x(cnt.size()*(M+1), 0.0);

                unsigned int offset = 0;
                for (unsigned int n=0; n<=N; n++)
                {
                    SpaceNodePDE sn;
                    sn.i = n; sn.x = n*hx;

                    if (v1x[n] == 1)
                    {
                        for (unsigned int m=0; m<=M; m++)
                        {
                            sn.j = m; sn.y = m*hy;

                            d[offset+m] = 2.0*p[m][n] - ht*f(sn, tn);

                            if (n==0)       d[offset+m] += ((a*a*ht)/(hx*hx))*(p[m][0]   - 2.0*p[m][1]   + p[m][2]);
                            if (n>0 && n<N) d[offset+m] += ((a*a*ht)/(hx*hx))*(p[m][n-1] - 2.0*p[m][n]   + p[m][n+1]);
                            if (n==N)       d[offset+m] += ((a*a*ht)/(hx*hx))*(p[m][N-2] - 2.0*p[m][N-1] + p[m][N]);

                            if (m == 0)
                            {
                                w[offset+0][offset+0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht - (2.0*a*a*lambda*ht)/(hy);
                                w[offset+0][offset+1] += -(2.0*a*a*ht)/(hy*hy);

                                d[offset+0] += ((2.0*a*a*ht)/hy)*g3(sn, tn);
                            }
                            else if (m == M)
                            {
                                w[offset+M][offset+(M-1)] += -(2.0*a*a*ht)/(hy*hy);
                                w[offset+M][offset+(M-0)] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht - (2.0*a*a*lambda*ht)/(hy);

                                d[offset+M] += ((2.0*a*a*ht)/(hy))*g4(sn, tn);
                            }
                            else
                            {
                                w[offset+m][offset+(m-1)] += -(a*a*ht)/(hy*hy);
                                w[offset+m][offset+(m+0)] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht;
                                w[offset+m][offset+(m+1)] += -(a*a*ht)/(hy*hy);
                            }

                            //------------------------------------- Adding delta part -------------------------------------//
                            for (unsigned int j=0; j<Lo; j++)
                            {
                                double _delta = delta(sn, setting.xi[j], j);
                                if (checkDelta(_delta))
                                {
                                    for (unsigned int s=0; s<controlNodes.size(); s++)
                                    {
                                        const ControlNode &cn = controlNodes[s];

                                        bool found = false;
                                        unsigned int offs = 0;
                                        for (unsigned int cs=0; cs<cnt.size(); cs++)
                                        {
                                            if (cn.n == cnt[cs])
                                            {
                                                found = true;
                                                offs = cs*(M+1);

                                                //unsigned int jinx = ons[s].n;
                                                unsigned int jiny = cn.m;
                                                w[offset+m][offs+jiny] += -ht*setting.k[cn.i][j] * _delta * cn.w;
                                            }
                                        }

                                        //                            if (found)
                                        //                            {
                                        //                                //printf("%d\n", on.m);
                                        //                                //d2[offset+n] += -ht*k[i][j]*z[i][j] * _delta * ons[j].w;

                                        //                                unsigned int jinx = ons[s].n;
                                        //                                unsigned int jiny = ons[s].m;

                                        //                                w2[offset+n][offs+jinx] += -ht*k[i][on.j] * _delta * ons[s].w;
                                        //                                printf("+ %4d %4d %4d %4d %4d %4d %f\n", i, sn2.i, sn2.j, on.m, on.n, offset+n, NAN);
                                        //                            }
                                        //                            else
                                        if (!found)
                                        {
                                            d[offset+m] += ht*setting.k[cn.i][j] * ph[cn.m][cn.n] * _delta * cn.w;
                                        }
                                    }
                                }
                            }
                        }
                        offset += M+1;
                    }
                }

                LinearEquation::GaussianElimination(w,d,x);

                offset = 0;
                for (unsigned int n=0; n<=N; n++)
                {
                    if (v1x[n] == 1)
                    {
                        for (unsigned int m=0; m<=M; m++)
                        {
                            ph[m][n] = x[offset+m];
                        }
                        offset += M+1;
                    }
                }

                w.clear();
                d.clear();
                x.clear();
                cnt.clear();
            }
        }
        //IPrinter::printMatrix(ph);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        {
            tn.i = l;
            tn.t = l*ht;
            {
                DoubleVector a1(N+1, 0.0);
                DoubleVector b1(N+1, 0.0);
                DoubleVector c1(N+1, 0.0);
                DoubleVector d1(N+1, 0.0);
                DoubleVector x1(N+1, 0.0);

                for (unsigned int m=0; m<=M; m++)
                {
                    SpaceNodePDE sn;
                    sn.j = m; sn.y = m*hy;

                    if (v1y[m] == 0)
                    {
                        for (unsigned int n=0; n<=N; n++)
                        {
                            sn.i = n; sn.x = n*hx;

                            d1[n] = 2.0*ph[m][n] - ht*f(sn, tn);

                            if (m==0)       d1[n] += ((a*a*ht)/(hy*hy))*(ph[0][n]   - 2.0*ph[1][n]   + ph[2][n]);
                            if (m>0 && m<M) d1[n] += ((a*a*ht)/(hy*hy))*(ph[m-1][n] - 2.0*ph[m][n]   + ph[m+1][n]);
                            if (m==M)       d1[n] += ((a*a*ht)/(hy*hy))*(ph[M-2][n] - 2.0*ph[M-1][n] + ph[M][n]);

                            if (n == 0)
                            {
                                b1[0] = 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht - (2.0*a*a*lambda*ht)/(hx);
                                c1[0] = -(2.0*a*a*ht)/(hx*hx);

                                d1[0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
                            }
                            else if (n == N)
                            {
                                a1[N] = -(2.0*a*a*ht)/(hx*hx);
                                b1[N] = 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht - (2.0*a*a*lambda*ht)/(hx);

                                d1[N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
                            }
                            else
                            {
                                a1[n] = -(a*a*ht)/(hx*hx);
                                b1[n] = 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                                c1[n] = -(a*a*ht)/(hx*hx);
                            }
                        }

                        tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), x1.data(), N+1);

                        for (unsigned int n=0; n<=N; n++) p[m][n] = x1[n];
                    }
                }

                a1.clear();
                b1.clear();
                c1.clear();
                d1.clear();
                x1.clear();
            }

            {
                std::vector<unsigned int> cnt;
                for (unsigned int m=0; m<=M; m++) if (v1y[m] != 0) cnt.push_back(m);

                DoubleMatrix w(cnt.size()*(N+1), cnt.size()*(N+1), 0.0);
                DoubleVector d(cnt.size()*(N+1), 0.0);
                DoubleVector x(cnt.size()*(N+1), 0.0);

                unsigned int offset = 0;
                for (unsigned int m=0; m<=M; m++)
                {
                    SpaceNodePDE sn2;
                    sn2.j = m; sn2.y = m*hy;

                    if (v1y[m] == 1)
                    {
                        for (unsigned int n=0; n<=N; n++)
                        {
                            sn2.i = n; sn2.x = n*hx;

                            d[offset+n] = 2.0*ph[m][n] - ht*f(sn2, tn);

                            if (m==0)       d[offset+n] += ((a*a*ht)/(hy*hy))*(ph[0][n]   - 2.0*ph[1][n]   + ph[2][n]);
                            if (m>0 && m<M) d[offset+n] += ((a*a*ht)/(hy*hy))*(ph[m-1][n] - 2.0*ph[m][n]   + ph[m+1][n]);
                            if (m==M)       d[offset+n] += ((a*a*ht)/(hy*hy))*(ph[M-2][n] - 2.0*ph[M-1][n] + ph[M][n]);

                            if (n == 0)
                            {
                                w[offset+0][offset+0] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht - (2.0*a*a*lambda*ht)/(hx);
                                w[offset+0][offset+1] += -(2.0*a*a*ht)/(hx*hx);

                                d[offset+0] += ((2.0*a*a*ht)/hx)*g1(sn2, tn);
                            }
                            else if (n == N)
                            {
                                w[offset+N][offset+(N-1)] += -(2.0*a*a*ht)/(hx*hx);
                                w[offset+N][offset+(N-0)] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht - (2.0*a*a*lambda*ht)/(hx);

                                d[offset+N] += ((2.0*a*a*ht)/(hx))*g2(sn2, tn);
                            }
                            else
                            {
                                w[offset+n][offset+(n-1)] += -(a*a*ht)/(hx*hx);
                                w[offset+n][offset+(n+0)] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                                w[offset+n][offset+(n+1)] += -(a*a*ht)/(hx*hx);
                            }

                            //------------------------------------ Adding delta part -------------------------------------//
                            for (unsigned int j=0; j<Lo; j++)
                            {
                                double _delta = delta(sn2, setting.xi[j], j);
                                if (checkDelta(_delta))
                                {
                                    for (unsigned int s=0; s<controlNodes.size(); s++)
                                    {
                                        const ControlNode &cn = controlNodes[s];

                                        bool found = false;
                                        unsigned int offs = 0;
                                        for (unsigned int cs=0; cs<cnt.size(); cs++)
                                        {
                                            if (cn.m == cnt[cs])
                                            {
                                                found = true;
                                                offs = cs*(N+1);

                                                unsigned int jinx = cn.n;
                                                //unsigned int jiny = ons[s].m;
                                                w[offset+n][offs+jinx] += -ht*setting.k[cn.i][j] * _delta * cn.w;
                                            }
                                        }

                                        //                            if (found)
                                        //                            {
                                        //                                //printf("%d\n", on.m);
                                        //                                //d2[offset+n] += -ht*k[i][j]*z[i][j] * _delta * ons[j].w;

                                        //                                unsigned int jinx = ons[s].n;
                                        //                                unsigned int jiny = ons[s].m;

                                        //                                w2[offset+n][offs+jinx] += -ht*k[i][on.j] * _delta * ons[s].w;
                                        //                                printf("+ %4d %4d %4d %4d %4d %4d %f\n", i, sn2.i, sn2.j, on.m, on.n, offset+n, NAN);
                                        //                            }
                                        //                            else
                                        if (!found)
                                        {
                                            d.at(offset+n) += ht*setting.k[cn.i][j] * p[cn.m][cn.n] * _delta * cn.w;
                                        }
                                    }
                                }
                            }
                        }
                        offset += N+1;
                    }
                }

                LinearEquation::GaussianElimination(w,d,x);

                offset = 0;
                for (unsigned int m=0; m<=M; m++)
                {
                    if (v1y[m] == 1)
                    {
                        for (unsigned int n=0; n<=N; n++)
                        {
                            p[m][n] = x[offset+n];
                        }
                        offset += N+1;
                    }
                }

                w.clear();
                d.clear();
                x.clear();
                cnt.clear();
            }
        }
        //IPrinter::printMatrix(p);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to x direction conditions -------------------------------------//
    }

    //--------------------------------------------------------------------------------------------//

    delete [] v1x;
    delete [] v1y;
    controlNodes.clear();
}

double IProblem2Backward2D::initial(const SpaceNodePDE &sn) const
{
    unsigned int i = sn.i;
    unsigned int j = sn.j;
    return -2.0 * mu[j][i] * (uT[j][i]-U[j][i]) + h(sn);
}

double IProblem2Backward2D::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double IProblem2Backward2D::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double y = sn.y;
    double t = tn.t;

    double res = 1.0 + 4.0*setting.a*setting.a - setting.lambda0*P(x,y,t);

    double W = 0.0;
    for (unsigned int j=0; j<setting.Lo; j++)
    {
        double _delta = delta(sn, setting.xi[j], 0);
        if (checkDelta(_delta))
        {
            double vi = 0.0;
            for (unsigned int i=0; i<setting.Lc; i++)
            {
                vi += setting.k[i][j] * P(setting.eta[i].x, setting.eta[i].y, t);
            }
            W += vi * _delta;
        }
    }
    res += W;

    return res;
}

bool IProblem2Backward2D::checkDelta(double _delta) const
{
    return _delta != 0.0;
}

double IProblem2Backward2D::delta(const SpaceNodePDE &sn, const SpaceNodePDE &xi, unsigned int j, unsigned int) const
{
    return delta4(sn, xi, j);
}

double IProblem2Backward2D::delta1(const SpaceNodePDE &sn, const SpaceNodePDE &xi, unsigned int) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();

    double res = 0.0;
    if ( sn.i == xi.i && sn.j == xi.j ) res = 1.0/(hx*hy);
    return res;
}

double IProblem2Backward2D::delta2(const SpaceNodePDE &sn, const SpaceNodePDE &xi, unsigned int) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();
    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double sigmaX = hx;
    double sigmaY = hy;

    unsigned int rx = (unsigned int)(round(xi.x*Nx));
    unsigned int ry = (unsigned int)(round(xi.y*Ny));

    double res = 0.0;
    if (rx-3 <= sn.i && sn.i <= rx+3 && ry-3 <= sn.j && sn.j <= ry+3)
    {
        res = (1.0/(2.0*M_PI*sigmaX*sigmaY))*exp(-0.5*(((sn.x-xi.x)*(sn.x-xi.x))/(sigmaX*sigmaX)+((sn.y-xi.y)*(sn.y-xi.y))/(sigmaY*sigmaY)));
    }

    return res;
}

double IProblem2Backward2D::delta3(const SpaceNodePDE &sn UNUSED_PARAM, const SpaceNodePDE &xi, unsigned int) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();

    double resX = 0.0;
    {
        double dhx = fabs(sn.x-xi.x);
        double hx3 = hx*hx*hx;
        double hx32 = (1.0/(2.0*hx3));
        double hx36 = (1.0/(6.0*hx3));

        if ( dhx <= hx )
        {
            resX = ((2.0*hx-dhx)*(hx-dhx)*(hx+dhx)) * hx32;
        }
        if ( hx < dhx && dhx <= 2.0*hx )
        {
            resX = ((2.0*hx-dhx)*(hx-dhx)*(3.0*hx-dhx)) * hx36;
        }
    }

    double resY = 0.0;
    {
        double dhy = fabs(sn.y-xi.y);
        double hy3 = hy*hy*hy;
        double hy32 = (1.0/(2.0*hy3));
        double hy36 = (1.0/(6.0*hy3));

        if ( dhy <= hy )
        {
            resY = ((2.0*hy-dhy)*(hy-dhy)*(hy+dhy)) * hy32;
        }
        if ( hy < dhy && dhy <= 2.0*hy )
        {
            resY = ((2.0*hy-dhy)*(hy-dhy)*(3.0*hy-dhy)) * hy36;
        }
    }

    return (resX*resY)/(hx*hy);
}

double IProblem2Backward2D::delta4(const SpaceNodePDE &sn UNUSED_PARAM, const SpaceNodePDE &xi, unsigned int) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();

    double resX = 0.0;
    {
        double dhx = fabs(sn.x-xi.x);
        if ( dhx <= hx )
        {
            resX = (hx-dhx) / hx;
        }
    }

    double resY = 0.0;
    {
        double dhy = fabs(sn.y-xi.y);
        if ( dhy <= hy )
        {
            resY = (hy-dhy) / hy;
        }
    }

    double res = (resX*resY)/(hx*hy);

    return res;
}

void IProblem2Backward2D::extendControlPoint(const SpaceNodePDE eta, std::vector<ControlNode> &ons, unsigned int i) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);

    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double hx = xd.step();
    double hy = yd.step();

    unsigned int rx = (unsigned int)(round(eta.x*Nx));
    unsigned int ry = (unsigned int)(round(eta.y*Ny));

    double hx3 = hx*hx*hx;
    double hx32 = (1.0/(2.0*hx3));
    double hx36 = (1.0/(6.0*hx3));

    double hy3 = hy*hy*hy;
    double hy32 = (1.0/(2.0*hy3));
    double hy36 = (1.0/(6.0*hy3));

    ControlNode cn;
    double dx = 0.0;
    double dy = 0.0;

    if (rx*hx <= eta.x && ry*hy <= eta.y ) // left bottom
    {
        cn.n = rx-1; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        /////////////////////////
        cn.n = rx+0; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);;
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        /////////////////////////
        cn.n = rx+1; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        ////////////////////////
        cn.n = rx+2; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+2; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+2; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+2; cn.m = ry+2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);
        return;
    }

    if ( rx*hx <= eta.x && ry*hy >= eta.y ) // left top
    {
        cn.n = rx-1; cn.m = ry-2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        /////////////////////////
        cn.n = rx+0; cn.m = ry-2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);;
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        /////////////////////////
        cn.n = rx+1; cn.m = ry-2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        ////////////////////////
        cn.n = rx+2; cn.m = ry-2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+2; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+2; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+2; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);
        return;
    }

    if ( rx*hx >= eta.x && ry*hy >= eta.y ) // right top
    {
        cn.n = rx-2; cn.m = ry-2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx-2; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-2; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-2; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        /////////////////////////
        cn.n = rx-1; cn.m = ry-2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);;
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        /////////////////////////
        cn.n = rx+0; cn.m = ry-2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);


        ////////////////////////
        cn.n = rx+1; cn.m = ry-2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);
        return;
    }

    if ( rx*hx >= eta.x && ry*hy <= eta.y ) // right bottom
    {
        cn.n = rx-2; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx-2; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-2; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-2; cn.m = ry+2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        /////////////////////////
        cn.n = rx-1; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx-1; cn.m = ry+2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);;
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        /////////////////////////
        cn.n = rx+0; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+0; cn.m = ry+2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);


        ////////////////////////
        cn.n = rx+1; cn.m = ry-1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+0; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+1; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(cn);

        cn.n = rx+1; cn.m = ry+2; cn.x = cn.n*hx; cn.y = cn.m*hy; cn.eta = eta; cn.i = i;
        dx = fabs(cn.x-eta.x);
        dy = fabs(cn.y-eta.y);
        cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(cn);
        return;
    }
}

double IProblem2Backward2D::g1(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return -setting.lambda*(sn.y*sn.y+tn.t);
}

double IProblem2Backward2D::g2(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 2.0 - setting.lambda*(1.0+sn.y*sn.y+tn.t);
}

double IProblem2Backward2D::g3(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return -setting.lambda*(sn.x*sn.x+tn.t);
}

double IProblem2Backward2D::g4(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 2.0 - setting.lambda*(sn.x*sn.x+1.0+tn.t);
}

double IProblem2Backward2D::h(const SpaceNodePDE &sn) const
{
    double x = sn.x; unsigned int i = sn.i;
    double y = sn.y; unsigned int j = sn.j;

    return (x*x + y*y + 1.0) + 2.0 * mu[j][i] * (uT[j][i] - U[j][i]);
}

double IProblem2Backward2D::P(double x, double y, double t) const
{
    return x*x + y*y + t;
}
