#include "iproblem2backward2d.h"

IProblem2Backward2D::~IProblem2Backward2D()
{}

void IProblem2Backward2D::calculateMVD(DoubleMatrix &p, vector<ExtendedSpaceNode2D> &info, bool use)
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

    unsigned int Lc = mParameter.Lc;
    unsigned int Lo = mParameter.Lo;

    if (use == true)
    {
        info.resize(Lc);
        for (unsigned int i=0; i<Lc; i++)
        {
            info[i].setSpaceNode(mParameter.eta[i]);
            info[i].id = i;
            info[i].extendWeights(xd, yd);
            info[i].extendLayers(L);
        }
    }

    p.clear();
    p.resize(M+1, N+1);

    DoubleMatrix ph(M+1, N+1);

    //--------------------------------------------------------------------------------------------//

    std::vector<ControlNode> controlNodes;
    for (unsigned int i=0; i<mParameter.Lc; i++) extendControlPoint(mParameter.eta[i], controlNodes, i);

    //std::vector<ControlDeltaNode> cndeltaNodes;
    //for (unsigned int i=0; i<setting.Lc; i++) extendContrlDeltaPoint(setting.eta[i], cndeltaNodes, i);

    unsigned int *v1y = new unsigned int[M+1]; for (unsigned int m=0; m<=M; m++) v1y[m] = 0;
    unsigned int *v1x = new unsigned int[N+1]; for (unsigned int n=0; n<=N; n++) v1x[n] = 0;

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            for (unsigned int j=0; j<Lo; j++)
            {
                double _delta = delta(sn, mParameter.xi[j], j);

                if (checkDelta(_delta))
                {
                    if ( v1y[m] == 0 ) v1y[m] = 1;
                    if ( v1x[n] == 0 ) v1x[n] = 1;
                }
            }
        }
    }

    //------------------------------------- initial conditions -------------------------------------//
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p[m][n] = initial(sn);
        }
    }
    layerInfo(p, L);
    //IPrinter::printMatrix(p[L]);
    //IPrinter::printSeperatorLine();
    //------------------------------------- initial conditions -------------------------------------//

    if (use ==true)
    {
        for (unsigned int i=0; i<Lc; i++)
        {
            ExtendedSpaceNode2D &pi = info[i];
            for (unsigned int r=0; r<4; r++)
            {
                for (unsigned int c=0; c<4; c++)
                {
                    unsigned int x = pi.wi[r][c].i;
                    unsigned int y = pi.wi[r][c].j;
                    pi.wi[r][c].u[L] = p[y][x];
                }
            }
        }
    }


    double a2_ht__hx2 = (a*a*ht)/(hx*hx);
    double a2_ht__hy2 = (a*a*ht)/(hy*hy);
    double lambda0_ht = lambda0*ht;
    double a2_lambda_ht__hy = (a*a*lambda*ht)/(hy);
    double a2_lambda_ht__hx = (a*a*lambda*ht)/(hx);

    std::vector<unsigned int> cntX; for (unsigned int n=0; n<=N; n++) if (v1x[n] != 0) cntX.push_back(n); unsigned int cntXSize = cntX.size();
    std::vector<unsigned int> cntY; for (unsigned int m=0; m<=M; m++) if (v1y[m] != 0) cntY.push_back(m); unsigned int cntYSize = cntY.size();

    double *a1Y = (double *) malloc(sizeof(double)*(M+1));
    double *b1Y = (double *) malloc(sizeof(double)*(M+1));
    double *c1Y = (double *) malloc(sizeof(double)*(M+1));
    double *d1Y = (double *) malloc(sizeof(double)*(M+1));
    double *x1Y = (double *) malloc(sizeof(double)*(M+1));

    double *a1X = (double *) malloc(sizeof(double)*(N+1));
    double *b1X = (double *) malloc(sizeof(double)*(N+1));
    double *c1X = (double *) malloc(sizeof(double)*(N+1));
    double *d1X = (double *) malloc(sizeof(double)*(N+1));
    double *x1X = (double *) malloc(sizeof(double)*(N+1));

    TimeNodePDE tn;
    for (unsigned int l=L-1; l!=UINT32_MAX; l--)
    {

        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        {
            tn.i = l;
            tn.t = l*ht + 0.5*ht;
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    if (v1x[n] == 0)
                    {
                        sn.i = n; sn.x = n*hx;
                        for (unsigned int m=0; m<=M; m++)
                        {
                            sn.j = m; sn.y = m*hy;

                            d1Y[m] = 2.0*p[m][n] - ht*f(sn, tn);

                            if (n==0)       d1Y[m] += a2_ht__hx2*(p[m][0]   - 2.0*p[m][1]   + p[m][2]);
                            if (n>0 && n<N) d1Y[m] += a2_ht__hx2*(p[m][n-1] - 2.0*p[m][n]   + p[m][n+1]);
                            if (n==N)       d1Y[m] += a2_ht__hx2*(p[m][N-2] - 2.0*p[m][N-1] + p[m][N]);

                            if (m == 0)
                            {
                                a1Y[0] = 0.0;
                                b1Y[0] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht - 2.0*a2_lambda_ht__hy;
                                c1Y[0] = -2.0*a2_ht__hy2;


                                d1Y[0] += ((2.0*a*a*ht)/(hy))*g3(sn, tn);
                            }
                            else if (m == M)
                            {
                                a1Y[M] = -2.0*a2_ht__hy2;
                                b1Y[M] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht - 2.0*a2_lambda_ht__hy;
                                c1Y[M] = 0.0;


                                d1Y[M] += ((2.0*a*a*ht)/(hy))*g4(sn, tn);
                            }
                            else
                            {
                                a1Y[m] = -a2_ht__hy2;
                                b1Y[m] = 2.0 + 2.0*a2_ht__hy2 + lambda0*ht;
                                c1Y[m] = -a2_ht__hy2;
                            }
                        }
                        tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
                        for (unsigned int m=0; m<=M; m++) ph[m][n] = x1Y[m];
                    }
                }
            }

            {
                double* a2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
                double* b2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
                double* c2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
                double* d2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
                double* x2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
                DoubleMatrix w2(cntXSize*(M+1), cntXSize*(M+1), 0.0);

                unsigned int offset = 0;
                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (v1x[n] == 1)
                    {
                        for (unsigned int m=0; m<=M; m++)
                        {
                            sn.j = m; sn.y = m*hy;

                            d2[offset+m] = 2.0*p[m][n] - ht*f(sn, tn);

                            if (n==0)       d2[offset+m] += a2_ht__hx2*(p[m][0]   - 2.0*p[m][1]   + p[m][2]);
                            if (n>0 && n<N) d2[offset+m] += a2_ht__hx2*(p[m][n-1] - 2.0*p[m][n]   + p[m][n+1]);
                            if (n==N)       d2[offset+m] += a2_ht__hx2*(p[m][N-2] - 2.0*p[m][N-1] + p[m][N]);

                            if (m == 0)
                            {
                                a2[offset+0] = 0.0;
                                b2[offset+0] = +2.0 + 2.0*a2_ht__hy2 + lambda0*ht - 2.0*a2_lambda_ht__hy;
                                c2[offset+0] = -2.0*a2_ht__hy2;

                                d2[offset+0] += ((2.0*a*a*ht)/hy)*g3(sn, tn);
                            }
                            else if (m == M)
                            {
                                a2[offset+M] = -2.0*a2_ht__hy2;
                                b2[offset+M] = 2.0 + 2.0*a2_ht__hy2 + lambda0*ht - 2.0*a2_lambda_ht__hy;
                                c2[offset+M] = 0.0;

                                d2[offset+M] += ((2.0*a*a*ht)/(hy))*g4(sn, tn);
                            }
                            else
                            {
                                a2[offset+m] = -a2_ht__hy2;
                                b2[offset+m] = +2.0 + 2.0*a2_ht__hy2 + lambda0*ht;
                                c2[offset+m] = -a2_ht__hy2;
                            }

                            //------------------------------------- Adding delta part -------------------------------------//
                            for (unsigned int j=0; j<Lo; j++)
                            {
                                double _delta = delta(sn, mParameter.xi[j], j);
                                if (checkDelta(_delta))
                                {
                                    for (unsigned int s=0; s<controlNodes.size(); s++)
                                    {
                                        const ControlNode &cn = controlNodes[s];

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntXSize; cs++)
                                        {
                                            if (cn.n == cntX[cs])
                                            {
                                                found = true;
                                                w2[offset+m][cs*(M+1)+cn.m] += -ht * mParameter.k[cn.i][j] * _delta * cn.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+m] += ht * mParameter.k[cn.i][j] * ph[cn.m][cn.n] * _delta * cn.w;
                                        }
                                    }
                                }
                            }
                        }
                        offset += M+1;
                    }
                }

                LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cntXSize*(M+1));

                offset = 0;
                for (unsigned int n=0; n<=N; n++)
                {
                    if (v1x[n] == 1)
                    {
                        for (unsigned int m=0; m<=M; m++)
                        {
                            ph[m][n] = x2[offset+m];
                        }
                        offset += M+1;
                    }
                }

                w2.clear();
                free(x2);
                free(d2);
                free(c2);
                free(b2);
                free(a2);
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
                for (unsigned int m=0; m<=M; m++)
                {
                    if (v1y[m] == 0)
                    {
                        sn.j = m; sn.y = m*hy;
                        for (unsigned int n=0; n<=N; n++)
                        {
                            sn.i = n; sn.x = n*hx;

                            d1X[n] = 2.0*ph[m][n] - ht*f(sn, tn);

                            if (m==0)       d1X[n] += a2_ht__hy2*(ph[0][n]   - 2.0*ph[1][n]   + ph[2][n]);
                            if (m>0 && m<M) d1X[n] += a2_ht__hy2*(ph[m-1][n] - 2.0*ph[m][n]   + ph[m+1][n]);
                            if (m==M)       d1X[n] += a2_ht__hy2*(ph[M-2][n] - 2.0*ph[M-1][n] + ph[M][n]);

                            if (n == 0)
                            {
                                a1X[0] = 0.0;
                                b1X[0] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht - 2.0*a2_lambda_ht__hx;
                                c1X[0] = -2.0*a2_ht__hx2;
                                d1X[0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
                            }
                            else if (n == N)
                            {
                                a1X[N] = -2.0*a2_ht__hx2;
                                b1X[N] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht - 2.0*a2_lambda_ht__hx;
                                c1X[N] = 0.0;
                                d1X[N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
                            }
                            else
                            {
                                a1X[n] = -a2_ht__hx2;
                                b1X[n] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
                                c1X[n] = -a2_ht__hx2;
                            }
                        }
                        tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                        for (unsigned int n=0; n<=N; n++) p[m][n] = x1X[n];
                    }
                }
            }

            {
                double* a2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
                double* b2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
                double* c2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
                double* d2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
                double* x2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
                DoubleMatrix w2(cntYSize*(N+1), cntYSize*(N+1), 0.0);

                unsigned int offset = 0;
                for (unsigned int m=0; m<=M; m++)
                {
                    if (v1y[m] == 1)
                    {
                        sn.j = m; sn.y = m*hy;
                        for (unsigned int n=0; n<=N; n++)
                        {
                            sn.i = n; sn.x = n*hx;

                            d2[offset+n] = 2.0*ph[m][n] - ht*f(sn, tn);

                            if (m==0)       d2[offset+n] += a2_ht__hy2*(ph[0][n]   - 2.0*ph[1][n]   + ph[2][n]);
                            if (m>0 && m<M) d2[offset+n] += a2_ht__hy2*(ph[m-1][n] - 2.0*ph[m][n]   + ph[m+1][n]);
                            if (m==M)       d2[offset+n] += a2_ht__hy2*(ph[M-2][n] - 2.0*ph[M-1][n] + ph[M][n]);

                            if (n == 0)
                            {
                                a2[offset+0] = 0.0;
                                b2[offset+0] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht - 2.0*a2_lambda_ht__hx;
                                c2[offset+0] = -2.0*a2_ht__hx2;
                                d2[offset+0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
                            }
                            else if (n == N)
                            {
                                a2[offset+N] = -2.0*a2_ht__hx2;
                                b2[offset+N] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht - 2.0*a2_lambda_ht__hx;
                                c2[offset+N] = 0.0;
                                d2[offset+N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
                            }
                            else
                            {
                                a2[offset+n] = -a2_ht__hx2;
                                b2[offset+n] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
                                c2[offset+n] = -a2_ht__hx2;
                            }

                            //------------------------------------ Adding delta part -------------------------------------//
                            for (unsigned int j=0; j<Lo; j++)
                            {
                                double _delta = delta(sn, mParameter.xi[j], j);
                                if (checkDelta(_delta))
                                {
                                    for (unsigned int s=0; s<controlNodes.size(); s++)
                                    {
                                        const ControlNode &cn = controlNodes[s];

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntYSize; cs++)
                                        {
                                            if (cn.m == cntY[cs])
                                            {
                                                found = true;
                                                w2[offset+n][cs*(N+1)+cn.n] += -ht * mParameter.k[cn.i][j] * _delta * cn.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+n] += ht * mParameter.k[cn.i][j] * p[cn.m][cn.n] * _delta * cn.w;
                                        }
                                    }
                                }
                            }
                        }
                        offset += N+1;
                    }
                }

                LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cntYSize*(N+1));

                offset = 0;
                for (unsigned int m=0; m<=M; m++)
                {
                    if (v1y[m] == 1)
                    {
                        for (unsigned int n=0; n<=N; n++)
                        {
                            p[m][n] = x2[offset+n];
                        }
                        offset += N+1;
                    }
                }

                w2.clear();
                free(x2);
                free(d2);
                free(c2);
                free(b2);
                free(a2);
            }
        }
        //IPrinter::printMatrix(p);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        if (use ==true)
        {
            for (unsigned int i=0; i<Lc; i++)
            {
                ExtendedSpaceNode2D &pi = info[i];
                for (unsigned int r=0; r<4; r++)
                {
                    for (unsigned int c=0; c<4; c++)
                    {
                        unsigned int x = pi.wi[r][c].i;
                        unsigned int y = pi.wi[r][c].j;
                        pi.wi[r][c].u[l] = p[y][x];
                    }
                }
            }
        }

        layerInfo(p, l);
    }

    free(x1X);
    free(d1X);
    free(c1X);
    free(b1X);
    free(a1X);

    free(x1Y);
    free(d1Y);
    free(c1Y);
    free(b1Y);
    free(a1Y);

    cntX.clear();
    cntY.clear();

    delete [] v1x;
    delete [] v1y;
    controlNodes.clear();
    ph.clear();
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

//void IProblem2Backward2D::calculateMVD(std::vector<DoubleMatrix> &p)
//{
//    Dimension xd = spaceDimension(Dimension::DimensionX);
//    Dimension yd = spaceDimension(Dimension::DimensionY);
//    Dimension td = timeDimension();

//    unsigned int N = xd.sizeN();
//    unsigned int M = yd.sizeN();
//    unsigned int L = td.sizeN();
//    double hx = xd.step();
//    double hy = yd.step();
//    double ht = td.step();

//    unsigned int Lc = setting.Lc;
//    unsigned int Lo = setting.Lo;

//    for (unsigned int l=0; l<p.size(); l++) p[l].clear();
//    p.clear();
//    p.resize(L+1);
//    DoubleMatrix ph(M+1, N+1);

//    //--------------------------------------------------------------------------------------------//

//    std::vector<ControlNode> controlNodes;
//    for (unsigned int i=0; i<setting.Lc; i++) extendControlPoint(setting.eta[i], controlNodes, i);

//    unsigned int *v1y = new unsigned int[M+1]; for (unsigned int m=0; m<=M; m++) v1y[m] = 0;
//    unsigned int *v1x = new unsigned int[N+1]; for (unsigned int n=0; n<=N; n++) v1x[n] = 0;

//    SpaceNodePDE sn;
//    for (unsigned int m=0; m<=M; m++)
//    {
//        sn.j = m; sn.y = m*hy;
//        for (unsigned int n=0; n<=N; n++)
//        {
//            sn.i = n; sn.x = n*hx;

//            for (unsigned int j=0; j<Lo; j++)
//            {
//                double _delta = delta(sn, setting.xi[j], j);

//                if (checkDelta(_delta))
//                {
//                    if ( v1y[m] == 0 ) v1y[m] = 1;
//                    if ( v1x[n] == 0 ) v1x[n] = 1;
//                }
//            }
//        }
//    }

//    //------------------------------------- initial conditions -------------------------------------//
//    p[L].resize(M+1,N+1);
//    for (unsigned int m=0; m<=M; m++)
//    {
//        sn.j = m; sn.y = m*hy;
//        for (unsigned int n=0; n<=N; n++)
//        {
//            sn.i = n; sn.x = n*hx;
//            p[L][m][n] = initial(sn);
//        }
//    }
//    //IPrinter::printMatrix(p[L]);
//    //IPrinter::printSeperatorLine();
//    //------------------------------------- initial conditions -------------------------------------//

//    double a2_ht__hx2 = (a*a*ht)/(hx*hx);
//    double a2_ht__hy2 = (a*a*ht)/(hy*hy);
//    double lambda0_ht = lambda0*ht;
//    double a2_lambda_ht__hy = (a*a*lambda*ht)/(hy);
//    double a2_lambda_ht__hx = (a*a*lambda*ht)/(hx);

//    std::vector<unsigned int> cntX; for (unsigned int n=0; n<=N; n++) if (v1x[n] != 0) cntX.push_back(n); unsigned int cntXSize = cntX.size();
//    std::vector<unsigned int> cntY; for (unsigned int m=0; m<=M; m++) if (v1y[m] != 0) cntY.push_back(m); unsigned int cntYSize = cntY.size();

//    double *a1Y = (double *) malloc(sizeof(double)*(M+1));
//    double *b1Y = (double *) malloc(sizeof(double)*(M+1));
//    double *c1Y = (double *) malloc(sizeof(double)*(M+1));
//    double *d1Y = (double *) malloc(sizeof(double)*(M+1));
//    double *x1Y = (double *) malloc(sizeof(double)*(M+1));

//    double *a1X = (double *) malloc(sizeof(double)*(N+1));
//    double *b1X = (double *) malloc(sizeof(double)*(N+1));
//    double *c1X = (double *) malloc(sizeof(double)*(N+1));
//    double *d1X = (double *) malloc(sizeof(double)*(N+1));
//    double *x1X = (double *) malloc(sizeof(double)*(N+1));

//    TimeNodePDE tn;
//    for (unsigned int l=L-1; l!=UINT32_MAX; l--)
//    {
//        p[l].resize(M+1, N+1);

//        //------------------------------------- approximatin to y direction conditions -------------------------------------//
//        {
//            tn.i = l;
//            tn.t = l*ht + 0.5*ht;
//            {
//                for (unsigned int n=0; n<=N; n++)
//                {
//                    if (v1x[n] == 0)
//                    {
//                        sn.i = n; sn.x = n*hx;
//                        for (unsigned int m=0; m<=M; m++)
//                        {
//                            sn.j = m; sn.y = m*hy;

//                            d1Y[m] = 2.0*p[l+1][m][n] - ht*f(sn, tn);

//                            if (n==0)       d1Y[m] += a2_ht__hx2*(p[l+1][m][0]   - 2.0*p[l+1][m][1]   + p[l+1][m][2]);
//                            if (n>0 && n<N) d1Y[m] += a2_ht__hx2*(p[l+1][m][n-1] - 2.0*p[l+1][m][n]   + p[l+1][m][n+1]);
//                            if (n==N)       d1Y[m] += a2_ht__hx2*(p[l+1][m][N-2] - 2.0*p[l+1][m][N-1] + p[l+1][m][N]);

//                            if (m == 0)
//                            {
//                                a1Y[0] = 0.0;
//                                b1Y[0] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht - 2.0*a2_lambda_ht__hy;
//                                c1Y[0] = -2.0*a2_ht__hy2;
//                                d1Y[0] += ((2.0*a*a*ht)/(hy))*g3(sn, tn);
//                            }
//                            else if (m == M)
//                            {
//                                a1Y[M] = -2.0*a2_ht__hy2;
//                                b1Y[M] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht - 2.0*a2_lambda_ht__hy;
//                                c1Y[M] = 0.0;
//                                d1Y[M] += ((2.0*a*a*ht)/(hy))*g4(sn, tn);
//                            }
//                            else
//                            {
//                                a1Y[m] = -a2_ht__hy2;
//                                b1Y[m] = 2.0 + 2.0*a2_ht__hy2 + lambda0*ht;
//                                c1Y[m] = -a2_ht__hy2;
//                            }
//                        }
//                        tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
//                        for (unsigned int m=0; m<=M; m++) ph[m][n] = x1Y[m];
//                    }
//                }
//            }

//            {

//                double* a2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
//                double* b2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
//                double* c2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
//                double* d2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
//                double* x2 = (double*) malloc(sizeof(double)*cntXSize*(M+1));
//                DoubleMatrix w2(cntXSize*(M+1), cntXSize*(M+1), 0.0);

//                unsigned int offset = 0;
//                for (unsigned int n=0; n<=N; n++)
//                {
//                    sn.i = n; sn.x = n*hx;

//                    if (v1x[n] == 1)
//                    {
//                        for (unsigned int m=0; m<=M; m++)
//                        {
//                            sn.j = m; sn.y = m*hy;

//                            d2[offset+m] = 2.0*p[l+1][m][n] - ht*f(sn, tn);

//                            if (n==0)       d2[offset+m] += a2_ht__hx2*(p[l+1][m][0]   - 2.0*p[l+1][m][1]   + p[l+1][m][2]);
//                            if (n>0 && n<N) d2[offset+m] += a2_ht__hx2*(p[l+1][m][n-1] - 2.0*p[l+1][m][n]   + p[l+1][m][n+1]);
//                            if (n==N)       d2[offset+m] += a2_ht__hx2*(p[l+1][m][N-2] - 2.0*p[l+1][m][N-1] + p[l+1][m][N]);

//                            if (m == 0)
//                            {
//                                a2[offset+0] = 0.0;
//                                b2[offset+0] = +2.0 + 2.0*a2_ht__hy2 + lambda0*ht - 2.0*a2_lambda_ht__hy;
//                                c2[offset+0] = -2.0*a2_ht__hy2;
//                                d2[offset+0] += ((2.0*a*a*ht)/hy)*g3(sn, tn);
//                            }
//                            else if (m == M)
//                            {
//                                a2[offset+M] = -2.0*a2_ht__hy2;
//                                b2[offset+M] = 2.0 + 2.0*a2_ht__hy2 + lambda0*ht - 2.0*a2_lambda_ht__hy;
//                                c2[offset+M] = 0.0;
//                                d2[offset+M] += ((2.0*a*a*ht)/(hy))*g4(sn, tn);
//                            }
//                            else
//                            {
//                                a2[offset+m] = -a2_ht__hy2;
//                                b2[offset+m] = +2.0 + 2.0*a2_ht__hy2 + lambda0*ht;
//                                c2[offset+m] = -a2_ht__hy2;
//                            }

//                            //------------------------------------- Adding delta part -------------------------------------//
//                            for (unsigned int j=0; j<Lo; j++)
//                            {
//                                double _delta = delta(sn, setting.xi[j], j);
//                                if (checkDelta(_delta))
//                                {
//                                    for (unsigned int s=0; s<controlNodes.size(); s++)
//                                    {
//                                        const ControlNode &cn = controlNodes[s];

//                                        bool found = false;
//                                        for (unsigned int cs=0; cs<cntXSize; cs++)
//                                        {
//                                            if (cn.n == cntX[cs])
//                                            {
//                                                found = true;
//                                                w2[offset+m][cs*(M+1)+cn.m] += -ht * setting.k[cn.i][j] * _delta * cn.w;
//                                            }
//                                        }

//                                        if (!found)
//                                        {
//                                            d2[offset+m] += ht * setting.k[cn.i][j] * ph[cn.m][cn.n] * _delta * cn.w;
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                        offset += M+1;
//                    }
//                }

//                LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cntXSize*(M+1));

//                offset = 0;
//                for (unsigned int n=0; n<=N; n++)
//                {
//                    if (v1x[n] == 1)
//                    {
//                        for (unsigned int m=0; m<=M; m++)
//                        {
//                            ph[m][n] = x2[offset+m];
//                        }
//                        offset += M+1;
//                    }
//                }

//                w2.clear();
//                free(x2);
//                free(d2);
//                free(c2);
//                free(b2);
//                free(a2);
//            }
//        }
//        //IPrinter::printMatrix(ph);
//        //IPrinter::printSeperatorLine();
//        //------------------------------------- approximatin to y direction conditions -------------------------------------//

//        //------------------------------------- approximatin to x direction conditions -------------------------------------//
//        {
//            tn.i = l;
//            tn.t = l*ht;
//            {
//                for (unsigned int m=0; m<=M; m++)
//                {
//                    if (v1y[m] == 0)
//                    {
//                        sn.j = m; sn.y = m*hy;
//                        for (unsigned int n=0; n<=N; n++)
//                        {
//                            sn.i = n; sn.x = n*hx;

//                            d1X[n] = 2.0*ph[m][n] - ht*f(sn, tn);

//                            if (m==0)       d1X[n] += a2_ht__hy2*(ph[0][n]   - 2.0*ph[1][n]   + ph[2][n]);
//                            if (m>0 && m<M) d1X[n] += a2_ht__hy2*(ph[m-1][n] - 2.0*ph[m][n]   + ph[m+1][n]);
//                            if (m==M)       d1X[n] += a2_ht__hy2*(ph[M-2][n] - 2.0*ph[M-1][n] + ph[M][n]);

//                            if (n == 0)
//                            {
//                                a1X[0] = 0.0;
//                                b1X[0] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht - 2.0*a2_lambda_ht__hx;
//                                c1X[0] = -2.0*a2_ht__hx2;
//                                d1X[0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
//                            }
//                            else if (n == N)
//                            {
//                                a1X[N] = -2.0*a2_ht__hx2;
//                                b1X[N] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht - 2.0*a2_lambda_ht__hx;
//                                c1X[N] = 0.0;
//                                d1X[N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
//                            }
//                            else
//                            {
//                                a1X[n] = -a2_ht__hx2;
//                                b1X[n] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
//                                c1X[n] = -a2_ht__hx2;
//                            }
//                        }
//                        tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
//                        for (unsigned int n=0; n<=N; n++) p[l][m][n] = x1X[n];
//                    }
//                }
//            }

//            {
//                double* a2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
//                double* b2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
//                double* c2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
//                double* d2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
//                double* x2 = (double*) malloc(sizeof(double)*cntYSize*(N+1));
//                DoubleMatrix w2(cntYSize*(N+1), cntYSize*(N+1), 0.0);

//                unsigned int offset = 0;
//                for (unsigned int m=0; m<=M; m++)
//                {
//                    if (v1y[m] == 1)
//                    {
//                        sn.j = m; sn.y = m*hy;
//                        for (unsigned int n=0; n<=N; n++)
//                        {
//                            sn.i = n; sn.x = n*hx;

//                            d2[offset+n] = 2.0*ph[m][n] - ht*f(sn, tn);

//                            if (m==0)       d2[offset+n] += a2_ht__hy2*(ph[0][n]   - 2.0*ph[1][n]   + ph[2][n]);
//                            if (m>0 && m<M) d2[offset+n] += a2_ht__hy2*(ph[m-1][n] - 2.0*ph[m][n]   + ph[m+1][n]);
//                            if (m==M)       d2[offset+n] += a2_ht__hy2*(ph[M-2][n] - 2.0*ph[M-1][n] + ph[M][n]);

//                            if (n == 0)
//                            {
//                                a2[offset+0] = 0.0;
//                                b2[offset+0] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht - 2.0*a2_lambda_ht__hx;
//                                c2[offset+0] = -2.0*a2_ht__hx2;
//                                d2[offset+0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
//                            }
//                            else if (n == N)
//                            {
//                                a2[offset+N] = -2.0*a2_ht__hx2;
//                                b2[offset+N] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht - 2.0*a2_lambda_ht__hx;
//                                c2[offset+N] = 0.0;
//                                d2[offset+N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
//                            }
//                            else
//                            {
//                                a2[offset+n] = -a2_ht__hx2;
//                                b2[offset+n] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
//                                c2[offset+n] = -a2_ht__hx2;
//                            }

//                            //------------------------------------ Adding delta part -------------------------------------//
//                            for (unsigned int j=0; j<Lo; j++)
//                            {
//                                double _delta = delta(sn, setting.xi[j], j);
//                                if (checkDelta(_delta))
//                                {
//                                    for (unsigned int s=0; s<controlNodes.size(); s++)
//                                    {
//                                        const ControlNode &cn = controlNodes[s];

//                                        bool found = false;
//                                        for (unsigned int cs=0; cs<cntYSize; cs++)
//                                        {
//                                            if (cn.m == cntY[cs])
//                                            {
//                                                found = true;
//                                                w2[offset+n][cs*(N+1)+cn.n] += -ht * setting.k[cn.i][j] * _delta * cn.w;
//                                            }
//                                        }

//                                        if (!found)
//                                        {
//                                            d2[offset+n] += ht * setting.k[cn.i][j] * p[l][cn.m][cn.n] * _delta * cn.w;
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                        offset += N+1;
//                    }
//                }

//                LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cntYSize*(N+1));

//                offset = 0;
//                for (unsigned int m=0; m<=M; m++)
//                {
//                    if (v1y[m] == 1)
//                    {
//                        for (unsigned int n=0; n<=N; n++)
//                        {
//                            p[l][m][n] = x2[offset+n];
//                        }
//                        offset += N+1;
//                    }
//                }

//                w2.clear();
//                free(x2);
//                free(d2);
//                free(c2);
//                free(b2);
//                free(a2);
//            }
//        }
//        //IPrinter::printMatrix(p);
//        //IPrinter::printSeperatorLine();
//        //------------------------------------- approximatin to x direction conditions -------------------------------------//
//    }

//    free(x1X);
//    free(d1X);
//    free(c1X);
//    free(b1X);
//    free(a1X);

//    free(x1Y);
//    free(d1Y);
//    free(c1Y);
//    free(b1Y);
//    free(a1Y);

//    cntX.clear();
//    cntY.clear();

//    delete [] v1x;
//    delete [] v1y;
//    controlNodes.clear();
//    ph.clear();
//}
