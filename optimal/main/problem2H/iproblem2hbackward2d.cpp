#include "iproblem2hbackward2d.h"

void IProblem2HBackward2D::calculateMVD(DoubleMatrix &p, vector<ExtendedSpaceNode2DH> &info, bool use) const
{
    Dimension dimX = spaceDimension(Dimension::DimensionX);
    Dimension dimY = spaceDimension(Dimension::DimensionY);
    Dimension time = timeDimension();

    unsigned int N = dimX.sizeN();
    unsigned int M = dimY.sizeN();
    unsigned int L = time.sizeN();

    double hx = dimX.step();
    double hy = dimY.step();
    double ht = time.step();

    double lambda = mParameter.lambda;
    double a = mParameter.a;

    DoubleMatrix p0(M+1, N+1);
    DoubleMatrix p1(M+1, N+1);
    DoubleMatrix ph(M+1, N+1);

    p.clear();
    p.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<IProblem2H2D::ExtendedSpacePointNode> obsDeltaNodes;
    for (unsigned int j=0; j<mParameter.No; j++) IProblem2H2D::distributeDelta(mParameter.xi[j], obsDeltaNodes, j, dimX, dimY);

    std::vector<IProblem2H2D::ExtendedSpacePointNode> cntPointNodes;
    for (unsigned int i=0; i<mParameter.Nc; i++) IProblem2H2D::distributeDelta(mParameter.eta[i], cntPointNodes, i, dimX, dimY);

    SpaceNodePDE sn;

    //--------------------------------------------------------------------------------------------//
    vector<unsigned int> rows0;
    vector<unsigned int> rows1;
    vector<unsigned int> rows2;
    for (unsigned int ny=0; ny<=M; ny++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<obsDeltaNodes.size(); i++)
        {
            const IProblem2H2D::ExtendedSpacePointNode &cdn = obsDeltaNodes.at(i);
            if (cdn.j == ny)
            {
                found1 = true;
                for (unsigned int j=0; j<cntPointNodes.size(); j++)
                {
                    const IProblem2H2D::ExtendedSpacePointNode &on = cntPointNodes.at(j);
                    if (on.j == ny)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == true  && found2 == true)  rows2.push_back(ny);
        if (found1 == true)                     rows1.push_back(ny);
        if (found1 == false && found2 == false) rows0.push_back(ny);
    }

    //--------------------------------------------------------------------------------------------//
    vector<unsigned int> cols0;
    vector<unsigned int> cols1;
    vector<unsigned int> cols2;
    for (unsigned int nx=0; nx<=N; nx++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<obsDeltaNodes.size(); i++)
        {
            const IProblem2H2D::ExtendedSpacePointNode &cdn = obsDeltaNodes.at(i);
            if (cdn.i == nx)
            {
                found1 = true;
                for (unsigned int j=0; j<cntPointNodes.size(); j++)
                {
                    const IProblem2H2D::ExtendedSpacePointNode &on = cntPointNodes.at(j);
                    if (on.i == nx)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == true  && found2 == true)  cols2.push_back(nx);
        if (found1 == true)                     cols1.push_back(nx);
        if (found1 == false && found2 == false) cols0.push_back(nx);
    }
    //--------------------------------------------------------------------------------------------//

    //-------------------------------------------- info --------------------------------------------//
    if (use == true)
    {
        info.resize(mParameter.Nc);
        for (unsigned int i=0; i<mParameter.Nc; i++)
        {
            ExtendedSpaceNode2DH &e = info[i];
            e.setSpaceNode(mParameter.eta[i]);
            e.id = i;
            e.extendWeights(dimX, dimY);
            e.extendLayers(L/2+1);
        }
    }
    //-------------------------------------------- info --------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p0[m][n] = initial1(sn);
        }
    }
    //IPrinter::printSeperatorLine();
    //IPrinter::printMatrix(p0);
    if (use == true) add2Info(p0, info, L/2);
    layerInfo(p0, L/2);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p1[m][n] = p0[m][n] + initial2(sn)*ht;

            double diff = 0.0;
            if (n==0 )     diff += a*a*(p0[m][n]-2.0*p0[m][n+1]+p0[m][n+2])/(hx*hx);
            else if (n==N) diff += a*a*(p0[m][n-2]-2.0*p0[m][n-1]+p0[m][n])/(hx*hx);
            else           diff += a*a*(p0[m][n-1]-2.0*p0[m][n]+p0[m][n+1])/(hx*hx);

            if (m==0 )     diff += a*a*(p0[m][n]-2.0*p0[m+1][n]+p0[m+2][n])/(hy*hy);
            else if (m==M) diff += a*a*(p0[m-2][n]-2.0*p0[m-1][n]+p0[m][n])/(hy*hy);
            else           diff += a*a*(p0[m-1][n]-2.0*p0[m][n]+p0[m+1][n])/(hy*hy);
            p1[m][n] += 4.0*ht*ht*0.5*diff;
        }
    }

    if (use == true) add2Info(p1, info, L/2-1);
    layerInfo(p1, L/2-1);
    //------------------------------------- initial conditions -------------------------------------//

    double *a1X = (double *) malloc(sizeof(double)*(N+1));
    double *b1X = (double *) malloc(sizeof(double)*(N+1));
    double *c1X = (double *) malloc(sizeof(double)*(N+1));
    double *d1X = (double *) malloc(sizeof(double)*(N+1));
    double *x1X = (double *) malloc(sizeof(double)*(N+1));

    double *a1Y = (double *) malloc(sizeof(double)*(M+1));
    double *b1Y = (double *) malloc(sizeof(double)*(M+1));
    double *c1Y = (double *) malloc(sizeof(double)*(M+1));
    double *d1Y = (double *) malloc(sizeof(double)*(M+1));
    double *x1Y = (double *) malloc(sizeof(double)*(M+1));

    //TimeNodePDE tn;

    for (unsigned int l1=2; l1<=L; l1+=2)
    {
        unsigned int l = L-l1;

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        //tn.i = l; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;

                d1X[n] = (2.0+2.0*lambda*ht)*p1[m][n] - (1.0+(lambda*ht)/2.0)*p0[m][n];

                if (m == 0)     d1X[n] += ((a*a*ht*ht)/(hy*hy))*(p1[0][n]   - 2.0*p1[1][n]   + p1[2][n]);
                if (m>0 && m<M) d1X[n] += ((a*a*ht*ht)/(hy*hy))*(p1[m-1][n] - 2.0*p1[m][n]   + p1[m+1][n]);
                if (m == M)     d1X[n] += ((a*a*ht*ht)/(hy*hy))*(p1[M-2][n] - 2.0*p1[M-1][n] + p1[M][n]);

                if (n==0)
                {
                    a1X[0] = 0.0;
                    b1X[0] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                    c1X[0] = -2.0*(a*a*ht*ht)/(hx*hx);
                    d1X[0] += 0.0;
                }
                else if(n == N)
                {
                    a1X[N] = -2.0*(a*a*ht*ht)/(hx*hx);
                    b1X[N] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                    c1X[N] = 0.0;
                    d1X[N] += 0.0;
                }
                else
                {
                    a1X[n] = -(a*a*ht*ht)/(hx*hx);
                    b1X[n] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                    c1X[n] = -(a*a*ht*ht)/(hx*hx);
                }
            }
            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
            for (unsigned int n=0; n<=N; n++) ph[m][n] = x1X[n];
        }

        if (rows2.size() == 0)
        {
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    d1X[n] = (2.0+2.0*lambda*ht)*p1[m][n] - (1.0+(lambda*ht)/2.0)*p0[m][n];

                    if (m == 0)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(p1[0][n]  -2.0*p1[1][n]  +p1[2][n]);
                    if (m>0 && m<M) d1X[n] = ((a*a*ht*ht)/(hy*hy))*(p1[m-1][n]-2.0*p1[m][n]  +p1[m+1][n]);
                    if (m == M)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(p1[M-2][n]-2.0*p1[M-1][n]+p1[M][n]);

                    if (n==0)
                    {
                        a1X[0] = 0.0;
                        b1X[0] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                        c1X[0] = -2.0*(a*a*ht*ht)/(hx*hx);
                        d1X[0] += 0.0;
                    }
                    else if(n==N)
                    {
                        a1X[N] = -2.0*(a*a*ht*ht)/(hx*hx);
                        b1X[N] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                        c1X[N] = 0.0;
                        d1X[N] += 0.0;
                    }
                    else
                    {
                        a1X[n] = -(a*a*ht*ht)/(hx*hx);
                        b1X[n] = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx)) + 3.0*(lambda*ht)/2.0;
                        c1X[n] = -(a*a*ht*ht)/(hx*hx);
                    }

                    for (unsigned int cni=0; cni<obsDeltaNodes.size(); cni++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &odn = obsDeltaNodes.at(cni);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int s=0; s<cntPointNodes.size(); s++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &cpn = cntPointNodes[s];
                                d1X[n] += ht*ht * mParameter.k[odn.id][cpn.id] * ph[cpn.j][cpn.i] * odn.w * (cpn.w * (hx*hy));
                            }

                            for (unsigned int j=0; j<mParameter.No; j++)
                            {
                                d1X[n] -= ht*ht * mParameter.k[odn.id][j] * mParameter.z[odn.id][j] * odn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                for (unsigned int n=0; n<=N; n++) ph[m][n] = x1X[n];
            }
        }
        else
        {

        }
        //--------------------------------------------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        //tn.i = l+1; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int col=0; col<cols0.size(); col++)
        {
            unsigned int n = cols0.at(col);
            sn.i = n; sn.x = n*hx;
            for (unsigned int m=0; m<=M; m++)
            {
                sn.j = m; sn.y = m*hy;

                d1Y[m] = (2.0+2.0*lambda*ht)*ph[m][n] - (1.0+(lambda*ht)/2.0)*p1[m][n];

                if (n==0)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(ph[m][0]   - 2.0*ph[m][1]   + ph[m][2]);
                if (n>0 && n<N) d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(ph[m][n-1] - 2.0*ph[m][n]   + ph[m][n+1]);
                if (n==N)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(ph[m][N-2] - 2.0*ph[m][N-1] + ph[m][N]);

                if (m == 0)
                {
                    a1Y[0] = 0.0;
                    b1Y[0] = +1.0 + 2.0*((a*a*ht*ht)/(hy*hy)) + 3.0*(lambda*ht)/2.0;
                    c1Y[0] = -2.0*(a*a*ht*ht)/(hy*hy);
                    d1Y[0] += 0.0;
                }
                else if (m == M)
                {
                    a1Y[M] = -2.0*(a*a*ht*ht)/(hy*hy);
                    b1Y[M] = +1.0 + 2.0*((a*a*ht*ht)/(hy*hy)) + 3.0*(lambda*ht)/2.0;
                    c1Y[M] = 0.0;
                    d1Y[M] += 0.0;
                }
                else
                {
                    a1Y[m] = -(a*a*ht*ht)/(hy*hy);
                    b1Y[m] = +1.0 + 2.0*(a*a*ht*ht)/(hy*hy) + 3.0*(lambda*ht)/2.0;
                    c1Y[m] = -(a*a*ht*ht)/(hy*hy);
                }
            }
            tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
            for (unsigned int m=0; m<=M; m++) p[m][n] = x1Y[m];
        }

        if (cols2.size() == 0)
        {
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    d1Y[m] = (2.0+2.0*lambda*ht)*ph[m][n] - (1.0+(lambda*ht)/2.0)*p1[m][n];

                    if (n==0)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(ph[m][0]   - 2.0*ph[m][1]   + ph[m][2]);
                    if (n>0 && n<N) d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(ph[m][n-1] - 2.0*ph[m][n]   + ph[m][n+1]);
                    if (n==N)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(ph[m][N-2] - 2.0*ph[m][N-1] + ph[m][N]);

                    if (m == 0)
                    {
                        a1Y[0] = 0.0;
                        b1Y[0] = +1.0 + 2.0*((a*a*ht*ht)/(hy*hy)) + 3.0*(lambda*ht)/2.0;
                        c1Y[0] = -2.0*(a*a*ht*ht)/(hy*hy);
                        d1Y[0] += 0.0;
                    }
                    else if (m == M)
                    {
                        a1Y[M] = -2.0*(a*a*ht*ht)/(hy*hy);
                        b1Y[M] = +1.0 + 2.0*((a*a*ht*ht)/(hy*hy)) + 3.0*(lambda*ht)/2.0;
                        c1Y[M] = 0.0;
                        d1Y[M] += 0.0;
                    }
                    else
                    {
                        a1Y[m] = -(a*a*ht*ht)/(hy*hy);
                        b1Y[m] = +1.0 + 2.0*(a*a*ht*ht)/(hy*hy) + 3.0*(lambda*ht)/2.0;
                        c1Y[m] = -(a*a*ht*ht)/(hy*hy);
                    }

                    for (unsigned int cni=0; cni<obsDeltaNodes.size(); cni++)
                    {
                        const IProblem2H2D::ExtendedSpacePointNode &odn = obsDeltaNodes.at(cni);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int s=0; s<cntPointNodes.size(); s++)
                            {
                                const IProblem2H2D::ExtendedSpacePointNode &cpn = cntPointNodes.at(s);
                                d1Y[m] += ht * mParameter.k[odn.id][cpn.id] * p[cpn.j][cpn.i] * odn.w * (cpn.w * (hx*hy));
                            }

                            for (unsigned int j=0; j<mParameter.No; j++)
                            {
                                d1Y[m] -= ht * mParameter.k[odn.id][j] * mParameter.z[odn.id][j] * odn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
                for (unsigned int m=0; m<=M; m++) p[m][n] = x1Y[m];
            }
        }
        else
        {

        }
        //--------------------------------------------------------------------------//

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p0[m][n] = ph[m][n];
                p1[m][n] = p[m][n];
            }
        }

        if (use == true) add2Info(p, info, l1/2);
        layerInfo(p, l1);
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

    rows0.clear();
    rows1.clear();
    rows2.clear();

    cols0.clear();
    cols1.clear();
    cols2.clear();

    obsDeltaNodes.clear();
    cntPointNodes.clear();

    //p3.clear();
    ph.clear();
    p1.clear();
    p0.clear();
}

void IProblem2HBackward2D::layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int layerNumber UNUSED_PARAM) const
{}

double IProblem2HBackward2D::initial1(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    //Dimension time = timeDimension();
    //double ht = time.step();
    return -2.0*UT1[sn.j][sn.i];
}

double IProblem2HBackward2D::initial2(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return 2.0*UT0[sn.j][sn.i]+mParameter.lambda*initial1(sn);
}

double IProblem2HBackward2D::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType) const
{
    return 0.0;
}

double IProblem2HBackward2D::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 0.0;
}

void IProblem2HBackward2D::add2Info(const DoubleMatrix &p, vector<ExtendedSpaceNode2DH> &info, unsigned int ln) const
{
    for (unsigned int i=0; i<mParameter.Nc; i++)
    {
        ExtendedSpaceNode2DH &pi = info[i];
        for (unsigned int r=0; r<4; r++)
        {
            for (unsigned int c=0; c<4; c++)
            {
                unsigned int x = pi.wi[r][c].i;
                unsigned int y = pi.wi[r][c].j;
                pi.wi[r][c].u[ln] = p[y][x];
            }
        }
    }
}

