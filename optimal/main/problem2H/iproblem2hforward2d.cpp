#include "iproblem2hforward2d.h"

void IProblem2HForward2D::calculateMVD(DoubleMatrix &u) const
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

    double lambda = 0.01;
    double a = 1.0;

    DoubleMatrix u0;
    u0.resize(M+1, N+1);

    DoubleMatrix u1;
    u1.resize(M+1, N+1);

    DoubleMatrix u2;
    u2.resize(M+1, N+1);

    u.clear();
    u.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<ObservationPointNode> observeNodes;
    for (unsigned int j=0; j<mParameter.No; j++) extendObservationPoint(mParameter.xi[j], observeNodes, j);

    std::vector<ControlDeltaNode> cndeltaNodes;
    for (unsigned int i=0; i<mParameter.Nc; i++) extendContrlDeltaPoint(mParameter.eta[i], cndeltaNodes, i);

    SpaceNodePDE sn;

    //--------------------------------------------------------------------------------------------//
    vector<unsigned int> rows0;
    vector<unsigned int> rows1;
    vector<unsigned int> rows2;
    for (unsigned int ny=0; ny<=M; ny++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cndeltaNodes.size(); i++)
        {
            const ControlDeltaNode &cdn = cndeltaNodes.at(i);
            if (cdn.j == ny)
            {
                found1 = true;
                for (unsigned int j=0; j<observeNodes.size(); j++)
                {
                    const ObservationPointNode &on = observeNodes.at(j);
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
        for (unsigned int i=0; i<cndeltaNodes.size(); i++)
        {
            const ControlDeltaNode &cdn = cndeltaNodes.at(i);
            if (cdn.i == nx)
            {
                found1 = true;
                for (unsigned int j=0; j<observeNodes.size(); j++)
                {
                    const ObservationPointNode &on = observeNodes.at(j);
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

    //------------------------------------- initial conditions -------------------------------------//
    for (unsigned int ny=0; ny<=M; ny++)
    {
        sn.j = ny; sn.y = ny*hy;
        for (unsigned int nx=0; nx<=N; nx++)
        {
            sn.i = nx; sn.x = nx*hx;
            u0[ny][nx] = initial1(sn);
            u1[ny][nx] = u0[ny][nx] + initial2(sn)*ht;
        }
    }

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

    TimeNodePDE tn;

    for (unsigned int nt=2; nt<=L; nt+2)
    {
        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        tn.i = nt; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;

                d1X[n] = (2.0+2.0*lambda*ht)*u1[m][n] - (1.0+(lambda*ht)/2.0)*u2[m][n];

                if (m == 0)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[0][n]  -2.0*u1[1][n]  +u1[2][n]);
                if (m>0 && m<M) d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[m-1][n]-2.0*u1[m][n]  +u1[m+1][n]);
                if (m == M)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[M-2][n]-2.0*u1[M-1][n]+u1[M][n]);

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
            }
            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
            for (unsigned int n=0; n<=N; n++) u2[m][n] = x1X[n];
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

                    d1X[n] = (2.0+2.0*lambda*ht)*u1[m][n] - (1.0+(lambda*ht)/2.0)*u2[m][n];

                    if (m == 0)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[0][n]  -2.0*u1[1][n]  +u1[2][n]);
                    if (m>0 && m<M) d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[m-1][n]-2.0*u1[m][n]  +u1[m+1][n]);
                    if (m == M)     d1X[n] = ((a*a*ht*ht)/(hy*hy))*(u1[M-2][n]-2.0*u1[M-1][n]+u1[M][n]);

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

                    for (unsigned int cni=0; cni<cndeltaNodes.size(); cni++)
                    {
                        const ControlDeltaNode &cdn = cndeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int s=0; s<observeNodes.size(); s++)
                            {
                                const ObservationPointNode &on = observeNodes[s];
                                d1X[n] += ht*ht * mParameter.k[cdn.id][on.id] * u[on.j][on.i] * cdn.w * on.w;
                            }

                            for (unsigned int j=0; j<mParameter.No; j++)
                            {
                                d1X[n] -= ht*ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                for (unsigned int n=0; n<=N; n++) u2[m][n] = x1X[n];
            }
        }
        else
        {

        }
        //--------------------------------------------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = nt+1; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int col=0; col<cols0.size(); col++)
        {
            unsigned int n = cols0[col];
            sn.i = n; sn.x = n*hx;
            for (unsigned int m=0; m<=M; m++)
            {
                sn.j = m; sn.y = m*hy;

                d1Y[m] = 2.0*u[m][n] + lambda0*theta*ht;
#ifdef USE_ADDITIONAL_FUNCTIONS
                d1Y[m] += ht*f(sn, tn);
#endif
                if (n==0)       d1Y[m] += a2_ht__hx2*(u.at(m,0)   - 2.0*u.at(m,1)   + u.at(m,2));
                if (n>0 && n<N) d1Y[m] += a2_ht__hx2*(u.at(m,n-1) - 2.0*u.at(m,n)   + u.at(m,n+1));
                if (n==N)       d1Y[m] += a2_ht__hx2*(u.at(m,N-2) - 2.0*u.at(m,N-1) + u.at(m,N));

                if (m == 0)
                {
                    a1Y[0] = 0.0;
                    b1Y[0] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                    c1Y[0] = -2.0*a2_ht__hy2;
                    d1Y[0] += 2.0*a2_lambda_ht__hy*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                    d1Y[0] += ((2.0*a*a*ht)/(hy))*g3(sn, tn);
#endif
                }
                else if (m == M)
                {
                    a1Y[M] = -2.0*a2_ht__hy2;
                    b1Y[M] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                    c1Y[M] = 0.0;
                    d1Y[M] += 2.0*a2_lambda_ht__hy*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                    d1Y[M] += ((2.0*a*a*ht)/(hy))*g4(sn, tn);
#endif
                }
                else
                {
                    a1Y[m] = -a2_ht__hy2;
                    b1Y[m] = 2.0 + 2.0*a2_ht__hy2 + lambda0_ht;
                    c1Y[m] = -a2_ht__hy2;
                }
            }
            tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
            for (unsigned int m=0; m<=M; m++) uh[m][n] = x1Y[m];
        }
        //--------------------------------------------------------------------------//
    }
}

void IProblem2HForward2D::layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const
{
    C_UNUSED(u);
    C_UNUSED(layerNumber);
}

double IProblem2HForward2D::initial1(const SpaceNodePDE &) const
{
    return 0.0;
}

double IProblem2HForward2D::initial2(const SpaceNodePDE &) const
{
    return 0.0;
}

double IProblem2HForward2D::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const
{
    return 0.0;
}

double IProblem2HForward2D::f(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double IProblem2HForward2D::a(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

void IProblem2HForward2D::extendObservationPoint(const SpacePoint &xi, std::vector<ObservationPointNode> &ons, unsigned int j) const
{
    double hx = spaceDimension(Dimension::DimensionX).step();
    double hy = spaceDimension(Dimension::DimensionY).step();

    std::vector<ExtendedSpacePointNode> extpoint;
    distributeDelta(xi, extpoint, j);

    for (unsigned int i=0; i<extpoint.size(); i++)
    {
        ExtendedSpacePointNode ep = extpoint.at(i);
        ObservationPointNode node;
        node.id = ep.id; node.w = ep.w * (hx*hy); node.pt = ep.pt;
        node.i = ep.i; node.x = ep.x;
        node.j = ep.j; node.y = ep.y;
        ons.push_back(node);
    }
    extpoint.clear();
}

void IProblem2HForward2D::extendContrlDeltaPoint(const SpacePoint &eta, std::vector<ControlDeltaNode> &cps, unsigned int id) const
{
    std::vector<ExtendedSpacePointNode> extpoint;
    distributeDelta(eta, extpoint, id);

    for (unsigned int i=0; i<extpoint.size(); i++)
    {
        ExtendedSpacePointNode ep = extpoint.at(i);
        ControlDeltaNode node;
        node.id = ep.id; node.w = ep.w; node.pt = ep.pt;
        node.i = ep.i; node.x = ep.x;
        node.j = ep.j; node.y = ep.y;
        cps.push_back(node);
    }
    extpoint.clear();
}

void IProblem2HForward2D::distributeDelta(const SpacePoint &pt, std::vector<ExtendedSpacePointNode> &nodes, unsigned int id) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();
    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double sigmaX = hx;
    double sigmaY = hy;

    unsigned int rx = (unsigned int)(round(pt.x*Nx));
    unsigned int ry = (unsigned int)(round(pt.y*Ny));

    unsigned int k=3;

    double sumX = 0.0;
    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        sumX += exp(-((n*hx-pt.x)*(n*hx-pt.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=ry-k; m<=ry+k; m++)
    {
        sumY += exp(-((m*hy-pt.y)*(m*hy-pt.y))/(2.0*sigmaY*sigmaY));
    }
    sumY *= hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/((2.0*M_PI)*sigma);

    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        for (unsigned int m=ry-k; m<=ry+k; m++)
        {
            ExtendedSpacePointNode node;
            node.i = n; node.x = n*hx; node.j = m; node.y = m*hy; node.pt = pt; node.id = id;
            node.w = factor*exp(-0.5*(((node.x-pt.x)*(node.x-pt.x))/(sigmaX*sigmaX)+((node.y-pt.y)*(node.y-pt.y))/(sigmaY*sigmaY)));
            nodes.push_back(node);
        }
    }
}
