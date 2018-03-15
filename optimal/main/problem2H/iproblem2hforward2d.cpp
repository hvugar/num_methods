#include "iproblem2hforward2d.h"

void IProblem2HForward2D::Main(int arg UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    IProblem2HForward2D p;
    p.addSpaceDimension(Dimension(0.01, 0, 100));
    p.addSpaceDimension(Dimension(0.01, 0, 100));
    p.setTimeDimension(Dimension(0.001, 0, 1000));

    IProblem2H::Parameter prm;
    prm.Ns = 2;
    prm.q.resize(prm.Ns);
    prm.q[0] = 1000.1; prm.q[1] = 2000.2;
    prm.theta.resize(prm.Ns);
    prm.Nc = 2;
    prm.No = 2;
    prm.eta.resize(prm.Nc);
    prm.eta[0].x = 0.3; prm.eta[0].y = 0.7;
    prm.eta[1].x = 0.7; prm.eta[1].y = 0.3;

    prm.xi.resize(prm.No);
    prm.xi[0].x = 0.4; prm.xi[0].y = 0.4;
    prm.xi[1].x = 0.6; prm.xi[1].y = 0.6;

    prm.k.resize(prm.Nc, prm.No, 0.1);
    prm.z.resize(prm.Nc, prm.No, 0.0);

    prm.a = 1.0;
    mParameter.lambda = 0.01;

    p.mParameter = prm;

    DoubleMatrix u;
    p.calculateMVD(u);
}

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

    double lambda = mParameter.lambda;
    double a = mParameter.a;

    DoubleMatrix u0;
    u0.resize(M+1, N+1);

    DoubleMatrix u1;
    u1.resize(M+1, N+1);

    DoubleMatrix u2;
    u2.resize(M+1, N+1);

    DoubleMatrix u3;
    u3.resize(M+1, N+1);

    u.clear();
    u.resize(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    std::vector<IProblem2H::ObservationPointNode> observeNodes;
    for (unsigned int j=0; j<mParameter.No; j++) extendObservationPoint(mParameter.xi[j], observeNodes, j);

    std::vector<IProblem2H::ControlDeltaNode> cndeltaNodes;
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
            const IProblem2H::ControlDeltaNode &cdn = cndeltaNodes.at(i);
            if (cdn.j == ny)
            {
                found1 = true;
                for (unsigned int j=0; j<observeNodes.size(); j++)
                {
                    const IProblem2H::ObservationPointNode &on = observeNodes.at(j);
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
            const IProblem2H::ControlDeltaNode &cdn = cndeltaNodes.at(i);
            if (cdn.i == nx)
            {
                found1 = true;
                for (unsigned int j=0; j<observeNodes.size(); j++)
                {
                    const IProblem2H::ObservationPointNode &on = observeNodes.at(j);
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
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u0[m][n] = initial1(sn);
            u1[m][n] = u0[m][n] + initial2(sn)*ht;

            if (m==20 && n==20) u1[m][n] += -(mParameter.q[0]*(1.0/ht))*((ht*ht)/2.0);
            if (m==80 && n==80) u1[m][n] += -(mParameter.q[1]*(1.0/ht))*((ht*ht)/2.0);
        }
    }
    IPrinter::printSeperatorLine("\n");
    IPrinter::printMatrix(12, 6, u1);
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

    TimeNodePDE tn;

    for (unsigned int l=2; l<=(L/2); l+=2)
    {
        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int row=0; row<rows0.size(); row++)
        {
            unsigned int m = rows0.at(row);
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;

                d1X[n] = (2.0+2.0*lambda*ht)*u1[m][n] - (1.0+(lambda*ht)/2.0)*u0[m][n];

                if (m == 0)     d1X[n] += ((a*a*ht*ht)/(hy*hy))*(u1[0][n]   - 2.0*u1[1][n]   + u1[2][n]);
                if (m>0 && m<M) d1X[n] += ((a*a*ht*ht)/(hy*hy))*(u1[m-1][n] - 2.0*u1[m][n]   + u1[m+1][n]);
                if (m == M)     d1X[n] += ((a*a*ht*ht)/(hy*hy))*(u1[M-2][n] - 2.0*u1[M-1][n] + u1[M][n]);

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

                    d1X[n] = (2.0+2.0*lambda*ht)*u1[m][n] - (1.0+(lambda*ht)/2.0)*u0[m][n];

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
                        const IProblem2H::ControlDeltaNode &cdn = cndeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int s=0; s<observeNodes.size(); s++)
                            {
                                const IProblem2H::ObservationPointNode &on = observeNodes[s];
                                d1X[n] += ht*ht * mParameter.k[cdn.id][on.id] * u2[on.j][on.i] * cdn.w * on.w;
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

        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        tn.i = l+1; tn.t = tn.i*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int col=0; col<cols0.size(); col++)
        {
            unsigned int n = cols0.at(col);
            sn.i = n; sn.x = n*hx;
            for (unsigned int m=0; m<=M; m++)
            {
                sn.j = m; sn.y = m*hy;

                d1Y[m] = (2.0+2.0*lambda*ht)*u2[m][n] - (1.0+(lambda*ht)/2.0)*u1[m][n];

                if (n==0)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(u2[m][0]   - 2.0*u2[m][1]   + u2[m][2]);
                if (n>0 && n<N) d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(u2[m][n-1] - 2.0*u2[m][n]   + u2[m][n+1]);
                if (n==N)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(u2[m][N-2] - 2.0*u2[m][N-1] + u2[m][N]);

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
            for (unsigned int m=0; m<=M; m++) u3[m][n] = x1Y[m];
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

                    d1Y[m] = (2.0+2.0*lambda*ht)*u2[m][n] - (1.0+(lambda*ht)/2.0)*u1[m][n];

                    if (n==0)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(u2[m][0]   - 2.0*u2[m][1]   + u2[m][2]);
                    if (n>0 && n<N) d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(u2[m][n-1] - 2.0*u2[m][n]   + u2[m][n+1]);
                    if (n==N)       d1Y[m] += ((a*a*ht*ht)/(hx*hx))*(u2[m][N-2] - 2.0*u2[m][N-1] + u2[m][N]);

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

                    for (unsigned int cni=0; cni<cndeltaNodes.size(); cni++)
                    {
                        const IProblem2H::ControlDeltaNode &cdn = cndeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int s=0; s<observeNodes.size(); s++)
                            {
                                const IProblem2H::ObservationPointNode &on = observeNodes.at(s);
                                d1Y[m] += ht * mParameter.k[cdn.id][on.id] * u3[on.j][on.i] * cdn.w * on.w;
                            }

                            for (unsigned int j=0; j<mParameter.No; j++)
                            {
                                d1Y[m] -= ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                }
                tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
                for (unsigned int m=0; m<=M; m++) u3[m][n] = x1Y[m];
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
                u0[m][n] = u2[m][n];
                u1[m][n] = u3[m][n];
            }
        }

        IPrinter::printSeperatorLine();
        printf(("%d\n"), l);
        IPrinter::printMatrix(12, 6, u1);
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            u[m][n] = u1[m][n];
        }
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

    observeNodes.clear();
    cndeltaNodes.clear();

    u3.clear();
    u2.clear();
    u1.clear();
    u0.clear();
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

void IProblem2HForward2D::extendObservationPoint(const SpacePoint &xi, std::vector<IProblem2H::ObservationPointNode> &ons, unsigned int j) const
{
    double hx = spaceDimension(Dimension::DimensionX).step();
    double hy = spaceDimension(Dimension::DimensionY).step();

    std::vector<IProblem2H::ExtendedSpacePointNode> extpoint;
    distributeDelta(xi, extpoint, j);

    for (unsigned int i=0; i<extpoint.size(); i++)
    {
        IProblem2H::ExtendedSpacePointNode ep = extpoint.at(i);
        IProblem2H::ObservationPointNode node;
        node.id = ep.id; node.w = ep.w * (hx*hy); node.pt = ep.pt;
        node.i = ep.i; node.x = ep.x;
        node.j = ep.j; node.y = ep.y;
        ons.push_back(node);
    }
    extpoint.clear();
}

void IProblem2HForward2D::extendContrlDeltaPoint(const SpacePoint &eta, std::vector<IProblem2H::ControlDeltaNode> &cps, unsigned int id) const
{
    std::vector<IProblem2H::ExtendedSpacePointNode> extpoint;
    distributeDelta(eta, extpoint, id);

    for (unsigned int i=0; i<extpoint.size(); i++)
    {
        IProblem2H::ExtendedSpacePointNode ep = extpoint.at(i);
        IProblem2H::ControlDeltaNode node;
        node.id = ep.id; node.w = ep.w; node.pt = ep.pt;
        node.i = ep.i; node.x = ep.x;
        node.j = ep.j; node.y = ep.y;
        cps.push_back(node);
    }
    extpoint.clear();
}

void IProblem2HForward2D::distributeDelta(const SpacePoint &pt, std::vector<IProblem2H::ExtendedSpacePointNode> &nodes, unsigned int id) const
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
            IProblem2H::ExtendedSpacePointNode node;
            node.i = n; node.x = n*hx; node.j = m; node.y = m*hy; node.pt = pt; node.id = id;
            node.w = factor*exp(-0.5*(((node.x-pt.x)*(node.x-pt.x))/(sigmaX*sigmaX)+((node.y-pt.y)*(node.y-pt.y))/(sigmaY*sigmaY)));
            nodes.push_back(node);
        }
    }
}
