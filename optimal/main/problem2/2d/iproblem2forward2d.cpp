#include "iproblem2forward2d.h"

IProblem2Forward2D::IProblem2Forward2D() {}

void IProblem2Forward2D::setSettings(P2Setting s)
{
    setting = s;
}

void IProblem2Forward2D::calculateMVD(std::vector<DoubleMatrix> &u) const
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

    for (unsigned int l=0; l<u.size(); l++) u[l].clear();
    u.clear();
    u.resize(L+1);
    DoubleMatrix uh(M+1, N+1);

    //--------------------------------------------------------------------------------------------//

    std::vector<ObservationNode> observeNodes;
    for (unsigned int j=0; j<setting.Lo; j++) extendObservationPoint(setting.xi[j], observeNodes, j);

    unsigned int *v1y = new unsigned int[M+1]; for (unsigned int m=0; m<=M; m++) v1y[m] = 0;
    unsigned int *v1x = new unsigned int[N+1]; for (unsigned int n=0; n<=N; n++) v1x[n] = 0;

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;

        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            for (unsigned int i=0; i<Lc; i++)
            {
                double _delta = delta(sn, setting.eta[i], i);

                if (checkDelta(_delta))
                {
                    if ( v1y[m] == 0 ) v1y[m] = 1;
                    if ( v1x[n] == 0 ) v1x[n] = 1;
                }
            }
        }
    }

    //------------------------------------- initial conditions -------------------------------------//
    u[0].resize(M+1,N+1);
    for (unsigned int m=0; m<=M; m++)
    {
        SpaceNodePDE sn;
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[0][m][n] = initial(sn);
        }
    }
    //IPrinter::printMatrix(u);
    //IPrinter::printSeperatorLine();
    //------------------------------------- initial conditions -------------------------------------//

    double a2_ht__hx2 = ((a*a*ht)/(hx*hx));
    double a2_ht__hy2 = ((a*a*ht)/(hy*hy));
    double lambda0_ht = lambda0*ht;
    double a2_lambda_ht__hy = (a*a*lambda*ht)/(hy);

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
    for (unsigned int l=1; l<=L; l++)
    {
        u[l].resize(M+1, N+1);

        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        {
            tn.i = l;
            tn.t = l*ht - 0.5*ht;
            clock_t t1 = clock();
            {
                SpaceNodePDE sn;
                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (v1x[n] == 0)
                    {
                        for (unsigned int m=0; m<=M; m++)
                        {
                            sn.j = m; sn.y = m*hy;

                            d1Y[m] = 2.0*u[l-1][m][n] + lambda0*theta*ht + ht*f(sn, tn);

                            if (n==0)       d1Y[m] += a2_ht__hx2*(u[l-1][m][0]   - 2.0*u[l-1][m][1]   + u[l-1][m][2]);
                            if (n>0 && n<N) d1Y[m] += a2_ht__hx2*(u[l-1][m][n-1] - 2.0*u[l-1][m][n]   + u[l-1][m][n+1]);
                            if (n==N)       d1Y[m] += a2_ht__hx2*(u[l-1][m][N-2] - 2.0*u[l-1][m][N-1] + u[l-1][m][N]);

                            if (m == 0)
                            {
                                b1Y[0] = 2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                                c1Y[0] = -2.0*a2_ht__hy2;
                                d1Y[0] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g3(sn, tn);
                            }
                            else if (m == M)
                            {
                                a1Y[M] = -2.0*a2_ht__hy2;
                                b1Y[M] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                                d1Y[M] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g4(sn, tn);
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
                }
            }
            clock_t t2 = clock();
            clock_t t5, t6;
            {
                DoubleMatrix w(cntXSize*(M+1), cntXSize*(M+1), 0.0);
                DoubleVector d(cntXSize*(M+1), 0.0);
                DoubleVector x(cntXSize*(M+1), 0.0);

                unsigned int offset = 0;
                for (unsigned int n=0; n<=N; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    if (v1x[n] == 1)
                    {
                        for (unsigned int m=0; m<=M; m++)
                        {
                            sn.j = m; sn.y = m*hy;

                            d[offset+m] = 2.0*u[l-1][m][n] + lambda0*theta*ht + ht*f(sn, tn);

                            if (n==0)       d[offset+m] += a2_ht__hx2*(u[l-1][m][0]   - 2.0*u[l-1][m][1]   + u[l-1][m][2]);
                            if (n>0 && n<N) d[offset+m] += a2_ht__hx2*(u[l-1][m][n-1] - 2.0*u[l-1][m][n]   + u[l-1][m][n+1]);
                            if (n==N)       d[offset+m] += a2_ht__hx2*(u[l-1][m][N-2] - 2.0*u[l-1][m][N-1] + u[l-1][m][N]);

                            if (m == 0)
                            {
                                w[offset+0][offset+0] += 2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                                w[offset+0][offset+1] += -2.0*a2_ht__hy2;
                                d[offset+0] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/hy)*g3(sn, tn);
                            }
                            else if (m == M)
                            {
                                w[offset+M][offset+(M-1)] += -2.0*a2_ht__hy2;
                                w[offset+M][offset+(M-0)] += 2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                                d[offset+M] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g4(sn, tn);
                            }
                            else
                            {
                                w[offset+m][offset+(m-1)] += -a2_ht__hy2;
                                w[offset+m][offset+(m+0)] += 2.0 + 2.0*a2_ht__hy2 + lambda0_ht;
                                w[offset+m][offset+(m+1)] += -a2_ht__hy2;
                            }

                            //------------------------------------- Adding delta part -------------------------------------//
                            for (unsigned int i=0; i<Lc; i++)
                            {
                                double _delta = delta(sn, setting.eta[i], i);
                                if (checkDelta(_delta))
                                {
                                    for (unsigned int s=0; s<observeNodes.size(); s++)
                                    {
                                        const ObservationNode &on = observeNodes[s];

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntXSize; cs++)
                                        {
                                            if (on.n == cntX[cs])
                                            {
                                                found = true;
                                                w[offset+m][cs*(M+1)+on.m] += -ht*setting.k[i][on.j] * _delta * on.w;
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
                                            d[offset+m] += ht*setting.k[i][on.j] * uh[on.m][on.n] * _delta * on.w;
                                        }
                                    }

                                    for (unsigned int j=0; j<Lo; j++)
                                    {
                                        d.at(offset+m) -= ht * setting.k[i][j] * setting.z[i][j] *_delta;
                                    }
                                }
                            }
                        }
                        offset += M+1;
                    }
                }

                t5 = clock();
                LinearEquation::GaussianElimination(w,d,x);
                t6 = clock();

                offset = 0;
                for (unsigned int n=0; n<=N; n++)
                {
                    if (v1x[n] == 1)
                    {
                        for (unsigned int m=0; m<=M; m++)
                        {
                            uh[m][n] = x[offset+m];
                        }
                        offset += M+1;
                    }
                }

                w.clear();
                d.clear();
                x.clear();
            }
            clock_t t3 = clock();
            printf("time: %8d %8d %f %8d %f %8d %f\n", l, t2-t1, ((float)(t2-t1))/CLOCKS_PER_SEC, t3-t2, (((float)(t3-t2))/CLOCKS_PER_SEC), t6-t5, (((float)(t6-t5))/CLOCKS_PER_SEC));
        }
        //IPrinter::printMatrix(uh);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        {
            tn.i = l;
            tn.t = l*ht;

            {
                SpaceNodePDE sn;
                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    if (v1y[m] == 0)
                    {
                        for (unsigned int n=0; n<=N; n++)
                        {
                            sn.i = n; sn.x = n*hx;

                            d1X[n] = 2.0*uh[m][n] + lambda0*theta*ht + ht*f(sn, tn);

                            if (m==0)       d1X[n] += a2_ht__hy2*(uh[0][n]   - 2.0*uh[1][n]   + uh[2][n]);
                            if (m>0 && m<M) d1X[n] += a2_ht__hy2*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]);
                            if (m==M)       d1X[n] += a2_ht__hy2*(uh[M-2][n] - 2.0*uh[M-1][n] + uh[M][n]);

                            if (n == 0)
                            {
                                b1X[0] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + (2.0*a*a*lambda*ht)/(hx);
                                c1X[0] = -2.0*a2_ht__hx2;
                                d1X[0] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/hx)*g1(sn, tn);
                            }
                            else if (n == N)
                            {
                                a1X[N] = -2.0*a2_ht__hx2;
                                b1X[N] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + (2.0*a*a*lambda*ht)/(hx);
                                d1X[N] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/(hx))*g2(sn, tn);
                            }
                            else
                            {
                                a1X[n] = -a2_ht__hx2;
                                b1X[n] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
                                c1X[n] = -a2_ht__hx2;
                            }
                        }
                        tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                        for (unsigned int n=0; n<=N; n++) u[l][m][n] = x1X[n];
                    }
                }

            }

            {
                DoubleMatrix w(cntYSize*(N+1), cntYSize*(N+1), 0.0);
                DoubleVector d(cntYSize*(N+1), 0.0);
                DoubleVector x(cntYSize*(N+1), 0.0);

                unsigned int offset = 0;
                for (unsigned int m=0; m<=M; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    if (v1y[m] == 1)
                    {
                        for (unsigned int n=0; n<=N; n++)
                        {
                            sn.i = n; sn.x = n*hx;

                            d[offset+n] = 2.0*uh[m][n] + lambda0*theta*ht + ht*f(sn, tn);

                            if (m==0)       d[offset+n] += a2_ht__hy2*(uh[0][n]   - 2.0*uh[1][n]   + uh[2][n]);
                            if (m>0 && m<M) d[offset+n] += a2_ht__hy2*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]);
                            if (m==M)       d[offset+n] += a2_ht__hy2*(uh[M-2][n] - 2.0*uh[M-1][n] + uh[M][n]);

                            if (n == 0)
                            {
                                w[offset+0][offset+0] += 2.0 + 2.0*a2_ht__hx2 + lambda0_ht + (2.0*a*a*lambda*ht)/(hx);
                                w[offset+0][offset+1] += -2.0*a2_ht__hx2;
                                d[offset+0] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/hx)*g1(sn, tn);
                            }
                            else if (n == N)
                            {
                                w[offset+N][offset+(N-1)] += -2.0*a2_ht__hx2;
                                w[offset+N][offset+(N-0)] += 2.0 + 2.0*a2_ht__hx2 + lambda0_ht + (2.0*a*a*lambda*ht)/(hx);
                                d[offset+N] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/(hx))*g2(sn, tn);
                            }
                            else
                            {
                                w[offset+n][offset+(n-1)] += -a2_ht__hx2;
                                w[offset+n][offset+(n+0)] += 2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
                                w[offset+n][offset+(n+1)] += -a2_ht__hx2;
                            }

                            //------------------------------------ Adding delta part -------------------------------------//
                            for (unsigned int i=0; i<Lc; i++)
                            {
                                double _delta = delta(sn, setting.eta[i], i);
                                if (checkDelta(_delta))
                                {
                                    for (unsigned int s=0; s<observeNodes.size(); s++)
                                    {
                                        const ObservationNode &on = observeNodes[s];

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntYSize; cs++)
                                        {
                                            if (on.m == cntY[cs])
                                            {
                                                found = true;
                                                w[offset+n][cs*(N+1)+on.n] += -ht*setting.k[i][on.j] * _delta * on.w;
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
                                            d[offset+n] += ht*setting.k[i][on.j] * u[l][on.m][on.n] * _delta * on.w;
                                        }
                                    }

                                    for (unsigned int j=0; j<Lo; j++)
                                    {
                                        d[offset+n] -= ht * setting.k[i][j] * setting.z[i][j] *_delta;
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
                            u[l][m][n] = x[offset+n];
                        }
                        offset += N+1;
                    }
                }

                w.clear();
                d.clear();
                x.clear();
            }
        }
        //IPrinter::printMatrix(u);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to x direction conditions -------------------------------------//
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
    observeNodes.clear();
    uh.clear();
}

double IProblem2Forward2D::initial(const SpaceNodePDE &sn) const
{
    double x = sn.x;
    double y = sn.y;
    return x*x + y*y;
}

double IProblem2Forward2D::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double IProblem2Forward2D::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double y = sn.y;
    double t = tn.t;

    double res = 1.0 - 4.0*setting.a*setting.a + setting.lambda0*(U(x,y,t) - setting.theta);

    double W = 0.0;
    for (unsigned int i=0; i<setting.Lc; i++)
    {
        double _delta = delta(sn, setting.eta[i], i, 0);
        if (checkDelta(_delta))
        {
            double vi = 0.0;
            for (unsigned int j=0; j<setting.Lo; j++)
            {
                vi += setting.k[i][j] * ( U(setting.xi[j].x, setting.xi[j].y, t) - setting.z[i][j]);
            }
            W += vi * _delta;
        }
    }
    res -= W;

    return res;
}

bool IProblem2Forward2D::checkDelta(double _delta) const
{
    return _delta != 0.0;
}

double IProblem2Forward2D::delta(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i UNUSED_PARAM, unsigned int source UNUSED_PARAM) const
{
    return delta4(sn, eta, i);
}

double IProblem2Forward2D::delta1(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();

    double res = 0.0;
    if ( sn.i == eta.i && sn.j == eta.j ) res = 1.0/(hx*hy);
    return res;
}

double IProblem2Forward2D::delta2(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();
    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double sigmaX = hx;
    double sigmaY = hy;

    unsigned int rx = (unsigned int)(round(eta.x*Nx));
    unsigned int ry = (unsigned int)(round(eta.y*Ny));

    double res = 0.0;
    if (rx-3 <= sn.i && sn.i <= rx+3 && ry-3 <= sn.j && sn.j <= ry+3)
    {
        res = (1.0/(2.0*M_PI*sigmaX*sigmaY))*exp(-0.5*(((sn.x-eta.x)*(sn.x-eta.x))/(sigmaX*sigmaX)+((sn.y-eta.y)*(sn.y-eta.y))/(sigmaY*sigmaY)));
    }

    return res;
}

double IProblem2Forward2D::delta3(const SpaceNodePDE &sn UNUSED_PARAM, const SpaceNodePDE &eta UNUSED_PARAM, unsigned int) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();

    double resX = 0.0;
    {
        double dhx = fabs(sn.x-eta.x);
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
        double dhy = fabs(sn.y-eta.y);
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

double IProblem2Forward2D::delta4(const SpaceNodePDE &sn UNUSED_PARAM, const SpaceNodePDE &eta UNUSED_PARAM, unsigned int) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();

    double resX = 0.0;
    {
        double dhx = fabs(sn.x-eta.x);
        if ( dhx <= hx )
        {
            resX = (hx-dhx) / hx;
        }
    }

    double resY = 0.0;
    {
        double dhy = fabs(sn.y-eta.y);
        if ( dhy <= hy )
        {
            resY = (hy-dhy) / hy;
        }
    }

    double res = (resX*resY)/(hx*hy);

    return res;
}

void IProblem2Forward2D::extendObservationPoint1(const SpaceNodePDE xi, std::vector<ObservationNode> &ons, unsigned int j) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);

    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double hx = xd.step();
    double hy = yd.step();

    unsigned int rx = (unsigned int)(round(xi.x*Nx));
    unsigned int ry = (unsigned int)(round(xi.y*Ny));

    ObservationNode on;
    if ( rx*hx <= xi.x && ry*hy <= xi.y ) // left bottom
    {
        on.n = rx+0; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx+1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx+0; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx+1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        return;
    }

    if ( rx*hx <= xi.x && ry*hy >= xi.y ) // left top
    {
        on.n = rx+0; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx+1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx+0; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx+1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        return;
    }

    if ( rx*hx >= xi.x && ry*hy <= xi.y ) // right bottom
    {
        on.n = rx+0; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx-1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx+0; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx-1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        return;
    }

    if ( rx*hx >= xi.x && ry*hy >= xi.y ) // right top
    {
        on.n = rx+0; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx-1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx+0; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        on.n = rx-1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy); ons.push_back(on);

        return;
    }
}

void IProblem2Forward2D::extendObservationPoint(const SpaceNodePDE xi, std::vector<ObservationNode> &ons, unsigned int j) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);

    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double hx = xd.step();
    double hy = yd.step();

    unsigned int rx = (unsigned int)(round(xi.x*Nx));
    unsigned int ry = (unsigned int)(round(xi.y*Ny));

    double hx3 = hx*hx*hx;
    double hx32 = (1.0/(2.0*hx3));
    double hx36 = (1.0/(6.0*hx3));

    double hy3 = hy*hy*hy;
    double hy32 = (1.0/(2.0*hy3));
    double hy36 = (1.0/(6.0*hy3));

    ObservationNode on;
    double dx = 0.0;
    double dy = 0.0;

    if ( rx*hx <= xi.x && ry*hy <= xi.y ) // left bottom
    {
        on.n = rx-1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        /////////////////////////
        on.n = rx+0; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);;
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        /////////////////////////
        on.n = rx+1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        ////////////////////////
        on.n = rx+2; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+2; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+2; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+2; on.m = ry+2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);
        return;
    }

    if ( rx*hx <= xi.x && ry*hy >= xi.y ) // left top
    {
        on.n = rx-1; on.m = ry-2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx-1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        /////////////////////////
        on.n = rx+0; on.m = ry-2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+0; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);;
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        /////////////////////////
        on.n = rx+1; on.m = ry-2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        ////////////////////////
        on.n = rx+2; on.m = ry-2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+2; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+2; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+2; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);
        return;
    }

    if ( rx*hx >= xi.x && ry*hy >= xi.y ) // right top
    {
        on.n = rx-2; on.m = ry-2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx-2; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-2; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-2; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        /////////////////////////
        on.n = rx-1; on.m = ry-2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx-1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);;
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        /////////////////////////
        on.n = rx+0; on.m = ry-2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+0; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);


        ////////////////////////
        on.n = rx+1; on.m = ry-2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);
        return;
    }

    if ( rx*hx >= xi.x && ry*hy <= xi.y ) // right bottom
    {
        on.n = rx-2; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx-2; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-2; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-2; on.m = ry+2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        /////////////////////////
        on.n = rx-1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx-1; on.m = ry+2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);;
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        /////////////////////////
        on.n = rx+0; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+0; on.m = ry+2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);


        ////////////////////////
        on.n = rx+1; on.m = ry-1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+0; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+1; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        ons.push_back(on);

        on.n = rx+1; on.m = ry+2; on.x = on.n*hx; on.y = on.m*hy; on.xi = xi; on.j = j;
        dx = fabs(on.x-xi.x);
        dy = fabs(on.y-xi.y);
        on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        ons.push_back(on);
        return;
    }
}

double IProblem2Forward2D::g1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    return setting.lambda*(y*y+t - setting.theta);
}

double IProblem2Forward2D::g2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    //return 2.0 + lambda*(U(1.0, y, t) - theta);
    return 2.0 + setting.lambda*(1.0 + y*y + t - setting.theta);
}

double IProblem2Forward2D::g3(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return setting.lambda*(x*x + t - setting.theta);
}

double IProblem2Forward2D::g4(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return 2.0 + setting.lambda*(1.0 + x*x + t - setting.theta);
}

double IProblem2Forward2D::U(double x, double y, double t) const
{
    return x*x + y*y + t;
}

