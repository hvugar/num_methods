#include "iproblem2backward2d.h"

IProblem2Backward2D::~IProblem2Backward2D()
{}

void IProblem2Backward2D::calculateMVD(DoubleMatrix &p, vector<ExtendedSpaceNode2D> &info, bool use)
{
    //puts("IProblem2Backward2D::calculateMVD...");
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
            ExtendedSpaceNode2D &e = info[i];
            e.setSpaceNode(mParameter.eta[i]);
            e.id = i;
            e.extendWeights(xd, yd, 4, 4);
            e.extendLayers(L+1);
        }
    }

    p.clear();
    p.resize(M+1, N+1);

    DoubleMatrix ph(M+1, N+1);

    //--------------------------------------------------------------------------------------------//

    std::vector<ControlNode> controlNodes;
    for (unsigned int i=0; i<mParameter.Lc; i++) extendControlPoint0(mParameter.eta[i], controlNodes, i);

#ifdef USE_B_VARIANT_2
    std::vector<ObservationDeltaNode> obdeltaNodes;
    for (unsigned int j=0; j<mParameter.Lo; j++) extendObservationDeltaPoint0(mParameter.xi[j], obdeltaNodes, j);
#endif

    unsigned int *v1y = new unsigned int[M+1]; for (unsigned int m=0; m<=M; m++) v1y[m] = 0;
    unsigned int *v1x = new unsigned int[N+1]; for (unsigned int n=0; n<=N; n++) v1x[n] = 0;

    SpaceNodePDE sn;

#ifdef USE_B_VARIANT_1
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
#endif

#ifdef USE_B_VARIANT_2
    for (unsigned int j=0; j<obdeltaNodes.size(); j++)
    {
        const ObservationDeltaNode &odn = obdeltaNodes.at(j);
        v1x[odn.i] = 1;
        v1y[odn.j] = 1;
    }
#endif

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

    if (use ==true)
    {
        for (unsigned int i=0; i<Lc; i++)
        {
            ExtendedSpaceNode2D &pi = info[i];
            for (unsigned int r=0; r<4; r++)
            {
                for (unsigned int c=0; c<4; c++)
                {
                    pi.wi[r][c].u[L] = p[pi.wi[r][c].j][pi.wi[r][c].i];
                }
            }
        }
    }

    layerInfo(p, L);
    //IPrinter::printMatrix(p[L]);
    //IPrinter::printSeperatorLine();
    //------------------------------------- initial conditions -------------------------------------//

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
        //puts("IProblem2Backward2D::calculateMVD.y-->...");
        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        {
            tn.i = l;
            tn.t = l*ht + 0.5*ht;
            //--------------------------------------------------------------------------//
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
            //--------------------------------------------------------------------------//
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
#ifdef USE_B_VARIANT_1
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
                                            if (cn.i == cntX[cs])
                                            {
                                                found = true;
                                                w2[offset+m][cs*(M+1)+cn.j] += -ht * mParameter.k[cn.id][j] * _delta * cn.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+m] += ht * mParameter.k[cn.id][j] * ph[cn.j][cn.i] * _delta * cn.w;
                                        }
                                    }

//                                    for (unsigned int i=0; i<mParameter.Lc; i++)
//                                    {
//                                        d2[offset+m] -= 2.0 * r * ht *  mParameter.k[i][odn.id] * penalty(i, tn) * odn.w;
//                                    }
                                }
                            }
#endif
                            //------------------------------------- Adding delta part -------------------------------------//
#ifdef USE_B_VARIANT_2
                            for (unsigned int onj=0; onj<obdeltaNodes.size(); onj++)
                            {
                                const ObservationDeltaNode &odn = obdeltaNodes.at(onj);
                                if (odn.i == sn.i && odn.j == sn.j)
                                {
                                    for (unsigned int s=0; s<controlNodes.size(); s++)
                                    {
                                        const ControlNode &cn = controlNodes[s];

                                        //if (cn.n == odn.n)
                                        //{
                                        //    w2[offset+m][odn.n*(M+1)+cn.m] += -ht * mParameter.k[cn.i][odn.j] * odn.w * cn.w;
                                        //}
                                        //else
                                        //{
                                        //    d2[offset+m] += ht * mParameter.k[cn.i][odn.j] * ph[cn.m][cn.n] * odn.w * cn.w;
                                        //}

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntXSize; cs++)
                                        {
                                            if (cn.i == cntX[cs])
                                            {
                                                found = true;
                                                w2[offset+m][cs*(M+1)+cn.j] += -ht * mParameter.k[cn.id][odn.id] * odn.w * cn.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+m] += ht * mParameter.k[cn.id][odn.id] * ph[cn.j][cn.i] * odn.w * cn.w;
                                        }
                                    }

                                    for (unsigned int i=0; i<Lc; i++)
                                    {
                                        d2[offset+m] -= 2.0 * r * ht *  mParameter.k[i][odn.id] * penalty(i, tn) * odn.w;
                                    }
                                }
                            }
#endif
                            //------------------------------------- Adding delta part -------------------------------------//
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
            //--------------------------------------------------------------------------//
        }
        //puts("IProblem2Backward2D::calculateMVD.y-->.");
        //IPrinter::printMatrix(ph);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        //puts("IProblem2Backward2D::calculateMVD.x-->...");
        {
            tn.i = l;
            tn.t = l*ht;
            //--------------------------------------------------------------------------//
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
            //--------------------------------------------------------------------------//
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
#ifdef USE_B_VARIANT_1
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
                                            if (cn.j == cntY[cs])
                                            {
                                                found = true;
                                                w2[offset+n][cs*(N+1)+cn.i] += -ht * mParameter.k[cn.id][j] * _delta * cn.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+n] += ht * mParameter.k[cn.id][j] * p[cn.j][cn.i] * _delta * cn.w;
                                        }
                                    }

                                    //
                                    //
                                    //
                                    //
                                }
                            }
#endif
                            //------------------------------------- Adding delta part -------------------------------------//
#ifdef USE_B_VARIANT_2
                            for (unsigned int odj=0; odj<obdeltaNodes.size(); odj++)
                            {
                                const ObservationDeltaNode &odn = obdeltaNodes.at(odj);
                                if (odn.i == sn.i && odn.j == sn.j)
                                {
                                    for (unsigned int s=0; s<controlNodes.size(); s++)
                                    {
                                        const ControlNode &cn = controlNodes[s];

                                        //if (cn.m == odn.m)
                                        //{
                                        //    w2[offset+n][odn.m*(N+1)+cn.n] += -ht * mParameter.k[cn.i][odn.j] * odn.w * cn.w;
                                        //}
                                        //else
                                        //{
                                        //    d2[offset+n] += ht * mParameter.k[cn.i][odn.j] * p[cn.m][cn.n] * odn.w * cn.w;
                                        //}

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntYSize; cs++)
                                        {
                                            if (cn.j == cntY[cs])
                                            {
                                                found = true;
                                                w2[offset+n][cs*(N+1)+cn.i] += -ht * mParameter.k[cn.id][odn.id] * cn.w * odn.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+n] += ht * mParameter.k[cn.id][odn.id] * p[cn.j][cn.i] * cn.w * odn.w;
                                        }
                                    }

                                    for (unsigned int i=0; i<Lc; i++)
                                    {
                                        d2[offset+n] -= 2.0 * r * ht *  mParameter.k[i][odn.id] * penalty(i, tn) * odn.w;
                                    }

                                }
                            }
#endif
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
            //--------------------------------------------------------------------------//
        }
        //puts("IProblem2Backward2D::calculateMVD.x-->.");
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
#ifdef USE_B_VARIANT_2
    obdeltaNodes.clear();
#endif
    ph.clear();
    //puts("IProblem2Backward2D::calculateMVD.");
}

void IProblem2Backward2D::setPenaltyCoefficient(double r)
{
    this->r = r;
}

double IProblem2Backward2D::penaltyCoefficient() const
{
    return r;
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

void IProblem2Backward2D::extendControlPoint0(const SpaceNodePDE &eta, std::vector<ControlNode> &cns, unsigned int i) const
{
#ifdef APPROX_B1_1
    extendControlPoint1(eta, cns, i);
#endif
#ifdef APPROX_B1_2
    extendControlPoint2(eta, cns, i);
#endif
#ifdef APPROX_B1_3
    extendControlPoint3(eta, cns, i);
#endif
}

void IProblem2Backward2D::extendControlPoint1(const SpaceNodePDE &eta, std::vector<ControlNode> &cns, unsigned int i) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);

    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double hx = xd.step();
    double hy = yd.step();

    unsigned int rx = (unsigned int)(floor(eta.x*Nx));
    unsigned int ry = (unsigned int)(floor(eta.y*Ny));

    ControlNode cn;

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    cn.w = ((hx-fabs(cn.x-eta.x))/hx)*((hy-fabs(cn.y-eta.y))/hy);
    cns.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    cn.w = ((hx-fabs(cn.x-eta.x))/hx)*((hy-fabs(cn.y-eta.y))/hy);
    cns.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    cn.w = ((hx-fabs(cn.x-eta.x))/hx)*((hy-fabs(cn.y-eta.y))/hy);
    cns.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    cn.w = ((hx-fabs(cn.x-eta.x))/hx)*((hy-fabs(cn.y-eta.y))/hy);
    cns.push_back(cn);
}

void IProblem2Backward2D::extendControlPoint2(const SpaceNodePDE &eta, std::vector<ControlNode> &cns, unsigned int i) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);

    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double hx = xd.step();
    double hy = yd.step();

    unsigned int rx = (unsigned int)(floor(eta.x*Nx));
    unsigned int ry = (unsigned int)(floor(eta.y*Ny));

    double hx3 = hx*hx*hx;
    double hx32 = (1.0/(2.0*hx3));
    double hx36 = (1.0/(6.0*hx3));

    double hy3 = hy*hy*hy;
    double hy32 = (1.0/(2.0*hy3));
    double hy36 = (1.0/(6.0*hy3));

    ControlNode cn;
    double dx = 0.0;
    double dy = 0.0;

    cn.i = rx-1; cn.x = cn.i*hx; cn.j = ry-1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    cns.push_back(cn);

    cn.i = rx-1; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    cns.push_back(cn);

    cn.i = rx-1; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    cns.push_back(cn);

    cn.i = rx-1; cn.x = cn.i*hx; cn.j = ry+2; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    cns.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry-1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    cns.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    cns.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    cns.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+2; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    cns.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry-1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    cns.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    cns.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    cns.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+2; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    cns.push_back(cn);

    cn.i = rx+2; cn.x = cn.i*hx; cn.j = ry-1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    cns.push_back(cn);

    cn.i = rx+2; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    cns.push_back(cn);

    cn.i = rx+2; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    cns.push_back(cn);

    cn.i = rx+2; cn.x = cn.i*hx; cn.j = ry+2; cn.y = cn.j*hy; cn.eta = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    cns.push_back(cn);
}

void IProblem2Backward2D::extendControlPoint3(const SpaceNodePDE &eta, std::vector<ControlNode> &ons, unsigned int i) const
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

    double factor = (1.0/(2.0*M_PI*sigmaX*sigmaY));

    unsigned int k=3;
    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        for (unsigned int m=ry-k; m<=ry+k; m++)
        {
            ControlNode cn;
            cn.i = n; cn.x = n*hx; cn.j = m; cn.y = m*hy; cn.eta = eta; cn.id = i;
            cn.w = factor*exp(-0.5*(((cn.x-eta.x)*(cn.x-eta.x))/(sigmaX*sigmaX)+((cn.y-eta.y)*(cn.y-eta.y))/(sigmaY*sigmaY)));
            ons.push_back(cn);
        }
    }
}

void IProblem2Backward2D::extendObservationDeltaPoint0(const SpaceNodePDE &xi, std::vector<ObservationDeltaNode> &ops, unsigned int j) const
{
#ifdef APPROX_BD_1
    extendObservationDeltaPoint1(xi, ops, j);
#endif
#ifdef APPROX_BD_2
    extendObservationDeltaPoint2(xi, ops, j);
#endif
#ifdef APPROX_BD_3
    extendObservationDeltaPoint3(xi, ops, j);
#endif
}

void IProblem2Backward2D::extendObservationDeltaPoint1(const SpaceNodePDE &xi, std::vector<ObservationDeltaNode> &ops, unsigned int j) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);

    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double hx = xd.step();
    double hy = yd.step();

    unsigned int rx = (unsigned int)(floor(xi.x*Nx));
    unsigned int ry = (unsigned int)(floor(xi.y*Ny));

    double factor = 1.0/(hx*hy);
    ObservationDeltaNode on;

    on.i = rx+0; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.xi = xi; on.id = j;
    on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy) * factor;
    ops.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.xi = xi; on.id = j;
    on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy) * factor;
    ops.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.xi = xi; on.id = j;
    on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy) * factor;
    ops.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.xi = xi; on.id = j;
    on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy) * factor;
    ops.push_back(on);
}

void IProblem2Backward2D::extendObservationDeltaPoint2(const SpaceNodePDE &xi, std::vector<ObservationDeltaNode> &ops, unsigned int j) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);

    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double hx = xd.step();
    double hy = yd.step();

    unsigned int rx = (unsigned int)(floor(xi.x*Nx));
    unsigned int ry = (unsigned int)(floor(xi.y*Ny));

    double hx3 = hx*hx*hx;
    double hx32 = (1.0/(2.0*hx3));
    double hx36 = (1.0/(6.0*hx3));

    double hy3 = hy*hy*hy;
    double hy32 = (1.0/(2.0*hy3));
    double hy36 = (1.0/(6.0*hy3));

    double factor = 1.0/(hx*hy);
    ObservationDeltaNode on;
    double dx = 0.0;
    double dy = 0.0;

    on.i = rx-1; on.x = on.i*hx; on.j = ry-1; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    ops.push_back(on);

    on.i = rx-1; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    ops.push_back(on);

    on.i = rx-1; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    ops.push_back(on);

    on.i = rx-1; on.x = on.i*hx; on.j = ry+2; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    ops.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry-1; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    ops.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    ops.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    ops.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry+2; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    ops.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry-1; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    ops.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    ops.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    ops.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+2; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    ops.push_back(on);

    on.i = rx+2; on.x = on.i*hx; on.j = ry-1; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    ops.push_back(on);

    on.i = rx+2; on.j = ry+0; on.x = on.i*hx; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    ops.push_back(on);

    on.i = rx+2; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    ops.push_back(on);

    on.i = rx+2; on.x = on.i*hx; on.j = ry+2; on.y = on.j*hy; on.xi = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    ops.push_back(on);
}

void IProblem2Backward2D::extendObservationDeltaPoint3(const SpaceNodePDE &xi, std::vector<ObservationDeltaNode> &ops, unsigned int j) const
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

    double factor = (1.0/(2.0*M_PI*sigmaX*sigmaY));

    unsigned int k=3;
    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        for (unsigned int m=ry-k; m<=ry+k; m++)
        {
            ObservationDeltaNode on;
            on.i = n; on.x = n*hx; on.j = m; on.y = m*hy; on.xi = xi; on.id = j;
            on.w = factor*exp(-0.5*(((on.x-xi.x)*(on.x-xi.x))/(sigmaX*sigmaX)+((on.y-xi.y)*(on.y-xi.y))/(sigmaY*sigmaY)));
            ops.push_back(on);
        }
    }
}
