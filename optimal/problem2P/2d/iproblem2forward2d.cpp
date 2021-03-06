#include "iproblem2forward2d.h"

IProblem2Forward2D::~IProblem2Forward2D()
{}

void IProblem2Forward2D::calculateMVD1(DoubleMatrix &u, vector<ExtendedSpaceNode2D> &info, bool use) const
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

    //unsigned int Lc = mParameter.Lc;
    unsigned int Lo = mParameter.Lo;

    if (use == true)
    {
        info.resize(Lo);
        for (unsigned int j=0; j<Lo; j++)
        {
            ExtendedSpaceNode2D &e = info[j];
            e.setSpaceNode(mParameter.xi[j]);
            e.id = j;
            e.extendWeights(xd, yd);
            e.extendLayers(L+1);
        }
    }

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix uh(M+1, N+1);

    //--------------------------------------------------------------------------------------------//

    std::vector<ObservationNode> observeNodes;
    for (unsigned int j=0; j<mParameter.Lo; j++) extendObservationPoint0(mParameter.xi[j], observeNodes, j);

#ifdef USE_F_VARIANT_2
    std::vector<ControlDeltaNode> cndeltaNodes;
    for (unsigned int i=0; i<mParameter.Lc; i++) extendContrlDeltaPoint0(mParameter.eta[i], cndeltaNodes, i);
#endif

    unsigned int *v1y = new unsigned int[M+1]; for (unsigned int m=0; m<=M; m++) v1y[m] = 0;
    unsigned int *v1x = new unsigned int[N+1]; for (unsigned int n=0; n<=N; n++) v1x[n] = 0;

    SpaceNodePDE sn;

#ifdef USE_F_VARIANT_1
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            for (unsigned int i=0; i<Lc; i++)
            {
                double _delta = delta(sn, mParameter.eta[i], i);

                if (checkDelta(_delta))
                {
                    if ( v1y[m] == 0 ) v1y[m] = 1;
                    if ( v1x[n] == 0 ) v1x[n] = 1;
                }
            }
        }
    }
#endif

#ifdef USE_F_VARIANT_2
    for (unsigned int i=0; i<cndeltaNodes.size(); i++)
    {
        const ControlDeltaNode &cdn = cndeltaNodes.at(i);

        v1x[cdn.i] = 1;
        v1y[cdn.j] = 1;
    }
#endif

    //------------------------------------- initial conditions -------------------------------------//
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[m][n] = initial(sn);
        }
    }

    if (use == true)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            ExtendedSpaceNode2D &pi = info[j];
            for (unsigned int r=0; r<4; r++)
            {
                for (unsigned int c=0; c<4; c++)
                {
                    pi.wi[r][c].u[0] = u[pi.wi[r][c].j][pi.wi[r][c].i];
                }
            }
        }
    }

    layerInfo(u, 0);
    //IPrinter::printMatrix(u[0]);
    //IPrinter::printSeperatorLine();
    //------------------------------------- initial conditions -------------------------------------//

    double a2_ht__hx2 = ((a*a*ht)/(hx*hx));
    double a2_ht__hy2 = ((a*a*ht)/(hy*hy));
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
    for (unsigned int l=1; l<=L; l++)
    {
        //puts("IProblem2Forward2D::calculateMVD.y-->...");
        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        {
            tn.i = l;
            tn.t = l*ht - 0.5*ht;

            {
                for (unsigned int n=0; n<=N; n++)
                {
                    if (v1x[n] == 0)
                    {
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

                            d2[offset+m] = 2.0*u[m][n] + lambda0*theta*ht + ht*f(sn, tn);

                            if (n==0)       d2[offset+m] += a2_ht__hx2*(u.at(m,0)   - 2.0*u.at(m,1)   + u.at(m,2));
                            if (n>0 && n<N) d2[offset+m] += a2_ht__hx2*(u.at(m,n-1) - 2.0*u.at(m,n)   + u.at(m,n+1));
                            if (n==N)       d2[offset+m] += a2_ht__hx2*(u.at(m,N-2) - 2.0*u.at(m,N-1) + u.at(m,N));

                            if (m == 0)
                            {
                                a2[offset+0] = 0.0;
                                b2[offset+0] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                                c2[offset+0] = -2.0*a2_ht__hy2;
                                d2[offset+0] += 2.0*a2_lambda_ht__hy*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                                d2[offset+0] += ((2.0*a*a*ht)/hy)*g3(sn, tn);
#endif
                            }
                            else if (m == M)
                            {
                                a2[offset+M] = -2.0*a2_ht__hy2;
                                b2[offset+M] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                                c2[offset+M] = 0.0;
                                d2[offset+M] += 2.0*a2_lambda_ht__hy*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                                d2[offset+M] += ((2.0*a*a*ht)/(hy))*g4(sn, tn);
#endif
                            }
                            else
                            {
                                a2[offset+m] = -a2_ht__hy2;
                                b2[offset+m] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht;
                                c2[offset+m] = -a2_ht__hy2;
                            }

                            //------------------------------------- Adding delta part -------------------------------------//
#ifdef USE_F_VARIANT_1
                            for (unsigned int i=0; i<Lc; i++)
                            {
                                double _delta = delta(sn, mParameter.eta[i], i);
                                if (checkDelta(_delta))
                                {
                                    for (unsigned int s=0; s<observeNodes.size(); s++)
                                    {
                                        const ObservationNode &on = observeNodes[s];

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntXSize; cs++)
                                        {
                                            if (on.i == cntX[cs])
                                            {
                                                found = true;
                                                w2[offset+m][cs*(M+1)+on.j] += -ht * mParameter.k[i][on.id] * _delta * on.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+m] += ht * mParameter.k[i][on.id] * uh[on.j][on.i] * _delta * on.w;
                                        }
                                    }

                                    for (unsigned int j=0; j<Lo; j++)
                                    {
                                        d2[offset+m] -= ht * mParameter.k[i][j] * mParameter.z[i][j] *_delta;
                                    }
                                }
                            }
#endif
                            //------------------------------------- Adding delta part -------------------------------------//
#ifdef USE_F_VARIANT_2
                            for (unsigned int cni=0; cni<cndeltaNodes.size(); cni++)
                            {
                                const ControlDeltaNode &cdn = cndeltaNodes.at(cni);
                                if (cdn.i == sn.i && cdn.j == sn.j)
                                {
                                    for (unsigned int s=0; s<observeNodes.size(); s++)
                                    {
                                        const ObservationNode &on = observeNodes.at(s);

                                        //if (on.n == cdn.n)
                                        //{
                                        //    printf("8887 %d %d %d %d %d %d %d\n", on.n, cdn.n, on.j, cdn.i, offset+m, cdn.n*(M+1)+on.m, cntXSize*(M+1));
                                        //    w2[offset+m][cdn.n*(M+1)+on.m] += -ht * mParameter.k[cdn.i][on.j] * cdn.w * on.w;
                                        //}
                                        //else
                                        //{
                                        //    d2[offset+m] += ht * mParameter.k[cdn.i][on.j] * uh[on.m][on.n] * cdn.w * on.w;
                                        //}

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntXSize; cs++)
                                        {
                                            if (on.i == cntX[cs])
                                            {
                                                found = true;
                                                w2[offset+m][cs*(M+1)+on.j] += -ht * mParameter.k[cdn.id][on.id] * cdn.w * on.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+m] += ht * mParameter.k[cdn.id][on.id] * uh[on.j][on.i] * cdn.w * on.w;
                                        }
                                    }

                                    for (unsigned int j=0; j<Lo; j++)
                                    {
                                        d2[offset+m] -= ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
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
                            uh[m][n] = x2[offset+m];
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
        //puts("IProblem2Forward2D::calculateMVD.y-->.");
        //IPrinter::printMatrix(uh);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        //puts("IProblem2Forward2D::calculateMVD.x-->...");
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

                            d1X[n] = 2.0*uh[m][n] + lambda0*theta*ht;
#ifdef USE_ADDITIONAL_FUNCTIONS
                            d1X[n] += ht*f(sn, tn);
#endif

                            if (m==0)       d1X[n] += a2_ht__hy2*(uh[0][n]   - 2.0*uh[1][n]   + uh[2][n]);
                            if (m>0 && m<M) d1X[n] += a2_ht__hy2*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]);
                            if (m==M)       d1X[n] += a2_ht__hy2*(uh[M-2][n] - 2.0*uh[M-1][n] + uh[M][n]);

                            if (n == 0)
                            {
                                a1X[0] = 0.0;
                                b1X[0] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                                c1X[0] = -2.0*a2_ht__hx2;
                                d1X[0] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                                d1X[0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
#endif
                            }
                            else if (n == N)
                            {
                                a1X[N] = -2.0*a2_ht__hx2;
                                b1X[N] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                                c1X[N] = 0.0;
                                d1X[N] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                                d1X[N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
#endif
                            }
                            else
                            {
                                a1X[n] = -a2_ht__hx2;
                                b1X[n] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
                                c1X[n] = -a2_ht__hx2;
                            }
                        }
                        tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                        for (unsigned int n=0; n<=N; n++) u[m][n] = x1X[n];
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

                            d2[offset+n] = 2.0*uh[m][n] + lambda0*theta*ht + ht*f(sn, tn);

                            if (m==0)       d2[offset+n] += a2_ht__hy2*(uh[0][n]   - 2.0*uh[1][n]   + uh[2][n]);
                            if (m>0 && m<M) d2[offset+n] += a2_ht__hy2*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]);
                            if (m==M)       d2[offset+n] += a2_ht__hy2*(uh[M-2][n] - 2.0*uh[M-1][n] + uh[M][n]);

                            if (n == 0)
                            {
                                a2[offset+0] = 0.0;
                                b2[offset+0] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                                c2[offset+0] = -2.0*a2_ht__hx2;
                                d2[offset+0] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                                d2[offset+0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
#endif
                            }
                            else if (n == N)
                            {
                                a2[offset+N] = -2.0*a2_ht__hx2;
                                b2[offset+N] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                                c2[offset+N] = 0.0;
                                d2[offset+N] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                                d2[offset+N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
#endif
                            }
                            else
                            {
                                a2[offset+n] = -a2_ht__hx2;
                                b2[offset+n] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
                                c2[offset+n] = -a2_ht__hx2;
                            }

                            //------------------------------------- Adding delta part -------------------------------------//
#ifdef USE_F_VARIANT_1
                            for (unsigned int i=0; i<Lc; i++)
                            {
                                double _delta = delta(sn, mParameter.eta[i], i);
                                if (checkDelta(_delta))
                                {
                                    for (unsigned int s=0; s<observeNodes.size(); s++)
                                    {
                                        const ObservationNode &on = observeNodes[s];

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntYSize; cs++)
                                        {
                                            if (on.j == cntY[cs])
                                            {
                                                found = true;
                                                w2[offset+n][cs*(N+1)+on.i] += -ht * mParameter.k[i][on.id] * _delta * on.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+n] += ht * mParameter.k[i][on.id] * u[on.j][on.i] * _delta * on.w;
                                        }
                                    }

                                    for (unsigned int j=0; j<Lo; j++)
                                    {
                                        d2[offset+n] -= ht * mParameter.k[i][j] * mParameter.z[i][j] *_delta;
                                    }
                                }
                            }
#endif
                            //------------------------------------- Adding delta part -------------------------------------//
#ifdef USE_F_VARIANT_2
                            for (unsigned int cni=0; cni<cndeltaNodes.size(); cni++)
                            {
                                const ControlDeltaNode &cdn = cndeltaNodes.at(cni);
                                if (cdn.i == sn.i && cdn.j == sn.j)
                                {
                                    //printf("%d %d\n", cdn.i, cdn.j);
                                    for (unsigned int s=0; s<observeNodes.size(); s++)
                                    {
                                        const ObservationNode &on = observeNodes[s];
                                        //printf("%d %d\n", on.i, on.j);

                                        //if (on.m == cdn.m)
                                        //{
                                        //    w2[offset+n][cdn.m*(N+1)+on.n] += -ht * mParameter.k[cdn.i][on.j] * cdn.w * on.w;
                                        //}
                                        //else
                                        //{
                                        //    d2[offset+n] += ht * mParameter.k[cdn.i][on.j] * u[on.m][on.n] * cdn.w * on.w;
                                        //}

                                        bool found = false;
                                        for (unsigned int cs=0; cs<cntYSize; cs++)
                                        {
                                            if (on.j == cntY[cs])
                                            {
                                                found = true;
                                                w2[offset+n][cs*(N+1)+on.i] += -ht * mParameter.k[cdn.id][on.id] * cdn.w * on.w;
                                            }
                                        }

                                        if (!found)
                                        {
                                            d2[offset+n] += ht * mParameter.k[cdn.id][on.id] * u[on.j][on.i] * cdn.w * on.w;
                                        }
                                    }

                                    for (unsigned int j=0; j<Lo; j++)
                                    {
                                        d2[offset+n] -= ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
                                    }

                                }
                            }
#endif
                        }
                        offset += N+1;
                    }
                }

                LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cntYSize*(N+1));

                if (l==L)
                {
                    FILE *file = fopen("data.txt", "w");
                    IPrinter::print(w2,w2.rows(),w2.cols(),14,10,file);
                    fclose(file);
                }

                offset = 0;
                for (unsigned int m=0; m<=M; m++)
                {
                    if (v1y[m] == 1)
                    {
                        for (unsigned int n=0; n<=N; n++)
                        {
                            u[m][n] = x2[offset+n];
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
        //puts("IProblem2Forward2D::calculateMVD.x-->.");

        //IPrinter::printMatrix(u[l]);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        //puts("IProblem2Forward2D::calculateMVD.use...");
        if (use == true)
        {
            for (unsigned int j=0; j<Lo; j++)
            {
                ExtendedSpaceNode2D &pi = info[j];
                for (unsigned int r=0; r<4; r++)
                {
                    for (unsigned int c=0; c<4; c++)
                    {
                        unsigned int x = pi.wi[r][c].i;
                        unsigned int y = pi.wi[r][c].j;
                        pi.wi[r][c].u[l] = u[y][x];
                    }
                }
            }
        }
        //puts("IProblem2Forward2D::calculateMVD.use.");

        layerInfo(u, l);
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
#ifdef USE_F_VARIANT_2
    cndeltaNodes.clear();
#endif
    uh.clear();
}

inline double rand_error()
{
    return 1.0;
}

void IProblem2Forward2D::calculateMVD(DoubleMatrix &u, vector<ExtendedSpaceNode2D> &info, bool use) const
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

    //unsigned int Lc = mParameter.Lc;
    unsigned int Lo = mParameter.Lo;

    if (use == true)
    {
        info.resize(Lo);
        for (unsigned int j=0; j<Lo; j++)
        {
            ExtendedSpaceNode2D &e = info[j];
            e.setSpaceNode(mParameter.xi[j]);
            e.id = j;
            e.extendWeights(xd, yd);
            e.extendLayers(L+1);
        }
    }

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix uh(M+1, N+1);

    //--------------------------------------------------------------------------------------------//

    std::vector<ObservationNode> observeNodes;
    for (unsigned int j=0; j<mParameter.Lo; j++) extendObservationPoint0(mParameter.xi[j], observeNodes, j);

    std::vector<ControlDeltaNode> cndeltaNodes;
    for (unsigned int i=0; i<mParameter.Lc; i++) extendContrlDeltaPoint0(mParameter.eta[i], cndeltaNodes, i);

    SpaceNodePDE sn;

    vector<unsigned int> rows0;
    vector<unsigned int> rows1;
    vector<unsigned int> rows2;
    for (unsigned int m=0; m<=M; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cndeltaNodes.size(); i++)
        {
            const ControlDeltaNode &cdn = cndeltaNodes.at(i);
            if (cdn.j == m)
            {
                for (unsigned int j=0; j<observeNodes.size(); j++)
                {
                    const ObservationNode &on = observeNodes.at(j);
                    if (on.j == m)
                    {
                        found2 = true;
                        break;
                    }
                }
                found1 = true;
                break;
            }
        }
        if (found1 == true  && found2 == true)  rows2.push_back(m);
        if (found1 == true)                     rows1.push_back(m);
        if (found1 == false && found2 == false) rows0.push_back(m);
    }

    vector<unsigned int> cols0;
    vector<unsigned int> cols1;
    vector<unsigned int> cols2;
    for (unsigned int n=0; n<=N; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cndeltaNodes.size(); i++)
        {
            const ControlDeltaNode &cdn = cndeltaNodes.at(i);
            if (cdn.i == n)
            {
                for (unsigned int j=0; j<observeNodes.size(); j++)
                {
                    const ObservationNode &on = observeNodes.at(j);
                    if (on.i == n)
                    {
                        found2 = true;
                        break;
                    }
                }
                found1 = true;
                break;
            }
        }
        if (found1 == true  && found2 == true)  cols2.push_back(n);
        if (found1 == true)                     cols1.push_back(n);
        if (found1 == false && found2 == false) cols0.push_back(n);
    }

    //------------------------------------- initial conditions -------------------------------------//
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[m][n] = initial(sn);
        }
    }

    if (use == true)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            ExtendedSpaceNode2D &pi = info[j];
            for (unsigned int r=0; r<4; r++)
            {
                for (unsigned int c=0; c<4; c++)
                {
                    pi.wi[r][c].u[0] = u[pi.wi[r][c].j][pi.wi[r][c].i];
                }
            }
        }
    }

    layerInfo(u, 0);
    //IPrinter::printMatrix(u[0]);
    //IPrinter::printSeperatorLine();
    //------------------------------------- initial conditions -------------------------------------//

    double a2_ht__hx2 = ((a*a*ht)/(hx*hx));
    double a2_ht__hy2 = ((a*a*ht)/(hy*hy));
    double lambda0_ht = lambda0*ht;
    double a2_lambda_ht__hy = (a*a*lambda*ht)/(hy);
    double a2_lambda_ht__hx = (a*a*lambda*ht)/(hx);

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
        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        {
            tn.i = l;
            tn.t = l*ht;
            {
                for (unsigned int row=0; row<rows0.size(); row++)
                {
                    unsigned int m = rows0.at(row);
                    sn.j = m; sn.y = m*hy;
                    for (unsigned int n=0; n<=N; n++)
                    {
                        sn.i = n; sn.x = n*hx;

                        d1X[n] = 2.0*u[m][n] + lambda0*theta*ht;
#ifdef USE_ADDITIONAL_FUNCTIONS
                        d1X[n] += ht*f(sn, tn);
#endif

                        if (m==0)       d1X[n] += a2_ht__hy2*(u[0][n]   - 2.0*u[1][n]   + u[2][n]);
                        if (m>0 && m<M) d1X[n] += a2_ht__hy2*(u[m-1][n] - 2.0*u[m][n]   + u[m+1][n]);
                        if (m==M)       d1X[n] += a2_ht__hy2*(u[M-2][n] - 2.0*u[M-1][n] + u[M][n]);

                        if (n == 0)
                        {
                            a1X[0] = 0.0;
                            b1X[0] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                            c1X[0] = -2.0*a2_ht__hx2;
                            d1X[0] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                            d1X[0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
#endif
                        }
                        else if (n == N)
                        {
                            a1X[N] = -2.0*a2_ht__hx2;
                            b1X[N] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                            c1X[N] = 0.0;
                            d1X[N] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                            d1X[N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
#endif
                        }
                        else
                        {
                            a1X[n] = -a2_ht__hx2;
                            b1X[n] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
                            c1X[n] = -a2_ht__hx2;
                        }
                    }

                    tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                    for (unsigned int n=0; n<=N; n++) uh[m][n] = x1X[n];
                }
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

                        d1X[n] = 2.0*u[m][n] + lambda0*theta*ht;
#ifdef USE_ADDITIONAL_FUNCTIONS
                        d1X[n] += ht*f(sn, tn);
#endif
                        if (m==0)       d1X[n] += a2_ht__hy2*(u[0][n]   - 2.0*u[1][n]   + u[2][n]);
                        if (m>0 && m<M) d1X[n] += a2_ht__hy2*(u[m-1][n] - 2.0*u[m][n]   + u[m+1][n]);
                        if (m==M)       d1X[n] += a2_ht__hy2*(u[M-2][n] - 2.0*u[M-1][n] + u[M][n]);

                        if (n == 0)
                        {
                            a1X[0] = 0.0;
                            b1X[0] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                            c1X[0] = -2.0*a2_ht__hx2;
                            d1X[0] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                            d1X[0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
#endif
                        }
                        else if (n == N)
                        {
                            a1X[N] = -2.0*a2_ht__hx2;
                            b1X[N] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                            c1X[N] = 0.0;
                            d1X[N] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                            d1X[N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
#endif
                        }
                        else
                        {
                            a1X[n] = -a2_ht__hx2;
                            b1X[n] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
                            c1X[n] = -a2_ht__hx2;
                        }

                        for (unsigned int cni=0; cni<cndeltaNodes.size(); cni++)
                        {
                            const ControlDeltaNode &cdn = cndeltaNodes.at(cni);
                            if (cdn.i == sn.i && cdn.j == sn.j)
                            {
                                for (unsigned int s=0; s<observeNodes.size(); s++)
                                {
                                    const ObservationNode &on = observeNodes[s];
                                    d1X[n] += ht * mParameter.k[cdn.id][on.id] * uh[on.j][on.i] * cdn.w * on.w*rand_error();
                                }

                                for (unsigned int j=0; j<Lo; j++)
                                {
                                    d1X[n] -= ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
                                }
                            }
                        }
                    }
                    tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
                    for (unsigned int n=0; n<=N; n++) uh[m][n] = x1X[n];
                }
            }
            else
            {
                double* a2 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
                double* b2 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
                double* c2 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
                double* d2 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
                double* x2 = (double*) malloc(sizeof(double)*rows1.size()*(N+1));
                DoubleMatrix w2(rows1.size()*(N+1), rows1.size()*(N+1), 0.0);

                unsigned int offset = 0;
                for (unsigned int row=0; row<rows1.size(); row++)
                {
                    unsigned int m = rows1.at(row);
                    sn.j = m; sn.y = m*hy;

                    for (unsigned int n=0; n<=N; n++)
                    {
                        sn.i = n; sn.x = n*hx;

                        d2[offset+n] = 2.0*u[m][n] + lambda0*theta*ht;
#ifdef USE_ADDITIONAL_FUNCTIONS
                        d2[offset+n] += ht*f(sn, tn);
#endif
                        if (m==0)       d2[offset+n] += a2_ht__hy2*(u[0][n]   - 2.0*u[1][n]   + u[2][n]);
                        if (m>0 && m<M) d2[offset+n] += a2_ht__hy2*(u[m-1][n] - 2.0*u[m][n]   + u[m+1][n]);
                        if (m==M)       d2[offset+n] += a2_ht__hy2*(u[M-2][n] - 2.0*u[M-1][n] + u[M][n]);

                        if (n == 0)
                        {
                            a2[offset+0] = 0.0;
                            b2[offset+0] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                            c2[offset+0] = -2.0*a2_ht__hx2;
                            d2[offset+0] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                            d2[offset+0] += ((2.0*a*a*ht)/hx)*g1(sn, tn);
#endif
                        }
                        else if (n == N)
                        {
                            a2[offset+N] = -2.0*a2_ht__hx2;
                            b2[offset+N] = +2.0 + 2.0*a2_ht__hx2 + lambda0_ht + 2.0*a2_lambda_ht__hx;
                            c2[offset+N] = 0.0;
                            d2[offset+N] += 2.0*a2_lambda_ht__hx*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                            d2[offset+N] += ((2.0*a*a*ht)/(hx))*g2(sn, tn);
#endif
                        }
                        else
                        {
                            a2[offset+n] = -a2_ht__hx2;
                            b2[offset+n] = 2.0 + 2.0*a2_ht__hx2 + lambda0_ht;
                            c2[offset+n] = -a2_ht__hx2;
                        }

                        //------------------------------------- Adding delta part -------------------------------------//
                        for (unsigned int cni=0; cni<cndeltaNodes.size(); cni++)
                        {
                            const ControlDeltaNode &cdn = cndeltaNodes.at(cni);
                            if (cdn.i == sn.i && cdn.j == sn.j)
                            {
                                //printf("%d %d\n", cdn.i, cdn.j);
                                for (unsigned int s=0; s<observeNodes.size(); s++)
                                {
                                    const ObservationNode &on = observeNodes[s];
                                    //printf("%d %d\n", on.i, on.j);

                                    //if (on.m == cdn.m)
                                    //{
                                    //    w2[offset+n][cdn.m*(N+1)+on.n] += -ht * mParameter.k[cdn.i][on.j] * cdn.w * on.w;
                                    //}
                                    //else
                                    //{
                                    //    d2[offset+n] += ht * mParameter.k[cdn.i][on.j] * u[on.m][on.n] * cdn.w * on.w;
                                    //}

                                    bool found = false;
                                    for (unsigned int cs=0; cs<rows1.size(); cs++)
                                    {
                                        if (on.j == rows1[cs])
                                        {
                                            found = true;
                                            w2[offset+n][cs*(N+1)+on.i] += -ht * mParameter.k[cdn.id][on.id] * cdn.w * on.w*rand_error();
                                        }
                                    }

                                    if (!found)
                                    {
                                        d2[offset+n] += ht * mParameter.k[cdn.id][on.id] * uh[on.j][on.i] * cdn.w * on.w*rand_error();
                                    }
                                }

                                for (unsigned int j=0; j<Lo; j++)
                                {
                                    d2[offset+n] -= ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
                                }

                            }
                        }
                    }
                    offset += N+1;
                }

                LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, rows1.size()*(N+1));

                offset = 0;
                for (unsigned int row=0; row<rows1.size(); row++)
                {
                    unsigned int m=rows1.at(row);
                    for (unsigned int n=0; n<=N; n++)
                    {
                        uh[m][n] = x2[offset+n];
                    }
                    offset += N+1;
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
        //IPrinter::printMatrix(u[l]);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        {
            tn.i = l;
            tn.t = l*ht - 0.5*ht;
            //--------------------------------------------------------------------------//
            {
                for (unsigned int col=0; col<cols0.size(); col++)
                {
                    unsigned int n = cols0[col];
                    sn.i = n; sn.x = n*hx;
                    for (unsigned int m=0; m<=M; m++)
                    {
                        sn.j = m; sn.y = m*hy;

                        d1Y[m] = 2.0*uh[m][n] + lambda0*theta*ht;
#ifdef USE_ADDITIONAL_FUNCTIONS
                        d1Y[m] += ht*f(sn, tn);
#endif
                        if (n==0)       d1Y[m] += a2_ht__hx2*(uh.at(m,0)   - 2.0*uh.at(m,1)   + uh.at(m,2));
                        if (n>0 && n<N) d1Y[m] += a2_ht__hx2*(uh.at(m,n-1) - 2.0*uh.at(m,n)   + uh.at(m,n+1));
                        if (n==N)       d1Y[m] += a2_ht__hx2*(uh.at(m,N-2) - 2.0*uh.at(m,N-1) + uh.at(m,N));

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
                    for (unsigned int m=0; m<=M; m++) u[m][n] = x1Y[m];
                }
            }
            //--------------------------------------------------------------------------//
            if (cols2.size() == 0)
            {
                for (unsigned int col=0; col<cols1.size(); col++)
                {
                    unsigned int n=cols1[col];
                    sn.i = n; sn.x = n*hx;
                    for (unsigned int m=0; m<=M; m++)
                    {
                        sn.j = m; sn.y = m*hy;

                        d1Y[m] = 2.0*uh[m][n] + lambda0*theta*ht;
#ifdef USE_ADDITIONAL_FUNCTIONS
                        d1Y[m] += ht*f(sn, tn);
#endif
                        if (n==0)       d1Y[m] += a2_ht__hx2*(uh.at(m,0)   - 2.0*uh.at(m,1)   + uh.at(m,2));
                        if (n>0 && n<N) d1Y[m] += a2_ht__hx2*(uh.at(m,n-1) - 2.0*uh.at(m,n)   + uh.at(m,n+1));
                        if (n==N)       d1Y[m] += a2_ht__hx2*(uh.at(m,N-2) - 2.0*uh.at(m,N-1) + uh.at(m,N));

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

                        for (unsigned int cni=0; cni<cndeltaNodes.size(); cni++)
                        {
                            const ControlDeltaNode &cdn = cndeltaNodes.at(cni);
                            if (cdn.i == sn.i && cdn.j == sn.j)
                            {
                                for (unsigned int s=0; s<observeNodes.size(); s++)
                                {
                                    const ObservationNode &on = observeNodes.at(s);
                                    d1Y[m] += ht * mParameter.k[cdn.id][on.id] * u[on.j][on.i] * cdn.w * on.w*rand_error();
                                }

                                for (unsigned int j=0; j<Lo; j++)
                                {
                                    d1Y[m] -= ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
                                }
                            }
                        }
                    }
                    tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
                    for (unsigned int m=0; m<=M; m++) u[m][n] = x1Y[m];
                }
            }
            else
            {
                double* a2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
                double* b2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
                double* c2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
                double* d2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
                double* x2 = (double*) malloc(sizeof(double)*cols1.size()*(M+1));
                DoubleMatrix w2(cols1.size()*(M+1), cols1.size()*(M+1), 0.0);

                unsigned int offset = 0;
                for (unsigned int col=0; col<cols1.size(); col++)
                {
                    unsigned int n = cols1.at(col);
                    sn.i = n; sn.x = n*hx;

                    for (unsigned int m=0; m<=M; m++)
                    {
                        sn.j = m; sn.y = m*hy;

                        d2[offset+m] = 2.0*uh[m][n] + lambda0*theta*ht;
#ifdef USE_ADDITIONAL_FUNCTIONS
                        d2[offset+m] += ht*f(sn, tn);
#endif
                        if (n==0)       d2[offset+m] += a2_ht__hx2*(uh.at(m,0)   - 2.0*uh.at(m,1)   + uh.at(m,2));
                        if (n>0 && n<N) d2[offset+m] += a2_ht__hx2*(uh.at(m,n-1) - 2.0*uh.at(m,n)   + uh.at(m,n+1));
                        if (n==N)       d2[offset+m] += a2_ht__hx2*(uh.at(m,N-2) - 2.0*uh.at(m,N-1) + uh.at(m,N));

                        if (m == 0)
                        {
                            a2[offset+0] = 0.0;
                            b2[offset+0] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                            c2[offset+0] = -2.0*a2_ht__hy2;
                            d2[offset+0] += 2.0*a2_lambda_ht__hy*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                            d2[offset+0] += ((2.0*a*a*ht)/hy)*g3(sn, tn);
#endif
                        }
                        else if (m == M)
                        {
                            a2[offset+M] = -2.0*a2_ht__hy2;
                            b2[offset+M] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht + 2.0*a2_lambda_ht__hy;
                            c2[offset+M] = 0.0;
                            d2[offset+M] += 2.0*a2_lambda_ht__hy*theta;
#ifdef USE_ADDITIONAL_FUNCTIONS
                            d2[offset+M] += ((2.0*a*a*ht)/(hy))*g4(sn, tn);
#endif
                        }
                        else
                        {
                            a2[offset+m] = -a2_ht__hy2;
                            b2[offset+m] = +2.0 + 2.0*a2_ht__hy2 + lambda0_ht;
                            c2[offset+m] = -a2_ht__hy2;
                        }

                        //------------------------------------- Adding delta part -------------------------------------//
                        for (unsigned int cni=0; cni<cndeltaNodes.size(); cni++)
                        {
                            const ControlDeltaNode &cdn = cndeltaNodes.at(cni);
                            if (cdn.i == sn.i && cdn.j == sn.j)
                            {
                                for (unsigned int s=0; s<observeNodes.size(); s++)
                                {
                                    const ObservationNode &on = observeNodes.at(s);

                                    //if (on.n == cdn.n)
                                    //{
                                    //    w2[offset+m][cdn.n*(M+1)+on.m] += -ht * mParameter.k[cdn.i][on.j] * cdn.w * on.w;
                                    //}
                                    //else
                                    //{
                                    //    d2[offset+m] += ht * mParameter.k[cdn.i][on.j] * uh[on.m][on.n] * cdn.w * on.w;
                                    //}

                                    bool found = false;
                                    for (unsigned int cs=0; cs<cols1.size(); cs++)
                                    {
                                        if (on.i == cols1[cs])
                                        {
                                            found = true;
                                            w2[offset+m][cs*(M+1)+on.j] += -ht * mParameter.k[cdn.id][on.id] * cdn.w * on.w*rand_error();
                                        }
                                    }

                                    if (!found)
                                    {
                                        d2[offset+m] += ht * mParameter.k[cdn.id][on.id] * u[on.j][on.i] * cdn.w * on.w*rand_error();
                                    }
                                }

                                for (unsigned int j=0; j<Lo; j++)
                                {
                                    d2[offset+m] -= ht * mParameter.k[cdn.id][j] * mParameter.z[cdn.id][j] * cdn.w;
                                }
                            }
                        }
                        //------------------------------------- Adding delta part -------------------------------------//
                    }
                    offset += M+1;
                }

                LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M+1));

                offset = 0;
                for (unsigned int col=0; col<cols1.size(); col++)
                {
                    unsigned int n=cols1.at(col);
                    for (unsigned int m=0; m<=M; m++)
                    {
                        u[m][n] = x2[offset+m];
                    }
                    offset += M+1;
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
        //IPrinter::printMatrix(uh);
        //IPrinter::printSeperatorLine();
        //------------------------------------- approximatin to y direction conditions -------------------------------------//



        if (use == true)
        {
            for (unsigned int j=0; j<Lo; j++)
            {
                ExtendedSpaceNode2D &pi = info[j];
                for (unsigned int r=0; r<4; r++)
                {
                    for (unsigned int c=0; c<4; c++)
                    {
                        unsigned int x = pi.wi[r][c].i;
                        unsigned int y = pi.wi[r][c].j;
                        pi.wi[r][c].u[l] = u[y][x];
                    }
                }
            }
        }

        layerInfo(u, l);

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
#ifdef USE_F_VARIANT_2
    cndeltaNodes.clear();
#endif
    uh.clear();
}

bool IProblem2Forward2D::checkDelta(double _delta) const
{
    return _delta != 0.0;
}

double IProblem2Forward2D::delta(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i UNUSED_PARAM, unsigned int source UNUSED_PARAM) const
{
    return delta3(sn, eta, i);
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
    double dhx = fabs(sn.x-eta.x);
    if ( dhx <= hx )
    {
        resX = (hx-dhx) / hx;
    }

    double resY = 0.0;
    double dhy = fabs(sn.y-eta.y);
    if ( dhy <= hy )
    {
        resY = (hy-dhy) / hy;
    }

    double res = (resX*resY)/(hx*hy);

    return res;
}

void IProblem2Forward2D::extendObservationPoint0(const SpaceNodePDE &xi, std::vector<ObservationNode> &ons, unsigned int j) const
{
#ifdef APPROX_F1_1
    extendObservationPoint1(xi, ons, j);
#endif
#ifdef APPROX_F1_2
    extendObservationPoint2(xi, ons, j);
#endif
#ifdef APPROX_F1_3
    extendObservationPoint3(xi, ons, j);
#endif
#ifdef APPROX_F1_4
    extendObservationPoint4(xi, ons, j);
#endif
}

void IProblem2Forward2D::extendObservationPoint1(const SpaceNodePDE &xi, std::vector<ObservationNode> &ons, unsigned int j) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);

    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double hx = xd.step();
    double hy = yd.step();

    unsigned int rx = (unsigned int)(floor(xi.x*Nx));
    unsigned int ry = (unsigned int)(floor(xi.y*Ny));

    ObservationNode on;

    on.i = rx+0; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.pt = xi; on.id = j;
    on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy);
    ons.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.pt = xi; on.id = j;
    on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy);
    ons.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.pt = xi; on.id = j;
    on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy);
    ons.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.pt = xi; on.id = j;
    on.w = ((hx-fabs(on.x-xi.x))/hx)*((hy-fabs(on.y-xi.y))/hy);
    ons.push_back(on);
}

void IProblem2Forward2D::extendObservationPoint2(const SpaceNodePDE &xi, std::vector<ObservationNode> &ons, unsigned int j) const
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

    ObservationNode on;
    double dx = 0.0;
    double dy = 0.0;

    on.i = rx-1; on.x = on.i*hx; on.j = ry-1; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    ons.push_back(on);

    on.i = rx-1; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    ons.push_back(on);

    on.i = rx-1; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    ons.push_back(on);

    on.i = rx-1; on.x = on.i*hx; on.j = ry+2; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    ons.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry-1; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    ons.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    ons.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    ons.push_back(on);

    on.i = rx+0; on.x = on.i*hx; on.j = ry+2; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    ons.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry-1; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    ons.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+0; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    ons.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    ons.push_back(on);

    on.i = rx+1; on.x = on.i*hx; on.j = ry+2; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    ons.push_back(on);

    on.i = rx+2; on.x = on.i*hx; on.j = ry-1; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    ons.push_back(on);

    on.i = rx+2; on.j = ry+0; on.x = on.i*hx; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    ons.push_back(on);

    on.i = rx+2; on.x = on.i*hx; on.j = ry+1; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
    ons.push_back(on);

    on.i = rx+2; on.x = on.i*hx; on.j = ry+2; on.y = on.j*hy; on.pt = xi; on.id = j;
    dx = fabs(on.x-xi.x); dy = fabs(on.y-xi.y); on.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    ons.push_back(on);
}

void IProblem2Forward2D::extendObservationPoint3(const SpaceNodePDE &xi, std::vector<ObservationNode> &ons, unsigned int j) const
{
    double hx = spaceDimension(Dimension::DimensionX).step();
    double hy = spaceDimension(Dimension::DimensionY).step();

    std::vector<ExtendedDeltaPoint> extpoint;
    distributeDelta(xi, extpoint, j);

    for (unsigned int i=0; i<extpoint.size(); i++)
    {
        ExtendedDeltaPoint ep = extpoint.at(i);
        ObservationNode node;
        node.id = ep.id; node.w = ep.w * (hx*hy); node.pt = ep.pt;
        node.i = ep.i; node.x = ep.x;
        node.j = ep.j; node.y = ep.y;
        ons.push_back(node);
    }
    extpoint.clear();

    //    Dimension xd = spaceDimension(Dimension::DimensionX);
    //    Dimension yd = spaceDimension(Dimension::DimensionY);
    //    double hx = xd.step();
    //    double hy = yd.step();
    //    unsigned int Nx = xd.sizeN();
    //    unsigned int Ny = yd.sizeN();

    //    double sigmaX = hx;
    //    double sigmaY = hy;

    //    unsigned int rx = (unsigned int)(round(xi.x*Nx));
    //    unsigned int ry = (unsigned int)(round(xi.y*Ny));

    //    double factor = (1.0/(2.0*M_PI*sigmaX*sigmaY));

    //    unsigned int k=5;
    //    for (unsigned int n=rx-k; n<=rx+k; n++)
    //    {
    //        for (unsigned int m=ry-k; m<=ry+k; m++)
    //        {
    //            ObservationNode on;
    //            on.i = n; on.x = n*hx; on.j = m; on.y = m*hy; on.pt = xi; on.id = j;
    //            on.w = factor*exp(-0.5*(((on.x-xi.x)*(on.x-xi.x))/(sigmaX*sigmaX)+((on.y-xi.y)*(on.y-xi.y))/(sigmaY*sigmaY))) * (hx*hy);
    //            ons.push_back(on);
    //        }
    //    }
}

void IProblem2Forward2D::extendContrlDeltaPoint0(const SpaceNodePDE &eta, std::vector<ControlDeltaNode> &cps, unsigned int i) const
{
#ifdef APPROX_FD_1
    extendContrlDeltaPoint1(eta, cps, i);
#endif
#ifdef APPROX_FD_2
    extendContrlDeltaPoint2(eta, cps, i);
#endif
#ifdef APPROX_FD_3
    extendContrlDeltaPoint3(eta, cps, i);
#endif
#ifdef APPROX_FD_4
    extendContrlDeltaPoint4(eta, cps, i);
#endif
}

void IProblem2Forward2D::extendContrlDeltaPoint1(const SpaceNodePDE &eta, std::vector<ControlDeltaNode> &cps, unsigned int i) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);

    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double hx = xd.step();
    double hy = yd.step();

    unsigned int rx = (unsigned int)(floor(eta.x*Nx));
    unsigned int ry = (unsigned int)(floor(eta.y*Ny));

    double factor = 1.0/(hx*hy);
    ControlDeltaNode cn;

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    cn.w = ((hx-fabs(cn.x-eta.x))/hx)*((hy-fabs(cn.y-eta.y))/hy) * factor;
    cps.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    cn.w = ((hx-fabs(cn.x-eta.x))/hx)*((hy-fabs(cn.y-eta.y))/hy) * factor;
    cps.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    cn.w = ((hx-fabs(cn.x-eta.x))/hx)*((hy-fabs(cn.y-eta.y))/hy) * factor;
    cps.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    cn.w = ((hx-fabs(cn.x-eta.x))/hx)*((hy-fabs(cn.y-eta.y))/hy) * factor;
    cps.push_back(cn);
}

void IProblem2Forward2D::extendContrlDeltaPoint2(const SpaceNodePDE &eta, std::vector<ControlDeltaNode> &cps, unsigned int i) const
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

    double factor = 1.0/(hx*hy);
    ControlDeltaNode cn;
    double dx = 0.0;
    double dy = 0.0;

    cn.i = rx-1; cn.x = cn.i*hx; cn.j = ry-1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    cps.push_back(cn);

    cn.i = rx-1; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    cps.push_back(cn);

    cn.i = rx-1; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    cps.push_back(cn);

    cn.i = rx-1; cn.x = cn.i*hx; cn.j = ry+2; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    cps.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry-1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    cps.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    cps.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    cps.push_back(cn);

    cn.i = rx+0; cn.x = cn.i*hx; cn.j = ry+2; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    cps.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry-1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    cps.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+0; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    cps.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    cps.push_back(cn);

    cn.i = rx+1; cn.x = cn.i*hx; cn.j = ry+2; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    cps.push_back(cn);

    cn.i = rx+2; cn.x = cn.i*hx; cn.j = ry-1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    cps.push_back(cn);

    cn.i = rx+2; cn.j = ry+0; cn.x = cn.i*hx; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    cps.push_back(cn);

    cn.i = rx+2; cn.x = cn.i*hx; cn.j = ry+1; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32)*factor;
    cps.push_back(cn);

    cn.i = rx+2; cn.x = cn.i*hx; cn.j = ry+2; cn.y = cn.j*hy; cn.pt = eta; cn.id = i;
    dx = fabs(cn.x-eta.x); dy = fabs(cn.y-eta.y); cn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36)*factor;
    cps.push_back(cn);
}

void IProblem2Forward2D::extendContrlDeltaPoint3(const SpaceNodePDE &eta, std::vector<ControlDeltaNode> &cps, unsigned int id) const
{
    std::vector<ExtendedDeltaPoint> extpoint;
    distributeDelta(eta, extpoint, id);

    for (unsigned int i=0; i<extpoint.size(); i++)
    {
        ExtendedDeltaPoint ep = extpoint.at(i);
        ControlDeltaNode node;
        node.id = ep.id; node.w = ep.w; node.pt = ep.pt;
        node.i = ep.i; node.x = ep.x;
        node.j = ep.j; node.y = ep.y;
        cps.push_back(node);
    }
    extpoint.clear();

    //    Dimension xd = spaceDimension(Dimension::DimensionX);
    //    Dimension yd = spaceDimension(Dimension::DimensionY);
    //    double hx = xd.step();
    //    double hy = yd.step();
    //    unsigned int Nx = xd.sizeN();
    //    unsigned int Ny = yd.sizeN();

    //    double sigmaX = hx;
    //    double sigmaY = hy;

    //    unsigned int rx = (unsigned int)(round(eta.x*Nx));
    //    unsigned int ry = (unsigned int)(round(eta.y*Ny));

    //    double factor = (1.0/(2.0*M_PI*sigmaX*sigmaY));

    //    unsigned int k=5;
    //    for (unsigned int n=rx-k; n<=rx+k; n++)
    //    {
    //        for (unsigned int m=ry-k; m<=ry+k; m++)
    //        {
    //            ControlDeltaNode cn;
    //            cn.i = n; cn.x = n*hx; cn.j = m; cn.y = m*hy; cn.pt = eta; cn.id = i;
    //            cn.w = factor*exp(-0.5*(((cn.x-eta.x)*(cn.x-eta.x))/(sigmaX*sigmaX)+((cn.y-eta.y)*(cn.y-eta.y))/(sigmaY*sigmaY)));
    //            cps.push_back(cn);
    //        }
    //    }
}

//--------------------------//

void IProblem2Forward2D::extendContrlDeltaPoint4(const SpaceNodePDE &eta, std::vector<ControlDeltaNode> &cps, unsigned int i) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();
    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double sigmaX = 3.0*hx;
    double sigmaY = 3.0*hy;

    unsigned int x0 = (unsigned int)(round(eta.x*Nx));
    unsigned int y0 = (unsigned int)(round(eta.y*Ny));

    double sumX = 0.0;
    for (unsigned int n=x0-3; n<=x0+3; n++)
    {
        sumX += exp(-((n*hx-eta.x)*(n*hx-eta.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=y0-3; m<=y0+3; m++)
    {
        sumY += exp(-((m*hy-eta.y)*(m*hy-eta.y))/(2.0*sigmaY*sigmaY));
    }
    sumY *= hy;

    double sigma = (sumX*sumY);// / (2.0*M_PI);
    double factor = 1.0/(/*(2.0*M_PI)*/sigma);

    for (unsigned int n=x0-3; n<=x0+3; n++)
    {
        for (unsigned int m=y0-3; m<=y0+3; m++)
        {
            ControlDeltaNode cn;
            cn.i = n; cn.x = n*hx; cn.j = m; cn.y = m*hy; cn.pt = eta; cn.id = i;
            cn.w = factor*exp(-((cn.x-eta.x)*(cn.x-eta.x))/(2.0*sigmaX*sigmaX)-((cn.y-eta.y)*(cn.y-eta.y))/(2.0*sigmaY*sigmaY));
            cps.push_back(cn);
        }
    }
}

void IProblem2Forward2D::extendObservationPoint4(const SpaceNodePDE &xi, std::vector<ObservationNode> &ons, unsigned int j) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();
    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double sigmaX = 3.0*hx;
    double sigmaY = 3.0*hy;

    unsigned int x0 = (unsigned int)(round(xi.x*Nx));
    unsigned int y0 = (unsigned int)(round(xi.y*Ny));

    double sumX = 0.0;
    for (unsigned int n=x0-3; n<=x0+3; n++)
    {
        sumX += exp(-((n*hx-xi.x)*(n*hx-xi.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=y0-3; m<=y0+3; m++)
    {
        sumY += exp(-((m*hy-xi.y)*(m*hy-xi.y))/(2.0*sigmaY*sigmaY));
    }
    sumY *= hy;

    double sigma = (sumX*sumY);// / (2.0*M_PI);
    double factor = 1.0/(/*(2.0*M_PI)*/sigma);

    for (unsigned int n=x0-3; n<=x0+3; n++)
    {
        for (unsigned int m=y0-3; m<=y0+3; m++)
        {
            ObservationNode on;
            on.i = n; on.x = n*hx; on.j = m; on.y = m*hy; on.pt = xi; on.id = j;
            on.w = factor*exp(-((on.x-xi.x)*(on.x-xi.x))/(2.0*sigmaX*sigmaX)-((on.y-xi.y)*(on.y-xi.y))/(2.0*sigmaY*sigmaY)) * (hx*hy);
            ons.push_back(on);
        }
    }
}
