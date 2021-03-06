#include "iloadedheatequation.h"
#include <float.h>

void ILoadedHeatEquation::getGridParameters(double &hx, unsigned int &N, double &ht, unsigned int &M,
                                            unsigned int &minN, unsigned int &maxN,
                                            unsigned int &minM, unsigned int &maxM) const
{
    Dimension time = timeDimension();
    Dimension dim1 = spaceDimensionX();

    ht = time.step();
    minM = time.min();
    maxM = time.max();
    M = maxM-minM;

    hx = dim1.step();
    minN = dim1.min();
    maxN = dim1.max();
    N = maxN-minN;
}

void ILoadedHeatEquation::transferringConditions2(DoubleVector &u)
{
    /* Grid parameters */
    double hx, ht;
    unsigned int N, M, minN, maxN, minM, maxM;
    getGridParameters(hx, N, ht, M, minN, maxN, minM, maxM);

//    Dimension time = mtimeDimension;
//    Dimension dim1 = mspaceDimension.at(0);

//    double ht = time.step();
//    unsigned int minM = time.minN();
//    unsigned int maxM = time.maxN();
//    unsigned int M = maxM-minM;

//    double hx = dim1.step();
//    unsigned int minN = dim1.minN();
//    unsigned int maxN = dim1.maxN();
//    unsigned int N = maxN-minN;

    /* initializing vector and coifficients */

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N+1));
    double *kb = (double*) malloc(sizeof(double)*(N+1));
    double *kc = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *ke = (double*) malloc(sizeof(double)*(N+1));

    /* first row */
    double *betta = (double *) malloc(sizeof(double)*(N+1));

    double a_a_ht_hx_hx         = (a*a*ht)/(hx*hx);
    double lambda1_a_a_ht_hx    = (a*a*lambda1*ht)/hx;
    double lambda2_a_a_ht_hx    = (a*a*lambda2*ht)/hx;
    double lambda2_a_a_ht_tt_hx = (a*a*lambda2*ht*theta)/hx;
    double lambda0_ht           = lambda0*ht;
    double lambda0_ht_tt        = lambda0*ht*theta;

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn, InitialCondition::InitialValue);
    }
    //layerInfo(u, 0);

    TimeNodePDE tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        double aa = 0.0;

        isn.i = 0+minN;
        isn.x = isn.i*hx;
        aa = a*a;

        double sum = 0.0;
        for (unsigned int s=0; s<L; s++) sum += params[s].k * params[s].z;

        /* preparing first row */
        ka[0] = 0.0;
        kb[0] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht + 2.0*lambda1_a_a_ht_hx;
        kc[0] = -2.0*a_a_ht_hx_hx;
        kd[0] = u[0] + lambda0_ht_tt + ht*f(isn,tn) - 2.0*lambda1_a_a_ht_hx*sum - 2.0*aa*(ht/hx)*g(tn);

        /* preparing intermadiate rows */

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            //aa = a*a;

            ka[n] = -a_a_ht_hx_hx;
            kb[n] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht;
            kc[n] = -a_a_ht_hx_hx;
            kd[n] = u[n] + lambda0_ht_tt + ht*f(isn,tn);
        }

        /* preparing last row */
        isn.i = N+minN;
        isn.x = isn.i*hx;
        //aa = a(isn,tn);

        ka[N] = -2.0*a_a_ht_hx_hx;
        kb[N] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht + 2.0*lambda2_a_a_ht_hx;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0_ht_tt + ht*f(isn,tn) + 2.0*lambda2_a_a_ht_tt_hx + 2.0*(aa*ht*h(tn))/hx;

        for (unsigned int n=0; n<=N; n++)
        {
            ke[n] = 0.0;
            for (unsigned int s=0; s<L; s++)
            {
                double diff = fabs(n*hx - params[s].e);
                if (diff <= hx) ke[n] += params[s].k * (1.0 - diff/hx);
             }
             ke[n] *= -2.0*lambda1_a_a_ht_hx;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        /* calculating... */

        for (unsigned int n=0; n<=N; n++) betta[n] = ke[n];
        betta[0] += kb[0];
        betta[1] += kc[0];
        double eta = kd[0];

        //IPrinter::printVector(betta,N+1, "betta ");

        unsigned int nonZeroCount = 0;
        unsigned int *nonZeroIndex0 = (unsigned int*) malloc(sizeof(unsigned int) * (N+1));
        for (unsigned int n=0; n<=N; n++)
        {
            if ( betta[n] != 0.0 )
            {
                nonZeroIndex0[nonZeroCount] = n;
                nonZeroCount++;
            }
        }
        unsigned int *nonZeroIndex = (unsigned int*) malloc(nonZeroCount*sizeof(unsigned int));
        for (unsigned int n=0; n<nonZeroCount; n++)
        {
            nonZeroIndex[n] = nonZeroIndex0[n];
        }
        free(nonZeroIndex0);

        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        DoubleMatrix ems(N+1, nonZeroCount+1);

        for (unsigned int n=0; n<nonZeroCount; n++)
        {
            ems[0][n] = betta[ nonZeroIndex[n] ];
        }
        ems[0][nonZeroCount] = eta;

        for (unsigned int n=1; n<=N-1; n++)
        {
            ems[n][0] = -ems[n-1][0]*(kb[n]/ka[n]) + ems[n-1][1];
            ems[n][1] = -ems[n-1][0]*(kc[n]/ka[n]);

            for (unsigned int s=2; s<nonZeroCount; s++)
            {
                if ( n+1 == nonZeroIndex[s] )
                {
                    ems[n][1] += ems[n-1][s];
                }

                if ( n+1 < nonZeroIndex[s] )
                {
                    ems[n][s] = ems[n-1][s];
                }
                else
                {
                    ems[n][s] = 0.0;
                }
            }
            ems[n][nonZeroCount] = -ems[n-1][0]*(kd[n]/ka[n]) + ems[n-1][nonZeroCount];
        }

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = ems[N-1][0]; M[0][1] = ems[N-1][1]; A[0] = ems[N-1][nonZeroCount] ;
        M[1][0] = ka[N];       M[1][1] = kb[N];       A[1] = kd[N];

        LinearEquation::GaussianElimination(M,A,x);

        for (unsigned int n=0; n<=N; n++) u[n] = 0.0;
        u[N-0] = x[1];
        u[N-1] = x[0];
        for (unsigned int n=N-2; n!=UINT32_MAX; n--)
        {
            double bb = 0.0;
            for (unsigned int s=2; s<nonZeroCount; s++)
            {
                bb += ems[n][s]*u[nonZeroIndex[s]];
            }
            u[n] = -(ems[n][1]*u[n+1] + bb - ems[n][nonZeroCount])/ems[n][0];
        }

        ems.clear();
        free(nonZeroIndex);

        //layerInfo(u, m);
    }

    free(betta);
    free(ke);
    free(rx);
    free(kd);
    free(kc);
    free(kb);
    free(ka);
}

void ILoadedHeatEquation::transferringConditions4(DoubleVector &u)
{
    /* Grid parameters */

    Dimension time = timeDimension();
    Dimension dim1 = spaceDimensionX();

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN-minN;

    /* initializing vector and coifficients */

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N+1));
    double *kb = (double*) malloc(sizeof(double)*(N+1));
    double *kc = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *ke = (double*) malloc(sizeof(double)*(N+1));

    /* first row */
    double *betta = (double *)malloc(sizeof(double)*(N+1));

    double a_a_ht_hx_hx         = (a*a*ht)/(hx*hx);
    double lambda1_a_a_ht_hx    = (a*a*lambda1*ht)/hx;
    double lambda2_a_a_ht_hx    = (a*a*lambda2*ht)/hx;
    double lambda2_a_a_ht_tt_hx = (a*a*lambda2*ht*theta)/hx;
    double lambda0_ht           = lambda0*ht;
    double lambda0_ht_tt        = lambda0*ht*theta;

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn, InitialCondition::InitialValue);
    }
    //layerInfo(u, 0);

    TimeNodePDE tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        double aa = 0.0;

        isn.i = 0+minN;
        isn.x = isn.i*hx;
        aa = a*a;

        double sum = 0.0;
        for (unsigned int s=0; s<L; s++) sum += params[s].k * params[s].z;

        /* preparing first row */
        ka[0] = 0.0;
        kb[0] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht + 2.0*lambda1_a_a_ht_hx;
        kc[0] = -2.0*a_a_ht_hx_hx;
        kd[0] = u[0] + lambda0_ht_tt + ht*f(isn,tn) - 2.0*lambda1_a_a_ht_hx*sum - 2.0*aa*(ht/hx)*g(tn);

        /* preparing intermadiate rows */

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            //aa = a*a;

            ka[n] = -a_a_ht_hx_hx;
            kb[n] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht;
            kc[n] = -a_a_ht_hx_hx;
            kd[n] = u[n] + lambda0_ht_tt + ht*f(isn,tn);
        }

        /* preparing last row */
        isn.i = N+minN;
        isn.x = isn.i*hx;
        //aa = a(isn,tn);

        ka[N] = -2.0*a_a_ht_hx_hx;
        kb[N] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht + 2.0*lambda2_a_a_ht_hx;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0_ht_tt + ht*f(isn,tn) + 2.0*lambda2_a_a_ht_tt_hx + 2.0*(aa*ht*h(tn))/hx;

        for (unsigned int n=0; n<=N; n++)
        {
            ke[n] = 0.0;
            for (unsigned int s=0; s<L; s++)
            {
                double diff = fabs(n*hx - params[s].e);
                if (diff <= hx) ke[n] += params[s].k * (1.0 - diff/hx);
             }
             ke[n] *= -2.0*lambda1_a_a_ht_hx;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        /* calculating... */

        for (unsigned int n=0; n<=N; n++) betta[n] = ke[n];
        betta[0] += kb[0];
        betta[1] += kc[0];
        double eta = kd[0];

        //IPrinter::printVector(betta,N+1, "betta ");

        unsigned int nonZeroCount = 0;
        unsigned int *nonZeroIndex0 = (unsigned int*) malloc(sizeof(unsigned int) * (N+1));
        for (unsigned int n=0; n<=N; n++)
        {
            if ( betta[n] != 0.0 )
            {
                nonZeroIndex0[nonZeroCount] = n;
                nonZeroCount++;
            }
        }
        unsigned int *nonZeroIndex = (unsigned int*) malloc(nonZeroCount*sizeof(unsigned int));
        for (unsigned int n=0; n<nonZeroCount; n++)
        {
            nonZeroIndex[n] = nonZeroIndex0[n];
        }
        free(nonZeroIndex0);

        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        DoubleMatrix ems(N+1, nonZeroCount+1);

        for (unsigned int n=0; n<nonZeroCount; n++)
        {
            ems[0][n] = betta[ nonZeroIndex[n] ];
        }
        ems[0][nonZeroCount] = eta;

        for (unsigned int n=1; n<=N-1; n++)
        {
            ems[n][0] = -ems[n-1][0]*(kb[n]/ka[n]) + ems[n-1][1];
            ems[n][1] = -ems[n-1][0]*(kc[n]/ka[n]);

            for (unsigned int s=2; s<nonZeroCount; s++)
            {
                if ( n+1 == nonZeroIndex[s] )
                {
                    ems[n][1] += ems[n-1][s];
                }

                if ( n+1 < nonZeroIndex[s] )
                {
                    ems[n][s] = ems[n-1][s];
                }
                else
                {
                    ems[n][s] = 0.0;
                }
            }
            ems[n][nonZeroCount] = -ems[n-1][0]*(kd[n]/ka[n]) + ems[n-1][nonZeroCount];
        }

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = ems[N-1][0]; M[0][1] = ems[N-1][1]; A[0] = ems[N-1][nonZeroCount] ;
        M[1][0] = ka[N];       M[1][1] = kb[N];       A[1] = kd[N];

        LinearEquation::GaussianElimination(M,A,x);

        for (unsigned int n=0; n<=N; n++) u[n] = 0.0;
        u[N-0] = x[1];
        u[N-1] = x[0];
        for (unsigned int n=N-2; n!=UINT32_MAX; n--)
        {
            double bb = 0.0;
            for (unsigned int s=2; s<nonZeroCount; s++)
            {
                bb += ems[n][s]*u[nonZeroIndex[s]];
            }
            u[n] = -(ems[n][1]*u[n+1] + bb - ems[n][nonZeroCount])/ems[n][0];
        }

        ems.clear();
        free(nonZeroIndex);

        //layerInfo(u, m);
    }

    free(betta);
    free(ke);
    free(rx);
    free(kd);
    free(kc);
    free(kb);
    free(ka);
}

void ILoadedHeatEquation::calculateM2(DoubleVector &u)
{
    /* Grid parameters */

    Dimension time = timeDimension();
    Dimension dim1 = spaceDimensionX();

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN-minN;

    /* initializing vector and coifficients */

    u.clear();
    u.resize(N+1);

    double *a1 = (double*) malloc(sizeof(double)*(N+1));
    double *b1 = (double*) malloc(sizeof(double)*(N+1));
    double *c1 = (double*) malloc(sizeof(double)*(N+1));
    double *d1 = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *e1 = (double*) malloc(sizeof(double)*(N+1));

    double a_a_ht_hx_hx         = (a*a*ht)/(hx*hx);
    double lambda1_a_a_ht_hx    = (a*a*lambda1*ht)/hx;
    double lambda2_a_a_ht_hx    = (a*a*lambda2*ht)/hx;
    double lambda2_a_a_ht_tt_hx = (a*a*lambda2*ht*theta)/hx;
    double lambda0_ht           = lambda0*ht;
    double lambda0_ht_tt        = lambda0*ht*theta;

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn, InitialCondition::InitialValue);
    }
    //layerInfo(u, 0);

    TimeNodePDE tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        double aa = 0.0;

        isn.i = 0+minN;
        isn.x = isn.i*hx;
        aa = a*a;

        double sum = 0.0;
        for (unsigned int s=0; s<L; s++) sum += params[s].k * params[s].z;

        /* preparing first row */
        a1[0] = 0.0;
        b1[0] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht + 2.0*lambda1_a_a_ht_hx;
        c1[0] = -2.0*a_a_ht_hx_hx;
        d1[0] = u[0] + lambda0_ht_tt + ht*f(isn,tn) - 2.0*lambda1_a_a_ht_hx*sum - 2.0*aa*(ht/hx)*g(tn);

        /* preparing intermadiate rows */

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            //aa = a*a;

            a1[n] = -a_a_ht_hx_hx;
            b1[n] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0_ht;
            c1[n] = -a_a_ht_hx_hx;
            d1[n] = u[n] + lambda0_ht_tt + ht*f(isn,tn);
        }

        /* preparing last row */
        isn.i = N+minN;
        isn.x = isn.i*hx;
        //aa = a(isn,tn);

        a1[N] = -2.0*a_a_ht_hx_hx;
        b1[N] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht + 2.0*lambda2_a_a_ht_hx;
        c1[N] = 0.0;
        d1[N] = u[N] + lambda0_ht_tt + ht*f(isn,tn) + 2.0*lambda2_a_a_ht_tt_hx + 2.0*aa*(ht/hx)*h(tn);

        for (unsigned int n=0; n<=N; n++)
        {
            e1[n] = 0.0;
            for (unsigned int s=0; s<L; s++)
            {
                double diff = fabs(n*hx - params[s].e);
                if (diff <= hx) e1[n] += params[s].k * (1.0 - diff/hx);
             }
             e1[n] *= -2.0*lambda1_a_a_ht_hx;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        /* calculating... */

        /* finding loaded points count which index >= 2 */
        unsigned int LPC = 0;
        for (unsigned int n=2; n<=N; n++) if ( fabs(e1[n]) > DBL_EPSILON ) LPC++;

        double * p = (double*) malloc(sizeof(double) * (N+1));
        double * q = (double*) malloc(sizeof(double) * (N+1));
        double **r = (double**) malloc(sizeof(double*) * LPC);
        for (unsigned int s=0; s<LPC; s++) r[s] = (double*) malloc(sizeof(double) * (N+1));

        unsigned int *kI = (unsigned int*) malloc(sizeof(unsigned int) * LPC);
        int kIndex = 0;
        for (unsigned int n=2; n<=N; n++)
        {
            if (fabs(e1[n]) > DBL_EPSILON)
            {
                kI[kIndex] = n;
                kIndex++;
            }
        }

        for (unsigned int n=0; n<=N; n++)
        {
            if (n==0)
            {
                b1[0] += e1[0];
                c1[0] += e1[1];

                p[0] = -c1[0]/b1[0];
                q[0] = +d1[0]/b1[0];

                for (unsigned int s=0; s<LPC; s++) r[s][0] = -e1[kI[s]]/b1[0];
            }
            else if (n==N)
            {
                double m = b1[n] + a1[n]*p[n-1];
                q[n] = (d1[n] - a1[n]*q[n-1])/m;
                p[n] = 0.0;
            }
            else
            {
                double m = b1[n] + a1[n]*p[n-1];
                double c = +c1[n];
                for (unsigned int s=0; s<LPC; s++)
                {
                    if (kI[s] == n+1) c += a1[n]*r[s][n-1];
                }

                p[n] = -c/m;
                q[n] = (d1[n] - a1[n]*q[n-1])/m;

                for (unsigned int s=0; s<LPC; s++)
                {
                    if (kI[s] <= n+1)
                    {
                        r[s][n] = 0.0;
                    }
                    else
                    {
                        r[s][n] = -(a1[n]*r[s][n-1])/m;
                    }
                }
            }
        }

        u[N] = q[N];

        for (unsigned int n=N; n!=0; n--)
        {
            u[n-1] = p[n-1]*u[n] + q[n-1];
            for (unsigned int s=0; s<LPC; s++)
            {
                u[n-1] += u[kI[s]] * r[s][n-1];
            }
        }

        free(kI);
        for (unsigned int s=0; s<LPC; s++) free(r[s]);
        free(r);
        free(p);
        free(q);

        //layerInfo(u, m);
    }

    free(rx);
    free(d1);
    free(c1);
    free(b1);
    free(a1);
}

