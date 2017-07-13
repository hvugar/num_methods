#include "iloadedheatequation.h"
#include <float.h>

void ILoadedHeatEquation::calculateM1(DoubleVector &u)
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.minN();
    unsigned int maxM = time.maxN();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.minN();
    unsigned int maxN = dim1.maxN();
    unsigned int N = maxN-minN;

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N+1));
    double *kb = (double*) malloc(sizeof(double)*(N+1));
    double *c1 = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));

    double a_a_ht_hx_hx = (a*a*ht)/(hx*hx);
    double lambda1_a_a_ht_hx = (lambda1*a*a*ht)/hx;
    double lambda2_a_a_ht_hx = (lambda2*a*a*ht)/hx;
    double lambda2_a_a_ht_tt_hx = (lambda2*a*a*ht*theta)/hx;
    double lambda0_ht = lambda0*ht;
    double lambda0_ht_tt = lambda0*ht*theta;

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn);
    }
    layerInfo(u, 0);

    TimeNode tn;
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

        ka[0] = 0.0;
        kb[0] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht + 2.0*lambda1_a_a_ht_hx;
        c1[0] = -2.0*a_a_ht_hx_hx;
        kd[0] = u[0] + lambda0_ht_tt + ht*f(isn,tn) - 2.0*lambda1_a_a_ht_hx*sum - 2.0*aa*(ht/hx)*g(tn);
        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            //aa = a*a;

            ka[n] = -a_a_ht_hx_hx;
            kb[n] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0_ht;
            c1[n] = -a_a_ht_hx_hx;
            kd[n] = u[n] + lambda0_ht_tt + ht*f(isn,tn);
        }
        isn.i = N+minN;
        isn.x = isn.i*hx;
        //aa = a(isn,tn);

        ka[N] = -2.0*a_a_ht_hx_hx;
        kb[N] = 1.0 + 2.0*a_a_ht_hx_hx + lambda0_ht + 2.0*lambda2_a_a_ht_hx;
        c1[N] = 0.0;
        kd[N] = u[N] + lambda0_ht_tt + ht*f(isn,tn) + 2.0*lambda2_a_a_ht_tt_hx + 2.0*aa*(ht/hx)*h(tn);

        double *betta = (double *)malloc(sizeof(double)*(N+1));
        for (unsigned int n=0; n<=N; n++) betta[n] = 0.0;
        double eta = kd[0];
        betta[0] = kb[0];
        betta[1] = c1[0];

        for (unsigned int s=0; s<L; s++)
        {
            betta[params[s].xi] = -2.0*lambda1*aa*(ht/hx)*params[s].k;
        }

        for (unsigned int n=1; n<=N-1; n++)
        {
            betta[n+0] -= betta[n-1]*(kb[n]/ka[n]);
            betta[n+1] -= betta[n-1]*(c1[n]/ka[n]);
            eta        -= betta[n-1]*(kd[n]/ka[n]);
        }

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = betta[N-1]; M[0][1] = betta[N]; A[0] = eta;
        M[1][0] = ka[N];      M[1][1] = kb[N];    A[1] = kd[N];

        GaussianElimination(M,A,x);

        u[N-0] = x[1];
        u[N-1] = x[0];
        double sum1 = 0.0;
        for (unsigned int n=N-1; n!=0; n--)
        {
            eta        += betta[n-1]*(kd[n]/ka[n]);
            betta[n+1] += betta[n-1]*(c1[n]/ka[n]);
            betta[n+0] += betta[n-1]*(kb[n]/ka[n]);

            //u[n-1] = eta; for (unsigned int i=n; i<=N; i++) u[n-1] -= betta[i]*u[i];
            u[n-1] = -betta[n+0]*u[n+0] - betta[n+1]*u[n+1] + eta - sum1;
            u[n-1] /= betta[n-1];

            if (n<=N-1)
            {
                sum1 += betta[n+1]*u[n+1];
                //printf("%d %18.6f %18.6f %18.6f %20.6f %20.6f\n", n, sum1, betta[n+1], u[n+1], betta[n+0], betta[n-1]);
            }
        }
        //printf("%d %18.6f %18.6f %18.6f %20.6f %20.6f\n", 0, sum1, betta[1], u[1], betta[0], NAN);

        free(betta);

        layerInfo(u, m);
        break;
    }

    free(ka);
    free(kb);
    free(c1);
    free(kd);
    free(rx);
}

void ILoadedHeatEquation::calculateM2(DoubleVector &u)
{
    /* Grid parameters */

    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.minN();
    unsigned int maxM = time.maxN();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.minN();
    unsigned int maxN = dim1.maxN();
    unsigned int N = maxN-minN;

    /* initializing vector and coifficients */

    u.clear();
    u.resize(N+1);

    double *a1 = (double*) malloc(sizeof(double)*(N+1));
    double *b1 = (double*) malloc(sizeof(double)*(N+1));
    double *c1 = (double*) malloc(sizeof(double)*(N+1));
    double *d1 = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));

    /* first row */
    double *e1 = (double*) malloc(sizeof(double)*(N+1));

    double a_a_ht_hx_hx         = (a*a*ht)/(hx*hx);
    double lambda1_a_a_ht_hx    = (lambda1*a*a*ht)/hx;
    double lambda2_a_a_ht_hx    = (lambda2*a*a*ht)/hx;
    double lambda2_a_a_ht_tt_hx = (lambda2*a*a*ht*theta)/hx;
    double lambda0_ht           = lambda0*ht;
    double lambda0_ht_tt        = lambda0*ht*theta;

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn);
    }
    layerInfo(u, 0);

    TimeNode tn;
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
                //double diff = fabs(n*hx - params[s].e);
                //if (diff <= hx)
                if (n == params[s].xi)
                {
                    e1[n] += params[s].k;// * (1.0 - diff/hx);
                }
             }
             e1[n] *= -2.0*lambda1_a_a_ht_hx;
        }

        /* calculating... */

        /* finding loaded points count which index >= 2 */
        unsigned int LPC = L;
        //for (unsigned int n=2; n<=N; n++) if ( fabs(e1[n]) > FLT_EPSILON ) LPC++;

        double* p = (double*) malloc(sizeof(double) * (N+1));
        double* q = (double*) malloc(sizeof(double) * (N+1));

        double** k = (double**) malloc(sizeof(double*) * (N+1));
        for (unsigned int n=0; n<=N; n++) k[n] = (double*) malloc(sizeof(double)*LPC);

        //unsigned int *kI = (unsigned int*) malloc(sizeof(unsigned int) * LPC);
        //int kIndex = 0;
        //for (unsigned int n=2; n<=N; n++)
        //{
        //    if (fabs(e1[n]) > FLT_EPSILON)
        //    {
        //        kI[kIndex] = n;
        //        kIndex++;
        //    }
        //}

        IPrinter::printVector(a1,N+1);
        IPrinter::printVector(b1,N+1);
        IPrinter::printVector(c1,N+1);
        IPrinter::printVector(d1,N+1);
        IPrinter::printSeperatorLine();

        printf("%18.10f\n", b1[0]);

        //printf("%4d c: %18.10f", 0, c1[0]); for (unsigned int s=0; s<LPC; s++) printf("k:%4d %18.10f", kI[s], e1[kI[s]]); printf("q:%18.10f\n", d1[0]);
        printf("%4d c: %18.10f", 0, c1[0]); for (unsigned int s=0; s<LPC; s++) printf("k:%4d %18.10f", params[s].xi, e1[params[s].xi]); printf("q:%18.10f\n", d1[0]);
        IPrinter::printSeperatorLine();

        for (unsigned int n=0; n<=N; n++)
        {
            if (n==0)
            {
                b1[0] += e1[0];
                c1[0] += e1[1];

                p[0] = -c1[0]/b1[0];
                q[0] = +d1[0]/b1[0];

                //for (unsigned int s=0; s<LPC; s++) k[0][s] = -e1[kI[s]]/b1[0];
                for (unsigned int s=0; s<LPC; s++) k[0][s] = -e1[params[s].xi]/b1[0];

                //printf("%4d p:%18.10f", 0, p[0]); for (unsigned int s=0; s<LPC; s++) printf("k:%4d %18.10f", kI[s], k[0][s]);  printf(" q:%18.10f\n", q[0]);
                printf("%4d p:%18.10f", 0, p[0]); for (unsigned int s=0; s<LPC; s++) printf("k:%4d %18.10f", params[s].xi, k[0][s]);  printf(" q:%18.10f\n", q[0]);
            }
            else if (n==N)
            {
                double m = b1[n] + a1[n]*p[n-1];
                q[n] = (d1[n] - a1[n]*q[n-1])/m;
                p[n] = 0.0;
                //printf("%4d p:%18.10f", n, p[n]); for (unsigned int s=0; s<LPC; s++) printf("k:%4d %18.10f", kI[s], k[n][s]); printf(" q:%18.10f\n", q[n]);
                printf("%4d p:%18.10f", n, p[n]); for (unsigned int s=0; s<LPC; s++) printf("k:%4d %18.10f", params[s].xi, k[n][s]); printf(" q:%18.10f\n", q[n]);
            }
            else
            {
                double m = b1[n] + a1[n]*p[n-1];
                double c = +c1[n];

                for (unsigned int s=0; s<LPC; s++)
                {
                    //if (kI[s] == n-1) c += a1[n]*k[n-1][s];
                    if (params[s].xi == n-1) c += a1[n]*k[n-1][s];
                }

                p[n] = -c/m;
                q[n] = (d1[n] - a1[n]*q[n-1])/m;

                for (unsigned int s=0; s<LPC; s++)
                {
                    //if (kI[s] <= n+1)
                    if (params[s].xi <= n+1)
                    {
                        k[n][s] = 0.0;
                    }
                    else
                    {
                        k[n][s] = -(a1[n]*k[n-1][s])/m;
                    }
                }
                //printf("%4d p:%18.10f", n, p[n]); for (unsigned int s=0; s<LPC; s++) printf("k:%4d %18.10f", kI[s], k[n][s]); printf(" q:%18.10f\n", q[n]);
                printf("%4d p:%18.10f", n, p[n]); for (unsigned int s=0; s<LPC; s++) printf("k:%4d %18.10f", params[s].xi, k[n][s]); printf(" q:%18.10f\n", q[n]);
            }
        }

        u[N] = q[N];

        for (unsigned int n=N; n!=0; n--)
        {
            u[n-1] = p[n-1]*u[n] + q[n-1];
            for (unsigned int s=0; s<LPC; s++)
            {
                //u[n-1] += u[kI[s]] * k[n-1][s];
            }
        }

        //free(kI);
        for (unsigned int n=0; n<=N; n++) free(k[n]);
        free(p);


        layerInfo(u, m);
        break;
    }

    free(a1);
    free(b1);
    free(c1);
    free(d1);
    free(rx);
}

void ILoadedHeatEquation::qovmaFirstRowM(double *e, double *a, double *b, double *c, double *d, double *x, unsigned int N) const
{
//    unsigned int L = 0;
//    for (unsigned int n=2; n<N; i++)
//    {
//        if (e[n]!=0.0) L++;
//    }

//    double **p = (double**) malloc(sizeof(double*) * (L+1));
//    for (unsigned int i=0; i<=L; i++)
//    {
//        p[i] = (double*) malloc(sizeof(double) * N);
//    }
//    double *q = (double*) malloc(sizeof(double) * N);

//    for (unsigned int i=0; i<N; i++)
//    {
//        if (i == 0)
//        {
//            p[0] = -e[1]/e[0];
//            for (unsigned int s=0; s<CNT; s++)
//            {
//                k[s][0] = -e[E[s]]/e[0];
//            }
//            q[0] = d[0]/e[0];
//            e[0] = 1.0;
//        }
//        else if (i == N-1)
//        {
//            q[i] = +(d[i]-a[i]*q[i-1])/(b[i]+a[i]*p[i-1]);
//            p[i] = 0.0;
//            for (unsigned int s=0; s<CNT; s++)
//            {
//                k[s][i] = 0.0;
//            }
//        }
//        else
//        {
//            double m = b[i]+a[i]*p[i-1];
//            p[i] = -c[i]/m;
//            q[i] = (d[i] - a[i]*q[i-1])/m;
//            for (unsigned int s=0; s<CNT; s++)
//            {
//                if (i<(E[s]-1))
//                {
//                    k[s][i] = -(a[i]*k[s][i-1])/m;
//                }
//                else
//                {
//                    k[s][i] = 0.0;
//                }
//                if (i==E[s]-1)
//                {
//                    p[i] += -(a[i]*k[s][i-1])/m;
//                }
//            }
//        }
//    }

//    for (unsigned int i=N-1; i != UINT_MAX; i--)
//    {
//        if (i==(N-1))
//        {
//            x[i] = q[i];
//        }
//        else
//        {
//            x[i] = q[i] + p[i]*x[i+1];
//            for (unsigned int s=0; s<CNT; s++)
//            {
//                if (i<=E[s]-1)
//                {
//                    x[i] = x[i] + k[s][i]*x[E[s]];
//                }
//            }
//        }
//    }

//    for (unsigned int i=0; i<=L; i++) free(p[i]);
//    free(p);
//    free(q);
}
