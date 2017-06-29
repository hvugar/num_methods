#include "iloadedheatequation.h"

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
    double *kc = (double*) malloc(sizeof(double)*(N+1));
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
        kb[0] = 1.0 + 2.0*a_a_ht_hx_hx + 2.0*lambda1_a_a_ht_hx + lambda0_ht;
        kc[0] = -2.0*a_a_ht_hx_hx;
        kd[0] = u[0] + lambda0_ht_tt - 2.0*lambda1_a_a_ht_hx*sum + ht*f(isn,tn) - 2.0*aa*(ht/hx)*g(tn);
        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            //aa = a*a;

            ka[n] = -a_a_ht_hx_hx;
            kb[n] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0_ht;
            kc[n] = -a_a_ht_hx_hx;
            kd[n] = u[n] + lambda0_ht_tt + ht*f(isn,tn);
        }
        isn.i = N+minN;
        isn.x = isn.i*hx;
        //aa = a(isn,tn);

        ka[N] = -2.0*a_a_ht_hx_hx;
        kb[N] = 1.0+2.0*a_a_ht_hx_hx+2.0*lambda2_a_a_ht_hx + lambda0_ht;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0_ht_tt + 2.0*lambda2_a_a_ht_tt_hx + ht*f(isn,tn) + 2.0*aa*(ht/hx)*h(tn);

//        IPrinter::printVector(ka,N+1,"a");
//        IPrinter::printVector(kb,N+1,"b");
//        IPrinter::printVector(kc,N+1,"c");
//        IPrinter::printVector(kd,N+1,"d");

        double *betta = (double *)malloc(sizeof(double)*(N+1));
        for (unsigned int n=0; n<=N; n++) betta[n] = 0.0;
        double eta = kd[0];
        betta[0] = kb[0];
        betta[1] = kc[0];

        for (unsigned int s=0; s<L; s++)
        {
            betta[params[s].xi] = -2.0*lambda1*aa*(ht/hx)*params[s].k;
        }

//        DoubleMatrix MM(N+1, N+1);
//        DoubleVector AA(N+1);
//        DoubleVector XX(N+1);

//        for (unsigned int n=0; n<=N; n++)
//        {
//            MM[0][n] = betta[n];
//            AA[n] = kd[n];
//        }
//        for (unsigned int n=1; n<=N; n++)
//        {
//            MM[n][n-1] = ka[n];
//            MM[n][n+0] = kb[n];
//            MM[n][n+1] = kc[n];
//        }

//        GaussianElimination(MM, AA, XX);

//        for (unsigned int n=0; n<=N; n++) u[n] = XX[n];

        for (unsigned int n=1; n<=N-1; n++)
        {
            betta[n+0] -= betta[n-1]*(kb[n]/ka[n]);
            betta[n+1] -= betta[n-1]*(kc[n]/ka[n]);
            eta        -= betta[n-1]*(kd[n]/ka[n]);

//            if (n==1)
//            {
//                printf("%.18f %.18f\n", betta[n], betta[n+1]);
//            }
        }

//        IPrinter::printVector(betta, N+1, "betta");

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = betta[N-1]; M[0][1] = betta[N]; A[0] = eta;
        M[1][0] = ka[N];      M[1][1] = kb[N];    A[1] = kd[N];

        GaussianElimination(M,A,x);

        //printf("Eta %18.10f\n", A[0]);
        //IPrinter::print(M,M.rows(),M.cols());

        u[N-0] = x[1];
        u[N-1] = x[0];
        double sum1 = 0.0;
        for (unsigned int n=N-1; n!=0; n--)
        {
            eta        += betta[n-1]*(kd[n]/ka[n]);
            betta[n+1] += betta[n-1]*(kc[n]/ka[n]);
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
    free(kc);
    free(kd);
    free(rx);
}

void ILoadedHeatEquation::calculateM2(DoubleVector &u)
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
    double *kc = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));
    double *de = (double*) malloc(sizeof(double)*(N+1));

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
        kb[0] = 1.0 + 2.0*a_a_ht_hx_hx + 2.0*lambda1_a_a_ht_hx + lambda0_ht;
        kc[0] = -2.0*a_a_ht_hx_hx;
        kd[0] = u[0] + lambda0_ht_tt - 2.0*lambda1_a_a_ht_hx*sum + ht*f(isn,tn) - 2.0*aa*(ht/hx)*g(tn);
        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            //aa = a*a;

            ka[n] = -a_a_ht_hx_hx;
            kb[n] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0_ht;
            kc[n] = -a_a_ht_hx_hx;
            kd[n] = u[n] + lambda0_ht_tt + ht*f(isn,tn);
        }
        isn.i = N+minN;
        isn.x = isn.i*hx;
        //aa = a(isn,tn);

        ka[N] = -2.0*a_a_ht_hx_hx;
        kb[N] = 1.0+2.0*a_a_ht_hx_hx+2.0*lambda2_a_a_ht_hx + lambda0_ht;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0_ht_tt + 2.0*lambda2_a_a_ht_tt_hx + ht*f(isn,tn) + 2.0*aa*(ht/hx)*h(tn);

        double *betta = (double *)malloc(sizeof(double)*(N+1));
        for (unsigned int n=0; n<=N; n++) betta[n] = 0.0;
        double eta = kd[0];
        betta[0] = kb[0];
        betta[1] = kc[0];

        for (unsigned int s=0; s<L; s++)
        {
            betta[params[s].xi] = -2.0*lambda1_a_a_ht_hx*params[s].k;
        }

        for (unsigned int n=1; n<=N-1; n++)
        {
            betta[n+0] -= betta[n-1]*(kb[n]/ka[n]);
            betta[n+1] -= betta[n-1]*(kc[n]/ka[n]);
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
            betta[n+1] += betta[n-1]*(kc[n]/ka[n]);
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
        //break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}
