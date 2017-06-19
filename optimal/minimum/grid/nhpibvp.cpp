#include "nhpibvp.h"

void NewtonHeatEquation::calculateGM1(DoubleVector &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

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

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn);
    }
    //layerInfo(u, 0);

    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        double aa = 0.0;

        isn.i = 0+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);

        ka[0] = 0.0;
        kb[0] = 1.0+aa*(ht/(hx*hx))+lambda1*aa*(ht/hx) + lambda0*ht;
        kc[0] = -aa*(ht/(hx*hx));
        kd[0] = u[0] + lambda0*ht*theta0(tn) + lambda1*aa*(ht/hx)*theta1(tn) + ht*f(isn,tn);
        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            aa = a(isn,tn);

            ka[n] = -aa*(ht/(hx*hx));
            kb[n] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
            kc[n] = -aa*(ht/(hx*hx));
            kd[n] = u[n] + lambda0*ht*theta0(tn) + ht*f(isn,tn);
        }
        isn.i = N+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);

        ka[N] = -aa*(ht/(hx*hx));
        kb[N] = 1.0+aa*(ht/(hx*hx))+lambda2*aa*(ht/hx) + lambda0*ht;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0*ht*theta0(tn) + lambda2*aa*(ht/hx)*theta2(tn) + ht*f(isn,tn);

        (*algorithm)(ka, kb, kc, kd, rx, N+1);

        for (unsigned int n=0; n<=N; n++) u[n] = rx[n];

        //layerInfo(u, m);
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void NewtonHeatEquation::calculateGM2(DoubleVector &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

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

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn);
    }
    //layerInfo(u, 0);

    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        double aa = 0.0;

        isn.i = 0+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);

        ka[0] = 0.0;
        kb[0] = 1.0+aa*(ht/(hx*hx))+lambda1*aa*(ht/hx) + lambda0*ht;
        kc[0] = -aa*(ht/(hx*hx));
        kd[0] = u[0] + lambda0*ht*theta0(tn) + lambda1*aa*(ht/hx)*theta1(tn) + ht*f(isn,tn);
        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            aa = a(isn,tn);

            ka[n] = -aa*(ht/(hx*hx));
            kb[n] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
            kc[n] = -aa*(ht/(hx*hx));
            kd[n] = u[n] + lambda0*ht*theta0(tn) + ht*f(isn,tn);
        }
        isn.i = N+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);

        ka[N] = -aa*(ht/(hx*hx));
        kb[N] = 1.0+aa*(ht/(hx*hx))+lambda2*aa*(ht/hx) + lambda0*ht;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0*ht*theta0(tn) + lambda2*aa*(ht/hx)*theta2(tn) + ht*f(isn,tn);

        //(*algorithm)(ka, kb, kc, kd, rx, N+1);

        double *betta = (double *)malloc(sizeof(double)*(N+1));
        double eta = kd[0];
        betta[0] = kb[0];
        betta[1] = kc[0];

        for (unsigned int n=1; n<=N-1; n++)
        {
            betta[n+0] = -betta[n-1]*(kb[n]/ka[n]) + betta[n];
            betta[n+1] = -betta[n-1]*(kc[n]/ka[n]);
            eta = eta - betta[n-1]*(kd[n]/ka[n]);
        }

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = betta[N-1]; M[0][1] = betta[N]; A[0] = eta;
        M[1][0] = ka[N];      M[1][1] = kb[N];    A[1] = kd[N];

        GaussianElimination(M,A,x);

        //printf("%18.10f %18.10f %18.10f\n", M[0][0], M[0][1], A[0]);
        //printf("%18.10f %18.10f %18.10f\n", M[1][0], M[1][1], A[1]);
        //IPrinter::printSeperatorLine();
        //printf("%18.10f %18.10f\n", x[0], x[1]);
        //IPrinter::printSeperatorLine();
        //printf("%18.10f %18.10f\n", M[0][0] * (hx*(N-1))*(m*ht) + M[0][1] * (hx*(N+0))*(m*ht), eta);
        //printf("%18.10f %18.10f\n", M[1][0] * (hx*(N-1))*(m*ht) + M[1][1] * (hx*(N+0))*(m*ht), kd[N]);
        //IPrinter::printSeperatorLine();

        //for (unsigned int n=0; n<=N; n++) u[n] = rx[n];
        u[N-0] = x[1];
        u[N-1] = x[0];
        for (unsigned int n=N-1; n!=1; n--)
        {
            //betta[n+1] = +betta[n-1]*(kc[n]/ka[n]);
            betta[n] = betta[n] + betta[n-1]*(kb[n]/ka[n]);
            eta      = eta      + betta[n-1]*(kd[n]/ka[n]);
            u[n-1] = (eta-betta[n]*u[n])/betta[n-1];
        }

        //IPrinter::printVector(18,10,u);

        //layerInfo(u, m);
        //break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}
