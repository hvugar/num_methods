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
