#include "bpibvp.h"

double BackwardParabolicIBVP::boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    C_UNUSED(boundary);
    return 0.0;
}

void BackwardParabolicIBVP::calculate(DoubleMatrix &psi, SweepMethodDirection direction)
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

    double h = ht/(hx*hx);

    psi.clear();
    psi.resize(M+1, N+1);

    double *ka = (double*) malloc(sizeof(double)*(N-1));
    double *kb = (double*) malloc(sizeof(double)*(N-1));
    double *kc = (double*) malloc(sizeof(double)*(N-1));
    double *kd = (double*) malloc(sizeof(double)*(N-1));
    double *rx = (double*) malloc(sizeof(double)*(N-1));

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        psi[0][n] = initial(isn);
    }

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    for (unsigned int m=M-1; m!=UINT32_MAX; m--)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            double alpha = -a(isn,tn)*h;
            double betta = 1.0 - 2.0*alpha;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = psi[m+1][n] - ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        psi[m][0] = boundary(lsn, tn, Left);
        psi[m][N] = boundary(rsn, tn, Right);

        kd[0]   -= a(lsn,tn) * h * psi[m][0];
        kd[N-2] -= a(rsn,tn) * h * psi[m][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) psi[m][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}
