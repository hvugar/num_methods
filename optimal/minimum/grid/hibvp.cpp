#include "hibvp.h"

void HyperbolicIBVP::gridMethod0(DoubleMatrix &u, SweepMethodDirection direction)
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

    double lambda = 0.25;
    double ht2 = ht*ht;
    double hx2 = hx*hx;
    double h = ht2/hx2;

    u.clear();
    u.resize(M+1, N+1);

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

        u[0][n] = initial1(isn);
        u[1][n] = u[0][n] + ht*initial2(isn);

        //u[2][n] = u[0][n] + 2.0*ht*initial2(isn);

        //TimeNode tn; tn.i = 0.0; tn.t = tn.i*ht;
        //u[2][n] = u[0][n] + ht*initial2(isn) + (2.0*a(isn,tn)+f(isn,tn))*((ht*ht)/2.0);
    }
    IPrinter::printSeperatorLine();
    IPrinter::printVector(14,10,u.row(0));
    IPrinter::printVector(14,10,u.row(1));
    IPrinter::printVector(14,10,u.row(2));
    IPrinter::printSeperatorLine();

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    TimeNode tn1;
    TimeNode tn2;
    for (unsigned int m=1; m<=M-1; m++)
    {
        tn.i = m+minM+1;
        tn.t = tn.i*ht;

        tn1.i = m+minM;
        tn1.t = tn1.i*ht;

        tn2.i = m+minM-1;
        tn2.t = tn2.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            //double a0 = a(isn,tn);
            //double a1 = a(isn,tn1);
            //double a2 = a(isn,tn2);

            double alpha = -a(isn,tn)*lambda*h;
            double betta = 1.0 - 2.0*alpha;

            double alpha1 = (1.0 - 2.0*lambda)*a(isn,tn1)*h;
            double alpha2 = lambda*a(isn,tn2)*h;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = (alpha1*(u[m][n+1]   - 2.0*u[m][n]   + u[m][n-1]))   + 2.0*u[m][n]
                    + (alpha2*(u[m-1][n+1] - 2.0*u[m-1][n] + u[m-1][n-1])) - u[m-1][n]
                    + (ht*ht)*f(isn, tn1);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m+1][0] = boundary(lsn, tn, Left);
        u[m+1][N] = boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * lambda * h * u[m+1][0];
        kd[N-2] += a(rsn,tn) * lambda * h * u[m+1][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m+1][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void HyperbolicIBVP::gridMethod1(DoubleMatrix &u, SweepMethodDirection direction)
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

    double lambda = 0.25;
    double h = (ht*ht)/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

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

        u[0][n] = initial1(isn);
        u[1][n] = u[0][n] + ht*initial2(isn);
    }

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    TimeNode tn1;
    TimeNode tn2;
    for (unsigned int m=2; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        tn1.i = m+minM-1;
        tn1.t = tn1.i*ht;

        tn2.i = m+minM-2;
        tn2.t = tn2.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            double a0 = a(isn,tn);
            double a1 = a(isn,tn1);
            double a2 = a(isn,tn2);

            double alpha = -a0*lambda*h;
            double betta = 1.0 - 2.0*alpha;

            double alpha1 = (1.0 - 2.0*lambda)*a1*h;
            double alpha2 = lambda*a2*h;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = (alpha1*(u[m-1][n+1] - 2.0*u[m-1][n] + u[m-1][n-1])) + 2.0*u[m-1][n]
                    + (alpha2*(u[m-2][n+1] - 2.0*u[m-2][n] + u[m-2][n-1])) - u[m-2][n]
                    + (ht*ht)*f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * lambda * h * u[m][0];
        kd[N-2] += a(rsn,tn) * lambda * h * u[m][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void HyperbolicIBVP::gridMethod2(DoubleMatrix &u, SweepMethodDirection direction)
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

    double lambda = 0.25;
    double h = (ht*ht)/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

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

        u[0][n] = initial1(isn);
        u[1][n] = u[0][n] + ht*initial2(isn);
    }

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    TimeNode tn1;
    TimeNode tn2;
    for (unsigned int m=2; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        tn1.i = m+minM-1;
        tn1.t = tn1.i*ht;

        tn2.i = m+minM-2;
        tn2.t = tn2.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            double a0 = a(isn,tn);
            double a1 = a(isn,tn1);
            double a2 = a(isn,tn2);

            double alpha = -a0*(1.0-2.0*lambda)*h;
            double betta = 1.0 - 2.0*alpha;

            double alpha1 = lambda*a1*h;
            double alpha2 = lambda*a2*h;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = (alpha1*(u[m-1][n+1] - 2.0*u[m-1][n] + u[m-1][n-1])) + 2.0*u[m-1][n]
                    + (alpha2*(u[m-2][n+1] - 2.0*u[m-2][n] + u[m-2][n-1])) - u[m-2][n]
                    + (ht*ht)*f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * (1.0-2.0*lambda) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * (1.0-2.0*lambda) * h * u[m][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}
