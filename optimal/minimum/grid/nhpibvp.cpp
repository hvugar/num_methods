#include "nhpibvp.h"

auto INewtonHeatEquation::layerInfo(DoubleVector &, unsigned int) const -> void {}
auto INewtonHeatEquation::layerInfo(DoubleMatrix &, unsigned int) const -> void {}

/**
 * @brief Shebeke usulu
 * @param u
 * @param direction
 */
void NewtonHeatEquation::calculateGM1(DoubleVector &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN-minN;

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N+1));
    double *kb = (double*) malloc(sizeof(double)*(N+1));
    double *kc = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn);
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
        aa = a(isn,tn);

        ka[0] = 0.0;
        kb[0] = 1.0+2.0*aa*(ht/(hx*hx))+2.0*lambda1*aa*(ht/hx) + lambda0*ht;
        kc[0] = -2.0*aa*(ht/(hx*hx));
        kd[0] = u[0] + lambda0*ht*theta0(tn) + 2.0*lambda1*aa*(ht/hx)*theta1(tn) + ht*f(isn,tn);
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

        ka[N] = -2.0*aa*(ht/(hx*hx));
        kb[N] = 1.0+2.0*aa*(ht/(hx*hx))+2.0*lambda2*aa*(ht/hx) + lambda0*ht;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0*ht*theta0(tn) + 2.0*lambda2*aa*(ht/hx)*theta2(tn) + ht*f(isn,tn);

        (*algorithm)(ka, kb, kc, kd, rx, N+1);

        for (unsigned int n=0; n<=N; n++) u[n] = rx[n];

        //layerInfo(u, m);
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

/**
 * @brief Soldan saga serti qovmaq
 * @param u
 * @param direction
 */
void NewtonHeatEquation::calculateGM2(DoubleVector &u, SweepMethodDirection direction)
{
    C_UNUSED(direction);
    //typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    //t_algorithm algorithm = &tomasAlgorithm;
    //if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    //if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN-minN;

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N+1));
    double *kb = (double*) malloc(sizeof(double)*(N+1));
    double *kc = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn);
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
        aa = a(isn,tn);

        ka[0] = 0.0;
        kb[0] = 1.0+2.0*aa*(ht/(hx*hx))+2.0*lambda1*aa*(ht/hx) + lambda0*ht;
        kc[0] = -2.0*aa*(ht/(hx*hx));
        kd[0] = u[0] + lambda0*ht*theta0(tn) + 2.0*lambda1*aa*(ht/hx)*theta1(tn) + ht*f(isn,tn);
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

        ka[N] = -2.0*aa*(ht/(hx*hx));
        kb[N] = 1.0+2.0*aa*(ht/(hx*hx))+2.0*lambda2*aa*(ht/hx) + lambda0*ht;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0*ht*theta0(tn) + 2.0* lambda2*aa*(ht/hx)*theta2(tn) + ht*f(isn,tn);

        //(*algorithm)(ka, kb, kc, kd, rx, N+1);

        double *betta = (double *)malloc(sizeof(double)*(N+1));
        for (unsigned int i=0; i<=N; i++) betta[i] = 0.0;
        double eta = kd[0];
        betta[0] = kb[0];
        betta[1] = kc[0];

        //IPrinter::printVector(betta,N+1);

        double betta_norm = 0.0;
        for (unsigned int i=0; i<=N; i++)
        {
            betta_norm += betta[i]*betta[i];
        }
        betta_norm += eta*eta;
        betta_norm = sqrt(betta_norm);
        //printf("norm: %f\n", betta_norm);

        for (unsigned int i=0; i<=N; i++)
        {
            betta[i] = betta[i]/betta_norm;
        }
        eta /= betta_norm;
        //.0IPrinter::printVector(betta,N+1);

        betta_norm = 0.0;
        for (unsigned int i=0; i<=N; i++)
        {
            betta_norm += betta[i]*betta[i];
        }
        betta_norm = sqrt(betta_norm);
        //printf("norm: %f\n", betta_norm);

        for (unsigned int n=1; n<=N-1; n++)
        {
            betta[n+0] = -betta[n-1]*(kb[n]/ka[n]) + betta[n];
            betta[n+1] = -betta[n-1]*(kc[n]/ka[n]);
            eta = eta - betta[n-1]*(kd[n]/ka[n]);

            double betta_norm = 0.0;
            for (unsigned int i=0; i<=N; i++)
            {
                betta_norm += betta[i]*betta[i];
            }
            betta_norm += eta*eta;
            betta_norm = sqrt(betta_norm);
            for (unsigned int i=0; i<=N; i++)
            {
                betta[i] = betta[i]/betta_norm;
            }
            eta /= betta_norm;

        }

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = betta[N-1]; M[0][1] = betta[N]; A[0] = eta;
        M[1][0] = ka[N];      M[1][1] = kb[N];    A[1] = kd[N];

        LinearEquation::GaussianElimination(M,A,x);

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
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

/**
 * @brief Kohne shebeke usulu
 * @param u
 * @param direction
 */
void NewtonHeatEquation::calculateGM3(DoubleVector &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN-minN;

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N-1));
    double *kb = (double*) malloc(sizeof(double)*(N-1));
    double *kc = (double*) malloc(sizeof(double)*(N-1));
    double *kd = (double*) malloc(sizeof(double)*(N-1));
    double *rx = (double*) malloc(sizeof(double)*(N-1));

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn);
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
        aa = a(isn,tn);
        double a0 = +1.0+lambda1*aa*(ht/hx) + lambda0*ht;
        double b0 = +aa*(ht/(hx*hx));
        double c0 = -aa*(ht/(hx*hx));
        double d0 = u[0] + lambda0*ht*theta0(tn) + lambda1*aa*(ht/hx)*theta1(tn) + ht*f(isn,tn);

        isn.i = 1+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);
        double a1 = -aa*(ht/(hx*hx));
        double b1 = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
        double c1 = -aa*(ht/(hx*hx));
        double d1 = u[1] + lambda0*ht*theta0(tn) + ht*f(isn,tn);

        ka[0] = 0.0;
        kb[0] = b1-b0*(a1/a0);
        kc[0] = c1-c0*(a1/a0);
        kd[0] = d1-d0*(a1/a0);

        for (unsigned int n=2; n<=N-2; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            aa = a(isn,tn);

            ka[n-1] = -aa*(ht/(hx*hx));
            kb[n-1] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
            kc[n-1] = -aa*(ht/(hx*hx));
            kd[n-1] = u[n] + lambda0*ht*theta0(tn) + ht*f(isn,tn);
        }

        isn.i = (N-1)+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);
        double aN1 = -aa*(ht/(hx*hx));
        double bN1 = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
        double cN1 = -aa*(ht/(hx*hx));
        double dN1 = u[N-1] + lambda0*ht*theta0(tn) + ht*f(isn,tn);

        isn.i = N+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);
        double aN = -aa*(ht/(hx*hx));
        double bN = +aa*(ht/(hx*hx));
        double cN = +1.0+lambda2*aa*(ht/hx) + lambda0*ht;
        double dN = u[N] + lambda0*ht*theta0(tn) + lambda2*aa*(ht/hx)*theta2(tn) + ht*f(isn,tn);

        ka[N-2] = aN1-aN*(cN1/cN);
        kb[N-2] = bN1-bN*(cN1/cN);
        kc[N-2] = 0.0;
        kd[N-2] = dN1-dN*(cN1/cN);

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++)
            u[n] = rx[n-1];
        u[0] = -(b0/a0)*u[1]  -(c0/a0)*u[2]  +d0/a0;
        u[N] = -(bN/cN)*u[N-1]-(aN/cN)*u[N-2]+dN/cN;

        //layerInfo(u, m);
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}
