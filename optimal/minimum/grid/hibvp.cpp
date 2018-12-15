#include "hibvp.h"

IHyperbolicIBVP::~IHyperbolicIBVP() {}

CCIHyperbolicIBVP::~CCIHyperbolicIBVP() {}

void CCIHyperbolicIBVP::calculateD11(DoubleVector &u, double a, double lambda) const
{
    const Dimension &dimx = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();
    const unsigned int N = static_cast<unsigned int>(dimx.size());
    const unsigned int M = static_cast<unsigned int>(time.size());
    const double hx = dimx.step();
    const double ht = time.step();

    const double alpha = -lambda*(a*a)*((ht*ht)/(hx*hx));
    const double betta = +1.0 + 2.0*lambda*(a*a)*((ht*ht)/(hx*hx));
    const double gamma = +(1.0-2.0*lambda)*(a*a)*((ht*ht)/(hx*hx));
    const double theta = +lambda*(a*a)*((ht*ht)/(hx*hx));
    const double ht_ht = ht*ht;

    u.clear();
    u.resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

    double *da = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *db = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dc = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dd = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<N-1; n++)
    {
        da[n] = dc[n] = alpha;
        db[n] = betta;
    }
    da[0] = 0.0; dc[N-2] = 0.0;

    SpaceNodePDE sn;
    TimeNodePDE tn0, tn1, tn2;

    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = n; sn.x = n*hx;
        u0[n] = initial1(sn);
    }
    layerInfo(u0, 0);

    tn0.i = 0; tn0.t = tn0.i*ht;
    tn1.i = 1; tn1.t = tn1.i*ht;
    sn.i = 0; sn.x = 0*hx;
    u1[0] = boundary(sn, tn1);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sn.i = n; sn.x = n*hx;
        u1[n] = u0[n] + ht*initial2(sn) + 0.5*ht*ht*(a*a*((u0[n-1]-2.0*u0[n]+u0[n+1]))/(hx*hx)+f(sn, tn0));
    }
    sn.i = N; sn.x = N*hx;
    u1[N] = boundary(sn, tn1);
    layerInfo(u1, 1);

    for (unsigned int m=2; m<=M; m++)
    {
        tn2.i = m+0; tn2.t = tn2.i*ht;
        tn1.i = m-1; tn1.t = tn1.i*ht;
        tn0.i = m-2; tn0.t = tn0.i*ht;

        sn.i = 0; sn.x = 0*hx;
        u[0] = boundary(sn, tn2);
        sn.i = N; sn.x = N*hx;
        u[N] = boundary(sn, tn2);
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;
            dd[n-1] = gamma*(u1[n-1]-2.0*u1[n]+u1[n+1])
                    + theta*(u0[n-1]-2.0*u0[n]+u0[n+1])
                    + 2.0*u1[n] - u0[n];
            //dd[n-1] += ht_ht*(lambda*f(sn, tn2) + (1.0-2.0*lambda)*f(sn, tn1) + lambda*f(sn, tn0));
            dd[n-1] += ht_ht*f(sn, tn1);
        }
        dd[0]   -= alpha*u[0];
        dd[N-2] -= alpha*u[N];

        tomasAlgorithm(da, db, dc, dd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[n] = rx[n-1];

        for (unsigned int n=0; n<=N; n++)
        {
            u0[n] = u1[n];
            u1[n] = u[n];
        }

        layerInfo(u, m);
    }

    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

void CCIHyperbolicIBVP::calculateD12(DoubleVector &u, double a) const
{
    const Dimension &dimx = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();
    const unsigned int N = static_cast<const unsigned int>(dimx.size());
    const unsigned int M = static_cast<const unsigned int>(time.size());
    const double hx = dimx.step();
    const double ht = time.step();

    const double alpha = -(a*a)*((ht*ht)/(hx*hx));
    const double betta = +(1.0 + 2.0*(a*a)*((ht*ht)/(hx*hx)));
    const double ht_ht = ht*ht;

    u.clear();
    u.resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

    double *da = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *db = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dc = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dd = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<N-1; n++)
    {
        da[n] = dc[n] = alpha;
        db[n] = betta;
    }
    da[0] = 0.0; dc[N-2] = 0.0;

    SpaceNodePDE sn;
    TimeNodePDE tn0, tn1;

    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = n; sn.x = n*hx;
        u0[n] = initial1(sn);
    }
    layerInfo(u0, 0);

    tn0.i = 0; tn0.t = tn0.i*ht;
    tn1.i = 1; tn1.t = tn1.i*ht;

    sn.i = 0; sn.x = 0*hx; u1[0] = boundary(sn, tn1);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sn.i = n; sn.x = n*hx;
        u1[n] = u0[n] + ht*initial2(sn) + 0.5*ht*ht*(a*a*((u0[n-1]-2.0*u0[n]+u0[n+1]))/(hx*hx)+f(sn, tn0));
    }
    sn.i = N; sn.x = N*hx; u1[N] = boundary(sn, tn1);
    layerInfo(u0, 1);

    TimeNodePDE tn;
    for (unsigned int m=2; m<=M; m++)
    {
        tn.i = m; tn.t = tn.i*ht;

        sn.i = 0; sn.x = 0*hx; u[0] = boundary(sn, tn);
        sn.i = N; sn.x = N*hx; u[N] = boundary(sn, tn);

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;
            dd[n-1] = 2.0*u1[n] - u0[n] + ht_ht*f(sn, tn);
        }
        dd[0]   -= alpha*u[0];
        dd[N-2] -= alpha*u[N];

        tomasAlgorithm(da, db, dc, dd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[n] = rx[n-1];

        for (unsigned int n=0; n<=N; n++)
        {
            u0[n] = u1[n];
            u1[n] = u[n];
        }
        layerInfo(u, m);
    }
}

void CCIHyperbolicIBVP::calculateD21(DoubleMatrix &u, double a, double lambda) const
{}

void HyperbolicIBVP::gridMethod(DoubleMatrix &u, SweepMethodDirection direction)
{
    C_UNUSED(u);
    C_UNUSED(direction);
}

void HyperbolicIBVP::gridMethod0(DoubleMatrix &u, SweepMethodDirection direction)
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
    SpaceNodePDE isn;
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

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
    TimeNodePDE tn1;
    TimeNodePDE tn2;
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
        u[m+1][0] = boundary(lsn, tn);
        u[m+1][N] = boundary(rsn, tn);

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
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
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
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;

        u[0][n] = initial1(isn);
        u[1][n] = u[0][n] + ht*initial2(isn);
    }

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
    TimeNodePDE tn1;
    TimeNodePDE tn2;
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
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

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
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
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
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;

        u[0][n] = initial1(isn);
        u[1][n] = u[0][n] + ht*initial2(isn);
    }

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
    TimeNodePDE tn1;
    TimeNodePDE tn2;
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
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

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
