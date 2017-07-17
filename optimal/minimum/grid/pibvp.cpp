#include "pibvp.h"

void ParabolicIBVP::gridMethod(DoubleVector &u, SweepMethodDirection direction) const
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

    u.clear();
    u.resize(N+1);

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
        u[n] = initial(isn);
    }
    layerInfo(u, 0);

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
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
            kd[n-1] = u[n] + ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[0] = boundary(lsn, tn, Left);
        u[N] = boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * h * u[0];
        kd[N-2] += a(rsn,tn) * h * u[N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[n] = rx[n-1];

        layerInfo(u, m);
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void ParabolicIBVP::gridMethod(DoubleMatrix &u, SweepMethodDirection direction) const
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
        u[0][n] = initial(isn);
    }

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
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
            kd[n-1] = u[m-1][n] + ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * h * u[m][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];

        IPrinter::printVector(u.row(m),"m1");
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void ParabolicIBVP::gridMethod1L(DoubleMatrix &u, SweepMethodDirection direction UNUSED_PARAM) const
{
    //typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    //t_algorithm algorithm = &tomasAlgorithm;
    //if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    //if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

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
        u[0][n] = initial(isn);
    }

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
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
            kd[n-1] = u[m-1][n] + ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * h * u[m][N];

        //printf("%d\n",N);
        //printf("a: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",ka[i]); printf("\n");
        //printf("b: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",kb[i]); printf("\n");
        //printf("c: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",kc[i]); printf("\n");
        //printf("d: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",kd[i]); printf("\n");
        //IPrinter::printSeperatorLine();
        //return;

        //(*algorithm)(ka, kb, kc, kd, rx, N-1);

#define NORMALIZE_1_1

        double *betta = (double *)malloc(sizeof(double)*(N-1));
        for (unsigned int i=0; i<=N-2; i++) betta[i] = 0.0;
        double eta = kd[0];
        betta[0] = kb[0];
        betta[1] = kc[0];

        double *norma = (double *)malloc(sizeof(double)*(N-1));

        //printf("--- %d %d %16.12f %16.12f %16.12f\n", m, 0, betta[0], betta[1], eta);

        //printf("0 %2d|", 0);
        //for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
        //printf("|%.10f\n", eta);

#ifdef NORMALIZE_1
        norma[0] = eta*eta;
        for (unsigned int i=0; i<=N-2; i++) norma[0] += betta[i]*betta[i]; norma[0] = sqrt(norma[0]);
        for (unsigned int i=0; i<=N-2; i++) betta[i] /= norma[0];
        eta /= norma[0];

        printf("1 %2d|", 0);
        for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
        printf("|%.10f|%.10f\n", eta,norma[0]);
#endif
        //printf("+++ %d %d %16.12f %16.12f %16.12f\n", m, 0, betta[0], betta[1], eta);

        //printf("---------------------------------\n");
        for (unsigned int n=1; n<=N-3; n++)
        {
            betta[n+0] = betta[n+0] - betta[n-1]*(kb[n]/ka[n]);
            betta[n+1] = betta[n+1] - betta[n-1]*(kc[n]/ka[n]);
            eta        = eta        - betta[n-1]*(kd[n]/ka[n]);

            //printf("0 %2d|", n);
            //for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
            //printf("|%.10f\n", eta);
#ifdef NORMALIZE_1
            norma[n] = eta*eta;
            for (unsigned int i=n; i<=N-2; i++) norma[n] += betta[i]*betta[i]; norma[n] = sqrt(norma[n]);
            for (unsigned int i=n; i<=N-2; i++) betta[i] /= norma[n];
            eta /= norma[n];

            printf("1 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
            printf("|%.10f|%.10f\n", eta,norma[n]);
#endif
        }
        //printf("---------------------------------\n");

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = betta[N-3]; M[0][1] = betta[N-2]; A[0] = eta;
        M[1][0] = ka[N-2];    M[1][1] = kb[N-2];    A[1] = kd[N-2];

        //IPrinter::print(M,2,2);

        GaussianElimination(M,A,x);

        //printf("%d %f %f\n", m, x[0], x[1]);

        //for (unsigned int n=0; n<=N; n++) u[n] = rx[n];
        u[m][N-1] = x[1];
        u[m][N-2] = x[0];
        for (unsigned int n=N-3; n!=0; n--)
        {
            //printf("0 %2d|", n);
            //for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
            //printf("|%.10f|%.10f\n", eta,norma[n]);
#ifdef NORMALIZE_1
            for (unsigned int i=n; i<=n+1; i++) betta[i] *= norma[n];
            eta *= norma[n];

            printf("1 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
            printf("|%.10f\n", eta);
#endif

            eta        = eta        + betta[n-1]*(kd[n]/ka[n]);
            betta[n+0] = betta[n+0] + betta[n-1]*(kb[n]/ka[n]);
            betta[n+1] = betta[n+1] + betta[n-1]*(kc[n]/ka[n]);
            //u[m][n]    = (eta - betta[n+1]*u[m][n+2] - betta[n+0]*u[m][n+1])/betta[n-1];
            u[m][n]    = (eta - betta[n+0]*u[m][n+1])/betta[n-1];

//            printf("2 %2d|", n);
//            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
//            printf("|%.10f\n", eta);

#ifdef NORMALIZE_1
            //norma[n] = eta*eta;
            //for (unsigned int i=n-1; i<=n; i++) norma[n] += betta[i]*betta[i]; norma[n] = sqrt(norma[n]);
#endif
            //printf("------ %d %.10f\n", n, u[m][n]);
        }

        IPrinter::printVector(16,12, u.row(m));
        //for (unsigned int n=0; n<=N-2; n++) u[m][n+1] = rx[n];
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void ParabolicIBVP::gridMethod1R(DoubleMatrix &u, SweepMethodDirection direction UNUSED_PARAM) const
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

    double h = ht/(hx*hx);

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
        u[0][n] = initial(isn);
    }

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
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
            kd[n-1] = u[m-1][n] + ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

         /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * h * u[m][N];

        //printf("%d\n",N);
        //printf("a: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",ka[i]); printf("\n");
        //printf("b: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",kb[i]); printf("\n");
        //printf("c: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",kc[i]); printf("\n");
        //printf("d: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",kd[i]); printf("\n");
        //IPrinter::printSeperatorLine();

        unsigned int N1 = N-2;
        double *betta = (double *)malloc(sizeof(double)*(N1+1));
        for (unsigned int i=0; i<=N1; i++) betta[i] = 0.0;

        betta[N1+0] = kb[N1];
        betta[N1-1] = ka[N1];
        double eta  = kd[N1];

        //printf("%18.10f %18.10f %18.10f\n", betta[N1-1], betta[N1-0], eta);

        //printf("---------------------------------\n");
        for (unsigned int n=N1; n>=2; n--)
        {
            betta[n-1] = betta[n-1] - betta[n]*(kb[n-1]/kc[n-1]);
            betta[n-2] = betta[n-2] - betta[n]*(ka[n-1]/kc[n-1]);
            eta        = eta        - betta[n]*(kd[n-1]/kc[n-1]);

            //printf("%18.10f %18.10f %18.10f\n", betta[n-2], betta[n-1], eta);
        }
        //printf("---------------------------------\n");

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = betta[0];   M[0][1] = betta[1];   A[0] = eta;
        M[1][0] = kb[0];      M[1][1] = kc[0];      A[1] = kd[0];

        //puts("---");
        //IPrinter::print(M,2,2);
        //puts("---");

        GaussianElimination(M,A,x);

        printf("%d %18.10f %18.10f\n", m, x[0], x[1]);
        //return;

        u[m][1] = x[0];
        u[m][2] = x[1];
        for (unsigned int n=2; n<=N1; n++)
        {
            //printf("0 %2d|", n);
            //for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
            //printf("|%.10f|%.10f\n", eta,norma[n]);
            eta        = eta        + betta[n]*(kd[n-1]/kc[n-1]);
            betta[n-2] = betta[n-2] + betta[n]*(ka[n-1]/kc[n-1]);
            betta[n-1] = betta[n-1] + betta[n]*(kb[n-1]/kc[n-1]);
            u[m][n+1]  = (eta - betta[n-1]*u[m][n])/betta[n];

            //printf("2 %2d|", n);
            //for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
            //printf("|%.10f\n", eta);
            //printf("------ %d %.10f\n", n, u[m][n]);
        }

        IPrinter::printVector(16,12, u.row(m));
        //for (unsigned int n=0; n<=N-2; n++) u[m][n+1] = rx[n];
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

#define real long double
void ParabolicIBVP::gridMethod11(DoubleMatrix &u, SweepMethodDirection direction UNUSED_PARAM) const
{
    //typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    //t_algorithm algorithm = &tomasAlgorithm;
    //if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    //if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    real ht = time.step();
    unsigned int minM = time.minN();
    unsigned int maxM = time.maxN();
    unsigned int M = maxM-minM;

    real hx = dim1.step();
    unsigned int minN = dim1.minN();
    unsigned int maxN = dim1.maxN();
    unsigned int N = maxN-minN;

    real h = ht/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    real *ka = (real*) malloc(sizeof(real)*(N-1));
    real *kb = (real*) malloc(sizeof(real)*(N-1));
    real *kc = (real*) malloc(sizeof(real)*(N-1));
    real *kd = (real*) malloc(sizeof(real)*(N-1));
    real *rx = (real*) malloc(sizeof(real)*(N-1));

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            real alpha = -a(isn,tn)*h;
            real betta = 1.0 - 2.0*alpha;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = u[m-1][n] + ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        IPrinter::printSeperatorLine();
        for (unsigned int i=0; i<=N-2; i++)
        {
            for (unsigned int j=0; j<=N-2; j++)
            {}
            if (i==0)
            printf("%14.10Lf %14.10Lf %14.10Lf %14.10Lf\n", ka[i], kb[i], kc[i], kd[i]);
        }
        IPrinter::printSeperatorLine();

        /* border conditions */
        u[m][0] = (real) boundary(lsn, tn, Left);
        u[m][N] = (real) boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * h * u[m][N];

        //(*algorithm)(ka, kb, kc, kd, rx, N-1);

#define NORMALIZE_11_1

        real *betta = (real *)malloc(sizeof(real)*(N-1));
        for (unsigned int i=0; i<=N-2; i++) betta[i] = 0.0;
        real eta = kd[0];
        betta[0] = kb[0];
        betta[1] = kc[0];

        real *norma = (real *)malloc(sizeof(real)*(N-1));

        //printf("--- %d %d %16.12f %16.12f %16.12f\n", m, 0, betta[0], betta[1], eta);

        printf("0 %2d|", 0);
        for (unsigned int i=0; i<=N-2; i++) printf("%16.12Lf ", betta[i]);
        printf("|%.10Lf\n", eta);

#ifdef NORMALIZE_11
        norma[0] = eta*eta;
        for (unsigned int i=0; i<=N-2; i++) norma[0] += betta[i]*betta[i]; norma[0] = sqrt(norma[0]);
        for (unsigned int i=0; i<=N-2; i++) betta[i] /= norma[0];
        eta /= norma[0];

        printf("1 %2d|", 0);
        for (unsigned int i=0; i<=N-2; i++) printf("%16.12Lf ", betta[i]);
        printf("|%.10Lf|%.10Lf\n", eta,norma[0]);
#endif
        //printf("+++ %d %d %16.12f %16.12f %16.12f\n", m, 0, betta[0], betta[1], eta);

        printf("---------------------------------\n");
        for (unsigned int n=1; n<=N-3; n++)
        {
            betta[n+0] = betta[n+0] - betta[n-1]*(kb[n]/ka[n]);
            betta[n+1] = betta[n+1] - betta[n-1]*(kc[n]/ka[n]);
            eta        = eta        - betta[n-1]*(kd[n]/ka[n]);

            printf("0 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12Lf ", betta[i]);
            printf("|%.10Lf\n", eta);
#ifdef NORMALIZE_11
            norma[n] = eta*eta;
            for (unsigned int i=n; i<=N-2; i++) norma[n] += betta[i]*betta[i]; norma[n] = sqrt(norma[n]);
            for (unsigned int i=n; i<=N-2; i++) betta[i] /= norma[n];
            eta /= norma[n];

            printf("1 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12Lf ", betta[i]);
            printf("|%.10Lf|%.10Lf\n", eta,norma[n]);
#endif
        }
        printf("---------------------------------\n");

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = betta[N-3]; M[0][1] = betta[N-2]; A[0] = eta;
        M[1][0] = ka[N-2];    M[1][1] = kb[N-2];    A[1] = kd[N-2];

        //IPrinter::print(M,2,2);

        GaussianElimination(M,A,x);

        //printf("%d %f %f\n", m, x[0], x[1]);

        //for (unsigned int n=0; n<=N; n++) u[n] = rx[n];
        u[m][N-1] = x[1];
        u[m][N-2] = x[0];
        for (unsigned int n=N-3; n!=0; n--)
        {
            printf("0 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12Lf ", betta[i]);
            printf("|%.10Lf|%.10Lf\n", eta,norma[n]);
#ifdef NORMALIZE_11
            for (unsigned int i=n; i<=n+1; i++) betta[i] *= norma[n];
            eta *= norma[n];

            printf("1 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12Lf ", betta[i]);
            printf("|%.10Lf\n", eta);
#endif

            eta        = eta        + betta[n-1]*(kd[n]/ka[n]);
            betta[n+0] = betta[n+0] + betta[n-1]*(kb[n]/ka[n]);
            betta[n+1] = betta[n+1] + betta[n-1]*(kc[n]/ka[n]);
            //u[m][n]    = (eta - betta[n+1]*u[m][n+2] - betta[n+0]*u[m][n+1])/betta[n-1];
            u[m][n]    = (eta - betta[n+0]*u[m][n+1])/betta[n-1];

//            printf("2 %2d|", n);
//            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
//            printf("|%.10f\n", eta);

#ifdef NORMALIZE_11
            //norma[n] = eta*eta;
            //for (unsigned int i=n-1; i<=n; i++) norma[n] += betta[i]*betta[i]; norma[n] = sqrt(norma[n]);
#endif
            //printf("------ %d %.10f\n", n, u[m][n]);
        }

        IPrinter::printVector(16,12, u.row(m), "\n   m1");
        //for (unsigned int n=0; n<=N-2; n++) u[m][n+1] = rx[n];
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void ParabolicIBVP::gridMethod2(DoubleMatrix &u, SweepMethodDirection direction UNUSED_PARAM) const
{
    //typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    //t_algorithm algorithm = &tomasAlgorithm;
    //if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    //if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

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
        u[0][n] = initial(isn);
    }

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
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
            kd[n-1] = u[m-1][n] + ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * h * u[m][N];

        //(*algorithm)(ka, kb, kc, kd, rx, N-1);

#define _NORMALIZE

        DoubleMatrix betta;
        betta.resize(N-1,N-1,0.0);
        DoubleVector eta(N-1,0.0);

//        double *betta = (double *)malloc(sizeof(double)*(N-1));
//        for (unsigned int i=0; i<=N-2; i++) betta[i] = 0.0;
        eta.at(0)     = kd[0];
        betta.at(0,0) = kb[0];
        betta.at(0,1) = kc[0];

        double *norma = (double *)malloc(sizeof(double)*(N-1));

        //printf("--- %d %d %16.12f %16.12f %16.12f\n", m, 0, betta[0], betta[1], eta);

        printf("0 %2d|", 0);
        for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta.at(0,i));
        printf("|%.10f\n", eta.at(0));

#ifdef _NORMALIZE
        norma[0] = eta.at(0)*eta.at(0);
        for (unsigned int i=0; i<=N-2; i++) norma[0] += betta.at(0,i)*betta.at(0,i); norma[0] = sqrt(norma[0]);
        for (unsigned int i=0; i<=N-2; i++) betta.at(0,i) /= norma[0];
        eta.at(0) /= norma[0];

        printf("1 %2d|", 0);
        for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta.at(0,i));
        printf("|%.10f|%.10f\n", eta.at(0), norma[0]);
#endif
        //printf("+++ %d %d %16.12f %16.12f %16.12f\n", m, 0, betta[0], betta[1], eta);

        printf("---------------------------------\n");
        for (unsigned int n=1; n<=N-3; n++)
        {
            betta.at(n,n+0) = betta.at(n-1,n+0) - betta.at(n-1,n-1)*(kb[n]/ka[n]);
            betta.at(n,n+1) = betta.at(n-1,n+1) - betta.at(n-1,n-1)*(kc[n]/ka[n]);
            eta.at(n)       = eta.at(n-1)       - betta.at(n-1,n-1)*(kd[n]/ka[n]);

            printf("0 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta.at(n,i));
            printf("|%.10f\n", eta.at(n));
#ifdef _NORMALIZE
            norma[n] = eta.at(n)*eta.at(n);
            for (unsigned int i=n; i<=N-2; i++) norma[n] += betta.at(n,i)*betta.at(n,i); norma[n] = sqrt(norma[n]);
            for (unsigned int i=n; i<=N-2; i++) betta.at(n,i) /= norma[n];
            eta.at(n) /= norma[n];

            printf("1 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta.at(n,i));
            printf("|%.10f|%.10f\n", eta.at(n), norma[n]);
#endif
        }
        printf("---------------------------------\n");

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = betta.at(N-3,N-3); M[0][1] = betta.at(N-3,N-2); A[0] = eta.at(N-3);
        M[1][0] = ka[N-2];           M[1][1] = kb[N-2];           A[1] = kd[N-2];

        //IPrinter::print(M,2,2);

        GaussianElimination(M,A,x);

        printf("%d %f %f\n", m, x[0], x[1]);

        //for (unsigned int n=0; n<=N; n++) u[n] = rx[n];
        u[m][N-1] = x[1];
        u[m][N-2] = x[0];
        for (unsigned int n=N-3; n!=0; n--)
        {
#ifdef _NORMALIZE
            printf("0 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta.at(n,i));
            printf("|%.10f|%.10f\n", eta.at(n), norma[n]);

            for (unsigned int i=n; i<=n+1; i++) betta.at(n,i) *= norma[n];
            eta.at(n) *= norma[n];

            printf("1 %2d|", n);
            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta.at(n,i));
            printf("|%.10f\n", eta.at(n));
#endif

            eta.at(n-1)       = eta.at(n)       + betta.at(n-1,n-1)*(kd[n]/ka[n]);
            betta.at(n-1,n+0) = betta.at(n,n+0) + betta.at(n-1,n-1)*(kb[n]/ka[n]);
            betta.at(n-1,n+1) = betta.at(n,n+1) + betta.at(n-1,n-1)*(kc[n]/ka[n]);
            //u[m][n]    = (eta - betta[n+1]*u[m][n+2] - betta[n+0]*u[m][n+1])/betta[n-1];
            u[m][n]    = (eta.at(n-1) - betta.at(n-1,n+0)*u[m][n+1])/betta.at(n-1,n-1);

//            printf("2 %2d|", n);
//            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
//            printf("|%.10f\n", eta);

#ifdef _NORMALIZE
            //norma[n] = eta*eta;
            //for (unsigned int i=n-1; i<=n; i++) norma[n] += betta[i]*betta[i]; norma[n] = sqrt(norma[n]);
#endif
            //printf("------ %d %.10f\n", n, u[m][n]);
        }

        IPrinter::printVector(u.row(m), "m1");
        //for (unsigned int n=0; n<=N-2; n++) u[m][n+1] = rx[n];
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void ParabolicIBVP::calculateN2L2RD(DoubleMatrix &u) const
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.minN();
    unsigned int maxM = time.maxN();
    unsigned int M = maxM - minM;

    double hx = dim1.step();
    unsigned int minN = dim1.minN();
    unsigned int maxN = dim1.maxN();
    unsigned int N = maxN - minN;

    const unsigned int k = 2;
    double h = ht/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    TimeNode tn;
    SpaceNode lsn;
    SpaceNode rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        /* n=1 */
        isn.i = minN+1;
        isn.x = isn.i*hx;

        double alpha = a(isn,tn)*h;
        A[0][0] = -2.0*alpha - 1.0;
        A[0][1] = alpha;
        b[0]    = -u[m-1][1] - alpha*u[m][0] - ht*f(isn,tn);

        A[0][1] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = b[0];

        for (unsigned int n=2; n<=N-k; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = a(isn,tn)*h;

            double g1 = alpha;
            double g2 = -2.0*alpha-1.0;
            double g3 = alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);

            g2 /= -g1;
            g3 /= -g1;
            fi /= +g1;
            g1  = 1.0;

            A[0][0] = A[0][1] + g2;
            A[0][1] = g3;
            b[0]    = b[0] - fi;
            \
            A[0][1] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems[n-1][0] = A[0][1];
            ems[n-1][1] = b[0];
        }

        isn.i = maxN-1;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[1][0] = alpha;
        A[1][1] = -2.0*alpha - 1.0;
        b[1]    = -u[m-1][N-1] - alpha*u[m][N] - ht*f(isn,tn);

        GaussianElimination(A, b, x);

        u[m][N-1] = x[1];
        u[m][N-2] = x[0];
        for (unsigned int n=N-(k+1); n>=1; n--)
        {
            u[m][n] = -ems[n-1][0]*u[m][n+1]+ems[n-1][1];
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void ParabolicIBVP::calculateN4L2RD(DoubleMatrix &u) const
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.minN();
    unsigned int maxM = time.maxN();
    unsigned int M = maxM - minM;

    double hx = dim1.step();
    unsigned int minN = dim1.minN();
    unsigned int maxN = dim1.maxN();
    unsigned int N = maxN - minN;

    const unsigned int k = 4;
    double h = ht/(24.0*hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    TimeNode tn;
    SpaceNode lsn;
    SpaceNode rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        /* n=1 */
        //        isn.i = minN+1;
        //        isn.x = isn.i*hx;
        //        double alpha = a(isn,tn)*h;
        //        A[0][0] = -40.0*alpha - 1.0;
        //        A[0][1] = +12.0*alpha;
        //        A[0][2] = +8.0*alpha;
        //        A[0][3] = -2.0*alpha;
        //        b[0]    = -u[m-1][1] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        /* n=4 */
        isn.i = minN+4;
        isn.x = isn.i*hx;
        double alpha = a(isn,tn)*h;
        A[0][0] = -112.0*alpha;
        A[0][1] = +228.0*alpha;
        A[0][2] = -208.0*alpha;
        A[0][3] = +70.0*alpha - 1.0;
        b[0]    = -u[m-1][4] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = A[0][2];
        ems[0][2] = A[0][3];
        ems[0][3] = b[0];

        for (unsigned int n=0; n<=N-(k+1)-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = a(isn,tn)*h;

            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            fi /= +g1;
            g1 = 1.0;

            A[0][0] = A[0][1] + g2;
            A[0][1] = A[0][2] + g3;
            A[0][2] = A[0][3] + g4;
            A[0][3] = g5;
            b[0]    = b[0] - fi;
            \
            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems[n][0] = A[0][1];
            ems[n][1] = A[0][2];
            ems[n][2] = A[0][3];
            ems[n][3] = b[0];
        }

        isn.i = maxN-3;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u[m-1][N-3] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        isn.i = maxN-2;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u[m-1][N-2] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        isn.i = maxN-1;
        isn.x = isn.i*hx;
        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u[m-1][N-1] - (+22.0*alpha)*u[m][N] - ht*f(isn,tn);

        GaussianElimination(A, b, x);

        u[m][N-1] = x[3];
        u[m][N-2] = x[2];
        u[m][N-3] = x[1];
        u[m][N-4] = x[0];
        for (unsigned int n=N-(k+1); n>=1; n--)
        {
            u[m][n] = -ems[n-1][0]*u[m][n+1]
                    -ems[n-1][1]*u[m][n+2]
                    -ems[n-1][2]*u[m][n+3]
                    +ems[n-1][3];
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void ParabolicIBVP::calculateN4L2RD_1(DoubleMatrix &u) const
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.minN();
    unsigned int maxM = time.maxN();
    unsigned int M = maxM - minM;

    double hx = dim1.step();
    unsigned int minN = dim1.minN();
    unsigned int maxN = dim1.maxN();
    unsigned int N = maxN - minN;

    const unsigned int k = 4;
    double h = ht/(24.0*hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    TimeNode tn;
    SpaceNode lsn;
    SpaceNode rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        /* n=1 */
        isn.i = minN+1;
        isn.x = isn.i*hx;
        double alpha = a(isn,tn)*h;
        A[0][0] = +40.0*alpha + 1.0;
        A[0][1] = -12.0*alpha;
        A[0][2] = -8.0*alpha;
        A[0][3] = +2.0*alpha;
        b[0]    = u[m-1][1] + 22.0*alpha*u[m][0] + ht*f(isn,tn);

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = A[0][2];
        ems[0][2] = A[0][3];
        ems[0][3] = b[0];

        for (unsigned int n=0; n<=N-(k+1)-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = a(isn,tn)*h;

            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            fi /= +g1;
            g1 = 1.0;

            A[0][0] = A[0][1] + g2;
            A[0][1] = A[0][2] + g3;
            A[0][2] = A[0][3] + g4;
            A[0][3] = g5;
            b[0]    = b[0] - fi;
            \
            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems[n][0] = A[0][1];
            ems[n][1] = A[0][2];
            ems[n][2] = A[0][3];
            ems[n][3] = b[0];
        }

        isn.i = maxN-3;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u[m-1][N-3] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        isn.i = maxN-2;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u[m-1][N-2] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        isn.i = maxN-1;
        isn.x = isn.i*hx;
        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u[m-1][N-1] - (+22.0*alpha)*u[m][N] - ht*f(isn,tn);

        GaussianElimination(A, b, x);

        u[m][N-1] = x[3];
        u[m][N-2] = x[2];
        u[m][N-3] = x[1];
        u[m][N-4] = x[0];
        for (unsigned int n=N-(k+1); n>=1; n--)
        {
            u[m][n] = -ems[n-1][0]*u[m][n+1]
                    -ems[n-1][1]*u[m][n+2]
                    -ems[n-1][2]*u[m][n+3]
                    +ems[n-1][3];
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void ParabolicIBVP::calculateMVD(DoubleMatrix &u) const
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);
    Dimension dim2 = mspaceDimension.at(1);

    double ht = time.step();
    double h1 = dim1.step();
    double h2 = dim2.step();
    unsigned int M = time.maxN();
    unsigned int N1 = dim1.maxN();
    unsigned int N2 = dim2.maxN();

    double a1 = 1.0;
    double a2 = 1.0;

    //cleaning matrix
    u.resize(N2+1, N1+1);

    DoubleMatrix uh(N2+1, N1+1);

    DoubleVector da1(N1-1);
    DoubleVector db1(N1-1);
    DoubleVector dc1(N1-1);
    DoubleVector dd1(N1-1);
    DoubleVector rx1(N1-1);

    DoubleVector da2(N2-1);
    DoubleVector db2(N2-1);
    DoubleVector dc2(N2-1);
    DoubleVector dd2(N2-1);
    DoubleVector rx2(N2-1);

    double x1_a = -(a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*ht)/(h1*h1);
    double x1_c = (a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*ht)/(h2*h2);
    double x2_c = (a1*ht)/(2.0*h1*h1);

    for (unsigned int j=0; j<=N2; j++)
    {
        SpaceNode sn;
        sn.j = j;
        sn.y = j*h2;
        for (unsigned int i=0; i<=N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            u[j][i] = initial(sn);
        }
    }

    TimeNode tn;
    SpaceNode sn;
    for (unsigned int k=1; k<=M; k++)
    {
        tn.i = 2*k-1;
        tn.t = 0.5*(2*k-1)*ht;
        // Approximation to x1 direction
        for (unsigned int j=1; j<N2; j++)
        {
            sn.j = j;
            sn.y = j*h2;
            for (unsigned int i=1; i<N1; i++)
            {
                sn.i = i;
                sn.x = i*h1;

                da1[i-1] = x1_a;
                db1[i-1] = x1_b;
                dc1[i-1] = x1_a;
                dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + (ht/2.0) * f(sn, tn);
            }

            da1[0]     = 0.0;
            dc1[N1-2]  = 0.0;

            SpaceNode sn0;
            sn0.i = 0;
            sn0.x = 0.0;
            sn0.j = j;
            sn0.y = j*h2;
            uh[j][0]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = N1;
            snN.x = N1*h1;
            snN.j = j;
            snN.y = j*h2;
            uh[j][N1] = boundary(snN, tn);

            dd1[0]    -= x1_a * uh[j][0];
            dd1[N1-2] -= x1_a * uh[j][N1];

            tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

            for (unsigned int i=1; i<N1; i++)
            {
                uh[j][i] = rx1[i-1];
            }
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            SpaceNode sn0;
            sn0.i = i;
            sn0.x = i*h1;
            sn0.j = 0;
            sn0.y = 0.0;
            uh[0][i]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = i;
            snN.x = i*h1;
            snN.j = N2;
            snN.y = N2*h2;
            uh[N2][i] = boundary(snN, tn);
        }

        tn.i = k;
        tn.t = k*ht;
        // Approximation to x2 direction
        for (unsigned int i=1; i<N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            for (unsigned int j=1; j<N2; j++)
            {
                sn.j = j;
                sn.y = j*h2;
                da2[j-1] = x2_a;
                db2[j-1] = x2_b;
                dc2[j-1] = x2_a;
                dd2[j-1] = x2_c*(uh[j][i-1] - 2.0*uh[j][i] + uh[j][i+1]) + uh[j][i] + (ht/2.0) * f(sn, tn);
            }
            da2[0]     = 0.0;
            dc2[N2-2]  = 0.0;

            SpaceNode sn0;
            sn0.i = i;
            sn0.x = i*h1;
            sn0.j = 0;
            sn0.y = 0.0;
            u[0][i]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = i;
            snN.x = i*h1;
            snN.j = N2;
            snN.y = N2*h2;
            u[N2][i] = boundary(snN, tn);

            dd2[0]    -= x2_a * u[0][i];
            dd2[N2-2] -= x2_a * u[N2][i];

            tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

            for (unsigned int j=1; j<N2; j++)
            {
                u[j][i] = rx2[j-1];
            }
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            SpaceNode sn0;
            sn0.i = 0;
            sn0.x = 0*h1;
            sn0.j = j;
            sn0.y = j*h2;
            u[j][0]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = N1;
            snN.x = N1*h1;
            snN.j = j;
            snN.y = j*h2;
            u[j][N1] = boundary(snN, tn);
        }
    }

    da1.clear();
    db1.clear();
    dc1.clear();
    dd1.clear();
    rx1.clear();

    da2.clear();
    db2.clear();
    dc2.clear();
    dd2.clear();
    rx2.clear();
}

void ParabolicIBVP::calculateMVD1(DoubleMatrix &u) const
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);
    Dimension dim2 = mspaceDimension.at(1);

    double ht = time.step();
    double h1 = dim1.step();
    double h2 = dim2.step();
    unsigned int M = time.maxN();
    unsigned int N1 = dim1.maxN();
    unsigned int N2 = dim2.maxN();

    double a1 = 1.0;
    double a2 = 1.0;

    //cleaning matrix
    u.resize(N2+1, N1+1);

    DoubleMatrix uh(N2+1, N1+1);

    double *ka1 = (double*) malloc(sizeof(double)*N1-1);
    double *kb1 = (double*) malloc(sizeof(double)*N1-1);
    double *kc1 = (double*) malloc(sizeof(double)*N1-1);
    double *kd1 = (double*) malloc(sizeof(double)*N1-1);
    double *rx1 = (double*) malloc(sizeof(double)*N1-1);

    double *ka2 = (double*) malloc(sizeof(double)*N2-1);
    double *kb2 = (double*) malloc(sizeof(double)*N2-1);
    double *kc2 = (double*) malloc(sizeof(double)*N2-1);
    double *kd2 = (double*) malloc(sizeof(double)*N2-1);
    double *rx2 = (double*) malloc(sizeof(double)*N2-1);

    double x1_a = -(a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*ht)/(h1*h1);
    double x1_c = (a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*ht)/(h2*h2);
    double x2_c = (a1*ht)/(2.0*h1*h1);

    for (unsigned int j=0; j<=N2; j++)
    {
        SpaceNode sn;
        sn.j = j;
        sn.y = j*h2;
        for (unsigned int i=0; i<=N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            u[j][i] = initial(sn);
        }
    }

    TimeNode tn;
    SpaceNode sn;
    for (unsigned int k=1; k<=M; k++)
    {
        tn.i = 2*k-1;
        tn.t = 0.5*(2*k-1)*ht;
        // Approximation to x1 direction
        for (unsigned int j=1; j<N2; j++)
        {
            sn.j = j;
            sn.y = j*h2;
            for (unsigned int i=1; i<N1; i++)
            {
                sn.i = i;
                sn.x = i*h1;

                ka1[i-1] = x1_a;
                kb1[i-1] = x1_b;
                kc1[i-1] = x1_a;
                kd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + (ht/2.0) * f(sn, tn);
            }

            ka1[0]     = 0.0;
            kc1[N1-2]  = 0.0;

            SpaceNode sn0;
            sn0.i = 0;
            sn0.x = 0.0;
            sn0.j = j;
            sn0.y = j*h2;
            uh[j][0]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = N1;
            snN.x = N1*h1;
            snN.j = j;
            snN.y = j*h2;
            uh[j][N1] = boundary(snN, tn);

            kd1[0]    -= x1_a * uh[j][0];
            kd1[N1-2] -= x1_a * uh[j][N1];

            tomasAlgorithm(ka1, kb1, kc1, kd1, rx1, N1-1);

            for (unsigned int i=1; i<N1; i++)
            {
                uh[j][i] = rx1[i-1];
            }
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            SpaceNode sn0;
            sn0.i = i;
            sn0.x = i*h1;
            sn0.j = 0;
            sn0.y = 0.0;
            uh[0][i]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = i;
            snN.x = i*h1;
            snN.j = N2;
            snN.y = N2*h2;
            uh[N2][i] = boundary(snN, tn);
        }

        tn.i = k;
        tn.t = k*ht;
        // Approximation to x2 direction
        for (unsigned int i=1; i<N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            for (unsigned int j=1; j<N2; j++)
            {
                sn.j = j;
                sn.y = j*h2;
                ka2[j-1] = x2_a;
                kb2[j-1] = x2_b;
                kc2[j-1] = x2_a;
                kd2[j-1] = x2_c*(uh[j][i-1] - 2.0*uh[j][i] + uh[j][i+1]) + uh[j][i] + (ht/2.0) * f(sn, tn);
            }
            ka2[0]     = 0.0;
            kc2[N2-2]  = 0.0;

            SpaceNode sn0;
            sn0.i = i;
            sn0.x = i*h1;
            sn0.j = 0;
            sn0.y = 0.0;
            u[0][i]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = i;
            snN.x = i*h1;
            snN.j = N2;
            snN.y = N2*h2;
            u[N2][i] = boundary(snN, tn);

            kd2[0]    -= x2_a * u[0][i];
            kd2[N2-2] -= x2_a * u[N2][i];

            tomasAlgorithm(ka2, kb2, kc2, kd2, rx2, N2-1);

            for (unsigned int j=1; j<N2; j++)
            {
                u[j][i] = rx2[j-1];
            }
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            SpaceNode sn0;
            sn0.i = 0;
            sn0.x = 0*h1;
            sn0.j = j;
            sn0.y = j*h2;
            u[j][0]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = N1;
            snN.x = N1*h1;
            snN.j = j;
            snN.y = j*h2;
            u[j][N1] = boundary(snN, tn);
        }
    }

    free(ka1);
    free(kb1);
    free(kc1);
    free(kd1);
    free(rx1);

    free(ka2);
    free(kb2);
    free(kc2);
    free(kd2);
    free(rx2);

    //    ka1.clear();
    //    kb1.clear();
    //    kc1.clear();
    //    kd1.clear();
    //    rx1.clear();

    //    ka2.clear();
    //    kb2.clear();
    //    kc2.clear();
    //    kd2.clear();
    //    rx2.clear();
}

void ParabolicIBVP::calculateMVD2(DoubleMatrix &u) const
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);
    Dimension dim2 = mspaceDimension.at(1);

    double ht = time.step();
    double h1 = dim1.step();
    double h2 = dim2.step();
    unsigned int M = time.maxN();
    unsigned int N1 = dim1.maxN();
    unsigned int N2 = dim2.maxN();

    double a1 = 1.0;
    double a2 = 1.0;

    double **v = (double**) malloc(sizeof(double*)*(N2+1));
    for (unsigned int i=0; i<=N2; i++) v[i] = (double*) malloc(sizeof(double)*(N1+1));

    double **w = (double**) malloc(sizeof(double*)*(N2+1));
    for (unsigned int i=0; i<=N2; i++) w[i] = (double*) malloc(sizeof(double)*(N1+1));

    double *ka1 = (double*) malloc(sizeof(double)*N1-1);
    double *kb1 = (double*) malloc(sizeof(double)*N1-1);
    double *kc1 = (double*) malloc(sizeof(double)*N1-1);
    double *kd1 = (double*) malloc(sizeof(double)*N1-1);
    double *rx1 = (double*) malloc(sizeof(double)*N1-1);

    double *ka2 = (double*) malloc(sizeof(double)*N2-1);
    double *kb2 = (double*) malloc(sizeof(double)*N2-1);
    double *kc2 = (double*) malloc(sizeof(double)*N2-1);
    double *kd2 = (double*) malloc(sizeof(double)*N2-1);
    double *rx2 = (double*) malloc(sizeof(double)*N2-1);

    double x1_a = -(a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*ht)/(h1*h1);
    double x1_c = (a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*ht)/(h2*h2);
    double x2_c = (a1*ht)/(2.0*h1*h1);

    for (unsigned int j=0; j<=N2; j++)
    {
        SpaceNode sn;
        sn.j = j;
        sn.y = j*h2;
        for (unsigned int i=0; i<=N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            v[j][i] = initial(sn);
        }
    }

    TimeNode tn;
    SpaceNode sn;
    TimeNode tn1;
    for (unsigned int k=1; k<=M; k++)
    {
        tn.i = 2*k-1;
        tn.t = (k-0.5)*ht;

        tn1.i = 2*k-1;
        tn1.t = (k-1.0)*ht;
        // Approximation to x1 direction
        for (unsigned int j=1; j<N2; j++)
        {
            sn.j = j;
            sn.y = j*h2;
            for (unsigned int i=1; i<N1; i++)
            {
                sn.i = i;
                sn.x = i*h1;

                ka1[i-1] = x1_a;
                kb1[i-1] = x1_b;
                kc1[i-1] = x1_a;
                kd1[i-1] = x1_c*(v[j-1][i] - 2.0*v[j][i] + v[j+1][i]) + v[j][i] + (ht/2.0) * /*f(sn, tn);*/(f(sn, tn)+f(sn, tn1))/2.0;
            }

            ka1[0]     = 0.0;
            kc1[N1-2]  = 0.0;

            SpaceNode sn0;
            sn0.i = 0;
            sn0.x = 0.0;
            sn0.j = j;
            sn0.y = j*h2;
            w[j][0]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = N1;
            snN.x = N1*h1;
            snN.j = j;
            snN.y = j*h2;
            w[j][N1] = boundary(snN, tn);

            kd1[0]    -= x1_a * w[j][0];
            kd1[N1-2] -= x1_a * w[j][N1];

            tomasAlgorithm(ka1, kb1, kc1, kd1, rx1, N1-1);

            for (unsigned int i=1; i<N1; i++)
            {
                w[j][i] = rx1[i-1];
            }
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            SpaceNode sn0;
            sn0.i = i;
            sn0.x = i*h1;
            sn0.j = 0;
            sn0.y = 0.0;
            w[0][i]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = i;
            snN.x = i*h1;
            snN.j = N2;
            snN.y = N2*h2;
            w[N2][i] = boundary(snN, tn);
        }

        tn.i = 2*k;
        tn.t = k*ht;
        // Approximation to x2 direction
        for (unsigned int i=1; i<N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            for (unsigned int j=1; j<N2; j++)
            {
                sn.j = j;
                sn.y = j*h2;
                ka2[j-1] = x2_a;
                kb2[j-1] = x2_b;
                kc2[j-1] = x2_a;
                kd2[j-1] = x2_c*(w[j][i-1] - 2.0*w[j][i] + w[j][i+1]) + w[j][i] + (ht/2.0) * f(sn, tn);
            }
            ka2[0]     = 0.0;
            kc2[N2-2]  = 0.0;

            SpaceNode sn0;
            sn0.i = i;
            sn0.x = i*h1;
            sn0.j = 0;
            sn0.y = 0.0;
            v[0][i]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = i;
            snN.x = i*h1;
            snN.j = N2;
            snN.y = N2*h2;
            v[N2][i] = boundary(snN, tn);

            kd2[0]    -= x2_a * v[0][i];
            kd2[N2-2] -= x2_a * v[N2][i];

            tomasAlgorithm(ka2, kb2, kc2, kd2, rx2, N2-1);

            for (unsigned int j=1; j<N2; j++)
            {
                v[j][i] = rx2[j-1];
            }
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            SpaceNode sn0;
            sn0.i = 0;
            sn0.x = 0*h1;
            sn0.j = j;
            sn0.y = j*h2;
            v[j][0]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = N1;
            snN.x = N1*h1;
            snN.j = j;
            snN.y = j*h2;
            v[j][N1] = boundary(snN, tn);
        }
    }

    for (unsigned int i=0; i<=N2; i++) free(w[i]);
    free(w);

    //cleaning matrix
    u.resize(N2+1, N1+1);
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            u[j][i] = v[j][i];
        }
    }

    for (unsigned int i=0; i<=N2; i++) free(v[i]);
    free(v);

    free(ka1);
    free(kb1);
    free(kc1);
    free(kd1);
    free(rx1);

    free(ka2);
    free(kb2);
    free(kc2);
    free(kd2);
    free(rx2);

    //    ka1.clear();
    //    kb1.clear();
    //    kc1.clear();
    //    kd1.clear();
    //    rx1.clear();

    //    ka2.clear();
    //    kb2.clear();
    //    kc2.clear();
    //    kd2.clear();
    //    rx2.clear();
}
