#include "pibvp.h"
#include <limits>
#include "../cmethods.h"
#include "../linearequation.h"

void CCIParabolicIBVP::calculate1(DoubleVector &u, double a) const
{
    const Dimension &dim1 = spaceDimension(Dimension::DimensionX);

    const double ht = mtimeDimension.step();
    const double hx = dim1.step();

    const int minM = mtimeDimension.min();
    const int maxM = mtimeDimension.max();
    const unsigned int difM = static_cast<const unsigned int>(maxM - minM);

    const int minN = dim1.min();
    const int maxN = dim1.max();
    const unsigned int N = static_cast<const unsigned int>(maxN - minN);

    const double h = (a*a*ht)/(hx*hx);

    u.clear();
    u.resize(N+1);

    double *ka = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *kb = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *kc = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *kd = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    /* initial condition */
    SpaceNodePDE sn;
    for (int n=minN; n<=maxM; n++)
    {
        sn.i = n; sn.x = sn.i*hx; u[n-minN] = initial(sn);
    }
    layerInfo(u, 0);

    SpaceNodePDE lsn; lsn.i = minN; lsn.x = minN*hx;
    SpaceNodePDE rsn; rsn.i = maxN; rsn.x = maxN*hx;

    TimeNodePDE tn;
    for (unsigned int m=1; m<=difM; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n+minN;
            sn.x = sn.i*hx;

            double alpha = -h;
            double betta = 1.0 - 2.0*alpha;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = u[n] + ht * f(sn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[0] = boundary(lsn, tn);
        u[N] = boundary(rsn, tn);

        kd[0]   += h * u[0];
        kd[N-2] += h * u[N];

        tomasAlgorithm(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[n] = rx[n-1];

        layerInfo(u, m);
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void CCIParabolicIBVP::calculate1(DoubleMatrix &u, double a) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionY);
    const Dimension &time = timeDimension();

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const unsigned int N = static_cast<unsigned int>(dimX.size());
    const unsigned int M = static_cast<unsigned int>(dimY.size());
    const unsigned int L = static_cast<unsigned int>(time.size());

    const double alpha = -(a*a*ht)/(hx*hx);
    const double betta = +1.0+2.0*(a*a*ht)/(hx*hx)+2.0*(a*a*ht)/(hy*hy);
    const double gamma = -(a*a*ht)/(hy*hy);

    const unsigned int N0 = N-1;
    const unsigned int M0 = M-1;

    double *da = static_cast<double*>(malloc(sizeof(double)*(M0*N0)));
    double *db = static_cast<double*>(malloc(sizeof(double)*(M0*N0)));
    double *dc = static_cast<double*>(malloc(sizeof(double)*(M0*N0)));
    double *dd = static_cast<double*>(malloc(sizeof(double)*(M0*N0)));
    double *de = static_cast<double*>(malloc(sizeof(double)*(M0*N0)));
    double *ff = static_cast<double*>(malloc(sizeof(double)*(M0*N0)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(M0*N0)));

    DoubleMatrix u0(M+1, N+1);
    u.clear();
    u.resize(M+1, N+1);
    TimeNodePDE tn;
    SpaceNodePDE sn;

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[m][n] = initial(sn);
        }
    }

    unsigned int inx = 0;
    for (unsigned int m=1; m<=M0; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=1; n<=N0; n++)
        {
            sn.i = n; sn.x = n*hx;

            da[inx] = alpha;
            db[inx] = betta;
            dc[inx] = alpha;
            dd[inx] = gamma;
            de[inx] = gamma;
            ff[inx] = ht*f(sn,tn) + u0[m][n];

            if (n==1)  { ff[inx] -= da[inx]*u[m][0]; da[inx] = 0.0; }
            if (n==N0) { ff[inx] -= dc[inx]*u[m][N]; dc[inx] = 0.0; }
            if (m==1)  { ff[inx] -= dd[inx]*u[0][n]; dd[inx] = 0.0; }
            if (m==M0) { ff[inx] -= de[inx]*u[M][n]; de[inx] = 0.0; }

            inx++;
        }
    }



    //puts("111");
    //IPrinter::printVector(da,M0*N0,nullptr,M0*N0,0,0,"data_a.txt");
    //IPrinter::printVector(db,M0*N0,nullptr,M0*N0,0,0,"data_b.txt");
    //IPrinter::printVector(dc,M0*N0,nullptr,M0*N0,0,0,"data_c.txt");
    //IPrinter::printVector(dd,M0*N0,nullptr,M0*N0,0,0,"data_d.txt");
    //IPrinter::printVector(de,M0*N0,nullptr,M0*N0,0,0,"data_e.txt");
    //IPrinter::printVector(ff,M0*N0,nullptr,M0*N0,0,0,"data_f.txt");
    //puts("222");

    free(rx);
    free(ff);
    free(de);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

void CCIParabolicIBVP::layerInfo(const DoubleMatrix &, unsigned int) const {}

void CCIParabolicIBVP::layerInfo(const DoubleVector &, unsigned int) const {}

void funcL(const double* a, const double *b, const double *c, const double *d, double *x, unsigned int N)
{
    double *e = (double*) malloc(sizeof(double)*N);
    for (unsigned int i=0; i<N; i++) e[i] = 0.0;

    double *e0 = (double*)malloc(sizeof(double)*N);
    double *e1 = (double*)malloc(sizeof(double)*N);
    double *e2 = (double*)malloc(sizeof(double)*N);
    for (unsigned int i=0; i<N; i++) e0[i] = e1[i] = e2[i] = 0.0;

    e0[0] = b[0];
    e1[0] = c[0];
    e2[0] = d[0];
    for (unsigned int n=1; n<=N-2; n++)
    {
        e0[n] = -e0[n-1]*(b[n]/a[n]) + e1[n-1];
        e1[n] = -e0[n-1]*(c[n]/a[n]);
        e2[n] = -e0[n-1]*(d[n]/a[n]) + e2[n-1];
    }

    DoubleMatrix M(2,2);
    DoubleVector A(2);
    DoubleVector y(2);

    M[0][0] = e0[N-2];   M[0][1] = e1[N-2];   A[0] = e2[N-2];
    M[1][0] = a[N-1];    M[1][1] = b[N-1];    A[1] = d[N-1];

    LinearEquation::GaussianElimination(M,A,y);

    x[N-1] = y[1];
    x[N-2] = y[0];
    for (unsigned int n=N-3; n!=UINT_MAX; n--)
    {
        x[n] = (e2[n] - e1[n]*x[n+1])/e0[n];
    }

    free(e2);
    free(e1);
    free(e0);
}

void ParabolicIBVP::gridMethod(DoubleVector &u, SweepMethodDirection direction) const
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

    double h = ht/(hx*hx);

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
    layerInfo(u, 0);

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
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
        //u[0] = boundary(lsn, tn, Left);
        //u[N] = boundary(rsn, tn, Right);

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
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
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
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
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
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

        kd[0]   += a(lsn,tn) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * h * u[m][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void ParabolicIBVP::gridMethod1L(DoubleMatrix &u, SweepMethodDirection direction UNUSED_PARAM) const
{
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

    double h = ht/(hx*hx);

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
        u[0][n] = initial(isn);
    }

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
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
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

        kd[0]   += a(lsn,tn) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * h * u[m][N];

        funcL(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void ParabolicIBVP::gridMethod1LT(DoubleMatrix &u, SweepMethodDirection direction UNUSED_PARAM) const
{
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

    double h = ht/(hx*hx);

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
        u[0][n] = initial(isn);
    }

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
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
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

        kd[0]   += a(lsn,tn) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * h * u[m][N];

        funcL(ka, kb, kc, kd, rx, N-1);
        for (unsigned int n=1; n<=N-1; n++)
        {
            u[m][n] = rx[n-1];
        }
        //        IPrinter::printVector(16,12, u.row(m));
        //        return;

        //        //printf("%d\n",N);
        //        //printf("a: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",ka[i]); printf("\n");
        //        //printf("b: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",kb[i]); printf("\n");
        //        //printf("c: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",kc[i]); printf("\n");
        //        //printf("d: "); for (unsigned int i=0;i<N-1;i++) printf("%14.10f",kd[i]); printf("\n");
        //        //IPrinter::printSeperatorLine();

        //        double *betta = (double *)malloc(sizeof(double)*(N-1));
        //        for (unsigned int i=0; i<=N-2; i++) betta[i] = 0.0;
        //        double eta = kd[0];
        //        betta[0] = kb[0];
        //        betta[1] = kc[0];

        //        for (unsigned int n=1; n<=N-3; n++)
        //        {
        //            printf("%4d %24.10f\t", n, betta[n+1]);

        //            double aaa = betta[n-1]/ka[n];
        //            betta[n+0] = betta[n+0] - aaa*kb[n];
        //            betta[n+1] = /*betta[n+1]*/ - aaa*kc[n];
        //            eta        = eta        - aaa*kd[n];

        //            printf("%4d %24.10f %24.10f\n", n, betta[n+1], aaa*kc[n]);
        //        }


        //        double *b0 = (double *)malloc(sizeof(double)*(N-1));
        //        double *b1 = (double *)malloc(sizeof(double)*(N-1));
        //        double *ee = (double *)malloc(sizeof(double)*(N-1));
        //        for (unsigned int i=0; i<=N-2; i++) b0[i] = b1[i] = ee[i] = 0.0;

        //        b0[0] = kb[0];
        //        b1[0] = kc[0];
        //        ee[0] = kd[0];

        //        for (unsigned int n=1; n<=N-3; n++)
        //        {
        //            b0[n] = -b0[n-1]*(kb[n]/ka[n]) + b1[n-1];
        //            b1[n] = -b0[n-1]*(kc[n]/ka[n]);
        //            ee[n] = -b0[n-1]*(kd[n]/ka[n]) + ee[n-1];
        //        }

        //        DoubleMatrix M(2,2);
        //        DoubleVector A(2);
        //        DoubleVector x(2);

        //        M[0][0] = b0[N-3];    M[0][1] = b1[N-3];    A[0] = ee[N-3];
        //        //M[0][0] = betta[N-3]; M[0][1] = betta[N-2]; A[0] = eta;
        //        M[1][0] = ka[N-2];    M[1][1] = kb[N-2];    A[1] = kd[N-2];

        //        //IPrinter::print(M,2,2);

        //        GaussianElimination(M,A,x);

        //        printf("%d %f %f\n", m, x[0], x[1]);
        //        return;
        //        puts("---");

        //        u[m][N-1] = x[1];
        //        u[m][N-2] = x[0];
        //        for (unsigned int n=N-3; n!=0; n--)
        //        {
        //            //printf("0 %2d|", n);
        //            //for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
        //            //printf("|%.10f|%.10f\n", eta,norma[n]);

        //            double aaa = betta[n-1]/ka[n];
        //            double ccc = aaa*kc[n];

        //            printf("%4d %24.10f\t", n, betta[n+1]);

        //            betta[n+0] = betta[n+0] + aaa*kb[n];
        //            betta[n+1] = betta[n+1] + ccc;
        //            eta        = eta        + aaa*kd[n];
        //            u[m][n]    = (eta - betta[n+0]*u[m][n+1])/betta[n-1];

        //            printf("%4d %24.10f %24.10f\n", n, betta[n+1], ccc);



        // //           printf("%18.10f %18.10f %18.10f %18.10f %18.10f %18.10f | %18.10f %18.10f %18.10f\n", eta, ee[n-1], betta[n+0], b1[n-1], betta[n-1], b0[n-1], betta[n+1], aaa*kc[n], betta[n+1] + aaa*kc[n]);

        //            //u[m][n] = (ee[n-1] - b1[n-1]*u[m][n+1])/b0[n-1];

        ////            printf("2 %2d|", n);
        ////            for (unsigned int i=0; i<=N-2; i++) printf("%16.12f ", betta[i]);
        ////            printf("|%.10f\n", eta);

        //            //printf("------ %d %.10f\n", n, u[m][n]);
        //        }

        //IPrinter::printVector(16,12, u.row(m));
        //for (unsigned int n=0; n<=N-2; n++) u[m][n+1] = rx[n];
        //break;
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
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
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
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
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
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

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

        LinearEquation::GaussianElimination(M,A,x);

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

void ParabolicIBVP::gridMethod2(DoubleMatrix &u, SweepMethodDirection direction UNUSED_PARAM) const
{
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

    double h = ht/(hx*hx);

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
        u[0][n] = initial(isn);
    }

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
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
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

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

        LinearEquation::GaussianElimination(M,A,x);

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
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM - minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
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
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    TimeNodePDE tn;
    SpaceNodePDE lsn;
    SpaceNodePDE rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

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

        LinearEquation::GaussianElimination(A, b, x);

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

#define SCHEME_2
#define __NORMALIZE__X
void ParabolicIBVP::calculateN4L2RD(DoubleMatrix &u) const
{
    /* get parameters */
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM - minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN - minN;

    const unsigned int k = 4;
    double h = ht/(24.0*hx*hx);

    /*****************************************/

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k+1);

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    TimeNodePDE tn;
    SpaceNodePDE lsn;
    SpaceNodePDE rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

        /* using 2nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = minN+1;
        isn.x = isn.i*hx;
        double alpha = a(isn,tn)*h;
        A[0][0] = -40.0*alpha - 1.0;
        A[0][1] = +12.0*alpha;
        A[0][2] = +8.0*alpha;
        A[0][3] = -2.0*alpha;
        b[0]    = -u[m-1][1] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        /* using 5th scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        //isn.i = minN+4;
        //isn.x = isn.i*hx;
        //double alpha = a(isn,tn)*h;
        //A[0][0] = -112.0*alpha;
        //A[0][1] = +228.0*alpha;
        //A[0][2] = -208.0*alpha;
        //A[0][3] = +70.0*alpha - 1.0;
        //b[0]    = -u[m-1][4] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        //double x1 = (1)*hx;
        //double x2 = (2)*hx;
        //double x3 = (3)*hx;
        //double x4 = (4)*hx;
        //printf("%4d %18.10f %18.10f\n", 0, b[0], A[0][0]*x1*x1*x1*tn.t+A[0][1]*x2*x2*x2*tn.t+A[0][2]*x3*x3*x3*tn.t+A[0][3]*x4*x4*x4*tn.t);

#ifdef __NORMALIZE__
        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;
#endif

        ems[0][0] = A[0][0];
        ems[0][1] = A[0][1];
        ems[0][2] = A[0][2];
        ems[0][3] = A[0][3];
        ems[0][4] = b[0];

#ifdef SCHEME_1
        unsigned int start = 1;
        unsigned int end   = N-k-1;
        unsigned int s     = 0;
#endif
#ifdef SCHEME_2
        unsigned int start = 2;
        unsigned int end   = N-k+0;
        unsigned int s     = 1;
#endif
#ifdef SCHEME_3
        unsigned int start = 3;
        unsigned int end   = N-k+1;
        unsigned int s     = 2;
#endif
#ifdef SCHEME_4
        unsigned int start = 4;
        unsigned int end   = N-k+2;
        unsigned int s     = 3;
#endif
#ifdef SCHEME_5
        unsigned int start = 5;
        unsigned int end   = N-k+3;
        unsigned int s     = 4;
#endif
        for (unsigned int n=start; n<=end; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = a(isn,tn)*h;

#ifdef SCHEME_1
            /* using 1nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_2
            /* using 2nd scheme, at point n=2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -40.0*alpha-1.0;
            double g3 = +12.0*alpha;
            double g4 = +8.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_3
            /* using 3nd scheme, at point n=3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +32.0*alpha;
            double g3 = -60.0*alpha-1.0;
            double g4 = +32.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_4
            /* using 4nd scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +8.0*alpha;
            double g3 = +12.0*alpha;
            double g4 = -40.0*alpha-1.0;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_5
            /* using 5nd scheme, at point n=5 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -112.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -208.0*alpha;
            double g5 = +70.0*alpha-1.0;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            fi /= +g1;
            g1 = 1.0;

            double A00 = A[0][0];
            A[0][0] = A[0][1] + g2*A00;
            A[0][1] = A[0][2] + g3*A00;
            A[0][2] = A[0][3] + g4*A00;
            A[0][3] = g5*A00;
            b[0]    = b[0] - fi*A00;

            //double x1 = (n+0)*hx;
            //double x2 = (n+1)*hx;
            //double x3 = (n+2)*hx;
            //double x4 = (n+3)*hx;
            //printf("%4d %18.10f %18.10f\n", n-s, b[0], A[0][0]*x1*x1*x1*tn.t+A[0][1]*x2*x2*x2*tn.t+A[0][2]*x3*x3*x3*tn.t+A[0][3]*x4*x4*x4*tn.t);

#ifdef __NORMALIZE__
            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;
#endif

            // printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);

            ems[n-s][0] = A[0][0];
            ems[n-s][1] = A[0][1];
            ems[n-s][2] = A[0][2];
            ems[n-s][3] = A[0][3];
            ems[n-s][4] = b[0];
        }

        /* using 2nd scheme, at point N-3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-3;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u[m-1][N-3] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 3rd scheme, at point N-2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-2;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u[m-1][N-2] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 4th scheme, at point N-1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-1;
        isn.x = isn.i*hx;
        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u[m-1][N-1] - (+22.0*alpha)*u[m][N] - ht*f(isn,tn);

        //double xN4 = (N-4)*hx;
        //double xN3 = (N-3)*hx;
        //double xN2 = (N-2)*hx;
        //double xN1 = (N-1)*hx;

        //puts("------------------------");
        //printf("0 %18.10f %18.10f\n", b[0], A[0][0]*xN4*xN4*xN4*tn.t+A[0][1]*xN3*xN3*xN3*tn.t+A[0][2]*xN2*xN2*xN2*tn.t+A[0][3]*xN1*xN1*xN1*tn.t);
        //printf("1 %18.10f %18.10f\n", b[1], A[1][0]*xN4*xN4*xN4*tn.t+A[1][1]*xN3*xN3*xN3*tn.t+A[1][2]*xN2*xN2*xN2*tn.t+A[1][3]*xN1*xN1*xN1*tn.t);
        //printf("2 %18.10f %18.10f\n", b[2], A[2][0]*xN4*xN4*xN4*tn.t+A[2][1]*xN3*xN3*xN3*tn.t+A[2][2]*xN2*xN2*xN2*tn.t+A[2][3]*xN1*xN1*xN1*tn.t);
        //printf("3 %18.10f %18.10f\n", b[3], A[3][0]*xN4*xN4*xN4*tn.t+A[3][1]*xN3*xN3*xN3*tn.t+A[3][2]*xN2*xN2*xN2*tn.t+A[3][3]*xN1*xN1*xN1*tn.t);
        //puts("------------------------");

        LinearEquation::GaussianElimination(A, b, x);

        //printf("x %14.10f %14.10f %14.10f %14.10f\n", x[0], x[1], x[2], x[3]);

        u[m][N-1] = x[3];
        u[m][N-2] = x[2];
        u[m][N-3] = x[1];
        u[m][N-4] = x[0];

        for (unsigned int n=end-s; n>=start-s; n--) u[m][n] = -(ems[n-1][1]*u[m][n+1] + ems[n-1][2]*u[m][n+2] + ems[n-1][3]*u[m][n+3] - ems[n-1][4]) / ems[n-1][0];

        //IPrinter::printVector(14, 10, u.row(m));
        //break;
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void ParabolicIBVP::calculateN4L2RDX(DoubleMatrix &u) const
{
    /* get parameters */
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM - minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN - minN;

    const unsigned int k = 4;
    double h = ht/(24.0*hx*hx);

    /*****************************************/

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    TimeNodePDE tn;
    SpaceNodePDE lsn;
    SpaceNodePDE rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

        /* using 2nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = minN+1;
        isn.x = isn.i*hx;
        double alpha = a(isn,tn)*h;
        A[0][0] = -40.0*alpha - 1.0;
        A[0][1] = +12.0*alpha;
        A[0][2] = +8.0*alpha;
        A[0][3] = -2.0*alpha;
        b[0]    = -u[m-1][1] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        /* using 5th scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        //isn.i = minN+4;
        //isn.x = isn.i*hx;
        //double alpha = a(isn,tn)*h;
        //A[0][0] = -112.0*alpha;
        //A[0][1] = +228.0*alpha;
        //A[0][2] = -208.0*alpha;
        //A[0][3] = +70.0*alpha - 1.0;
        //b[0]    = -u[m-1][4] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        //double x1 = (1)*hx;
        //double x2 = (2)*hx;
        //double x3 = (3)*hx;
        //double x4 = (4)*hx;
        //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*tn.t+A[0][1]*x2*x2*tn.t+A[0][2]*x3*x3*tn.t+A[0][3]*x4*x4*tn.t);

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = A[0][2];
        ems[0][2] = A[0][3];
        ems[0][3] = b[0];

#ifdef SCHEME_1
        unsigned int start = 1;
        unsigned int end   = N-k-1;
#endif
#ifdef SCHEME_2
        unsigned int start = 2;
        unsigned int end   = N-k+0;
#endif
#ifdef SCHEME_3
        unsigned int start = 3;
        unsigned int end   = N-k+1;
#endif
#ifdef SCHEME_4
        unsigned int start = 4;
        unsigned int end   = N-k+2;
#endif
#ifdef SCHEME_5
        unsigned int start = 5;
        unsigned int end   = N-k+3;
#endif
        for (unsigned int n=start; n<=end; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = a(isn,tn)*h;

#ifdef SCHEME_1
            /* using 1nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_2
            /* using 2nd scheme, at point n=2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -40.0*alpha-1.0;
            double g3 = +12.0*alpha;
            double g4 = +8.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_3
            /* using 3nd scheme, at point n=3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +32.0*alpha;
            double g3 = -60.0*alpha-1.0;
            double g4 = +32.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_4
            /* using 4nd scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +8.0*alpha;
            double g3 = +12.0*alpha;
            double g4 = -40.0*alpha-1.0;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_5
            /* using 5nd scheme, at point n=5 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -112.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -208.0*alpha;
            double g5 = +70.0*alpha-1.0;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif

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

            //double x1 = (n+1)*hx;
            //double x2 = (n+2)*hx;
            //double x3 = (n+3)*hx;
            //double x4 = (n+4)*hx;
            //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);
            \
            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            // printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);

#ifdef SCHEME_1
            ems[n+0][0] = A[0][1]; ems[n+0][1] = A[0][2]; ems[n+0][2] = A[0][3]; ems[n+0][3] = b[0];
#endif
#ifdef SCHEME_2
            ems[n-1][0] = A[0][1]; ems[n-1][1] = A[0][2]; ems[n-1][2] = A[0][3]; ems[n-1][3] = b[0];
#endif
#ifdef SCHEME_3
            ems[n-2][0] = A[0][1]; ems[n-2][1] = A[0][2]; ems[n-2][2] = A[0][3]; ems[n-2][3] = b[0];
#endif
#ifdef SCHEME_4
            ems[n-3][0] = A[0][1]; ems[n-3][1] = A[0][2]; ems[n-3][2] = A[0][3]; ems[n-3][3] = b[0];
#endif
#ifdef SCHEME_5
            ems[n-4][0] = A[0][1]; ems[n-4][1] = A[0][2]; ems[n-4][2] = A[0][3]; ems[n-4][3] = b[0];
#endif
        }

        /* using 2nd scheme, at point N-3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-3;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u[m-1][N-3] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 3rd scheme, at point N-2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-2;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u[m-1][N-2] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 4th scheme, at point N-1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-1;
        isn.x = isn.i*hx;
        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u[m-1][N-1] - (+22.0*alpha)*u[m][N] - ht*f(isn,tn);

        //double xN4 = (N-4)*hx;
        //double xN3 = (N-3)*hx;
        //double xN2 = (N-2)*hx;
        //double xN1 = (N-1)*hx;

        //puts("------------------------");
        //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*xN4*xN4*tn.t+A[0][1]*xN3*xN3*tn.t+A[0][2]*xN2*xN2*tn.t+A[0][3]*xN1*xN1*tn.t);
        //printf("1 %14.10f %14.10f\n", b[1], A[1][0]*xN4*xN4*tn.t+A[1][1]*xN3*xN3*tn.t+A[1][2]*xN2*xN2*tn.t+A[1][3]*xN1*xN1*tn.t);
        //printf("2 %14.10f %14.10f\n", b[2], A[2][0]*xN4*xN4*tn.t+A[2][1]*xN3*xN3*tn.t+A[2][2]*xN2*xN2*tn.t+A[2][3]*xN1*xN1*tn.t);
        //printf("3 %14.10f %14.10f\n", b[3], A[3][0]*xN4*xN4*tn.t+A[3][1]*xN3*xN3*tn.t+A[3][2]*xN2*xN2*tn.t+A[3][3]*xN1*xN1*tn.t);
        //puts("------------------------");

        LinearEquation::GaussianElimination(A, b, x);

        printf("x %14.10f %14.10f %14.10f %14.10f\n", x[0], x[1], x[2], x[3]);

        u[m][N-1] = x[3];
        u[m][N-2] = x[2];
        u[m][N-3] = x[1];
        u[m][N-4] = x[0];

#ifdef SCHEME_1
        for (unsigned int n=N-k-1; n>=1; n--) u[m][n+0] = -ems[n-1][0]*u[m][n+1] - ems[n-1][1]*u[m][n+2] - ems[n-1][2]*u[m][n+3] + ems[n-1][3];
#endif
#ifdef SCHEME_2
        for (unsigned int n=N-k+0; n>=2; n--) u[m][n-1] = -ems[n-2][0]*u[m][n+0] - ems[n-2][1]*u[m][n+1] - ems[n-2][2]*u[m][n+2] + ems[n-2][3];
#endif
#ifdef SCHEME_3
        for (unsigned int n=N-k+1; n>=3; n--) u[m][n-2] = -ems[n-3][0]*u[m][n-1] - ems[n-3][1]*u[m][n+0] - ems[n-3][2]*u[m][n+1] + ems[n-3][3];
#endif
#ifdef SCHEME_4
        for (unsigned int n=N-k+2; n>=4; n--) u[m][n-3] = -ems[n-4][0]*u[m][n-2] - ems[n-4][1]*u[m][n-1] - ems[n-4][2]*u[m][n+0] + ems[n-4][3];
#endif
#ifdef SCHEME_5
        for (unsigned int n=N-k+3; n>=5; n--) u[m][n-4] = -ems[n-5][0]*u[m][n-3] - ems[n-5][1]*u[m][n-2] - ems[n-5][2]*u[m][n-1] + ems[n-5][3];
#endif
        IPrinter::printVector(18, 10, u.row(m));
        break;
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void ParabolicIBVP::calculateN6L2RD(DoubleMatrix &u) const
{
    /* get parameters */
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM - minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN - minN;

    const unsigned int k = 6;
    double h = ht/(24.0*hx*hx);

    /*****************************************/

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);
    DoubleVector emk(N-k, k);

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    TimeNodePDE tn;
    SpaceNodePDE lsn;
    SpaceNodePDE rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

        /* using 2nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = minN+1;
        isn.x = isn.i*hx;
        double alpha = a(isn,tn)*h;
        A[0][0] = -40.0*alpha - 1.0;
        A[0][1] = +12.0*alpha;
        A[0][2] = +8.0*alpha;
        A[0][3] = -2.0*alpha;
        b[0]    = -u[m-1][1] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        /* using 5th scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        //isn.i = minN+4;
        //isn.x = isn.i*hx;
        //double alpha = a(isn,tn)*h;
        //A[0][0] = -112.0*alpha;
        //A[0][1] = +228.0*alpha;
        //A[0][2] = -208.0*alpha;
        //A[0][3] = +70.0*alpha - 1.0;
        //b[0]    = -u[m-1][4] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        //double x1 = (1)*hx;
        //double x2 = (2)*hx;
        //double x3 = (3)*hx;
        //double x4 = (4)*hx;
        //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*tn.t+A[0][1]*x2*x2*tn.t+A[0][2]*x3*x3*tn.t+A[0][3]*x4*x4*tn.t);

        //A[0][1] /= A[0][0];
        //A[0][2] /= A[0][0];
        //A[0][3] /= A[0][0];
        //b[0]    /= A[0][0];
        //A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = A[0][2];
        ems[0][2] = A[0][3];
        ems[0][3] = b[0];
        emk[0] = A[0][0];

#define SCHEME_11

#ifdef SCHEME_1
        unsigned int start = 1;
        unsigned int end   = N-k-1;
#endif
#ifdef SCHEME_2
        unsigned int start = 2;
        unsigned int end   = N-k+0;
#endif
#ifdef SCHEME_3
        unsigned int start = 3;
        unsigned int end   = N-k+1;
#endif
#ifdef SCHEME_4
        unsigned int start = 4;
        unsigned int end   = N-k+2;
#endif
#ifdef SCHEME_5
        unsigned int start = 5;
        unsigned int end   = N-k+3;
#endif
        for (unsigned int n=start; n<=end; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = a(isn,tn)*h;

#ifdef SCHEME_1
            /* using 1nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_2
            /* using 2nd scheme, at point n=2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -40.0*alpha-1.0;
            double g3 = +12.0*alpha;
            double g4 = +8.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_3
            /* using 3nd scheme, at point n=3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +32.0*alpha;
            double g3 = -60.0*alpha-1.0;
            double g4 = +32.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_4
            /* using 4nd scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +8.0*alpha;
            double g3 = +12.0*alpha;
            double g4 = -40.0*alpha-1.0;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_5
            /* using 5nd scheme, at point n=5 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -112.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -208.0*alpha;
            double g5 = +70.0*alpha-1.0;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            fi /= +g1;
            g1 = 1.0;

            double a00 = A[0][0];
            A[0][0] = A[0][1] + g2*a00;
            A[0][1] = A[0][2] + g3*a00;
            A[0][2] = A[0][3] + g4*a00;
            A[0][3] = g5*a00;
            b[0]    = b[0] - fi*a00;

            //double x1 = (n+1)*hx;
            //double x2 = (n+2)*hx;
            //double x3 = (n+3)*hx;
            //double x4 = (n+4)*hx;
            //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);
            \
            //A[0][1] /= A[0][0];
            //A[0][2] /= A[0][0];
            //A[0][3] /= A[0][0];
            //b[0]    /= A[0][0];
            //A[0][0] = 1.0;

            // printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);

#ifdef SCHEME_1
            ems[n+0][0] = A[0][1]; ems[n+0][1] = A[0][2]; ems[n+0][2] = A[0][3]; ems[n+0][3] = b[0]; emk[n+0] = A[0][0];
#endif
#ifdef SCHEME_2
            ems[n-1][0] = A[0][1]; ems[n-1][1] = A[0][2]; ems[n-1][2] = A[0][3]; ems[n-1][3] = b[0];
#endif
#ifdef SCHEME_3
            ems[n-2][0] = A[0][1]; ems[n-2][1] = A[0][2]; ems[n-2][2] = A[0][3]; ems[n-2][3] = b[0];
#endif
#ifdef SCHEME_4
            ems[n-3][0] = A[0][1]; ems[n-3][1] = A[0][2]; ems[n-3][2] = A[0][3]; ems[n-3][3] = b[0];
#endif
#ifdef SCHEME_5
            ems[n-4][0] = A[0][1]; ems[n-4][1] = A[0][2]; ems[n-4][2] = A[0][3]; ems[n-4][3] = b[0];
#endif
        }

        /* using 2nd scheme, at point N-3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-3;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u[m-1][N-3] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 3rd scheme, at point N-2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-2;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u[m-1][N-2] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 4th scheme, at point N-1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-1;
        isn.x = isn.i*hx;
        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u[m-1][N-1] - (+22.0*alpha)*u[m][N] - ht*f(isn,tn);

        //double xN4 = (N-4)*hx;
        //double xN3 = (N-3)*hx;
        //double xN2 = (N-2)*hx;
        //double xN1 = (N-1)*hx;

        //puts("------------------------");
        //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*xN4*xN4*tn.t+A[0][1]*xN3*xN3*tn.t+A[0][2]*xN2*xN2*tn.t+A[0][3]*xN1*xN1*tn.t);
        //printf("1 %14.10f %14.10f\n", b[1], A[1][0]*xN4*xN4*tn.t+A[1][1]*xN3*xN3*tn.t+A[1][2]*xN2*xN2*tn.t+A[1][3]*xN1*xN1*tn.t);
        //printf("2 %14.10f %14.10f\n", b[2], A[2][0]*xN4*xN4*tn.t+A[2][1]*xN3*xN3*tn.t+A[2][2]*xN2*xN2*tn.t+A[2][3]*xN1*xN1*tn.t);
        //printf("3 %14.10f %14.10f\n", b[3], A[3][0]*xN4*xN4*tn.t+A[3][1]*xN3*xN3*tn.t+A[3][2]*xN2*xN2*tn.t+A[3][3]*xN1*xN1*tn.t);
        //puts("------------------------");

        LinearEquation::GaussianElimination(A, b, x);

        //printf("x %14.10f %14.10f %14.10f %14.10f\n", x[0], x[1], x[2], x[3]);

        u[m][N-1] = x[3];
        u[m][N-2] = x[2];
        u[m][N-3] = x[1];
        u[m][N-4] = x[0];

#ifdef SCHEME_1
        for (unsigned int n=N-k-1; n>=1; n--) u[m][n+0] = (-ems[n-1][0]*u[m][n+1] - ems[n-1][1]*u[m][n+2] - ems[n-1][2]*u[m][n+3] + ems[n-1][3])/emk[n-1];
#endif
#ifdef SCHEME_2
        for (unsigned int n=N-k+0; n>=2; n--) u[m][n-1] = -ems[n-2][0]*u[m][n+0] - ems[n-2][1]*u[m][n+1] - ems[n-2][2]*u[m][n+2] + ems[n-2][3];
#endif
#ifdef SCHEME_3
        for (unsigned int n=N-k+1; n>=3; n--) u[m][n-2] = -ems[n-3][0]*u[m][n-1] - ems[n-3][1]*u[m][n+0] - ems[n-3][2]*u[m][n+1] + ems[n-3][3];
#endif
#ifdef SCHEME_4
        for (unsigned int n=N-k+2; n>=4; n--) u[m][n-3] = -ems[n-4][0]*u[m][n-2] - ems[n-4][1]*u[m][n-1] - ems[n-4][2]*u[m][n+0] + ems[n-4][3];
#endif
#ifdef SCHEME_5
        for (unsigned int n=N-k+3; n>=5; n--) u[m][n-4] = -ems[n-5][0]*u[m][n-3] - ems[n-5][1]*u[m][n-2] - ems[n-5][2]*u[m][n-1] + ems[n-5][3];
#endif
        //IPrinter::printVector(18, 10, u.row(m));
        //break;
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
    unsigned int M = time.max();
    unsigned int N1 = dim1.max();
    unsigned int N2 = dim2.max();

    double a1 = 1.0;
    double a2 = 1.0;

    //cleaning matrix
    u.clear();
    u.resize(N2+1, N1+1);

    DoubleMatrix v(N2+1, N1+1);

    double* da1 = (double*) malloc(sizeof(double)*(N1-1));
    double* db1 = (double*) malloc(sizeof(double)*(N1-1));
    double* dc1 = (double*) malloc(sizeof(double)*(N1-1));
    double* dd1 = (double*) malloc(sizeof(double)*(N1-1));
    double* rx1 = (double*) malloc(sizeof(double)*(N1-1));

    double* da2 = (double*) malloc(sizeof(double)*(N2-1));
    double* db2 = (double*) malloc(sizeof(double)*(N2-1));
    double* dc2 = (double*) malloc(sizeof(double)*(N2-1));
    double* dd2 = (double*) malloc(sizeof(double)*(N2-1));
    double* rx2 = (double*) malloc(sizeof(double)*(N2-1));

    double x1_a = -(a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*ht)/(h1*h1);
    double x1_c = (a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*ht)/(h2*h2);
    double x2_c = (a1*ht)/(2.0*h1*h1);

    for (unsigned int j=0; j<=N2; j++)
    {
        SpaceNodePDE sn;
        sn.j = j;
        sn.y = j*h2;
        for (unsigned int i=0; i<=N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            u[j][i] = initial(sn);
        }
    }

    TimeNodePDE tn;
    SpaceNodePDE sn;
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
                dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + 0.5*ht * f(sn, tn);
            }

            da1[0]     = 0.0;
            dc1[N1-2]  = 0.0;

            SpaceNodePDE sn0; sn0.i = 0;  sn0.x = 0.0;   sn0.j = j; sn0.y = j*h2;
            SpaceNodePDE snN; snN.i = N1; snN.x = N1*h1; snN.j = j; snN.y = j*h2;
            v[j][0]  = boundary(sn0, tn);
            v[j][N1] = boundary(snN, tn);

            dd1[0]    -= x1_a * v[j][0];
            dd1[N1-2] -= x1_a * v[j][N1];

            tomasAlgorithm(da1, db1, dc1, dd1, rx1, N1-1);

            for (unsigned int i=1; i<N1; i++) v[j][i] = rx1[i-1];
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            SpaceNodePDE sn0; sn0.i = i; sn0.x = i*h1; sn0.j = 0;  sn0.y = 0.0;
            SpaceNodePDE snN; snN.i = i; snN.x = i*h1; snN.j = N2; snN.y = N2*h2;
            v[0][i]  = boundary(sn0, tn);
            v[N2][i] = boundary(snN, tn);
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

                da2[j-1] = x2_a;
                db2[j-1] = x2_b;
                dc2[j-1] = x2_a;
                dd2[j-1] = x2_c*(v[j][i-1] - 2.0*v[j][i] + v[j][i+1]) + v[j][i] + 0.5*ht * f(sn, tn);
            }
            da2[0]     = 0.0;
            dc2[N2-2]  = 0.0;

            SpaceNodePDE sn0; sn0.i = i; sn0.x = i*h1; sn0.j = 0;  sn0.y = 0.0;
            SpaceNodePDE snN; snN.i = i; snN.x = i*h1; snN.j = N2; snN.y = N2*h2;
            u[0][i]  = boundary(sn0, tn);
            u[N2][i] = boundary(snN, tn);

            dd2[0]    -= x2_a * u[0][i];
            dd2[N2-2] -= x2_a * u[N2][i];

            tomasAlgorithm(da2, db2, dc2, dd2, rx2, N2-1);

            for (unsigned int j=1; j<N2; j++) u[j][i] = rx2[j-1];
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            SpaceNodePDE sn0; sn0.i = 0;  sn0.x = 0*h1;  sn0.j = j; sn0.y = j*h2;
            SpaceNodePDE snN; snN.i = N1; snN.x = N1*h1; snN.j = j; snN.y = j*h2;
            u[j][0]  = boundary(sn0, tn);
            u[j][N1] = boundary(snN, tn);
        }
    }

    free(rx2);
    free(dd2);
    free(dc2);
    free(db2);
    free(da2);

    free(rx1);
    free(dd1);
    free(dc1);
    free(db1);
    free(da1);
}

void ParabolicIBVP::calculateMVD_TEST(DoubleMatrix &u) const
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);
    Dimension dim2 = mspaceDimension.at(1);

    double ht = time.step();
    double h1 = dim1.step();
    double h2 = dim2.step();
    unsigned int M = time.max();
    unsigned int N1 = dim1.max();
    unsigned int N2 = dim2.max();

    double a1 = 1.0;
    double a2 = 1.0;

    //cleaning matrix
    u.clear();
    u.resize(N2+1, N1+1);

    DoubleMatrix v(N2+1, N1+1);

    double* da1 = (double*) malloc(sizeof(double)*(N1-1));
    double* db1 = (double*) malloc(sizeof(double)*(N1-1));
    double* dc1 = (double*) malloc(sizeof(double)*(N1-1));
    double* dd1 = (double*) malloc(sizeof(double)*(N1-1));
    double* rx1 = (double*) malloc(sizeof(double)*(N1-1));

    double* da2 = (double*) malloc(sizeof(double)*(N2-1));
    double* db2 = (double*) malloc(sizeof(double)*(N2-1));
    double* dc2 = (double*) malloc(sizeof(double)*(N2-1));
    double* dd2 = (double*) malloc(sizeof(double)*(N2-1));
    double* rx2 = (double*) malloc(sizeof(double)*(N2-1));

    double x1_a = -(a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*ht)/(h1*h1);
    double x1_c = (a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*ht)/(h2*h2);
    double x2_c = (a1*ht)/(2.0*h1*h1);

    for (unsigned int j=0; j<=N2; j++)
    {
        SpaceNodePDE sn;
        sn.j = j;
        sn.y = j*h2;
        for (unsigned int i=0; i<=N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            u[j][i] = initial(sn);
        }
    }

    TimeNodePDE tn, tn0;
    SpaceNodePDE sn;
    for (unsigned int k=1; k<=M; k++)
    {
        tn.i = k;    tn0.i = k-1;
        tn.t = k*ht; tn0.t = (k-1)*ht;

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
                //dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + 0.5*ht * f(sn, tn);
                //dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + 0.5*ht * f(sn, tn0);
                dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + 0.5*ht * (f(sn, tn)+f(sn, tn0))*0.5;
            }

            da1[0]     = 0.0;
            dc1[N1-2]  = 0.0;

            SpaceNodePDE sn0; sn0.i = 0;  sn0.x = 0.0;   sn0.j = j; sn0.y = j*h2;
            SpaceNodePDE snN; snN.i = N1; snN.x = N1*h1; snN.j = j; snN.y = j*h2;
            //v[j][0]  = boundary(sn0, tn);                          v[j][N1] = boundary(snN, tn);
            //v[j][0]  = boundary(sn0, tn0);                         v[j][N1] = boundary(snN, tn0);
            v[j][0]  = (boundary(sn0, tn)+boundary(sn0, tn0))*0.5; v[j][N1] = (boundary(snN, tn)+boundary(snN, tn0))*0.5;

            dd1[0]    -= x1_a * v[j][0];
            dd1[N1-2] -= x1_a * v[j][N1];

            tomasAlgorithm(da1, db1, dc1, dd1, rx1, N1-1);

            for (unsigned int i=1; i<N1; i++) v[j][i] = rx1[i-1];
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            SpaceNodePDE sn0; sn0.i = i; sn0.x = i*h1; sn0.j = 0;  sn0.y = 0.0;
            SpaceNodePDE snN; snN.i = i; snN.x = i*h1; snN.j = N2; snN.y = N2*h2;
            //v[0][i]  = boundary(sn0, tn);                          v[N2][i] = boundary(snN, tn);
            //v[0][i]  = boundary(sn0, tn0);                         v[N2][i] = boundary(snN, tn0);
            v[0][i]  = (boundary(sn0, tn)+boundary(sn0, tn0))*0.5; v[N2][i] = (boundary(snN, tn)+boundary(snN, tn0))*0.5;
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
                dd2[j-1] = x2_c*(v[j][i-1] - 2.0*v[j][i] + v[j][i+1]) + v[j][i] + 0.5*ht * f(sn, tn);
            }
            da2[0]     = 0.0;
            dc2[N2-2]  = 0.0;

            SpaceNodePDE sn0; sn0.i = i; sn0.x = i*h1; sn0.j = 0;  sn0.y = 0.0;
            SpaceNodePDE snN; snN.i = i; snN.x = i*h1; snN.j = N2; snN.y = N2*h2;
            u[0][i]  = boundary(sn0, tn);
            u[N2][i] = boundary(snN, tn);

            dd2[0]    -= x2_a * u[0][i];
            dd2[N2-2] -= x2_a * u[N2][i];

            tomasAlgorithm(da2, db2, dc2, dd2, rx2, N2-1);

            for (unsigned int j=1; j<N2; j++) u[j][i] = rx2[j-1];
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            SpaceNodePDE sn0; sn0.i = 0;  sn0.x = 0*h1;  sn0.j = j; sn0.y = j*h2;
            SpaceNodePDE snN; snN.i = N1; snN.x = N1*h1; snN.j = j; snN.y = j*h2;
            u[j][0]  = boundary(sn0, tn);
            u[j][N1] = boundary(snN, tn);
        }
    }

    free(rx2);
    free(dd2);
    free(dc2);
    free(db2);
    free(da2);

    free(rx1);
    free(dd1);
    free(dc1);
    free(db1);
    free(da1);
}
