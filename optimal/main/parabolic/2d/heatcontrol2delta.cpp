#include "heatcontrol2delta.h"

// Optimal points [0.50,0.80] [0.70,0.20] [0.20,0.30]
// Working points [0.60,0.70] [0.65,0.25] [0.25,0.35] epsilon1 0.0001 epsilon2 0.0001 epsilon3: 0.0001. min:1.0 0.0001
// Working points [0.60,0.70] [0.60,0.30] [0.30,0.40] epsilon1 0.0001 epsilon2 0.0001 epsilon3: 0.0001. min:1.0 0.0001

void HeatControl2Delta::main()
{
    HeatControl2Delta hc(100, 100, 100);

    hc.O.resize(2*hc.L);
    hc.O[0] = 0.50;
    hc.O[1] = 0.80;
    hc.O[2] = 0.70;
    hc.O[3] = 0.20;
    hc.O[4] = 0.20;
    hc.O[5] = 0.30;

    hc.initialize();

    DoubleVector x(2*hc.L + (hc.M+1)*hc.L);
    x[0] = 0.60; x[1] = 0.70; x[2] = 0.65; x[3] = 0.25; x[4] = 0.25; x[5] = 0.35;

    for (unsigned int k=0; k<=hc.M; k++)
    {
        x[2*hc.L + 0*(hc.M+1) + k] = 1.0;//hc.v1(k*hc.ht);
        x[2*hc.L + 1*(hc.M+1) + k] = 1.0;//hc.v2(k*hc.ht);
        x[2*hc.L + 2*(hc.M+1) + k] = 1.0;//hc.v3(k*hc.ht);
    }

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setGradient(&hc);
    g2.setEpsilon1(0.0001);
    g2.setEpsilon2(0.0001);
    g2.setEpsilon3(0.0001);
    g2.setR1MinimizeEpsilon(1.0, 0.0001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(x);

    DoubleVector gr1(x.size());
    IGradient::Gradient(&hc, 0.00001, x, gr1);
    gr1.L2Normalize();

    DoubleVector gr2(x.size());
    hc.gradient(x, gr2);
    gr2.L2Normalize();

    printf("J[%d]: %.16f\n", 0, hc.fx(x));
    printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", hc.O[0], hc.O[1], hc.O[2], hc.O[3], hc.O[4], hc.O[5]);
    printf("e1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
    printf("gr1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", gr1[0], gr1[1], gr1[2], gr1[3], gr1[4], gr1[5]);
    printf("gr2: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", gr2[0], gr2[1], gr2[2], gr2[3], gr2[4], gr2[5]);
}

HeatControl2Delta::HeatControl2Delta(unsigned int M, unsigned int N2, unsigned int N1)
{
    alpha = 1.0;

    this->M  = M;
    this->N2 = N2;
    this->N1 = N1;
    this->L  = 3;

    t0 = 0.0;
    t1 = 1.0;

    x10 = 0.0;
    x11 = 1.0;

    x20 = 0.0;
    x21 = 1.0;

    a1 = a2 = 1.0;

    ht = (t1-t0)/M;
    h1 = (x11-x10)/N1;
    h2 = (x21-x20)/N2;

    double sgm1 = 3.0*h1;
    double sgm2 = 3.0*h2;
    gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    gause_b = 2.0*sgm1*sgm2;
}

double HeatControl2Delta::fx(const DoubleVector& x)
{
    px = &x;
    DoubleMatrix u;
    IParabolicEquation2D::caluclateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

    double sum = 0.0;
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            double betta = 1.0;
            if (i==0 || i==N1 || j==0 || j==N2) betta = 0.5;
            if (i==0  && j==0)  betta = 0.25;
            if (i==0  && j==N2) betta = 0.25;
            if (i==N1 && j==0)  betta = 0.25;
            if (i==N1 && j==N2) betta = 0.25;
            sum += betta*(u[j][i] - U[j][i])*(u[j][i] - U[j][i]);
        }
    }
    sum = (h1*h2)*sum;

    double nrm = 0.0;
    //nrm = norm(x);
    return sum + alpha*nrm;
}

double HeatControl2Delta::norm(const DoubleVector& x) const
{
    double nrm = 0.0;
    for (unsigned int k=0; k<=M; k++)
    {
        double betta = 1.0;
        if (k==0 || k==M) betta = 0.5;
        nrm += betta*(x[2*L+0*(M+1)+k] - v1(k*ht))*(x[2*L+0*(M+1)+k] - v1(k*ht));
        nrm += betta*(x[2*L+1*(M+1)+k] - v2(k*ht))*(x[2*L+1*(M+1)+k] - v2(k*ht));
        nrm += betta*(x[2*L+2*(M+1)+k] - v3(k*ht))*(x[2*L+2*(M+1)+k] - v3(k*ht));
    }
    nrm = ht * nrm;
    if (nrm < 0.00001) const_cast<HeatControl2Delta*>(this)->alpha = 0.0;
    return nrm;
}

void HeatControl2Delta::gradient(const DoubleVector& x, DoubleVector& g)
{
    px = &x;
    DoubleMatrix u;
    IParabolicEquation2D::caluclateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

    pu = &u;
    DoubleCube psi;
    IBackwardParabolicEquation2D::caluclateMVD(psi, h1, h2, ht, N1, N2, M, a1, a2);

    for (unsigned int i=0; i<g.size(); i++) g[i] = 0.0;

    for (unsigned int k=M; k!=(unsigned int)0-1; k--)
    {
        calculateGX(x, psi[k], g, k);
        //calculateGF(x, psi[k], g, k);
    }

//    for (unsigned int k=0; k<=M; k++)
//    {
//        unsigned int i1 = (unsigned int)round(E[0]/h1);
//        unsigned int j1 = (unsigned int)round(E[1]/h2);
//        g[0*(M+1)+k] = -psi[k][j1][i1] + 2.0*alpha*(v[0*(M+1)+k] - v1(k*ht));

//        unsigned int i2 = (unsigned int)round(E[2]/h1);
//        unsigned int j2 = (unsigned int)round(E[3]/h2);
//        g[1*(M+1)+k] = -psi[k][j2][i2] + 2.0*alpha*(v[1*(M+1)+k] - v2(k*ht));

//        unsigned int i3 = (unsigned int)round(E[4]/h1);
//        unsigned int j3 = (unsigned int)round(E[5]/h2);
//        g[2*(M+1)+k] = -psi[k][j3][i3] + 2.0*alpha*(v[2*(M+1)+k] - v3(k*ht));
//    }

    psi.clear();
//    IGradient::Gradient(this, 0.0001, x, g);
}

void HeatControl2Delta::calculateGX(const DoubleVector& x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k)
{
    double psiX1;
    double psiX2;
    if (k==0 || k==M)
    {
        psiDerivative(psiX1, psiX2, x[0], x[1], psi);
        g[0] = g[0] + v1(k*ht) * psiX1;
        g[1] = g[1] + v1(k*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], psi);
        g[2] = g[2] + v2(k*ht) * psiX1;
        g[3] = g[3] + v2(k*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[4], x[5], psi);
        g[4] = g[4] + v3(k*ht) * psiX1;
        g[5] = g[5] + v3(k*ht) * psiX2;
    }
    else
    {
        psiDerivative(psiX1, psiX2, x[0], x[1], psi);
        g[0] = g[0] + 2.0*v1(k*ht) * psiX1;
        g[1] = g[1] + 2.0*v1(k*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], psi);
        g[2] = g[2] + 2.0*v2(k*ht) * psiX1;
        g[3] = g[3] + 2.0*v2(k*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[4], x[5], psi);
        g[4] = g[4] + 2.0*v3(k*ht) * psiX1;
        g[5] = g[5] + 2.0*v3(k*ht) * psiX2;
    }

    if (k==0)
    {
        g[0] = -(ht/2.0)*g[0];
        g[1] = -(ht/2.0)*g[1];
        g[2] = -(ht/2.0)*g[2];
        g[3] = -(ht/2.0)*g[3];
        g[4] = -(ht/2.0)*g[4];
        g[5] = -(ht/2.0)*g[5];
    }
}

void HeatControl2Delta::psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi)
{
    unsigned int i = (unsigned int)round(x1/h1);
    unsigned int j = (unsigned int)round(x2/h2);

    if (i==0) psiX1  = (psi[j][i+1] - psi[j][i])/h1;
    else if (i==N1) psiX1 = (psi[j][i] - psi[j][i-1])/h1;
    else psiX1 = (psi[j][i+1] - psi[j][i-1])/(2.0*h1);

    if (j==0) psiX2 = (psi[j+1][i] - psi[j][i])/h2;
    else if (j==N2) psiX2 = (psi[j][i] - psi[j-1][i])/h2;
    else psiX2 = (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
}

void HeatControl2Delta::calculateGF(const DoubleVector &x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k)
{
    C_UNUSED(k);
    for (unsigned k=0; k<=M; k++)
    {
        if (alpha < 1.0)
        {
            g[2*L+0*(M+1)+k] = 0.0;
            g[2*L+1*(M+1)+k] = 0.0;
            g[2*L+2*(M+1)+k] = 0.0;
        }
        else
        {
            double p[L];
            p[0]=p[1]=p[2]=0.0;
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    if (fabs(x[0]-i*h1)<h1 && fabs(x[1]-j*h2)<h2) p[0] += psi[j][i]*((h1-fabs(x[0]-i*h1))/h1)*((h2-fabs(x[1]-j*h2))/h2);
                    if (fabs(x[2]-i*h1)<h1 && fabs(x[3]-j*h2)<h2) p[1] += psi[j][i]*((h1-fabs(x[2]-i*h1))/h1)*((h2-fabs(x[3]-j*h2))/h2);
                    if (fabs(x[4]-i*h1)<h1 && fabs(x[5]-j*h2)<h2) p[2] += psi[j][i]*((h1-fabs(x[4]-i*h1))/h1)*((h2-fabs(x[5]-j*h2))/h2);
                }
            }
            g[2*L+0*(M+1)+k] = -p[0] + 2.0*alpha*(x[2*L+0*(M+1)+k] - v1(k*ht));
            g[2*L+1*(M+1)+k] = -p[1] + 2.0*alpha*(x[2*L+1*(M+1)+k] - v2(k*ht));
            g[2*L+2*(M+1)+k] = -p[2] + 2.0*alpha*(x[2*L+2*(M+1)+k] - v3(k*ht));
        }
    }
}

void HeatControl2Delta::calculateG2(const DoubleVector &x, DoubleVector& g1)
{
    C_UNUSED(x);
    C_UNUSED(g1);

    //    double h = 0.01;
    //    DoubleVector E(2*L);
    //    DoubleVector g(2*L);
    //    double f0 = fx(x);

    //    for (unsigned int i=0; i<x.size(); i++)
    //    {
    //        E = x;
    //        E[i] += h;
    //        double f1 = fx(E);
    //        g[i] = (f1-f0)/h;
    //    }

    //    printf("e2: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", x[0], x[1], x[2], x[3], x[4], x[5]);
    //    printf("g2: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", g[0], g[1], g[2], g[3], g[4], g[5]);
}

double HeatControl2Delta::initial(unsigned int i, unsigned int j) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    return u(x1, x2, t0);
}

double HeatControl2Delta::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    double x2 = j*h2;
    double t  = 0.5*k*ht;
    return u(x10, x2, t);
}

double HeatControl2Delta::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    //double t  = 0.5*k*ht;
    double k1 = (k%2==0 ? k/2 : (k+1)/2);
    const DoubleVector &x = *px;

    double _v1 = x[2*L+0*(M+1)+k1];
    double _v2 = x[2*L+1*(M+1)+k1];
    double _v3 = x[2*L+2*(M+1)+k1];

    double sum = 0.0;
    sum += _v1 * gause_a * exp(-((x1-x[0])*(x1-x[0]) + (x2-x[1])*(x2-x[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-x[2])*(x1-x[2]) + (x2-x[3])*(x2-x[3]))/gause_b);
    sum += _v3 * gause_a * exp(-((x1-x[4])*(x1-x[4]) + (x2-x[5])*(x2-x[5]))/gause_b);

    //    if (fabs(x1-e[0])<=h1 && fabs(x2-e[1])<=h2)
    //    {
    //        sum += g1(t) * ((h1-fabs(x1-e[0]))/(h1*h1))*((h2-fabs(x2-e[1]))/(h2*h2));
    //    }
    //    if (fabs(x1-e[2])<=h1 && fabs(x2-e[3])<=h2)
    //    {
    //        sum += g2(t) * ((h1-fabs(x1-e[2]))/(h1*h1))*((h2-fabs(x2-e[3]))/(h2*h2));
    //    }
    //    if (fabs(x1-e[4])<=h1 && fabs(x2-e[5])<=h2)
    //    {
    //        sum += g3(t) * ((h1-fabs(x1-e[4]))/(h1*h1))*((h2-fabs  v c  v (x2-e[5]))/(h2*h2));
    //    }
    return sum;
}

double HeatControl2Delta::binitial(unsigned int i, unsigned int j) const
{
    return -2.0*((*pu)[j][i] - U[j][i]);
}

double HeatControl2Delta::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HeatControl2Delta::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

void HeatControl2Delta::initialize()
{
    // cleaning U
    for (unsigned int j=0; j<U.size(); j++) U[j].clear();
    U.clear();

    DoubleVector x = O;
    // initializing U
    U.resize(N2+1);
    for (unsigned int j=0; j<=N2; j++) U[j].resize(N1+1);

#ifdef POWER_OPTIMIZE
    x.resize( 2*L + (M+1)*L );
#else
    x.resize( 2*L );
#endif

#ifdef POWER_OPTIMIZE
    for (unsigned int k=0; k<=M; k++)
    {
        x[2*L + 0*(M+1) + k] = g1(k*ht);
        x[2*L + 1*(M+1) + k] = g2(k*ht);
        x[2*L + 2*(M+1) + k] = g3(k*ht);
    }
#endif

    px = &x;
    IParabolicEquation2D::caluclateMVD(U, h1, h2, ht, N1, N2, M, a1, a2);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    IPrinter::printMatrix(U, 10, 10);
    printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

    write("optimal.txt", U);

    FILE* f = fopen("optimal1.txt", "w");
    IPrinter::printMatrix(U, N1, N2, NULL, f);
    fclose(f);
}

void HeatControl2Delta::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double alpha, RnFunction* fn) const
{
    C_UNUSED(alpha);

    HeatControl2Delta *hc = dynamic_cast<HeatControl2Delta*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(x));
    //printf("Norm: %.16f Alpha: %.16f %.16f\n", hc->norm(x), hc->alpha, alpha);
    printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", O[0], O[1], O[2], O[3], O[4], O[5]);
    printf("e1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
    printf("g1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", g[0], g[1], g[2], g[3], g[4], g[5]);

#ifdef POWER_OPTIMIZE
    IPrinter::printVector(x, "f1:", 10, 2*L+0*(M+1), 2*L+0*(M+1) + M);
    IPrinter::printVector(x, "f2:", 10, 2*L+1*(M+1), 2*L+1*(M+1) + M);
    IPrinter::printVector(x, "f3:", 10, 2*L+2*(M+1), 2*L+2*(M+1) + M);
#endif

    //    DoubleMatrix u;
    //    hc->calculateU(x, u);
    //    char filename1[100];
    //    char filename2[100];
    //    int count1 = sprintf(filename1, "optimal%d.txt", i);
    //    filename1[count1] = '\0';
    //    int count2 = sprintf(filename2, "optimal%d.png", i);
    //    filename2[count2] = '\0';
    //    hc->write(filename2, u);
    //    system("imager.exe -w 101 -h 101 -i optimal1.txt -o optimal3.png");

    //    DoubleVector f(hc->M+1);
    //    for (unsigned int k=0; k<=M; k++) f[k] = x[2*hc->L+0*(M+1)+k];
    //    Printer::printVector(f, 10, "g1");
    //    for (unsigned int k=0; k<=M; k++) f[k] = x[2*hc->L+1*(M+1)+k];
    //    Printer::printVector(f, 10, "g2");
    //    for (unsigned int k=0; k<=M; k++) f[k] = x[2*hc->L+2*(M+1)+k];
    //    Printer::printVector(f, 10, "g3");

    //    DoubleVector fg(hc->M+1);
    //    for (unsigned int k=0; k<=M; k++) fg[k] = g[2*hc->L+0*(M+1)+k];
    //    Printer::printVector(fg, 10, "fg1");
    //    for (unsigned int k=0; k<=M; k++) fg[k] = g[2*hc->L+1*(M+1)+k];
    //    Printer::printVector(fg, 10, "fg2");
    //    for (unsigned int k=0; k<=M; k++) fg[k] = g[2*hc->L+2*(M+1)+k];
    //    Printer::printVector(fg, 10, "fg3");
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

    //    hc->calculateU(e, hc->uT);
    //    char buffer [12];
    //    int n=sprintf (buffer, "file%d.txt", i);
    //    hc->write(buffer, hc->uT);
}

void HeatControl2Delta::project(DoubleVector &e, int index)
{
    if (index<6)
    {
        if (e[index] > 1.0) e[index] = 1.0;
        if (e[index] < 0.0) e[index] = 0.0;
    }
    //    for (unsigned int i=0; i<e.size(); i++)
    //    {
    //        if (e[i]>1.0) e[i]=1.0;
    //        if (e[i]<0.0) e[i]=0.0;
    //    }
}

void HeatControl2Delta::write(const char *fileName, const DoubleMatrix& m)
{
    FILE* f = fopen(fileName, "w");
    for (unsigned int j=0; j<m.size(); j++)
    {
        for (unsigned int i=0; i<m[j].size(); i++)
        {
            if (i==0)
                fprintf(f, "%.10f", m[j][i]);
            else
                fprintf(f, " %.10f", m[j][i]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}
