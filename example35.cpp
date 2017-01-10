#include "example35.h"

void Parabolic1DControl35::main(int argc, char ** argv)
{
    /* Function */
    Parabolic1DControl35 pc;

    DoubleVector v0((pc.M+1)*pc.L);
    for (unsigned int j=0; j<=pc.M; j++)
    {
        double t = j*pc.ht;
        v0[0*(pc.M+1)+j] = sin(t);
        v0[1*(pc.M+1)+j] = cos(t);
    }

    IPrinter::printVector(v0, "v1:", 10 , 0*(pc.M+1), 0*(pc.M+1)+pc.M, stdout);
    IPrinter::printVector(v0, "v2:", 10 , 1*(pc.M+1), 1*(pc.M+1)+pc.M, stdout);
	puts("-------------------------------------------------------------------");

	DoubleVector ga(v0.size());
    pc.gradient(v0, ga);
    ga.L2Normalize();
    IPrinter::printVector(ga, "ga1:", 10 , 0*(pc.M+1), 0*(pc.M+1)+pc.M, stdout);
    IPrinter::printVector(ga, "ga2:", 10 , 1*(pc.M+1), 1*(pc.M+1)+pc.M, stdout);
	puts("-------------------------------------------------------------------");
	
    DoubleVector gn(v0.size());
    IGradient::Gradient(&pc, pc.h, v0, gn);
    gn.L2Normalize();
    IPrinter::printVector(gn, "gn1:", 10 , 0*(pc.M+1), 0*(pc.M+1)+pc.M, stdout);
    IPrinter::printVector(gn, "gn2:", 10 , 1*(pc.M+1), 1*(pc.M+1)+pc.M, stdout);
	puts("-------------------------------------------------------------------");
	
	IPrinter::printVector(v0, "v1:", 10 , 0*(pc.M+1), 0*(pc.M+1)+pc.M, stdout);
    IPrinter::printVector(v0, "v2:", 10 , 1*(pc.M+1), 1*(pc.M+1)+pc.M, stdout);
	printf("J[%d]: %.16f\n", 0, pc.fx(v0));
	
    /* Minimization */
    ConjugateGradient g2;
    g2.setGradient(&pc);
    g2.setFunction(&pc);
    g2.setEpsilon1(0.00001);
    g2.setEpsilon2(0.00001);
    g2.setEpsilon3(0.00001);
    g2.setR1MinimizeEpsilon(0.1, 0.00001);
    g2.setPrinter(&pc);
    g2.setNormalize(true);
    g2.calculate(v0);

	//puts("-------------------------------------------------------------------");
    //IPrinter::printVector(v0, "v1:", 10 , 0*(pc.M+1), 0*(pc.M+1)+pc.M, stdout);
    //IPrinter::printVector(v0, "v2:", 10 , 1*(pc.M+1), 1*(pc.M+1)+pc.M, stdout);
}

Parabolic1DControl35::Parabolic1DControl35()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    hx = 0.01;
    ht = 0.01;
    N = 100;
    M = 100;
    L = 2;
    a = 1.0;
    h = 0.01;

    e.resize(L);
    e[0] = 0.25;
    e[1] = 0.65;

    DoubleVector v((M+1)*L);
    for (unsigned int j=0; j<=M; j++)
    {
        double t = j*ht;
        v[0*(M+1)+j] = v1(t);
        v[1*(M+1)+j] = v2(t);
    }
    pv = &v;
	//hx=0.001 ht=0.001
    IParabolicEquation::calculateU(U, hx, ht, N, M, a);
    IPrinter::printVector(U, "U: ");
    puts("-------------------------------------------------------------------");
}

double Parabolic1DControl35::fx(const DoubleVector &v)
{
    pv = &v;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    double sum = 0.0;
    double alpha;
    for (unsigned int i=0; i<=N; i++)
    {
        alpha = 1.0;
        if (i==0 || i==N) alpha = 0.5;
        sum += alpha*(u[i] - U[i])*(u[i] - U[i]);
    }
    sum = hx*sum;

    double norm = 0.0;
    for (unsigned int j=0; j<=M; j++)
    {
        double betta = 1.0;
        if (j==0 || j==M) alpha = 0.5;
        norm += betta*(v[0*(M+1)+j] - v1(j*ht))*(v[0*(M+1)+j] - v1(j*ht));
        norm += betta*(v[1*(M+1)+j] - v2(j*ht))*(v[1*(M+1)+j] - v2(j*ht));
    }
    norm = ht*norm;

    return sum + norm;
}

void Parabolic1DControl35::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    pu = &u;
    DoubleMatrix psi;
    IBackwardParabolicEquation::calculateU(psi, hx, ht, N, M, a);

    for (unsigned int j=0; j<=M; j++)
    {
        unsigned int i1 = (unsigned int)(e[0]/hx);
        unsigned int i2 = (unsigned int)(e[1]/hx);

        g[0*(M+1)+j] = -psi[j][i1] + 2.0*(v[0*(M+1)+j]-v1(j*ht));
        g[1*(M+1)+j] = -psi[j][i2] + 2.0*(v[1*(M+1)+j]-v2(j*ht));
    }
}

void Parabolic1DControl35::print(unsigned int i, const DoubleVector &v, const DoubleVector &gr, double alpha, RnFunction *fn) const
{
	puts("-------------------------------------------------------------------");
	IPrinter::printVector(v, "v1:", 10 , 0*(M+1), 0*(M+1)+M, stdout);
    IPrinter::printVector(v, "v2:", 10 , 1*(M+1), 1*(M+1)+M, stdout);
    printf("J[%d]: %.16f\n", i, fn->fx(v));
}

double Parabolic1DControl35::initial(unsigned int i) const
{
    double x = i*hx;
    return x+cos(x);
}

double Parabolic1DControl35::boundary(Boundary type, unsigned int j) const
{
     double t = j*ht;
     if (type == Left) return t*t+1.0;
     if (type == Right) return t*t+1.0+cos(1.0);
     return 0.0;
}

double Parabolic1DControl35::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    double sum = 0.0;
    double sgm = 3.0*hx;
    double gause_a = 1.0/(sqrt(2.0*M_PI)*sgm);
    double gause_b = 2.0*sgm*sgm;
    sum += v1(t) * gause_a * exp(-((x-e[0])*(x-e[0]))/gause_b);
    sum += v2(t) * gause_a * exp(-((x-e[1])*(x-e[1]))/gause_b);

//    if (i==e[0]/hx) sum += (1/hx)*v1(t);
//    if (i==e[1]/hx) sum += (1/hx)*v2(t);
    return sum;
}

double Parabolic1DControl35::binitial(unsigned int i) const
{
    return -2.0 * ((*pu)[i] - U[i]);
}

double Parabolic1DControl35::bboundary(Boundary type, unsigned int j) const
{
    return 0.0;
}

double Parabolic1DControl35::bf(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double Parabolic1DControl35::v1(double t) const
{
    return 10.0*t;
}

double Parabolic1DControl35::v2(double t) const
{
    return 12.0*t;
}

