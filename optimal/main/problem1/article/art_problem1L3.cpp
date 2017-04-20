#include "art_problem1L3.h"

void ArticleProblem1L3::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ArticleProblem1L3 a3;
}

ArticleProblem1L3::ArticleProblem1L3()
{
    L = 3;

    alpha0 = 1.0;
    alpha1 = 0.0001;
    alpha2 = 0.0001;
    alpha3 = 0.0001;

    lambda0 = 0.0001;
    lambda1 = 1000.0;
    lambda2 = 0.0010;

    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    N = 100;
    M = 100;
    hx = 0.01;
    ht = 0.01;

    fi = 0.0;
    tt = 2.0;
    vfi << fi;
    vtt << tt;

    V.clear();
    V.resize(N+1);
    for (unsigned int n=0; n<=N; n++) V[n] = 5.0;

    k0 << -5.20 << -2.80 << -8.70;
    z0 << +10.0 << +5.30 << +4.32;
    e0 << +0.25 << +0.50 << +0.75;

    R = 0.0;
    DD = 1;
    vmin = 0.0;
    vmax = 0.0;
    d0 = 0.0;
    d1 = 0.0;

    DoubleVector y0;
    y0 << +5.20 << +2.80 << +8.70;
    y0 << -10.0 << -5.30 << -4.32;
    y0 << +0.52 << +0.05 << +0.57;

    optimize(y0);
}

void ArticleProblem1L3::optimize(DoubleVector &y0) const
{
    ArticleProblem1L3* p = const_cast<ArticleProblem1L3*>(this);
    ConjugateGradient g;
    g.setFunction(p);
    g.setGradient(p);
    g.setPrinter(p);
    g.setProjection(p);
    g.setEpsilon1(0.000001);//0.00000001
    g.setEpsilon2(0.000001);//0.00000001
    g.setEpsilon3(0.000001);//0.00000001
    g.setR1MinimizeEpsilon(10.0, 0.00001); //0.00000001
    g.setNormalize(true);
    g.showEndMessage(true);
    g.setResetIteration(false);
    g.calculate(y0);
}

void ArticleProblem1L3::print(unsigned int i, const DoubleVector &y, const DoubleVector &g, double r, GradientMethod::MethodResult) const
{
    ArticleProblem1L3 *pm = const_cast<ArticleProblem1L3*>(this);
    DoubleVector k,z,e;
    getParameters(k,z,e,y);

    IPrinter::printSeperatorLine();

    DoubleVector ng(y.size());
    IGradient::Gradient(pm, 0.001, y, ng);

    DoubleMatrix u;
    pm->py = &y;
    pm->calculateU(u, k, z, e);

    double v = vf(M,k,z,e,u);

    DoubleVector ag = g;

    DoubleVector nag = ag;
    DoubleVector nng = ng;

    nag.L2Normalize();
    nng.L2Normalize();

    printf("J[%d]: %.10f v: %.10f\n", i, r, v);
    printf("k: %14.10f %14.10f %14.10f\n", k[0], k[1], k[2]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ag[0], ag[1], ag[2], nag[0], nag[1], nag[2]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ng[0], ng[1], ng[2], nng[0], nng[1], nng[2]);
    printf("z: %14.10f %14.10f %14.10f\n", z[0], z[1], z[2]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ag[3], ag[4], ag[5], nag[3], nag[4], nag[5]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ng[3], ng[4], ng[5], nng[3], nng[4], nng[5]);
    printf("e: %14.10f %14.10f %14.10f\n", e[0], e[1], e[2]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ag[6], ag[7], ag[8], nag[6], nag[7], nag[8]);
    printf("a: %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f\n", ng[6], ng[7], ng[8], nng[6], nng[7], nng[8]);
    IPrinter::printVector(14,10,u.row(u.rows()-1),"u: ");
}

void ArticleProblem1L3::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
{
    if (i == 6) { if (x.at(6) < 0.05) x.at(6) = 0.05; if (x.at(6) > 0.95) x.at(6) = 0.95; }
    if (i == 7) { if (x.at(7) < 0.05) x.at(7) = 0.05; if (x.at(7) > 0.95) x.at(7) = 0.95; }
    if (i == 8) { if (x.at(8) < 0.05) x.at(8) = 0.05; if (x.at(8) > 0.95) x.at(8) = 0.95; }

}
