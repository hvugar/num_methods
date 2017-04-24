#include "art_problem1L3.h"

void ArticleProblem1L3::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ArticleProblem1L3 a3;
}

ArticleProblem1L3::ArticleProblem1L3()
{
    L = 3;

    alpha0 = 1.0;
    alpha1 = 1.0;
    alpha2 = 1.0;
    alpha3 = 1.0;

    a = 1.0;
    lambda0 = 0.001;
    lambda1 = 1000.0;
    lambda2 = 0.0010;

    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    N = 1000;
    M = 1000;
    hx = 0.001;
    ht = 0.001;

    fi = 1.0;
    tt = 2.0;
    vfi << fi;
    vtt << tt;

    k0 << -5.20 << -2.80 << -8.70;
    z0 << +10.0 << +5.30 << +4.32;
    e0 << +0.25 << +0.50 << +0.75;
//    k0 << +0.00 << -0.00 << -0.00;
//    z0 << +0.00 << +0.00 << +0.00;
//    e0 << +0.00 << +0.00 << +0.00;

    R = 1.0;
    DD = 1;
    vmin = 0.0;
    vmax = 80.0;
    d0 = (vmax+vmin)/2.0;
    d1 = (vmax-vmin)/2.0;

    DoubleVector k,z,e;
    k << -5.20 << -2.80 << -8.70;
    z << +10.0 << +5.30 << +4.32;
    e << +0.25 << +0.50 << +0.75;
    DoubleMatrix u;
    calculateU(u, k,z,e);
    V.clear();
    V.resize(N+1);
    for (unsigned int n=0; n<=N; n++) V[n] = 5.895;
//    V = u.row(M);
    IPrinter::printVector(V);

    DoubleVector y0;
    y0 << -1.20 << -1.80 << -0.70;
    y0 << -10.0 << -5.30 << -4.32;
    y0 << +0.25 << +0.50 << +0.75;

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

void ArticleProblem1L3::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
{
    if (i == 6) { if (x.at(6) < 0.005) x.at(6) = 0.005; if (x.at(6) > 0.095) x.at(6) = 0.095; }
    if (i == 7) { if (x.at(7) < 0.005) x.at(7) = 0.005; if (x.at(7) > 0.095) x.at(7) = 0.095; }
    if (i == 8) { if (x.at(8) < 0.005) x.at(8) = 0.005; if (x.at(8) > 0.095) x.at(8) = 0.095; }

}
