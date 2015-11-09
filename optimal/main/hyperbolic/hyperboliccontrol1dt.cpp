#include "hyperboliccontrol1dt.h"
#include "hyperboliccontrol1d2.h"

HyperbolicControl1DT::HyperbolicControl1DT()
{
}

HyperbolicControl1DT::~HyperbolicControl1DT()
{
}

double HyperbolicControl1DT::fx(double t)
{
    printf("t: %f\n", t);
    DoubleVector v;
    HyperbolicControl1D2 hc;
    hc.t1 = t;
    hc.doSettings();

    v.resize(2*(hc.M+hc.DM+1));
    hc.initialize();
    for (unsigned int j=0; j<=hc.M+hc.DM; j++)
    {
        double t = j*hc.ht;
        v[j] = 2.0*t;
        v[(hc.M+hc.DM+1)+j] = 2.0*t + 2.0;
    }
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.00001);
    g2.setEpsilon2(0.00001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.01, 0.00001);
    //g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(false);
    g2.calculate(v);

    DoubleMatrix u;
    hc.calculateU(v, u);
    Printer::printVector(u[hc.M]);

    //DoubleVector v1(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v1[j] = v[j];
    //DoubleVector v2(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v2[j] = v[hc.M+hc.DM+1+j];
    //Printer::printVector(v1);
    //Printer::printVector(v2);

    double rf = hc.fx(v);
    printf("************************\n%f %.16f\n", t, rf);
    return hc.fx(v);
}

void HyperbolicControl1DT::main()
{
    HyperbolicControl1DT hc;

    double alpha0 = 0.6;
    double a,b,alpha;
    R1Minimize::StranghLineSearch(alpha0, 0.2, a, b, &hc);
    R1Minimize::GoldenSectionSearch(a, b, alpha, &hc, 0.0001);
    printf("alpha: %f\n", alpha);
}

