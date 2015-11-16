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
    puts("------------------------------------------------------------------------------------------------------");
    printf("T: %.8f\n", t);

    HyperbolicControl1D2 hc(0.0, t);

    DoubleVector v;
    v.resize(2*(hc.M+hc.DM+1));

    printf("M: %d\n", hc.M);

    for (unsigned int j=0; j<=(hc.M+hc.DM); j++)
    {
        double t = j*hc.ht;
        v[j] = t*t+1.5;
        v[(hc.M+hc.DM+1)+j] = t*t+2.5;
    }

    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.000000001);
    g2.setEpsilon2(0.000000001);
    //g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.2, 0.000000001);
    //g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(v);

    DoubleMatrix u;
    hc.calculateU(v, u);

    DoubleVector v1(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v1[j] = v[j];
    DoubleVector v2(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v2[j] = v[hc.M+hc.DM+1+j];
    Printer::printVector(v1, 10, "v1:\t");
    Printer::printVector(v2, 10, "v2:\t");

    puts("++++++++++");
    FILE* f = fopen("d:/u1000.txt", "a");
    fprintf(f, "------------------------------------------------------------\n");
    fprintf(f, "T: %f\n", t);
    for (unsigned int j=hc.M; j<=hc.M+hc.DM; j++)
    {
        printf("u[%d]:\t", j);
        Printer::printVector(u[j]);

        for (unsigned int i=0; i<=hc.N; i++)
        {
            fprintf(f, "%.8f ", u[j][i]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
//    Printer::printVector(u[hc.M], 10, "u[M]:\t");
//    Printer::printVector(u[hc.M+hc.DM], 10, "u[M+DM]:");

    double rf = hc.fx(v);
    printf("Integral: %.16f\n", rf);
    return rf;
}

void HyperbolicControl1DT::main()
{
    HyperbolicControl1DT hc;

    double t0 = 0.8;
    double a,b,t;
    R1Minimize::StranghLineSearch(t0, 0.2, a, b, &hc);
    R1Minimize::GoldenSectionSearch(a, b, t, &hc, 0.000001);
    printf("T: %f\n", t);
}

