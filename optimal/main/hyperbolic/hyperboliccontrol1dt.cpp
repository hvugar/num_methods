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
    printf("M: %d N: %d DM: %d ht: %.8f hx: %.8f\n", hc.M, hc.N, hc.DM, hc.ht, hc.hx);

    DoubleVector v;
    v.resize(2*(hc.M+hc.DM+1));

    //printf("M: %d\n", hc.M);

    for (unsigned int j=0; j<=(hc.M+hc.DM); j++)
    {
        //double t = j*hc.ht;
        v[j] = hc.UT;///t*t+1.5;
        v[(hc.M+hc.DM+1)+j] = hc.UT;//t*t+2.5;
    }

    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.000000001);
    g2.setEpsilon2(0.000000001);
    //g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.1, 0.000000001);
    //g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(v);

    DoubleMatrix u;
    hc.calculateU(v, u);

    //DoubleVector v1(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v1[j] = v[j];
    //DoubleVector v2(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v2[j] = v[hc.M+hc.DM+1+j];
    //DoubleVector v1(hc.DM+1); for (unsigned j=hc.M; j<=hc.M+hc.DM; j++) v1[j] = v[j];
    //DoubleVector v2(hc.DM+1); for (unsigned j=hc.M; j<=hc.M+hc.DM; j++) v2[j] = v[hc.M+hc.DM+1+j];
    //Printer::printVector(v1, 10, "v1:\t");
    //Printer::printVector(v2, 10, "v2:\t");
    double rf = hc.fx(v);

    //puts("++++++++++");
    //FILE* f = fopen("vugar1.txt", "a");
    //fprintf(f, "------------------------------------------------------------\n");
    //fprintf(f, "T: %.8f Integral: %.16f\n", t, rf);
    for (unsigned int j=hc.M; j<=hc.M+hc.DM; j++)
    {
        printf("u[%d]:\t", j);
        Printer::printVector(u[j]);

        //fprintf(f, "u[%d]:\t", j);
        //for (unsigned int i=0; i<=hc.N; i++)
        //{
        //    fprintf(f, "%.8f ", u[j][i]);
        //}
        //fprintf(f, "\n");
    }
    //fclose(f);
    //Printer::printVector(u[hc.M], 10, "u[M]:\t");
    //Printer::printVector(u[hc.M+hc.DM], 10, "u[M+DM]:");

    printf("T: %.8f Integral: %.16f M: %d\n", t, rf, hc.M);
    printf("%.8f %.16f\n", t, rf);
    return rf;
}

void HyperbolicControl1DT::main()
{
    HyperbolicControl1DT hc;
    //hc.fx(1.0);

    for (double t=0.1; t<1.11; t+=0.1)
        hc.fx(t);

//    double t0 = 0.8;
//    double a,b,t;
//    R1Minimize::StranghLineSearch(t0, 0.2, a, b, &hc);
//    printf("a: %f b: %f T: %f\n", a, b, t);
//    R1Minimize::GoldenSectionSearch(a, b, t, &hc, 0.000001);
////    printf("a: %f b: %f T: %f\n", a, b, t);
//    printf("T: %f Integral: %.16f\n", t, hc.fx(t));

//    for (double t=0.8; t<2.0; t+=0.01)
    {
       // hc.fx(1.0);
    }

}

