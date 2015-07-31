#include "samplecontrol.h"
#include "cfunction1.h"
#include "cfunction2.h"


SampleControl::SampleControl() : ConjugateGradient()
{}

SampleControl::~SampleControl()
{
}

void SampleControl::calculateGradient()
{
    CFunction1* func = dynamic_cast<CFunction1*>(function());

    func->gradientJ(grad_step, m_g, m_x);

    printf("J[%d] = %.10f\n", iterationCount, func->fx(m_x));
}

void SampleControl::print()
{
    //    SteepestDescentGradient::print();
}

void SampleControl::Main()
{
    CFunction1* f = new CFunction1(0.0, 1.0, 0.01);
    DoubleVector u(f->n);
    for (int i=0; i<f->n; i++) u[i] = 0.00001;//3*f->t[i];//*f->t[i];

    SampleControl sc;
    sc.setFunction(f);
    sc.setX(u);
    sc.setEpsilon(0.01);
    sc.setGradientStep(0.000001);
    sc.setR1MinimizeEpsilon(0.01, 0.000001);
    sc.calculate();

    /* Function */
    CFunction1 c(0.0, 1.0, 0.01);

    /* initial point */
    DoubleVector u0(c.n);

    for (int i=0; i<c.n; i++) u0[i] = 0.00001;
    /* Minimization */
    ConjugateGradient g1;
    g1.setFunction(&c);
    g1.setEpsilon(0.01);
    g1.setGradientStep(0.000001);
    g1.setR1MinimizeEpsilon(0.01, 0.000001);
    g1.setX(u0);
    g1.setPrinter(new CFunction1Printer);
    g1.calculate();

    puts("-----------------------------------------------------------------");
    for (int i=0; i<c.n; i++) u0[i] = 0.00001;
    /* Minimization */
    SteepestDescentGradient g2;
    g2.setFunction(&c);
    g2.setEpsilon(0.01);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.01, 0.000001);
    g2.setX(u0);
    g2.setPrinter(new CFunction1Printer);
    g2.calculate();
}
