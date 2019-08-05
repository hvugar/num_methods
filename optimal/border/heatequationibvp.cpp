#include "heatequationibvp.h"

void HeatEquationIBVP::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationIBVP w;
    w.setThermalDiffusivity(1.0);
    w.setThermalConductivity(0.0);
    w.setTimeDimension(Dimension(0.005, 0, 200));
    w.setSpaceDimensionX(Dimension(0.01, 0, 100));
    w.setSpaceDimensionY(Dimension(0.01, 0, 100));

    Benchmark bm;
    bm.tick();
    //w.implicit_calculate_D1V1();
    w.implicit_calculate_D2V1CN();
    //IPrinter::printSeperatorLine();
    bm.tock();
    bm.printDuration();
}

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    TimeNodePDE tn; tn.t = 0.0;
    return U(sn, tn);
}

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
    return U(sn, tn)*(condition.alpha()/condition.gamma());
}

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a = thermalDiffusivity();
    const double c = thermalConductivity();

    //return 1.0 - 6.0*a*a*sn.x + c * U(sn,tn);
    //return 2.0*tn.t - 6.0*a*a*sn.x + c * U(sn,tn);
    //return 3.0*tn.t*tn.t - 6.0*a*a*sn.x + c * U(sn,tn);
    //return 4.0*tn.t*tn.t*tn.t - 6.0*a*a*sn.x + c * U(sn,tn);
    //return 5.0*tn.t*tn.t*tn.t*tn.t - 6.0*a*a*sn.x + c * U(sn,tn);

    return 1.0 - 6.0*a*a*(sn.x+sn.y) + c * U(sn,tn);
    //return 2.0*tn.t - 6.0*a*a*(sn.x+sn.y) + c * U(sn,tn);
    //return 1.0 - 6.0*a*a*(sn.x+sn.y) + c * (sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t);
}

double HeatEquationIBVP::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    //return sn.x*sn.x*sn.x + tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t*tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t*tn.t*tn.t;

    return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t;
    //return (sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t)*(condition.alpha()/condition.gamma());
}

void HeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
    //if (tn.i==200 || tn.i==199 || tn.i==198 || tn.i==197 || tn.i==196 ||
    //    tn.i==4   || tn.i==3   || tn.i==2   || tn.i==1 || tn.i==0)
    {
        IPrinter::printVector(u);
        //IPrinter::printSeperatorLine();
    }

    if (tn.i==200)
    {
        double norm = 0.0;
        double max = 0.0;
        TimeNodePDE tn; tn.t = 1.0;
        SpaceNodePDE sn;
        for (unsigned int i=0; i<=100; i++)
        {
            sn.x = i*0.01;
            double k = 1.0; if (i==0 || i== 100) k = 0.5;
            norm += 0.01*k*(u[i]-U(sn, tn))*(u[i]-U(sn, tn));

            if (max < fabs(u[i]-U(sn, tn))) max = fabs(u[i]-U(sn, tn));
        }
        printf("norm: %.10f max: %.10f\n", sqrt(norm), max);
    }

    return;

}

void HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    if (tn.i==400 || tn.i==399 || tn.i==398 || /*tn.i==397 || tn.i==396 ||*/
        tn.i==4   || tn.i==3   || tn.i==2   || tn.i==1 || tn.i==0)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
    }

    if (tn.i==400)
    {
        double norm = 0.0;
        double max = 0.0;
        TimeNodePDE tn; tn.t = 1.0;
        SpaceNodePDE sn;
        for (unsigned int j=0; j<=100; j++)
        {
        for (unsigned int i=0; i<=100; i++)
        {
            sn.x = i*0.01;
            sn.y = j*0.01;
            double k1 = 1.0; if (j==0 || j== 100) k1 = 0.5;
            double k2 = 1.0; if (i==0 || i== 100) k2 = 0.5;
            norm += 0.01*0.01*k1*k2*(u[j][i]-U(sn, tn))*(u[j][i]-U(sn, tn));

            if (max < fabs(u[j][i]-U(sn, tn))) max = fabs(u[j][i]-U(sn, tn));
        }
        }
        printf("norm: %.10f max: %.10f\n", sqrt(norm), max);
    }

    return;
}

//---------------------------------------------------------------------------------------------//

void FinalHeatEquationIBVP::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

//    FinalHeatEquationIBVP w;
//    w.setThermalDiffusivity(1.0);
//    w.setTimeDimension(Dimension(0.005, 0, 200));
//    w.setSpaceDimensionX(Dimension(0.01, 0, 100));
//    w.setSpaceDimensionY(Dimension(0.01, 0, 100));

//    Benchmark bm;
//    bm.tick();
//    IPrinter::printSeperatorLine();
//    bm.tock();
//    bm.printDuration();
}

double FinalHeatEquationIBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    return 0.0;
}

double FinalHeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    return 0.0;
}

double FinalHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a = thermalDiffusivity();
    return 0.0;
}

void FinalHeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
}

void FinalHeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
}

