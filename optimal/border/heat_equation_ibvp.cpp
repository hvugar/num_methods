#include "heat_equation_ibvp.h"

void HeatEquationIBVP::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationIBVP w;
    w.setThermalDiffusivity(1.0);
    w.setThermalConductivity(0.1);
    w.setTimeDimension(Dimension(0.01, 0, 100));
    w.setSpaceDimensionX(Dimension(0.01, 0, 100));
    w.setSpaceDimensionY(Dimension(0.01, 0, 100));

    Benchmark bm;
    bm.tick();
    w.implicit_calculate_D1V1CN();
//    w.implicit_calculate_D2V1CN();
    //w.explicit_calculate_D2V1();
    //IPrinter::printSeperatorLine();
    bm.tock();
    bm.printDuration();
}

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    TimeNodePDE tn; tn.t = 0.0;
    return U(sn, tn);

    //return 0.0;
}

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    if (sn.i == 0)
    {
        //condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        //return U(sn, tn)*(condition.alpha()/condition.gamma());

        //condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +1.0);
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();

        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +1.0, +1.0, +1.0);
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();

        //double lambda = 0.1; condition = BoundaryConditionPDE(BoundaryCondition::Robin, -lambda, +1.0, -lambda); return 2.0;
        //condition = BoundaryConditionPDE(BoundaryCondition::Neumann, 0.0, +1.0, 0.0); return 0.0;
    }
    else
    {
        //condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        //return U(sn, tn)*(condition.alpha()/condition.gamma());

        //condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +1.0);
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();

        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +1.0, +1.0, +1.0);
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();

        //condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +0.0); return 0.0;
        //double lambda = 0.1; condition = BoundaryConditionPDE(BoundaryCondition::Robin, lambda, +1.0, lambda); return +1.0;
    }
}

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double td = thermalDiffusivity();
    const double tc = thermalConductivity();

    //return 1.0 - 0.0*td + tc*U(sn,tn);
    //return 2.0*tn.t - 0.0*td + tc*U(sn,tn);
    //return 3.0*tn.t*tn.t - 0.0*td + tc*U(sn,tn);

    //return 1.0 - 2.0*td + tc*U(sn,tn);
    return 2.0*tn.t - 2.0*td + tc*U(sn,tn);
    //return 3.0*tn.t*tn.t - 2.0*td + tc*U(sn,tn);

    //return 1.0 - 6.0*td*sn.x + tc*U(sn,tn);
    //return 2.0*tn.t - 6.0*td*sn.x + tc*U(sn,tn);
    //return 3.0*tn.t*tn.t - 6.0*td*sn.x + tc*U(sn,tn);

    //return 4.0*tn.t*tn.t*tn.t - 6.0*td*sn.x + tc*U(sn,tn);
    //return 5.0*tn.t*tn.t*tn.t*tn.t - 6.0*td*sn.x + tc*U(sn,tn);

    //return 1.0 - 6.0*td*(sn.x+sn.y) + tc*U(sn,tn);
    //return 2.0*tn.t - 6.0*td*(sn.x+sn.y) + tc*U(sn,tn);
    //return 3.0*tn.t*tn.t - 6.0*td*(sn.x+sn.y) + tc*U(sn,tn);

    //return 0.0;
}

double HeatEquationIBVP::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    //return sn.x + tn.t;
    //return sn.x + tn.t*tn.t;
    //return sn.x + tn.t*tn.t*tn.t;

    //return sn.x*sn.x + tn.t;
    return sn.x*sn.x + tn.t*tn.t;
    //return sn.x*sn.x + tn.t*tn.t*tn.t;

    //return sn.x*sn.x*sn.x + tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;

    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t*tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t*tn.t*tn.t;

    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t*tn.t;
}

void HeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
    //if (tn.i==200 || tn.i==199 || tn.i==198 || tn.i==197 || tn.i==196 ||
    //    tn.i==4   || tn.i==3   || tn.i==2   || tn.i==1 || tn.i==0)
    {
        IPrinter::printVector(u);
    }

    if (tn.i==100)
    {
        double L2Norm = 0.0;
        double EuNorm = 0.0;
        double L1Norm = 0.0;
        TimeNodePDE tn; tn.t = 1.0;
        SpaceNodePDE sn;
        for (unsigned int i=0; i<=100; i++)
        {
            sn.i = static_cast<int>(i);
            sn.x = sn.i*0.01;

            double k = 1.0; if (i==0 || i== 100) k = 0.5;
            L2Norm += 0.01*k*(u[i]-U(sn, tn))*(u[i]-U(sn, tn));

            EuNorm += (u[i]-U(sn, tn))*(u[i]-U(sn, tn));

            if (L1Norm < fabs(u[i]-U(sn, tn))) L1Norm = fabs(u[i]-U(sn, tn));
        }
        printf("L2Norm: %.10f EuNorm: %.10f L1Norm: %.10f\n", sqrt(L2Norm), sqrt(EuNorm), L1Norm);
    }

    return;

}

void HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    if (tn.i==40000 || /*tn.i==199 || tn.i==198 || tn.i==397 || tn.i==396 ||*/
        tn.i==4   || tn.i==3   || tn.i==2   || tn.i==1 || tn.i==0)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
    }

    if (tn.i==40000)
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
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    FinalHeatEquationIBVP w;
    w.setThermalDiffusivity(-1.0);
    w.setThermalConductivity(0.1);
    w.setTimeDimension(Dimension(0.01, 0, 100));
    w.setSpaceDimensionX(Dimension(0.01, 0, 100));
    w.setSpaceDimensionY(Dimension(0.01, 0, 100));

    Benchmark bm;
    bm.tick();
    w.implicit_calculate_D1V1CN();
//    w.implicit_calculate_D2V1CN();
    //w.explicit_calculate_D2V1();
    //IPrinter::printSeperatorLine();
    bm.tock();
    bm.printDuration();

}

double FinalHeatEquationIBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    TimeNodePDE tn; tn.t = 1.0;
    return U(sn, tn);
}

double FinalHeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    if (sn.i == 0)
    {
        //condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        //return U(sn, tn)*(condition.alpha()/condition.gamma());

        //condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +1.0);
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();

        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +1.0, +1.0, +1.0);
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();

        //double lambda = 0.1; condition = BoundaryConditionPDE(BoundaryCondition::Robin, -lambda, +1.0, -lambda); return 2.0;
        //condition = BoundaryConditionPDE(BoundaryCondition::Neumann, 0.0, +1.0, 0.0); return 0.0;
    }
    else
    {
        //condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        //return U(sn, tn)*(condition.alpha()/condition.gamma());

        //condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +1.0);
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();

        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +1.0, +1.0, +1.0);
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
        //return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();

        //condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +0.0); return 0.0;
        //double lambda = 0.1; condition = BoundaryConditionPDE(BoundaryCondition::Robin, lambda, +1.0, lambda); return +1.0;
    }
}

double FinalHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double td = thermalDiffusivity();
    const double uk = 0.0;
    const double tc = thermalConductivity();

    //return 1.0 - 0.0*td - 1.0*uk - tc*U(sn,tn);
    //return 2.0*tn.t - 0.0*td - tc*U(sn,tn);
    //return 3.0*tn.t*tn.t - 0.0*td - tc*U(sn,tn);

    //return 1.0 - 2.0*td - tc*U(sn,tn);
    return 2.0*tn.t - 2.0*td - tc*U(sn,tn);
    //return 3.0*tn.t*tn.t - 2.0*td - tc*U(sn,tn);

    //return 1.0 - 6.0*td*sn.x + tc*U(sn,tn);
    //return 2.0*tn.t - 6.0*td*sn.x + tc*U(sn,tn);
    //return 3.0*tn.t*tn.t - 6.0*td*sn.x + tc*U(sn,tn);

    //return 4.0*tn.t*tn.t*tn.t - 6.0*td*sn.x + tc*U(sn,tn);
    //return 5.0*tn.t*tn.t*tn.t*tn.t - 6.0*td*sn.x + tc*U(sn,tn);

    //return 1.0 - 6.0*td*(sn.x+sn.y) + tc*U(sn,tn);
    //return 2.0*tn.t - 6.0*td*(sn.x+sn.y) + tc*U(sn,tn);
    //return 3.0*tn.t*tn.t - 6.0*td*(sn.x+sn.y) + tc*U(sn,tn);
}

double FinalHeatEquationIBVP::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    //return sn.x + tn.t;
    //return sn.x + tn.t*tn.t;
    //return sn.x + tn.t*tn.t*tn.t;

    //return sn.x*sn.x + tn.t;
    return sn.x*sn.x + tn.t*tn.t;
    //return sn.x*sn.x + tn.t*tn.t*tn.t;

    //return sn.x*sn.x*sn.x + tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;

    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t*tn.t;
    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t*tn.t*tn.t;

    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t*tn.t;
}

void FinalHeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    //if (tn.i==200 || tn.i==199 || tn.i==198 || tn.i==197 || tn.i==196 ||
    //    tn.i==4   || tn.i==3   || tn.i==2   || tn.i==1 || tn.i==0)
    {
        IPrinter::printVector(u);
    }
}

void FinalHeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
}
