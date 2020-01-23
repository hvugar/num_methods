#include "heat_equation_ibvp.h"

#define __DIRICHLET_
//#define __NEUMANN__
//#define __ROBIN__
//#define x3_t2
#define x1_y1_t1

void HeatEquationIBVP::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationIBVP h;
    h.setTimeDimension(Dimension(0.01, 0, 100));
    h.setSpaceDimensionX(Dimension(0.010, 100, 200));
    h.setSpaceDimensionY(Dimension(0.005, 200, 400));

    h.setThermalDiffusivity(1.2);
//    h.setThermalConductivity(-0.6);
//    h.setThermalConvection(-0.8);

    //    Benchmark bm;
    //    bm.tick();
    //    w.implicit_calculate_D1V1();
    //    bm.tock();
    //    bm.printDuration();

    //    bm.tick();
    //    w.implicit_calculate_D1V1_1();
    //    bm.tock();
    //    bm.printDuration();

    Benchmark bm;
    bm.tick();
    h.implicit_calculate_D2V1_2();
    bm.tock();
    bm.printDuration();

//    bm.tick();
//    h.implicit_calculate_D2V1();
//    bm.tock();
//    bm.printDuration();
}

HeatEquationIBVP::HeatEquationIBVP() : IHeatEquationIBVP()
{}

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition) const
{
    TimeNodePDE tn; tn.t = 0.0;
    return U(sn, tn);
}

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
#if defined (x1_y1_t1) || defined ( x2_y2_t2 ) || defined ( x2_y2_t1 )
    //condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +2.0, +0.0, +1.0);
    //return U(sn, tn)*(condition.alpha()/condition.gamma());

    condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
    return (condition.alpha()*U(sn,tn)+condition.beta()*Un(sn,tn))/condition.gamma();
#else


    const int min = spaceDimensionX().min();
    const int max = spaceDimensionX().max();

    if (sn.i == min)
    {
#if defined(__DIRICHLET_)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return U(sn, tn)*(condition.alpha()/condition.gamma());
#endif

#if defined(__NEUMANN__)
        condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +1.0);
#if defined(x1_t1) || defined(x1_t2) || defined(x1_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
#endif
#if defined(x2_t1) || defined(x2_t2) || defined(x2_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
#endif
#if defined(x3_t1) || defined(x3_t2) || defined(x3_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();
#endif
#endif

#if defined(__ROBIN__)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
#if defined(x1_t1) || defined(x1_t2) || defined(x1_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
#endif
#if defined(x2_t1) || defined(x2_t2) || defined(x2_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
#endif
#if defined(x3_t1) || defined(x3_t2) || defined(x3_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();
#endif
#endif
    }

    if (sn.i == max)
    {
#if defined(__DIRICHLET_)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return U(sn, tn)*(condition.alpha()/condition.gamma());
#endif

#if defined(__NEUMANN__)
#if defined(x1_t1) || defined(x1_t2) || defined(x1_t3)
        condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +1.0);
        return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
#endif
#if defined(x2_t1) || defined(x2_t2) || defined(x2_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
#endif
#if defined(x3_t1) || defined(x3_t2) || defined(x3_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();
#endif
#endif

#if defined(__ROBIN__)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
#if defined(x1_t1) || defined(x1_t2) || defined(x1_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
#endif
#if defined(x2_t1) || defined(x2_t2) || defined(x2_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
#endif
#if defined(x3_t1) || defined(x3_t2) || defined(x3_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();
#endif
#endif
    }
#endif
}

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a1 = thermalDiffusivity();
    const double a2 = thermalDiffusivity();
    const double b1 = thermalConductivity();
    const double b2 = thermalConductivity();
    const double c = thermalConvection();

#if defined( x1_t1 )
    return 1.0 - 0.0*td + tc*U(sn,tn) - 1.0*tv;
#elif defined( x1_t2 )
    return 2.0*tn.t - 0.0*td + tc*U(sn,tn) - 1.0*tv;
#elif defined( x1_t3 )
    return 3.0*tn.t*tn.t - 0.0*td + tc*U(sn,tn) - 1.0*tv;
#endif

#if defined( x2_t1 )
    return 1.0 - 2.0*td - 2.0*cv*sn.x - tc*U(sn,tn);
#elif defined( x2_t2 )
    return 2.0*tn.t - 2.0*a - 2.0*b*sn.x - c*U(sn,tn);
#elif defined( x2_t3 )
    return 3.0*tn.t*tn.t - 2.0*td - 2.0*cv*sn.x - tc*U(sn,tn);
#endif

#if defined( x3_t1 )
    return 1.0 - 6.0*td*sn.x - 3.0*cv*sn.x*sn.x - tc*U(sn,tn);
#elif defined( x3_t2 )
    return 2.0*tn.t - 6.0*a*sn.x - 3.0*b*sn.x*sn.x - c*U(sn,tn);
#elif defined( x3_t3 )
    return 3.0*tn.t*tn.t - 6.0*td*sn.x + tc*U(sn,tn) - 3.0*tv*sn.x*sn.x;
#endif

#if defined (x1_y1_t1) || defined( x2_y2_t2 ) || defined( x2_y2_t1 )
    return Ut(sn,tn) - (a1*Uxx(sn,tn) + a2*Uyy(sn,tn) + b1*Ux(sn,tn) + b2*Uy(sn,tn) + c*U(sn,tn));
#endif
}

double HeatEquationIBVP::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x1_y1_t1 )
    return sn.x + sn.y + tn.t;
#endif
#if defined( x2_y2_t1 )
    return sn.x*sn.x + sn.y*sn.y + tn.t;
#endif
#if defined( x2_y2_t2 )
    return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
#endif
}

double HeatEquationIBVP::Ut(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x1_y1_t1 ) || defined( x2_y2_t1 )
    return 1.0;
#endif
#if defined( x2_y2_t2 )
    return 2.0*tn.t;
#endif
}

double HeatEquationIBVP::Uxx(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x1_y1_t1 )
    return 0.0;
#endif
#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
    return 2.0;
#endif
}

double HeatEquationIBVP::Uyy(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x1_y1_t1 )
    return 0.0;
#endif
#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
    return 2.0;
#endif
}

double HeatEquationIBVP::Ux(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x1_y1_t1 )
    return 1.0;
#endif
#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
    return 2.0*sn.x;
#endif
}

double HeatEquationIBVP::Uy(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x1_y1_t1 )
    return 1.0;
#endif
#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
    return 2.0*sn.y;
#endif
}

double HeatEquationIBVP::Un(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

    const int xmin = spaceDimensionX().min();
    const int xmax = spaceDimensionX().max();
    const int ymin = spaceDimensionY().min();
    const int ymax = spaceDimensionY().max();

#if defined( x1_y1_t1 )
    //printf("%d %d %f %f\n",sn.i, sn.j, sn.x, sn.y);
    if (tn.i%2==1)
    {
        if (sn.i == xmin) return -1.0;
        if (sn.i == xmax) return +1.0;
    }
    if (tn.i%2==0)
    {
        if (sn.j == ymin) return -1.0;
        if (sn.j == ymax) return +1.0;
    }
#endif

#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
    if (tn.i%2==1)
    {
        if (sn.i == xmin) return -2.0*sn.x;
        if (sn.i == xmax) return +2.0*sn.x;
    }
    if (tn.i%2==0)
    {
        if (sn.j == ymin) return -2.0*sn.y;
        if (sn.j == ymax) return +2.0*sn.y;
    }
#endif
}

void HeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    int min = spaceDimensionX().min();
    int max = spaceDimensionX().max();

    if (tn.i==0 || tn.i==1 || tn.i==99 || tn.i==100) IPrinter::printVector(u);

    if (tn.i==1001)
    {
        double L2Norm = 0.0;
        double EuNorm = 0.0;
        double L1Norm = 0.0;
        TimeNodePDE tn; tn.t = 1.0;
        SpaceNodePDE sn;
        unsigned int n = 0;
        for (int i=min; i<=max; i++, n++)
        {
            sn.i = static_cast<int>(i);
            sn.x = sn.i*0.01;

            double k = 1.0; if (i==min || i== max) k = 0.5;
            L2Norm += 0.01*k*(u[n]-U(sn, tn))*(u[n]-U(sn, tn));

            EuNorm += (u[n]-U(sn, tn))*(u[n]-U(sn, tn));

            if (L1Norm < fabs(u[n]-U(sn, tn))) L1Norm = fabs(u[n]-U(sn, tn));
        }
        printf("L2Norm: %.10f EuNorm: %.10f L1Norm: %.10f\n", sqrt(L2Norm), sqrt(EuNorm), L1Norm);
    }
}

void HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    if (tn.i==0 || tn.i==1 || tn.i==2 || tn.i==3 || tn.i==4 || tn.i==5 || tn.i==6 || tn.i==199 || tn.i==200)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
        if (tn.i%2==0) IPrinter::printSeperatorLine();
    }
    return;

    //    IPrinter::printMatrix(u);
    //    IPrinter::printSeperatorLine();
    return;

    if (tn.i==200 || /*tn.i==199 || tn.i==198 || tn.i==397 || tn.i==396 ||*/
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

void HeatEquationFBVP::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationFBVP w;
    w.setTimeDimension(Dimension(0.01, 0, 100));
    w.setSpaceDimensionX(Dimension(0.01, 100, 200));
    //w.setSpaceDimensionY(Dimension(0.01, 100, 200));

    w.setThermalDiffusivity(-1.2);
    w.setThermalConductivity(+0.6);
    w.setThermalConvection(+0.8);

    Benchmark bm;
    bm.tick();
    w.implicit_calculate_D1V1();
    bm.tock();
    bm.printDuration();

    bm.tick();
    w.implicit_calculate_D1V1_1();
    bm.tock();
    bm.printDuration();
}

HeatEquationFBVP::HeatEquationFBVP() : IHeatEquationFBVP()
{}

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition) const
{
    TimeNodePDE tn; tn.t = 1.0;
    return U(sn, tn);
}

double HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    const int min = spaceDimensionX().min();
    const int max = spaceDimensionX().max();

    if (sn.i == min)
    {
#if defined(__DIRICHLET_)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return U(sn, tn)*(condition.alpha()/condition.gamma());
#endif

#if defined(__NEUMANN__)
        condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +1.0);
#if defined(x1_t1) || defined(x1_t2) || defined(x1_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
#endif
#if defined(x2_t1) || defined(x2_t2) || defined(x2_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
#endif
#if defined(x3_t1) || defined(x3_t2) || defined(x3_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();
#endif
#endif

#if defined(__ROBIN__)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
#if defined(x1_t1) || defined(x1_t2) || defined(x1_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
#endif
#if defined(x2_t1) || defined(x2_t2) || defined(x2_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
#endif
#if defined(x3_t1) || defined(x3_t2) || defined(x3_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();
#endif
#endif
    }

    if (sn.i == max)
    {
#if defined(__DIRICHLET_)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return U(sn, tn)*(condition.alpha()/condition.gamma());
#endif

#if defined(__NEUMANN__)
        condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +1.0);
#if defined(x1_t1) || defined(x1_t2) || defined(x1_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
#endif
#if defined(x2_t1) || defined(x2_t2) || defined(x2_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
#endif
#if defined(x3_t1) || defined(x3_t2) || defined(x3_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();
#endif
#endif

#if defined(__ROBIN__)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
#if defined(x1_t1) || defined(x1_t2) || defined(x1_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(1.0))/condition.gamma();
#endif
#if defined(x2_t1) || defined(x2_t2) || defined(x2_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(2.0*sn.x))/condition.gamma();
#endif
#if defined(x3_t1) || defined(x3_t2) || defined(x3_t3)
        return (condition.alpha()*U(sn, tn)+condition.beta()*(3.0*sn.x*sn.x))/condition.gamma();
#endif
#endif
    }
}

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();

#if defined( x1_t1 )
    return 1.0 + 0.0*a + 1.0*b - tc*U(sn,tn);
#elif defined( x1_t2 )
    return 2.0*tn.t + 0.0*td - tc*U(sn,tn) + 1.0*tv;
#elif defined( x1_t3 )
    return 3.0*tn.t*tn.t + 0.0*td - tc*U(sn,tn) + 1.0*tv;
#endif

#if defined( x2_t1 )
    return 1.0 + 2.0*td - tc*U(sn,tn) + 2.0*tv*sn.x;
#elif defined( x2_t2 )
    return 2.0*tn.t - 2.0*a - 2.0*b*sn.x - c*U(sn,tn);
#elif defined( x2_t3 )
    return 3.0*tn.t*tn.t + 2.0*td - tc*U(sn,tn) + 2.0*tv*sn.x;
#endif

#if defined( x3_t1 )
    return 1.0 + 6.0*td*sn.x - tc*U(sn,tn) + 3.0*tv*sn.x*sn.x;
#elif defined( x3_t2 )
    return 2.0*tn.t - 6.0*a*sn.x - 3.0*b*sn.x*sn.x - c*U(sn,tn);
#elif defined( x3_t3 )
    return 3.0*tn.t*tn.t + 6.0*td*sn.x - tc*U(sn,tn) + 3.0*tv*sn.x*sn.x;
#endif

    return 0.0;
}

double HeatEquationFBVP::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

    return 0.0;
}

void HeatEquationFBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    if (tn.i==0 || tn.i==1 || tn.i==99 || tn.i==100) IPrinter::printVector(u);
}

void HeatEquationFBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
}
