#include "heat_equation_ibvp.h"

#define __DIRICHLET_
//#define __NEUMANN__
//#define __ROBIN__
//#define x3_t2
//#define x2_y2_t1

#define __DIMENSION_2__
#define __QUADRATIC__

#if defined(__QUADRATIC__)
#define x1
#define y1
#define t1
#endif

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
    h.setSpaceDimensionY(Dimension(0.005, 400, 600));

    h.setThermalDiffusivity(1.2);
    h.setThermalConductivity(-0.6);
    h.setThermalConvection(-0.8);

    Benchmark bm;
    bm.tick();
    //    h.implicit_calculate_D1V1();
    h.implicit_calculate_D2V1();
    bm.tock();
    bm.printDuration();
}

HeatEquationIBVP::HeatEquationIBVP() : IHeatEquationIBVP()
{}

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition) const
{
    TimeNodePDE tn; tn.t = 0.0;
    return ::u_fx(this, sn, tn);
}

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
#if defined(__DIMENSION_2__)
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +2.0, +0.0, +1.0);
    return u_fx(this, sn, tn)*(condition.alpha()/condition.gamma());

    //condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
    //return (condition.alpha()*::u_fx(this, sn, tn)
    //        +condition.beta()*::u_fx(this, sn, tn, -1, 3, 3))/condition.gamma();
#else


    const int min = spaceDimensionX().min();
    const int max = spaceDimensionX().max();

    if (sn.i == min)
    {
#if defined(__DIRICHLET_)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return ::u_fx(sn, tn)*(condition.alpha()/condition.gamma());
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
        return u_fx(sn, tn)*(condition.alpha()/condition.gamma());
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
#if defined(__DIMENSION_1__)
    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();
    return ::u_fx(sn,tn,1,true,false,false)
            -::u_fx(sn,tn,2,false,true,false)*a
            -::u_fx(sn,tn,1,false,true,false)*b
            -::u_fx(sn,tn,0,false,true,false)*c;
#endif

#if defined(__DIMENSION_2__)
    const double a1 = thermalDiffusivity();
    const double a2 = thermalDiffusivity();
    const double b1 = thermalConductivity();
    const double b2 = thermalConductivity();
    const double c = thermalConvection();
    return ::u_fx(this, sn,tn,1,0,0)
            -::u_fx(this, sn,tn,-1,+2,-1)*a1
            -::u_fx(this, sn,tn,-1,-1,+2)*a2
            -::u_fx(this, sn,tn,-1,+1,-1)*b1
            -::u_fx(this, sn,tn,-1,-2,+1)*b2
            -::u_fx(this, sn,tn,+0,+0,+0)*c;
#endif


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
    return 2.0*tn.t - 6.0*a1*sn.x - 3.0*b1*sn.x*sn.x - c*U(sn,tn);
#elif defined( x3_t3 )
    return 3.0*tn.t*tn.t + 6.0*td*sn.x - tc*U(sn,tn) + 3.0*tv*sn.x*sn.x;
#endif

#if defined (x1_y1_t1) || defined( x2_y2_t2 ) || defined( x2_y2_t1 )
    return Ut(sn,tn) - (a1*Uxx(sn,tn) + a2*Uyy(sn,tn) + b1*Ux(sn,tn) + b2*Uy(sn,tn) + c*U(sn,tn));
#endif
}

//double HeatEquationIBVP::Ut(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
//{
//#if defined( x1_t1 )
//    return sn.x + tn.t;
//#elif defined( x1_t2 )
//    return sn.x + tn.t*tn.t;
//#elif defined( x1_t3 )
//    return sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x2_t1 )
//    return sn.x*sn.x + tn.t;
//#elif defined( x2_t2 )
//    return sn.x*sn.x + tn.t*tn.t;
//#elif defined( x2_t3 )
//    return sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x3_t1 )
//    return sn.x*sn.x*sn.x + tn.t;
//#elif defined( x3_t2 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t;
//#elif defined( x3_t3 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x1_y1_t1 ) || defined( x2_y2_t1 )
//    return 1.0;
//#endif
//#if defined( x2_y2_t2 )
//    return 2.0*tn.t;
//#endif
//}

//double HeatEquationIBVP::Uxx(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
//{
//#if defined( x1_t1 )
//    return sn.x + tn.t;
//#elif defined( x1_t2 )
//    return sn.x + tn.t*tn.t;
//#elif defined( x1_t3 )
//    return sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x2_t1 )
//    return sn.x*sn.x + tn.t;
//#elif defined( x2_t2 )
//    return sn.x*sn.x + tn.t*tn.t;
//#elif defined( x2_t3 )
//    return sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x3_t1 )
//    return sn.x*sn.x*sn.x + tn.t;
//#elif defined( x3_t2 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t;
//#elif defined( x3_t3 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x1_y1_t1 )
//    return 0.0;
//#endif
//#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
//    return 2.0;
//#endif
//}

//double HeatEquationIBVP::Uyy(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
//{
//#if defined( x1_t1 )
//    return sn.x + tn.t;
//#elif defined( x1_t2 )
//    return sn.x + tn.t*tn.t;
//#elif defined( x1_t3 )
//    return sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x2_t1 )
//    return sn.x*sn.x + tn.t;
//#elif defined( x2_t2 )
//    return sn.x*sn.x + tn.t*tn.t;
//#elif defined( x2_t3 )
//    return sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x3_t1 )
//    return sn.x*sn.x*sn.x + tn.t;
//#elif defined( x3_t2 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t;
//#elif defined( x3_t3 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x1_y1_t1 )
//    return 0.0;
//#endif
//#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
//    return 2.0;
//#endif
//}

//double HeatEquationIBVP::Ux(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
//{
//#if defined( x1_t1 )
//    return sn.x + tn.t;
//#elif defined( x1_t2 )
//    return sn.x + tn.t*tn.t;
//#elif defined( x1_t3 )
//    return sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x2_t1 )
//    return sn.x*sn.x + tn.t;
//#elif defined( x2_t2 )
//    return sn.x*sn.x + tn.t*tn.t;
//#elif defined( x2_t3 )
//    return sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x3_t1 )
//    return sn.x*sn.x*sn.x + tn.t;
//#elif defined( x3_t2 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t;
//#elif defined( x3_t3 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x1_y1_t1 )
//    return 1.0;
//#endif
//#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
//    return 2.0*sn.x;
//#endif
//}

//double HeatEquationIBVP::Uy(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
//{
//#if defined( x1_t1 )
//    return sn.x + tn.t;
//#elif defined( x1_t2 )
//    return sn.x + tn.t*tn.t;
//#elif defined( x1_t3 )
//    return sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x2_t1 )
//    return sn.x*sn.x + tn.t;
//#elif defined( x2_t2 )
//    return sn.x*sn.x + tn.t*tn.t;
//#elif defined( x2_t3 )
//    return sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x3_t1 )
//    return sn.x*sn.x*sn.x + tn.t;
//#elif defined( x3_t2 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t;
//#elif defined( x3_t3 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x1_y1_t1 )
//    return 1.0;
//#endif
//#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
//    return 2.0*sn.y;
//#endif
//}

//double HeatEquationIBVP::Un(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
//{
//#if defined( x1_t1 )
//    return sn.x + tn.t;
//#elif defined( x1_t2 )
//    return sn.x + tn.t*tn.t;
//#elif defined( x1_t3 )
//    return sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x2_t1 )
//    return sn.x*sn.x + tn.t;
//#elif defined( x2_t2 )
//    return sn.x*sn.x + tn.t*tn.t;
//#elif defined( x2_t3 )
//    return sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//#if defined( x3_t1 )
//    return sn.x*sn.x*sn.x + tn.t;
//#elif defined( x3_t2 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t;
//#elif defined( x3_t3 )
//    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
//#endif

//    const int xmin = spaceDimensionX().min();
//    const int xmax = spaceDimensionX().max();
//    const int ymin = spaceDimensionY().min();
//    const int ymax = spaceDimensionY().max();

//#if defined( x1_y1_t1 )
//    //printf("%d %d %f %f\n",sn.i, sn.j, sn.x, sn.y);

//    if (tn.i%2==1)
//    {
//        if (sn.i == xmin) return -1.0;
//        if (sn.i == xmax) return +1.0;
//        if (sn.j == ymin) return -1.0;
//        if (sn.j == ymax) return +1.0;
//    }

//    if (tn.i%2==0)
//    {
//        if (sn.j == ymin) return -1.0;
//        if (sn.j == ymax) return +1.0;
//    }
//#endif

//#if defined( x2_y2_t1 ) || defined( x2_y2_t2 )
//    if (tn.i%2==1)
//    {
//        if (sn.i == xmin && (sn.j != ymin && sn.j != ymax)) return -2.0*sn.x;
//        if (sn.i == xmax && (sn.j != ymin && sn.j != ymax)) return +2.0*sn.x;
//        if (sn.j == ymin && (sn.j != ymin || sn.j != ymax)) return -2.0*sn.y;
//        if (sn.j == ymax && (sn.j != ymin || sn.j != ymax)) return +2.0*sn.y;
//    }
//    if (tn.i%2==0)
//    {
//        if (sn.j == ymin && (sn.i != xmin && sn.i != xmax)) return -2.0*sn.y;
//        if (sn.j == ymax && (sn.i != xmin && sn.i != xmax)) return +2.0*sn.y;
//        if (sn.i == xmin && (sn.i != xmin || sn.i != xmax)) return -2.0*sn.x;
//        if (sn.i == xmax && (sn.i != xmin || sn.i != xmax)) return +2.0*sn.x;
//    }
//#endif

//    throw std::exception();
//}

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
            L2Norm += 0.01*k*(u[n]-::u_fx(this, sn, tn))*(u[n]-::u_fx(this, sn, tn));

            EuNorm += (u[n]-u_fx(this, sn, tn))*(u[n]-u_fx(this, sn, tn));

            if (L1Norm < fabs(u[n]-::u_fx(this, sn, tn))) L1Norm = fabs(u[n]-::u_fx(this, sn, tn));
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
                norm += 0.01*0.01*k1*k2*(u[j][i]-::u_fx(this,sn, tn))*(u[j][i]-::u_fx(this,sn, tn));

                if (max < fabs(u[j][i]-::u_fx(this, sn, tn))) max = fabs(u[j][i]-::u_fx(this, sn, tn));
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

    HeatEquationFBVP h;
    h.setTimeDimension(Dimension(0.01, 0, 100));
    h.setSpaceDimensionX(Dimension(0.010, 100, 200));
    h.setSpaceDimensionY(Dimension(0.005, 400, 600));

    h.setThermalDiffusivity(-1.2);
    h.setThermalConductivity(-0.6);
    h.setThermalConvection(-0.8);

    Benchmark bm;
    bm.tick();
    h.implicit_calculate_D1V1();
    //h.implicit_calculate_D2V1();
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

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a1 = thermalDiffusivity();
    const double a2 = thermalDiffusivity();
    const double b1 = thermalConductivity();
    const double b2 = thermalConductivity();
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
    return 2.0*tn.t - 6.0*a1*sn.x - 3.0*b1*sn.x*sn.x - c*U(sn,tn);
#elif defined( x3_t3 )
    return 3.0*tn.t*tn.t + 6.0*td*sn.x - tc*U(sn,tn) + 3.0*tv*sn.x*sn.x;
#endif

#if defined (x1_y1_t1) || defined( x2_y2_t2 ) || defined( x2_y2_t1 )
    return Ut(sn,tn) - (a1*Uxx(sn,tn) + a2*Uyy(sn,tn) + b1*Ux(sn,tn) + b2*Uy(sn,tn) + c*U(sn,tn));
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

#if defined( x1_y1_t1 )
    return sn.x + sn.y + tn.t;
#endif
#if defined( x2_y2_t1 )
    return sn.x*sn.x + sn.y*sn.y + tn.t;
#endif
#if defined( x2_y2_t2 )
    return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
#endif

    return 0.0;
}

double HeatEquationFBVP::Ut(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
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

    return 0.0;
}

double HeatEquationFBVP::Uxx(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
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

    return 0.0;
}

double HeatEquationFBVP::Uyy(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
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

    return 0.0;
}

double HeatEquationFBVP::Ux(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
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

    return 0.0;
}

double HeatEquationFBVP::Uy(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
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

    return 0.0;
}

double HeatEquationFBVP::Un(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
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
        if (sn.j == ymin) return -1.0;
        if (sn.j == ymax) return +1.0;
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
        if (sn.i == xmin && (sn.j != ymin && sn.j != ymax)) return -2.0*sn.x;
        if (sn.i == xmax && (sn.j != ymin && sn.j != ymax)) return +2.0*sn.x;
        if (sn.j == ymin && (sn.j != ymin || sn.j != ymax)) return -2.0*sn.y;
        if (sn.j == ymax && (sn.j != ymin || sn.j != ymax)) return +2.0*sn.y;
    }
    if (tn.i%2==0)
    {
        if (sn.j == ymin && (sn.i != xmin && sn.i != xmax)) return -2.0*sn.y;
        if (sn.j == ymax && (sn.i != xmin && sn.i != xmax)) return +2.0*sn.y;
        if (sn.i == xmin && (sn.i != xmin || sn.i != xmax)) return -2.0*sn.x;
        if (sn.i == xmax && (sn.i != xmin || sn.i != xmax)) return +2.0*sn.x;
    }
#endif

    throw std::exception();
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

    //if (/*tn.i==0 || tn.i==1 || tn.i==2 || tn.i==3 || tn.i==4 || tn.i==5 || tn.i==6 ||*/ tn.i==198 || tn.i==199 || tn.i==200)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
        if (tn.i%2==0) IPrinter::printSeperatorLine();
    }
    return;
}

double u_fx(const IParabolicIBVP *p, const SpaceNodePDE &sn, const TimeNodePDE &tn, int dt, int dx, int dy)
{
    double res = 0.0;

#if defined(__QUADRATIC__)

#if defined(x1)
    if (dx == 0) res += sn.x;
    if (dx == 1) res += 1.0;
    if (dx == 2) res += 0.0;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2 == 1)
        {
            if (sn.i == xmin && (sn.j != ymin && sn.j != ymax)) res += -1.0;
            if (sn.i == xmax && (sn.j != ymin && sn.j != ymax)) res += +1.0;
            if (sn.j == ymin && (sn.j != ymin || sn.j != ymax)) res += -1.0;
            if (sn.j == ymax && (sn.j != ymin || sn.j != ymax)) res += +1.0;
        }
        if (tn.i%2==0)
        {
            if (sn.j == ymin && (sn.i != xmin && sn.i != xmax)) res += -1.0;
            if (sn.j == ymax && (sn.i != xmin && sn.i != xmax)) res += +1.0;
            if (sn.i == xmin && (sn.i != xmin || sn.i != xmax)) res += -1.0;
            if (sn.i == xmax && (sn.i != xmin || sn.i != xmax)) res += +1.0;
        }
    }
#endif

#if defined(x2)
    if (x)
    {
        if (derivativeX == 0) res += sn.x*sn.x;
        if (derivativeX == 1) res += 2.0*sn.x;
        if (derivativeX == 2) res += 2.0;
        if (derivativeX == 3)
        {
            const int xmin = p->spaceDimensionX().min();
            const int xmax = p->spaceDimensionX().max();
            const int ymin = p->spaceDimensionY().min();
            const int ymax = p->spaceDimensionY().max();

            if (tn.i%2 == 1)
            {
                if (sn.i == xmin && (sn.j != ymin && sn.j != ymax)) return -2.0*sn.x;
                if (sn.i == xmax && (sn.j != ymin && sn.j != ymax)) return +2.0*sn.x;
                if (sn.j == ymin && (sn.j != ymin || sn.j != ymax)) return -2.0*sn.y;
                if (sn.j == ymax && (sn.j != ymin || sn.j != ymax)) return +2.0*sn.y;
            }
            if (tn.i%2==0)
            {
                if (sn.j == ymin && (sn.i != xmin && sn.i != xmax)) return -2.0*sn.y;
                if (sn.j == ymax && (sn.i != xmin && sn.i != xmax)) return +2.0*sn.y;
                if (sn.i == xmin && (sn.i != xmin || sn.i != xmax)) return -2.0*sn.x;
                if (sn.i == xmax && (sn.i != xmin || sn.i != xmax)) return +2.0*sn.x;
            }
        }
    }
#endif

#if defined(x3)
    if (x)
    {
        if (derivativeX == 0) res += sn.x*sn.x*sn.x;
        if (derivativeX == 1) res += 3.0*sn.x*sn.x;
        if (derivativeX == 2) res += 6.0*sn.x;
    }
#endif

#if defined(y1)
    if (dy == 0) res += sn.y;
    if (dt == 1) res += 1.0;
    if (dy == 2) res += 0.0;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2 == 1)
        {
            if (sn.i == xmin && (sn.j != ymin && sn.j != ymax)) res += -1.0;
            if (sn.i == xmax && (sn.j != ymin && sn.j != ymax)) res += +1.0;
            if (sn.j == ymin && (sn.j != ymin || sn.j != ymax)) res += -1.0;
            if (sn.j == ymax && (sn.j != ymin || sn.j != ymax)) res += +1.0;
        }
        if (tn.i%2==0)
        {
            if (sn.j == ymin && (sn.i != xmin && sn.i != xmax)) res += -1.0;
            if (sn.j == ymax && (sn.i != xmin && sn.i != xmax)) res += +1.0;
            if (sn.i == xmin && (sn.i != xmin || sn.i != xmax)) res += -1.0;
            if (sn.i == xmax && (sn.i != xmin || sn.i != xmax)) res += +1.0;
        }
    }
#endif
#if defined(y2)
    if (y)
    {
        if (derivativeY == 0) res += sn.y*sn.y;
        if (derivativeY == 1) res += 2.0*sn.y;
        if (derivativeY == 2) res += 2.0;
    }
#endif
#if defined(y3)
    if (y)
    {
        if (derivativeY == 0) res += sn.y*sn.y*sn.y;
        if (derivativeY == 1) res += 3.0*sn.y*sn.y;
        if (derivativeY == 2) res += 6.0*sn.y;
    }
#endif
#if defined(t1)
    if (dt == 0) res += tn.t;
    if (dt == 1) res += 1.0;
    if (dt == 2) res += 0.0;
#endif
#if defined(t2)
    if (time)
    {
        if (derivativeT == 0) res += tn.t*tn.t;
        if (derivativeT == 1) res += 2.0*tn.t;
        if (derivativeT == 2) res += 2.0;
    }
#endif
#if defined(t3)
    if (time)
    {
        if (derivativeT == 0) res += tn.t*tn.t*tn.t;
        if (derivativeT == 1) res += 3.0*tn.t*tn.t;
        if (derivativeT == 2) res += 6.0*tn.t;
    }
#endif
#endif

    return res;
}


