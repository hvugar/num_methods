#include "heat_equation_ibvp.h"

#define DIMENSION_1
#define QUADRATIC

#define ROBIN

#if defined(QUADRATIC)
#define x2
#define y2
#define t1
#else

#endif

void HeatEquationIBVP::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationIBVP h;
    h.setTimeDimension(Dimension(0.000025, 0, 40000));
    h.setSpaceDimensionX(Dimension(0.010, 100, 200));
#ifdef __DIMENSION_2__
    h.setSpaceDimensionY(Dimension(0.005, 400, 600));
#endif

    //h.setThermalDiffusivity(1.2);
    h.setThermalConductivity(-0.8);
    //h.setThermalConvection(-0.6);

    Benchmark bm;
    bm.tick();
#ifdef DIMENSION_1
    //h.implicit_calculate_D1V1();
    h.explicit_calculate_D1V1();
#endif
#ifdef DIMENSION_2
    h.implicit_calculate_D2V1();
#endif
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
#if defined(DIMENSION_1)
#if defined(DIRICHLET)
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +2.0, +0.0, +1.0);
    return ::u_fx(this, sn, tn)*(condition.alpha()/condition.gamma());
#endif
#if defined(ROBIN)
    condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
    return (condition.alpha()*::u_fx(this, sn, tn)+condition.beta()*::u_fx(this, sn, tn, -1, 1, -1))/condition.gamma();
#endif
#endif

#if defined(DIMENSION_2)
    #if defined(DIRICHLET)
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +2.0, +0.0, +1.0);
    return ::u_fx(this, sn, tn)*(condition.alpha()/condition.gamma());
#endif
#if defined(ROBIN)
    condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
    return (condition.alpha()*::u_fx(this, sn, tn)+condition.beta()*::u_fx(this, sn, tn, -1, 3, 3))/condition.gamma();
#endif
#endif
}

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined(DIMENSION_1)
    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();
    return ::u_fx(this,sn,tn,+1,-1,-1)-::u_fx(this,sn,tn,-1,+2,-1)*a-::u_fx(this,sn,tn,-1,+1,-1)*b-::u_fx(this,sn,tn,0,0,0)*c;
#endif

#if defined(DIMENSION_2)
    const double a1 = thermalDiffusivity();
    const double a2 = thermalDiffusivity();
    const double b1 = thermalConductivity();
    const double b2 = thermalConductivity();
    const double c = thermalConvection();
    return ::u_fx(this, sn,tn,+1,-1,-1) - ::u_fx(this, sn,tn,-1,+2,-1)*a1 - ::u_fx(this, sn,tn,-1,-1,+2)*a2 - ::u_fx(this, sn,tn,-1,+1,-1)*b1 - ::u_fx(this, sn,tn,-1,-1,+1)*b2 - ::u_fx(this, sn,tn,+0,+0,+0)*c;
#endif
}

void HeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
    //IPrinter::printVector(u);
    if (tn.i==0 || tn.i==1 || tn.i==timeDimension().max()-1 || tn.i==timeDimension().max()) IPrinter::printVector(u);
    return;

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

    if (tn.i==0 || tn.i==1 || tn.i==2 || tn.i==198 || tn.i==199 || tn.i==200)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
        if (tn.i%2==0) IPrinter::printSeperatorLine();
    }
    return;

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
    h.setTimeDimension(Dimension(0.000025, 0, 40000));
    h.setSpaceDimensionX(Dimension(0.010, 100, 200));
#ifdef DIMENSION_2
    h.setSpaceDimensionY(Dimension(0.005, 400, 600));
#endif

    h.setThermalDiffusivity(-1.2);
    //h.setThermalConductivity(+0.8);
    //h.setThermalConvection(+0.6);

    Benchmark bm;
    bm.tick();
#ifdef DIMENSION_1
    //h.implicit_calculate_D1V1();
    h.explicit_calculate_D1V1();
#endif
#ifdef DIMENSION_2
    h.implicit_calculate_D2V1();
#endif
    bm.tock();
    bm.printDuration();
}

HeatEquationFBVP::HeatEquationFBVP() : IHeatEquationFBVP()
{}

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition) const
{
    TimeNodePDE tn; tn.t = 1.0;
    return ::p_fx(this, sn, tn);
}

double HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
#if defined(DIMENSION_1)
#if defined(DIRICHLET)
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +2.0, +0.0, +1.0);
    return ::p_fx(this, sn, tn)*(condition.alpha()/condition.gamma());
#endif
#if defined(ROBIN)
    condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
    return (condition.alpha()*::p_fx(this, sn, tn)+condition.beta()*::p_fx(this, sn, tn, -1, 1, -1))/condition.gamma();
#endif
#endif

#if defined(DIMENSION_2)
#if defined(DIRICHLET)
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +2.0, +0.0, +1.0);
    return ::p_fx(this, sn, tn)*(condition.alpha()/condition.gamma());
#endif
#if defined(ROBIN)
    condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
    return (condition.alpha()*::p_fx(this, sn, tn)+condition.beta()*::p_fx(this, sn, tn, -1, 3, 3))/condition.gamma();
#endif
#endif
}

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined(DIMENSION_1)
    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();
    return ::p_fx(this,sn,tn,+1,-1,-1)-::p_fx(this,sn,tn,-1,+2,-1)*a-::p_fx(this,sn,tn,-1,+1,-1)*b-::p_fx(this,sn,tn,0,0,0)*c;
#endif

#if defined(DIMENSION_2)
    const double a1 = thermalDiffusivity();
    const double a2 = thermalDiffusivity();
    const double b1 = thermalConductivity();
    const double b2 = thermalConductivity();
    const double c = thermalConvection();
    return ::p_fx(this, sn,tn,+1,-1,-1) - ::p_fx(this, sn,tn,-1,+2,-1)*a1 - ::p_fx(this, sn,tn,-1,-1,+2)*a2 - ::p_fx(this, sn,tn,-1,+1,-1)*b1 - ::p_fx(this, sn,tn,-1,-1,+1)*b2 - ::p_fx(this, sn,tn,+0,+0,+0)*c;
#endif
}

void HeatEquationFBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);;
    //IPrinter::printVector(u);
    if (tn.i==0 || tn.i==1 || tn.i==timeDimension().max()-1 || tn.i==timeDimension().max()) IPrinter::printVector(u);
}

void HeatEquationFBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    if (tn.i==0 || tn.i==1 || tn.i==2 || tn.i==198 || tn.i==199 || tn.i==200)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
        if (tn.i%2==0) IPrinter::printSeperatorLine();
    }
}

double u_fx(const IParabolicIBVP *p, const SpaceNodePDE &sn, const TimeNodePDE &tn, int dt, int dx, int dy)
{
    double res = 0.0;

#if defined(QUADRATIC)

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

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -1.0;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +1.0;
        }
        else
        {
            if (sn.i == xmin) res += -1.0;
            if (sn.i == xmax) res += +1.0;
        }
    }
#endif

#if defined(x2)
    if (dx == 0) res += sn.x*sn.x;
    if (dx == 1) res += 2.0*sn.x;
    if (dx == 2) res += 2.0;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -2.0*sn.x;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +2.0*sn.x;
        }
        else
        {
            if (sn.i == xmin) res += -2.0*sn.x;
            if (sn.i == xmax) res += +2.0*sn.x;
        }
    }
#endif

#if defined(x3)
    if (dx == 0) res += sn.x*sn.x*sn.x;
    if (dx == 1) res += 3.0*sn.x*sn.x;
    if (dx == 2) res += 6.0*sn.x;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -3.0*sn.x*sn.x;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +3.0*sn.x*sn.x;
        }
        else
        {
            if (sn.i == xmin) res += -3.0*sn.x*sn.x;
            if (sn.i == xmax) res += +3.0*sn.x*sn.x;
        }
    }
#endif

#if defined(y1)
    if (dy == 0) res += sn.y;
    if (dy == 1) res += 1.0;
    if (dy == 2) res += 0.0;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -1.0;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +1.0;
        }
        else
        {
            if (sn.j == ymin) res += -1.0;
            if (sn.j == ymax) res += +1.0;
        }
    }
#endif
#if defined(y2)
    if (dy == 0) res += sn.y*sn.y;
    if (dy == 1) res += 2.0*sn.y;
    if (dy == 2) res += 2.0;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -2.0*sn.y;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +2.0*sn.y;
        }
        else
        {
            if (sn.j == ymin) res += -2.0*sn.y;
            if (sn.j == ymax) res += +2.0*sn.y;
        }
    }
#endif
#if defined(y3)
    if (dy == 0) res += sn.y*sn.y*sn.y;
    if (dy == 1) res += 3.0*sn.y*sn.y;
    if (dy == 2) res += 6.0*sn.y;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -3.0*sn.y*sn.y;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +3.0*sn.y*sn.y;
        }
        else
        {
            if (sn.j == ymin) res += -3.0*sn.y*sn.y;
            if (sn.j == ymax) res += +3.0*sn.y*sn.y;
        }
    }
#endif
#if defined(t1)
    if (dt == 0) res += tn.t;
    if (dt == 1) res += 1.0;
    if (dt == 2) res += 0.0;
#endif
#if defined(t2)
    if (dt == 0) res += tn.t*tn.t;
    if (dt == 1) res += 2.0*tn.t;
    if (dt == 2) res += 2.0;
#endif
#if defined(t3)
    if (dt == 0) res += tn.t*tn.t*tn.t;
    if (dt == 1) res += 3.0*tn.t*tn.t;
    if (dt == 2) res += 6.0*tn.t;
#endif
#endif

    return res;
}

double p_fx(const IParabolicFBVP *p, const SpaceNodePDE &sn, const TimeNodePDE &tn, int dt, int dx, int dy)
{
    double res = 0.0;

#if defined(QUADRATIC)

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

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -1.0;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +1.0;
        }
        else
        {
            if (sn.i == xmin) res += -1.0;
            if (sn.i == xmax) res += +1.0;
        }
    }
#endif

#if defined(x2)
    if (dx == 0) res += sn.x*sn.x;
    if (dx == 1) res += 2.0*sn.x;
    if (dx == 2) res += 2.0;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -2.0*sn.x;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +2.0*sn.x;
        }
        else
        {
            if (sn.i == xmin) res += -2.0*sn.x;
            if (sn.i == xmax) res += +2.0*sn.x;
        }
    }
#endif

#if defined(x3)
    if (dx == 0) res += sn.x*sn.x*sn.x;
    if (dx == 1) res += 3.0*sn.x*sn.x;
    if (dx == 2) res += 6.0*sn.x;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -3.0*sn.x*sn.x;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +3.0*sn.x*sn.x;
        }
        else
        {
            if (sn.i == xmin) res += -3.0*sn.x*sn.x;
            if (sn.i == xmax) res += +3.0*sn.x*sn.x;
        }
    }
#endif

#if defined(y1)
    if (dy == 0) res += sn.y;
    if (dy == 1) res += 1.0;
    if (dy == 2) res += 0.0;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -1.0;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +1.0;
        }
        else
        {
            if (sn.j == ymin) res += -1.0;
            if (sn.j == ymax) res += +1.0;
        }
    }
#endif
#if defined(y2)
    if (dy == 0) res += sn.y*sn.y;
    if (dy == 1) res += 2.0*sn.y;
    if (dy == 2) res += 2.0;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -2.0*sn.y;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +2.0*sn.y;
        }
        else
        {
            if (sn.j == ymin) res += -2.0*sn.y;
            if (sn.j == ymax) res += +2.0*sn.y;
        }
    }
#endif
#if defined(y3)
    if (dy == 0) res += sn.y*sn.y*sn.y;
    if (dy == 1) res += 3.0*sn.y*sn.y;
    if (dy == 2) res += 6.0*sn.y;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -3.0*sn.y*sn.y;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +3.0*sn.y*sn.y;
        }
        else
        {
            if (sn.j == ymin) res += -3.0*sn.y*sn.y;
            if (sn.j == ymax) res += +3.0*sn.y*sn.y;
        }
    }
#endif
#if defined(t1)
    if (dt == 0) res += tn.t;
    if (dt == 1) res += 1.0;
    if (dt == 2) res += 0.0;
#endif
#if defined(t2)
    if (dt == 0) res += tn.t*tn.t;
    if (dt == 1) res += 2.0*tn.t;
    if (dt == 2) res += 2.0;
#endif
#if defined(t3)
    if (dt == 0) res += tn.t*tn.t*tn.t;
    if (dt == 1) res += 3.0*tn.t*tn.t;
    if (dt == 2) res += 6.0*tn.t;
#endif
#endif

    return res;
}
