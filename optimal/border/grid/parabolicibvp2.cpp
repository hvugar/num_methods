#include "parabolicibvp2.h"

void ParabolicIBVP2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ParabolicIBVP2 p;
    p.setTimeDimension(Dimension(0.1, 0, 10));
    p.addSpaceDimension(Dimension(0.1, 0, 10));
    p.addSpaceDimension(Dimension(0.1, 0, 10));

    {
        Dimension time = p.timeDimension();
        Dimension dimx = p.spaceDimension(Dimension::DimensionX);
        Dimension dimy = p.spaceDimension(Dimension::DimensionY);

        double ht = time.step();
        double hx = dimx.step();
        double hy = dimy.step();

        unsigned int M  = time.sizeN();
        unsigned int Nx = dimx.sizeN();
        unsigned int Ny = dimy.sizeN();

        DoubleMatrix u(Ny+1, Nx+1);
        clock_t t = clock();
        TimeNodePDE tn; tn.i = 10; tn.t = 1.0;
        for (unsigned int n2=0; n2<=Ny; n2++)
        {
            for (unsigned int n1=0; n1<=Nx; n1++)
            {
                SpaceNodePDE sn; sn.i = n1; sn.x = n1*hx; sn.j = n2; sn.y = n2*hy;
                u[n2][n1] = p.U(sn,tn);
            }
        }
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }

    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateMVD(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    return;
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateMVD1(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateMVD2(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
}

double ParabolicIBVP2::initial(const SpaceNodePDE &sn) const
{
    C_UNUSED(sn);

#ifdef SAMPLE_0
    return 0.0;
#endif

#ifdef SAMPLE_1
    return sn.x*sn.x + sn.y*sn.y;
#endif

#ifdef SAMPLE_2
    return 0.0;
#endif

#ifdef SAMPLE_3
    return sn.x*sn.x * sn.y*sn.y;
#endif

#ifdef SAMPLE_4
    return sn.x*sn.x * sn.y*sn.y;
#endif

#ifdef SAMPLE_5
    return sn.x * sn.y;
#endif

    return NAN;
}

double ParabolicIBVP2::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType) const
{

#ifdef SAMPLE_0
    return (sn.x*sn.x + sn.y*sn.y)*tn.t;
#endif

#ifdef SAMPLE_1
    return sn.x*sn.x + sn.y*sn.y + tn.t;
#endif

#ifdef SAMPLE_2
    return sn.x*sn.x * sn.y*sn.y * tn.t;
#endif

#ifdef SAMPLE_3
    return sn.x*sn.x * sn.y*sn.y;
#endif

#ifdef SAMPLE_4
    return sn.x*sn.x * sn.y*sn.y + tn.t;
#endif

#ifdef SAMPLE_5
    return sn.x * sn.y + tn.t;
#endif

    return NAN;
}

double ParabolicIBVP2::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#ifdef SAMPLE_0
    return (sn.x*sn.x + sn.y*sn.y) - 4.0*a(sn, tn)*tn.t;
#endif

#ifdef SAMPLE_1
    return 1.0 - 4.0*a(sn,tn);
#endif

#ifdef SAMPLE_2
    return sn.x*sn.x * sn.y*sn.y - 2.0*a(sn,tn)*tn.t*(sn.x*sn.x + sn.y*sn.y);
#endif

#ifdef SAMPLE_3
    return -2.0*a(sn,tn)*(sn.x*sn.x + sn.y*sn.y);
#endif

#ifdef SAMPLE_4
    return sn.x*sn.x * sn.y*sn.y - 2.0*a(sn,tn)*(sn.x*sn.x + sn.y*sn.y);
#endif

#ifdef SAMPLE_5
    return 1.0;
#endif

    return NAN;
}

double ParabolicIBVP2::a(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 1.0;
}

double ParabolicIBVP2::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#ifdef SAMPLE_0
    return (sn.x*sn.x + sn.y*sn.y) * tn.t;
#endif

#ifdef SAMPLE_1
    return (sn.x*sn.x + sn.y*sn.y) + tn.t;
#endif

#ifdef SAMPLE_2
    return sn.x*sn.x * sn.y*sn.y * tn.t;
#endif

#ifdef SAMPLE_3
    return sn.x*sn.x * sn.y*sn.y;
#endif

#ifdef SAMPLE_4
    return sn.x*sn.x * sn.y*sn.y + tn.t;
#endif

#ifdef SAMPLE_5
    return sn.x * sn.y + tn.t;
#endif

    return NAN;
}
