#include "heatequationibvp1.h"

HeatEquationIBVP1::HeatEquationIBVP1()
{}

HeatEquationIBVP1::~HeatEquationIBVP1() {}

//double HeatEquationIBVP1::initial(const SpaceNodePDE &) const
//{
//    return 0.0;
//}

double HeatEquationIBVP1::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const
{
    return 0.0;
}

//double HeatEquationIBVP1::f(const SpaceNodePDE &, const TimeNodePDE &) const
//{
//    return 0.0;
//}

void HeatEquationIBVP1::layerInfo(const DoubleVector &u, unsigned int)
{
    //    IPrinter::printVector(u);
}

//void HeatEquationIBVP1::layerInfo(const DoubleMatrix &, unsigned int)
//{ }

void HeatEquationIBVP1::gridMethod1(DoubleVector &u, double a)
{
    Dimension time = mtimeDimension;
    Dimension dimX = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int M = time.sizeN();

    double hx = dimX.step();
    unsigned int N =dimX.sizeN();

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N+1));
    double *kb = (double*) malloc(sizeof(double)*(N+1));
    double *kc = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));

    /* initial condition */
    for (unsigned int n=0; n<=N; n++) u[n] = 0.0;
    layerInfo(u, 0);

    double q = 1.0;
    for (unsigned int m=1; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            ka[n] = -(a*a*ht)/(hx*hx);
            kb[n] = 1.0 + 2.0*(a*a*ht)/(hx*hx);
            kc[n] = -(a*a*ht)/(hx*hx);
            kd[n] = u[n];

            if (m==1)
            {
                kd[n] += ht * q * (1.0/sqrt(2.0*M_PI*hx*hx)) * exp(-((n*hx-0.5)*(n*hx-0.5))/(2.0*hx*hx)) * (1.0/ht);
            }
        }

        ka[0] = 0.0;
        kc[0] = -2.0*(a*a*ht)/(hx*hx);

        ka[N] = -2.0*(a*a*ht)/(hx*hx);
        kc[N] = 0.0;

        tomasAlgorithmR2L(ka, kb, kc, kd, rx, N+1);

        for (unsigned int n=0; n<=N; n++) u[n] = rx[n];
        layerInfo(u, m);
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}
