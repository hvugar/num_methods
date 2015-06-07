#include "prjgradient.h"

ProjectionGradient::ProjectionGradient() : FastProximalGradient(), a(-INFINITY), b(INFINITY)
{
}

ProjectionGradient::~ProjectionGradient()
{

}

void ProjectionGradient::calculate()
{
    mcount = 0;
    do
    {
        /* calculating function gradient at current point */
        calcGradient();

        /* if norm of gradinet at current point is less than epsilon break. no minimize */
        if (gradientNorm() < epsilon())
            break;

        mcount++;

        double gn = 0.0;
        for (unsigned int i=0; i<mg.size(); i++) gn = gn + mg[i]*mg[i];
        gn = sqrt(gn);
        for (unsigned int i=0; i<mg.size(); i++) mg[i] = mg[i] / gn;

        /* R1 minimization in direct of antigradient */
        malpha = minimize();

        print();

        for (unsigned int i=0; i<mx.size(); i++)
        {
            mx[i] = mx[i] - malpha * mg[i];
            if (mx[i] < a) { mx[i] = a; }
            if (mx[i] > b) { mx[i] = b; }
        }

        /* calculating distance previous and new point */
    } while (distance() > epsilon());
}

double ProjectionGradient::fx(double alpha)
{
    unsigned int n = mx.size();
    std::vector<double> x2(n);
    for (unsigned int i=0; i<n; i++)
    {
        x2[i] = mx[i] - alpha * mg[i];
        if ( x2[i] < a ) x2[i] = a;
        if ( x2[i] > b ) x2[i] = b;
    }
    return mfn->fx(x2);
}

