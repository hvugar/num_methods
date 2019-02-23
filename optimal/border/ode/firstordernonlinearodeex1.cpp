#include "firstordernonlinearodeex1.h"

double FirstOrderNonLinearODEEx1::f(double x, double y, unsigned int k) const
{
    double sigma = 0.01;
//    return  3.0*x*x;// + 4.0*(1.0/(sqrt(2.0*M_PI*sigma*sigma))) * exp(((x-0.5)*(x-0.5))/(-2.0*sigma*sigma));
//    return  4.0*x*x*x*exp(pow(x,4.0));// + 4.0*(1.0/(sqrt(2.0*M_PI*sigma*sigma))) * exp(((x-0.8)*(x-0.8))/(-2.0*sigma*sigma));

    double w = 0.0;
    if (k==80) w = 100.0;
    return  4.0*x*x*x*exp(pow(x,4.0)) + 4.0*w;
}
