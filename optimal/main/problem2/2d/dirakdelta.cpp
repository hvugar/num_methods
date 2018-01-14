#include "dirakdelta.h"
#include <stdio.h>

void DirakDelta::payla(const SpaceNodePDE &point, const Dimension &dimension, std::vector<SpaceNodePDE> &nodes)
{
    unsigned int rx = (unsigned int)(round(point.x*dimension.sizeN()));
    printf("%d\n", rx);

    double hx = dimension.step();
    double sigmaX = 3.0*hx;


    double sigma = 0.0;
//    sigma += 0.5*exp(-((((rx-2)*hx-point.x)*((rx-2)*hx-point.x))/(2.0*sigmaX*sigmaX)));
//    for (unsigned int i=rx-1; i<=rx+1; i++)
//    {
//        sigma += exp(-(((i*hx-point.x)*(i*hx-point.x))/(2.0*sigmaX*sigmaX)));
//    }
//    sigma += 0.5*exp(-((((rx-2)*hx-point.x)*((rx-2)*hx-point.x))/(2.0*sigmaX*sigmaX)));
//    sigma *= (1.0/(sqrt(2.0*M_PI)))*hx;

    for (unsigned int i=rx-2; i<=rx+2; i++)
    {
        sigma += exp(-(((i*hx-point.x)*(i*hx-point.x))/(2.0*sigmaX*sigmaX)));
    }
    sigma *= (1.0/(sqrt(2.0*M_PI)))*hx;

    double factor = (1.0/(sqrt(2.0*M_PI)*sigma));
    for (unsigned int i=rx-2; i<=rx+2; i++)
    {
        SpaceNodePDE node;
        node.i = i;
        node.x = i*hx;
        node.y = factor*exp(-(((node.x-point.x)*(node.x-point.x))/(2.0*sigmaX*sigmaX)));
        printf("-- %d %f %f\n", node.i, node.x, node.y);
        nodes.push_back(node);
    }
}

