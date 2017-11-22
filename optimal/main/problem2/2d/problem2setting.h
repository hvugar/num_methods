#ifndef PROBLEM2SETTING_H
#define PROBLEM2SETTING_H

#include <vector>
#include <vector2d.h>
#include <matrix2d.h>
#include <grid/grid.h>

using namespace std;

struct P2Setting
{
    double a;
    double lambda0;
    double lambda;
    double theta;

    unsigned int Lc;
    unsigned int Lo;

    DoubleMatrix k;
    DoubleMatrix z;

    vector<SpaceNodePDE> xi;
    vector<SpaceNodePDE> eta;

    void toVector(DoubleVector &prms) const;
    void fromVector(const DoubleVector &prms);
};

#endif // PROBLEM2SETTING_H
