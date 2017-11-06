#ifndef IPROBLEM2BACKWARD2D_H
#define IPROBLEM2BACKWARD2D_H

#include <grid/pibvp.h>
#include <vector>
#include <printer.h>

using namespace std;

class IProblem2Backward2D : public IParabolicIBVP
{
public:
    IProblem2Backward2D(double a, double lambda0, double lambda, double theta, unsigned int Lc, unsigned int Lo);

    void calculateMVD(DoubleMatrix &p);

    DoubleMatrix U;
    DoubleMatrix mu;
    DoubleMatrix uT;

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    double delta(const SpaceNodePDE &sn, unsigned int j, unsigned int source) const;

    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double h(const SpaceNodePDE &sn) const;

private:
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
};

#endif // IPROBLEM2BACKWARD2D_H
