#ifndef PROBLEM22DIPARABOLICIBVP_H
#define PROBLEM22DIPARABOLICIBVP_H

#include <grid/pibvp.h>
#include <vector>
#include <printer.h>
#include <time.h>

using namespace std;

struct P2Setting
{
    unsigned int Lc;
    unsigned int Lo;

    DoubleMatrix k;
    DoubleMatrix z;

    vector<SpaceNodePDE> xi;
    vector<SpaceNodePDE> eta;

    void toVector(DoubleVector &prms) const;
    void fromVector(const DoubleVector &prms);
};

struct WISpaceNodePDE : public SpaceNodePDE
{
    double w;
    double *u;
};

struct ExtendedSpaceNode
{
    unsigned int id;
    SpaceNodePDE xi;
    WISpaceNodePDE **wi;
    unsigned int rows;
    unsigned int cols;
    void init(unsigned int rows, unsigned int cols, unsigned int M);
    void clear();
};

class IProblem22DPIBVP : public IParabolicIBVP
{
public:
    //    struct ExSpaceNodePDE : public SpaceNodePDE
    //    {
    //        double w;
    //        unsigned int index;

    //        double *timeValues;
    //        unsigned int timeN;
    //    };

    //    struct ControlPoint : public SpaceNodePDE
    //    {
    //        vector<ExSpaceNodePDE> extNodes;
    //    };

    //    struct ObservationPoint : public SpaceNodePDE
    //    {
    //        vector<ExSpaceNodePDE> extNodes;
    //    };

    //    struct Parameters
    //    {
    //        unsigned int Lc;
    //        unsigned int Lo;
    //        DoubleMatrix k;
    //        DoubleMatrix z;
    //        vector<ControlPoint> eta;
    //        vector<ObservationPoint> xi;

    //        void toVector(DoubleVector &prms) const;
    //        void fromVector(const DoubleVector &prms, unsigned int Lc, unsigned int Lo);
    //    };

public:
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    //    void setParametr(const Parameters& p);
    //    const Parameters& parametrs() const;

    //    void extendDeltaControlPointsToGrid();
    //    void extendDeltaControlPointToGrid1(ControlPoint &cp, unsigned int index);

    //    void extendObservationPointToGrid();
    //    void extendObservationPointToGrid1(ObservationPoint &op, unsigned int index);

    double a;
    double lambda0;
    double lambda;
    double theta;

protected:
    //    Parameters mParameters;
};

#endif // PROBLEM22DIPARABOLICIBVP_H
