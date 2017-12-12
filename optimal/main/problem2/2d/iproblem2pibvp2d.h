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

class ExtendedSpaceNode2D : public SpaceNodePDE
{
public:
    explicit ExtendedSpaceNode2D();

    void setSpaceNode(const SpaceNodePDE& sn);

    void extendWeights(const Dimension &dimX, const Dimension &dimY, unsigned int rows = 4, unsigned int cols = 4);
    void clearWeights();

    void extendLayers(unsigned int layerNumber);
    void clearLayers();

    double value(double x, double y, unsigned int layer) const;

    unsigned int id;
    unsigned int rows;
    unsigned int cols;
    WISpaceNodePDE **wi;
    unsigned int layerNumber;
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
