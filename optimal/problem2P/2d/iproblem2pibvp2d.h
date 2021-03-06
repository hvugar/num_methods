#ifndef PROBLEM22DIPARABOLICIBVP_H
#define PROBLEM22DIPARABOLICIBVP_H

#include <grid/pibvp.h>
#include <vector>
#include <printer.h>
#include <time.h>

using namespace std;
#define USE_F_VARIANT_2
#define USE_B_VARIANT_2

#define APPROX_F1_3
#define APPROX_FD_3

#define APPROX_B1_3
#define APPROX_BD_3

//#define USE_ADDITIONAL_FUNCTIONS

struct ExtendedGridNode
{
    SpaceNodePDE point;
    unsigned int id;
    unsigned int i;
    unsigned int j;
    double x;
    double y;
    double w;
};

void extendPointToGridNodes(const SpaceNodePDE& point, int id, vector<ExtendedGridNode> nodes, const Dimension &dimensionX, const Dimension &dimensionY);

struct ExtendedDeltaPoint
{
    SpaceNodePDE pt;
    unsigned int id;
    unsigned int i;
    unsigned int j;
    double x;
    double y;
    double w;
};

struct ObservationNode : public ExtendedDeltaPoint {};

struct ObservationDeltaNode : public ExtendedDeltaPoint {};

struct ControlNode : public ExtendedDeltaPoint {};

struct ControlDeltaNode : public ExtendedDeltaPoint {};

struct Parameter
{
    Parameter(unsigned int Lc=0, unsigned int Lo=0);
    Parameter(const DoubleVector& prmtrs, unsigned int Lc, unsigned int Lo);

    unsigned int Lc;
    unsigned int Lo;

    DoubleMatrix k;
    DoubleMatrix z;

    vector<SpaceNodePDE> xi;
    vector<SpaceNodePDE> eta;

    //void toVector(DoubleVector &prms) const;
    //void fromVector(const DoubleVector &prms);
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
    virtual ~ExtendedSpaceNode2D();

    void setSpaceNode(const SpaceNodePDE& sn);

    void extendWeights(const Dimension &dimX, const Dimension &dimY, unsigned int rows = 4, unsigned int cols = 4);
    void clearWeights();

    void extendLayers(unsigned int layerNumber);
    void clearLayers();

    double value(unsigned int layer) const;
    double valueDx(unsigned int layer) const;
    double valueDy(unsigned int layer) const;

    unsigned int id;
    unsigned int rows;
    unsigned int cols;
    WISpaceNodePDE **wi;
    unsigned int layerNumber;

private:
    double value(double x, double y, unsigned int layer) const;
};

class IProblem22DPIBVP : public IParabolicIBVP
{
public:
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    void setEquationParameters(double a, double lambda0, double lambda);
    void setIntTemperature(double fi);
    void setEnvTemperature(double theta);

    void setParameter(const Parameter &parameter);
    const Parameter& parameter() const;

protected:
    double a;
    double lambda0;
    double lambda;
    double theta;
    double fi;
    Parameter mParameter;

    void distributeDelta(const SpaceNodePDE &pt, std::vector<ExtendedDeltaPoint> &nodes, unsigned int id) const;
};

#endif // PROBLEM22DIPARABOLICIBVP_H
