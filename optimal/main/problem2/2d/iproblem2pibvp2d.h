#ifndef PROBLEM22DIPARABOLICIBVP_H
#define PROBLEM22DIPARABOLICIBVP_H

#include <grid/pibvp.h>
#include <vector>
#include <printer.h>
#include <time.h>

using namespace std;

struct ObservationNode
{
    SpaceNodePDE xi;
    unsigned int j;
    unsigned int n;
    unsigned int m;
    double x;
    double y;
    double w;
};

struct ControlNode
{
    SpaceNodePDE eta;
    unsigned int i;
    unsigned int n;
    unsigned int m;
    double x;
    double y;
    double w;
};

struct ControlDeltaNode
{
    SpaceNodePDE eta;
    unsigned int i;
    unsigned int n;
    unsigned int m;
    double x;
    double y;
    double w;
};

struct ObservationDeltaNode
{
    SpaceNodePDE xi;
    unsigned int j;
    unsigned int n;
    unsigned int m;
    double x;
    double y;
    double w;
};


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

    double value(unsigned int layer) const;
    double valueDx(unsigned int layer) const;
    double valueDy(unsigned int layer) const;

    //double valueDx(double x, double y, unsigned int layer) const;
    //double valueDy(double x, double y, unsigned int layer) const;

    //double valueDxN(unsigned int layer, double h) const;
    //double valueDyN(unsigned int layer, double h) const;

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

    void setEquationParameters(double a, double lambda0, double lambda, double theta);

    void setParamter(const Parameter &parameter);
    const Parameter& parameter() const;

protected:
    double a;
    double lambda0;
    double lambda;
    double theta;
    Parameter mParameter;
};

#endif // PROBLEM22DIPARABOLICIBVP_H
