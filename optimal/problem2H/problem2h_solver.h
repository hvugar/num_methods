#ifndef PROBLEM2H_SOLVER_H
#define PROBLEM2H_SOLVER_H

#include "problem2h_global.h"
#include <stdio.h>
#include <grid/hibvp.h>
#include <deltagrid.h>

struct SpacePointInfo : public SpacePoint
{
    std::vector<double> vl;
    std::vector<double> dx;
    std::vector<double> dy;
    DeltaGrid2D deltaGrid;
};

class Problem2HCommon
{
public:
    unsigned int Nt;
    std::vector<double> q;
    std::vector<SpacePoint> theta;
    std::vector<DeltaGrid2D> thetaDeltaGrid;
    DoubleMatrix initialMatrix;

    unsigned int Nc;
    unsigned int No;
    std::vector<SpacePointInfo> eta;
    std::vector<SpacePointInfo> ksi;
    DoubleMatrix layerMatrix;

    DoubleMatrix k;
    DoubleMatrix z;
};

class PROBLEM2HSHARED_EXPORT Problem2HForward : virtual public IWaveEquationIBVP, virtual public Problem2HCommon
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

public:
    void initInitialConditionMatrix(uint32_t Nt, const std::vector<double> &q, const std::vector<SpacePoint> &theta);
    void clearInitialConditionMatrix();

    void initControlMeasurementDeltaGrid(std::vector<SpacePoint> &eta, std::vector<SpacePoint> &ksi);
    void prepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn);


    void saveToTextF(const DoubleMatrix &, const TimeNodePDE &) const;
};

class PROBLEM2HSHARED_EXPORT Problem2HBackward : virtual public IConjugateWaveEquationIBVP, virtual public Problem2HCommon
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;
};

class PROBLEM2HSHARED_EXPORT Problem2HFunctional : virtual public Problem2HCommon {};

class PROBLEM2HSHARED_EXPORT Problem2HSolver : virtual public Problem2HForward//, virtual public Problem2HBackward
{
public:
    static void Main(int argc, char** argv);

    //Problem2HForward& forward();
    //Problem2HBackward& backward();
};

#endif // PROBLEM2H_SOLVER1_H
