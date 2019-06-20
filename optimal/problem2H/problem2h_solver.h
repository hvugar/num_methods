#ifndef PROBLEM2H_SOLVER_H
#define PROBLEM2H_SOLVER_H

#include "problem2h_global.h"
#include <stdio.h>
#include <grid/hibvp.h>
#include <deltagrid.h>

struct SpacePointInfo
{
    std::vector<double> vl;
    std::vector<double> dx;
    std::vector<double> dy;
    unsigned int length;
};

class Problem2HCommon
{
public:
    unsigned int Nq;
    std::vector<double> q;
    std::vector<DeltaGrid2D> zta;

    unsigned int Nc;
    unsigned int No;
    std::vector<DeltaGrid2D> eta;
    std::vector<SpacePointInfo> p_info;
    std::vector<DeltaGrid2D> ksi;
    std::vector<SpacePointInfo> u_info;

    DoubleMatrix k;
    DoubleMatrix z;

    DoubleMatrix f_initialMatrix;
    DoubleMatrix f_layerMatrix;
    DoubleMatrix b_initialMatrix;
    DoubleMatrix b_layerMatrix;

    unsigned int Nt;
    std::vector<TimeNodePDE> times;

    void initMatrixes(const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension, unsigned int Nq, unsigned int Nc, unsigned No);
    void clearMatrixes();
};

class PROBLEM2HSHARED_EXPORT Problem2HForward : virtual public IWaveEquationIBVP, virtual public Problem2HCommon
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

public:
    void initInitialConditionMatrix(const std::vector<SpacePoint> &zta, const std::vector<double> &q);
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
