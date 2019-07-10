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
    //unsigned int Nq;
    //std::vector<double> q;
    //std::vector<DeltaGrid2D> zta;

    unsigned int Nc;
    unsigned int No;
    std::vector<DeltaGrid2D> eta;
    std::vector<SpacePointInfo> p_info;
    std::vector<DeltaGrid2D> ksi;
    std::vector<SpacePointInfo> u_info;

    DoubleMatrix k;
    DoubleMatrix z;

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
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

public:
    void setInitialConditionMatrix(const SpacePoint *zta, const double* q, unsigned int Nq);
    void clrInitialConditionMatrix();

    void initControlMeasurementDeltaGrid(unsigned int Nc, unsigned int No);
    void prepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn);

    void saveToTextF(const DoubleMatrix &, const TimeNodePDE &) const;

private:
    DoubleMatrix f_initialMatrix;
    DoubleMatrix f_crLayerMatrix;

    DeltaGrid2D *_deltaGridControl;
    DeltaGrid2D *_deltaGridMeasurement;
};

class PROBLEM2HSHARED_EXPORT Problem2HBackward : virtual public IConjugateWaveEquationIBVP, virtual public Problem2HCommon
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;
};

class PROBLEM2HSHARED_EXPORT Problem2HFunctional : virtual public Problem2HCommon {};

class PROBLEM2HSHARED_EXPORT Problem2HSolver : virtual public Problem2HForward,
                                               virtual public Problem2HBackward,
                                               virtual public Problem2HFunctional
{
public:
    static void Main(int argc, char** argv);

    InitialBoundaryValueProblemPDE& ibvp() { return *(dynamic_cast<InitialBoundaryValueProblemPDE*>(this)); }
    Problem2HForward& fw()   { return *(dynamic_cast<Problem2HForward*>(this)); }
    Problem2HBackward& bw()  { return *(dynamic_cast<Problem2HBackward*>(this)); }
};

#endif // PROBLEM2H_SOLVER1_H
