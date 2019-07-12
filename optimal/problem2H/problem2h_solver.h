#ifndef PROBLEM2H_SOLVER_H
#define PROBLEM2H_SOLVER_H

#include "problem2h_global.h"
#include <stdio.h>
#include <grid/hibvp.h>
#include <deltagrid.h>

namespace Problem2H {

};

struct InitialPulse
{
    SpacePoint point;
    double blow;
    //DeltaGrid2D pointDeltaGrid;
};

struct PROBLEM2HSHARED_EXPORT SpacePointInfo
{
    std::vector<double> vl;
    std::vector<double> dx;
    std::vector<double> dy;
    unsigned int length;
};

class PROBLEM2HSHARED_EXPORT Problem2HCommon
{
public:
    virtual ~Problem2HCommon();

    virtual void setParameterCounts(unsigned int Nc, unsigned No, const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension);
    virtual void setOptimizedParameters(const DoubleMatrix &k, const DoubleMatrix &z, const std::vector<SpacePoint> &eta, const std::vector<SpacePoint> &ksi);
    //virtual void setOptimizedParameters(const DoubleVector &x);
    virtual void distributeControlDeltaGrid();
    virtual void distributeMeasurementDeltaGrid();

    unsigned int Nt;
    std::vector<TimeNodePDE> times;

protected:
    InitialPulse *initialPulses;
    unsigned int initialPulsesCount;

    unsigned int Nc;
    unsigned int No;
    std::vector<SpacePoint> eta;
    std::vector<SpacePoint> ksi;
    DoubleMatrix k;
    DoubleMatrix z;

    DeltaGrid2D *_deltaGridControl;
    SpacePointInfo *p_info;
    DeltaGrid2D *_deltaGridMeasurement;
    SpacePointInfo *u_info;

    std::vector<DoubleMatrix> u_list;
};

class PROBLEM2HSHARED_EXPORT Problem2HWaveEquationIBVP : virtual public IWaveEquationIBVP, virtual public Problem2HCommon
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;
public:
    virtual void setSpaceDimensions(const Dimension &dimensionX, const Dimension &dimensionY);
    void setInitialConditionMatrix(InitialPulse *initialPulses, unsigned int initialPulsesCount);
    void clrInitialConditionMatrix();
private:
    void layerInfoPrepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn);
    void layerInfoSave2TextFile(const DoubleMatrix &, const TimeNodePDE &) const;
    DoubleMatrix f_initialMatrix;
    DoubleMatrix f_crLayerMatrix;
};

class PROBLEM2HSHARED_EXPORT Problem2HConjugateWaveEquationIBVP : virtual public IConjugateWaveEquationIBVP, virtual public Problem2HCommon
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;
public:
    virtual void setSpaceDimensions(const Dimension &dimensionX, const Dimension &dimensionY);
private:
    void layerInfoPrepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn);
    void layerInfoSave2TextFile(const DoubleMatrix &u, const TimeNodePDE & tn) const;
    DoubleMatrix b_crLayerMatrix;
};

class PROBLEM2HSHARED_EXPORT Problem2HSolver : virtual protected Problem2HWaveEquationIBVP,
                                               virtual protected Problem2HConjugateWaveEquationIBVP,
                                               virtual protected RnFunction
{
public:
    static void Main(int argc, char* argv[]);
    static void checkGradient3(const Problem2HSolver &prob);

    void setDimensions(const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension);
    void setEquationParameters(double waveSpeed, double waveDissipation);

    Problem2HWaveEquationIBVP& fw()   { return *(dynamic_cast<Problem2HWaveEquationIBVP*>(this)); }
    Problem2HConjugateWaveEquationIBVP& bw()  { return *(dynamic_cast<Problem2HConjugateWaveEquationIBVP*>(this)); }

    void OptimalParameterFromVector(const DoubleVector &x);
    void OptimalParameterToVector(DoubleVector &x) const;

protected:
    virtual double fx(const DoubleVector &x) const;
    virtual void gradient(const DoubleVector &, DoubleVector &) const;

    double integral(const std::vector<DoubleMatrix> &vu) const;
    double integralU(const DoubleMatrix &u) const;

    unsigned int L;
    unsigned int D;

private:

};

#endif // PROBLEM2H_SOLVER1_H
