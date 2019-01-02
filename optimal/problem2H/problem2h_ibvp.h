#ifndef PROBLEM2H_IBVP_H
#define PROBLEM2H_IBVP_H

#include "problem2h_common.h"

class  Problem2HNDirichletForward1;
class  Problem2HNDirichletBackward1;

class PROBLEM2HSHARED_EXPORT Problem2HNDirichletForward1 : public CdIHyperbolicIBVP
{
public:
    Problem2HNDirichletForward1();

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

    std::vector<DoubleMatrix> vu;
    spif_vectorH u_info;
    unsigned int LD = 50;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    void setParameters(const EquationParameterH &equationParameter, const OptimizeParameterH &optimizeParameter, unsigned int LD);
    void clear();

private:
    std::vector<DeltaGrid2D> tetaGridList;
    std::vector<DeltaGrid2D> msrmGridList;
    std::vector<DeltaGrid2D> cntrGridList;

    EquationParameterH mEquationParameter;
    OptimizeParameterH mOptimizeParameter;

    DoubleMatrix ixv;
    DoubleMatrix fxv;

    double noise = 0.0;

    friend class Problem2HNDirichletBackward1;
};

//**********************************************************************************************************//

class PROBLEM2HSHARED_EXPORT Problem2HNDirichletBackward1 : public ConjugateCdIHyperbolicIBVP
{
public:
    Problem2HNDirichletBackward1(const Problem2HNDirichletForward1& fw);

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    const Problem2HNDirichletForward1 &_fw;
    DoubleMatrix ixv;
    DoubleMatrix fxv;    
};

//**********************************************************************************************************//

#endif // PROBLEM2H_IBVP_H
