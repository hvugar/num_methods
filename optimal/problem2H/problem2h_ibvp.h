#ifndef PROBLEM2H_IBVP_H
#define PROBLEM2H_IBVP_H

#include "problem2h_common.h"

class PROBLEM2HSHARED_EXPORT Problem2HNDirichletForward1 : public CC1IHyperbolicIBVP
{
public:
    const EquationParameterH *e_prm = nullptr;
    const OptimizeParameterH *o_prm = nullptr;
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

    std::vector<DoubleMatrix> vu;
    unsigned int LD = 50;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    void setEquationParameters(const EquationParameterH &e_prm, const OptimizeParameterH &o_prm, unsigned int N, double hx, unsigned int M, double hy);

private:
    std::vector<DeltaGrid> tetaGridList;
    std::vector<DeltaGrid> msrmGridList;
    std::vector<DeltaGrid> cntrGridList;

    DoubleMatrix ixv;
    DoubleMatrix fxv;
};

//**********************************************************************************************************//

class PROBLEM2HSHARED_EXPORT Problem2HNDirichletBackward1 : public ConjugateCC1IHyperbolicIBVP
{
public:
    const EquationParameterH *e_prm = nullptr;
    const OptimizeParameterH *o_prm = nullptr;
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    void setEquationParameters(const EquationParameterH &e_prm, const OptimizeParameterH &o_prm, unsigned int N, double hx, unsigned int M, double hy);

private:
    std::vector<DeltaGrid> msrmGridList;
    std::vector<DeltaGrid> cntrGridList;

    DoubleMatrix ixv;
    DoubleMatrix fxv;

    std::vector<DoubleMatrix> vu;
    unsigned int LD = 50;
};

//**********************************************************************************************************//

#endif // PROBLEM2H_IBVP_H
