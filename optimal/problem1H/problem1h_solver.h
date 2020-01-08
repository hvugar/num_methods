#ifndef PROBLEM1H_SOLVER_H
#define PROBLEM1H_SOLVER_H

#include "problem1h_global.h"
#include <function.h>
#include <gradient.h>
#include <printer.h>
#include <grid/hibvp.h>

#include <iostream>
#include <functional>
#include <memory>

//#define _LEFT_BORDER_DIRICHLET
#define _LEFT_BORDER_ROBIN
#define _RIGHT_BORDER_DIRICHLET
//#define _RIGHT_BORDER_ROBIN

using namespace std;
using namespace std::placeholders;

namespace h1p
{

class ProblemSolver;

struct EquationParameters
{
    double waveSpeed = 1.0;
    double waveDissipation = 0.0;
    double unknownB = 0.0;
    double restoration = 0.0;

    double lambda1 = 1.0;
    double lambda2 = 1.0;

    double *v0;
    double *v1;
};

class PROBLEM1HSHARED_EXPORT WaveEquationIBVP : public IWaveEquationIBVP
{
public:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

    ProblemSolver *solver;

protected:
    virtual double weight() const { return 0.25; }

public:
    virtual const Dimension& timeDimension() const { return _timeDimension; }
    virtual const Dimension& spaceDimensionX() const { return _spaceDimensionX; }
    virtual const Dimension& spaceDimensionY() const { return _spaceDimensionY; }
    virtual const Dimension& spaceDimensionZ() const { return _spaceDimensionZ; }

    void setTimeDimension(const Dimension &timeDimension) { this->_timeDimension = timeDimension; }
    void setSpaceDimensionX(const Dimension &spaceDimensionX) { this->_spaceDimensionX = spaceDimensionX; }
    void setSpaceDimensionY(const Dimension &spaceDimensionY) { this->_spaceDimensionY = spaceDimensionY; }

    bool showLayers = false;

private:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;
    Dimension _spaceDimensionZ;
};

class PROBLEM1HSHARED_EXPORT WaveEquationFBVP : public IWaveEquationFBVP
{
public:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

    ProblemSolver *solver;

protected:
    virtual double weight() const { return 0.25; }

public:
    virtual const Dimension& timeDimension() const { return _timeDimension; }
    virtual const Dimension& spaceDimensionX() const { return _spaceDimensionX; }
    virtual const Dimension& spaceDimensionY() const { return _spaceDimensionY; }
    virtual const Dimension& spaceDimensionZ() const { return _spaceDimensionZ; }

    void setTimeDimension(const Dimension &timeDimension) { this->_timeDimension = timeDimension; }
    void setSpaceDimensionX(const Dimension &spaceDimensionX) { this->_spaceDimensionX = spaceDimensionX; }
    void setSpaceDimensionY(const Dimension &spaceDimensionY) { this->_spaceDimensionY = spaceDimensionY; }

    bool showLayers = false;

private:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;
    Dimension _spaceDimensionZ;
};


class PROBLEM1HSHARED_EXPORT ProblemSolver : public RnFunction, public IGradient
{
public:
    static void Main(int argc, char* argv[]);

    ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX);

    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;
    virtual double fx(const DoubleVector &x) const;

    void setTimeDimension(const Dimension &timeDimension);
    void setSpaceDimensionX(const Dimension &spaceDimensionX);

    virtual double frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void frw_layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

    virtual double bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void bcw_layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

protected:
    WaveEquationIBVP forward;
    WaveEquationFBVP backward;

    EquationParameters params;
    Dimension _timeDimension;
    Dimension _spaceDimensionX;

    DoubleVector U1;
    DoubleVector V1;
    DoubleVector U2;
    DoubleVector V2;
    DoubleVector p0;
    DoubleVector p0x;

    ProblemSolver *const_this;

    friend class WaveEquationIBVP;
    friend class WaveEquationFBVP;
};

};

#endif // PROBLEM1H_SOLVER_H
