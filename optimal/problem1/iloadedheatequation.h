#ifndef ILOADEDHEATEQUATION_H
#define ILOADEDHEATEQUATION_H

#include <grid/pibvp.h>

class ILoadedHeatEquation : public IParabolicIBVP
{
public:
    struct Parameter
    {
        double k;
        double z;
        double e;
        //unsigned int xi;
    };

    double a;
    double lambda0;
    double lambda1;
    double lambda2;
    double theta;

    //void LoadMatrixParameters(double *ka, double *b, double *c, double *d, double *e, unsigned int N, unsigned int m);

    /* Transferring conditions */
    void transferringConditions2(DoubleVector &u);
    void transferringConditions4(DoubleVector &u);

    /* Classic method */
    void calculateM2(DoubleVector &u);

    unsigned int L;
    std::vector<Parameter> params;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

public:
    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const {}
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}

    //virtual double a(const SpaceNode &sn, const TimeNode &tn) const;
    virtual double g(const TimeNodePDE &tn) const = 0;
    virtual double h(const TimeNodePDE &tn) const = 0;

private:
    void getGridParameters(double &hx, unsigned int &N, double &ht, unsigned int &M,
                           unsigned int &minN, unsigned int &maxN,
                           unsigned int &minM, unsigned int &maxM) const;

    //void qovmaFirstRowM(double *e, double *a, double *b, double *c, double *d, double *x, unsigned int N) const;
};

#endif // ILOADEDHEATEQUATION_H
