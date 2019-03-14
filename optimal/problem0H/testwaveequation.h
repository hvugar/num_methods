#ifndef TESTWAVEEQUATION_H
#define TESTWAVEEQUATION_H

#include <grid/hibvp.h>

#ifdef USE_IMAGING
#include <QPixmap>
#include <QGuiApplication>
#include <imaging.h>
#endif

class TestWaveEquation : public CdIHyperbolicIBVP
{
public:
    void static Main(int agrc, char** argv);
    double U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

protected:
    double _waveSpeed;
    double _dissipation;

    double alpha1;
    double alpha2;
    double alpha3;

    DoubleMatrix u1;
    DoubleMatrix u2;

    SpacePoint p1;

    auto saveToImage(const DoubleMatrix &U, unsigned int ln) const -> void;
};

#endif // TESTWAVEEQUATION_H
