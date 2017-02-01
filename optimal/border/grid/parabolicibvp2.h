#ifndef PARABOLICIBVP2_H
#define PARABOLICIBVP2_H

#include <grid/pibvp.h>
#include <time.h>

class MINIMUMSHARED_EXPORT ParabolicIBVP2 : public ParabolicIBVP
{
public:
    ParabolicIBVP2();

protected:
    virtual double initial(unsigned int n UNUSED_PARAM) const { return 0.0; }
    virtual double boundary(unsigned int m UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const { return 0.0; }
    virtual double f(unsigned int n UNUSED_PARAM, unsigned int m UNUSED_PARAM) const  { return 0.0; }
    virtual double a(unsigned int n UNUSED_PARAM, unsigned int m UNUSED_PARAM) const  { return 0.0; }

    virtual double initial(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const;
    virtual double boundary(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const;
    virtual double f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const;
    virtual double a(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const;

    virtual double initial(const SpaceNode &sn) const;
    virtual double boundary(const SpaceNode &sn, const TimeNode& tn) const;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const;
    virtual double a(const SpaceNode &sn, const TimeNode &tn) const;

public:
    void static Main(int argc, char* argv[]);

private:
    Dimension time;
    Dimension dim1;
    Dimension dim2;
};

#endif // PARABOLICIBVP2_H
