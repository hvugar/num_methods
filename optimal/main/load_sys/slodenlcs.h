#ifndef SYSTEMLINEARODENONLOCALCONTIONS_H
#define SYSTEMLINEARODENONLOCALCONTIONS_H

#define SAMPLE_3

#include <load_sys/islodenlcs.h>

class SystemLinearODENonLocalContions : public ISystemLinearODENonLocalContions
{
public:
    SystemLinearODENonLocalContions(const ODEGrid &grid);
    static void Main(int agrc, char *argv[]);

    void initialize();
    double x(double t, int i) const;

public:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const;
};

#endif // SYSTEMLINEARODENONLOCALCONTIONS_H
