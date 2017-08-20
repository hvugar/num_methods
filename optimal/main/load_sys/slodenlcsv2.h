#ifndef SYSTEMLINEARODENONLOCALCONTIONSV2_H
#define SYSTEMLINEARODENONLOCALCONTIONSV2_H

#include <load_sys/islodenlcsv2.h>

class SystemLinearODENonLocalContionsV2 : public ISystemLinearODENonLocalContionsV2
{
public:
    static void Main(int agrc, char *argv[]);

public:
    SystemLinearODENonLocalContionsV2();

    virtual double A(TimeNode node, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(TimeNode node, unsigned int s, unsigned int row = 0, unsigned int col = 0) const;
    virtual double C(TimeNode node, unsigned int row = 0) const;

    virtual void calculateForward();

    double X(double t, unsigned int row) const;
private:
    void fillMatrix(DoubleMatrix &M, DoubleVector &P, unsigned int st, unsigned int n, unsigned int k1, unsigned int k2,
                    unsigned int row);
};

#endif // SYSTEMLINEARODENONLOCALCONTIONSV2_H
