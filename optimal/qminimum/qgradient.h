#ifndef QGRADIENT_H
#define QGRADIENT_H

#include "qminimum_global.h"
#include <stdlib.h>
#include <string.h>

typedef double (*R1Function)(double);
typedef double (*R2Function)(double x, double y);
typedef double (*RnFunction)(double *x, int n);

#include <QVector>

class QMINIMUMSHARED_EXPORT QGradient
{

public:
    QGradient(RnFunction f, double *x0, int n);
    virtual ~QGradient();

    void calculate();
    double minimize(double);

    double* getX() const;

private:
    RnFunction fn;
    R1Function f1;
    double *x;
    int n;
};

#endif // QGRADIENT_H
