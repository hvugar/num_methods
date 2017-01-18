#ifndef GRIDMETHOD_H
#define GRIDMETHOD_H

#include "global.h"
#include "vector2d.h"
#include "matrix2d.h"
#include "cmethods.h"

class MINIMUMSHARED_EXPORT Grid1DSetting
{
public:
    Grid1DSetting();
    virtual ~Grid1DSetting() {}

public:
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double x0;
    double x1;
    double t0;
    double t1;
};

class MINIMUMSHARED_EXPORT GridMethod
{
public:
    GridMethod();
    virtual ~GridMethod() {}

    void setGridSetting(const Grid1DSetting &setting);

protected:
    Grid1DSetting setting;
};

#endif // GRIDMETHOD_H
