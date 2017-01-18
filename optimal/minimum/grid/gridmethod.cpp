#include "gridmethod.h"

Grid1DSetting::Grid1DSetting()
{
    hx = 0.01;
    ht = 0.01;
    N  = M  = 100;
    x1 = x0 = 0.0;
    t1 = t0 = 1.0;
}

GridMethod::GridMethod()
{}

void GridMethod::setGridSetting(const Grid1DSetting &setting)
{
    this->setting = setting;
}
