#include "gradient.h"

GradientMethod::GradientMethod() : m_fn(NULL), printer(NULL)
{
    m_epsilon = 0.0;
    grad_step = 0.0;
    min_step = 0.0;
    min_epsilon = 0.0;
    iterationCount = 0;
}

GradientMethod::~GradientMethod()
{
}

void GradientMethod::setFunction(RnFunction *f)
{
    m_fn = f;
}

RnFunction* GradientMethod::function() const
{
    return m_fn;
}

double GradientMethod::epsilon() const
{
    return m_epsilon;
}

void GradientMethod::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void GradientMethod::setR1MinimizeEpsilon(double step, double epsilon)
{
    min_step = step;
    min_epsilon = epsilon;
}

void GradientMethod::setGradientStep(double step)
{
    grad_step = step;
}

int GradientMethod::count() const
{
    return iterationCount;
}

void GradientMethod::setPrinter(GrPrinter *printer)
{
    this->printer = printer;
}
