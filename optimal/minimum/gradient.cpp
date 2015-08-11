#include "gradient.h"

GradientMethod::GradientMethod() : m_fn(NULL), printer(NULL), projection(NULL)
{
    m_epsilon1 = 0.1;
    m_epsilon1 = 0.1;
    grad_step = 0.01;
    min_step = 0.1;
    min_epsilon = 0.01;
    iterationCount = 0;
    normalize = true;
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

double GradientMethod::epsilon1() const
{
    return m_epsilon1;
}

void GradientMethod::setEpsilon1(double epsilon)
{
    m_epsilon1 = epsilon;
}

double GradientMethod::epsilon2() const
{
    return m_epsilon2;
}

void GradientMethod::setEpsilon2(double epsilon)
{
    m_epsilon2 = epsilon;
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

void GradientMethod::setPrinter(Printer *printer)
{
    this->printer = printer;
}

void GradientMethod::setProjection(Projection *proj)
{
    this->projection = proj;
}

void GradientMethod::setNormalize(bool normalize)
{
    this->normalize = normalize;
}
