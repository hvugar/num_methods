#include "gradient.h"
#include "function.h"
#include <math.h>

GradientMethod::GradientMethod() : m_fn(NULL), m_printer(NULL), m_projection(NULL)
{
    m_epsilon1 = m_epsilon2 = m_epsilon3 = 0.1;
    min_step = 0.1;
    min_epsilon = 0.01;
    iterationCount = 0;
    m_normalize = true;
    mshowEndMessage = true;
}

GradientMethod::~GradientMethod()
{
}

void GradientMethod::setFunction(RnFunction *fn)
{
    m_fn = fn;
}

RnFunction* GradientMethod::function() const
{
    return m_fn;
}

IGradient* GradientMethod::gradient() const
{
    return m_gr;
}

void GradientMethod::setGradient(IGradient *gr)
{
    m_gr = gr;
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

void GradientMethod::setEpsilon3(double epsilon)
{
    m_epsilon3 = epsilon;
}

double GradientMethod::epsilon3() const
{
    return m_epsilon3;
}

void GradientMethod::setR1MinimizeEpsilon(double step, double epsilon)
{
    min_step = step;
    min_epsilon = epsilon;
}

int GradientMethod::count() const
{
    return iterationCount;
}

void GradientMethod::setPrinter(IPrinter *printer)
{
    this->m_printer = printer;
}

void GradientMethod::setProjection(IProjection *projection)
{
    this->m_projection = projection;
}

void GradientMethod::setNormalize(bool normalize)
{
    this->m_normalize = normalize;
}

void GradientMethod::showEndMessage(bool showEndMessage)
{
    mshowEndMessage = showEndMessage;
}
