#include "gradient.h"
#include "function.h"
#include <math.h>

GradientMethod::GradientMethod() : m_fn(NULL), m_gr(NULL), m_printer(NULL), m_printer_gr(NULL), m_projection(NULL),
    m_epsilon1(0.1), m_epsilon2(0.1), m_epsilon3(0.1), min_step(0.1), min_epsilon(0.01),
    m_iteration_count(0), m_normalize(true), m_show_end_message(true)
{}

GradientMethod::~GradientMethod() {}

/**
 * @brief Objective function
 * @param f
 */
void GradientMethod::setFunction(RnFunction *fn)
{
    m_fn = fn;
}

/**
 * @brief Objective function
 * @return
 */
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

/**
 * @brief Epsilon for gradient norm
 * @return
 */
double GradientMethod::epsilon1() const
{
    return m_epsilon1;
}

/**
 * @brief Epsilon for gradient norm
 * @param epsilon
 */
void GradientMethod::setEpsilon1(double epsilon)
{
    m_epsilon1 = epsilon;
}

/**
 * @brief Epsilon for distance between points
 * @return
 */
double GradientMethod::epsilon2() const
{
    return m_epsilon2;
}

/**
 * @brief Epsilon for distance between points
 * @param epsilon
 */
void GradientMethod::setEpsilon2(double epsilon)
{
    m_epsilon2 = epsilon;
}

/**
 * @brief Epsilon for distance between points
 * @param epsilon
 */
void GradientMethod::setEpsilon3(double epsilon)
{
    m_epsilon3 = epsilon;
}

/**
 * @brief Epsilon for distance between points
 * @return
 */
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
    return m_iteration_count;
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
    m_show_end_message = showEndMessage;
}

