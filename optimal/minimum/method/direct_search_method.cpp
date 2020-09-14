#include "direct_search_method.h"
#include "vectornormalizer.h"
#include "printer.h"
#include "function.h"
#include <math.h>
#include <float.h>
#include <stdexcept>

/*****************************************************************************************************/

DirectSearchMethod::DirectSearchMethod() : m_fn(nullptr)
{}

DirectSearchMethod::~DirectSearchMethod() {}

/**
 * @brief Objective function
 * @param f
 */
void DirectSearchMethod::setFunction(RnFunction *fn)
{
    m_fn = fn;
}

/**
 * @brief Objective function
 * @return
 */
RnFunction* DirectSearchMethod::function() const
{
    return m_fn;
}
