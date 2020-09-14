#ifndef DIRECT_SEARCH_METHOD_H
#define DIRECT_SEARCH_METHOD_H

#include "global.h"
#include "vector2d.h"
#include "function.h"
#include "exceptions.h"

class IPrinter;
class IProjection;
class IGradientPrinter;
class IVectorNormalizer;

/**
 * @brief The Abstract Gradient Method class
 */
class MINIMUMSHARED_EXPORT DirectSearchMethod
{
public:
    DirectSearchMethod();
    virtual ~DirectSearchMethod();

    virtual void calculate(DoubleVector &x) = 0;

    virtual RnFunction* function() const;
    virtual void setFunction(RnFunction *function);

protected:
    RnFunction *m_fn;
};

#endif // DIRECT_SEARCH_METHOD_H
