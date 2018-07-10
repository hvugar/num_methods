#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include "global.h"
#include <stdexcept>

class exceptions
{
public:
    exceptions();
    ~exceptions();
};


class invalid_step_exception : public std::invalid_argument
{

};

class FunctionException : public std::exception
{
public:
    MINIMUMSHARED_EXPORT FunctionException() throw() {}
    MINIMUMSHARED_EXPORT virtual ~FunctionException() throw();
    MINIMUMSHARED_EXPORT virtual const char* what() const throw();
};

#endif // EXCEPTIONS_H
