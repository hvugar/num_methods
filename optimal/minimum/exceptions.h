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

class MINIMUMSHARED_EXPORT FunctionException : public std::exception
{
public:
    FunctionException() throw() {}
    virtual ~FunctionException() throw();
    virtual const char* what() const throw();
};

#endif // EXCEPTIONS_H
