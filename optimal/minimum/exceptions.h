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

class MINIMUMSHARED_EXPORT MatrixException : public std::exception
{
public:
    MatrixException(int messageType);
    virtual const char* what() const noexcept;
private:
    int messageType;
};


#endif // EXCEPTIONS_H
