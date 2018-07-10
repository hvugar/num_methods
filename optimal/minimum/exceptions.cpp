#include "exceptions.h"


FunctionException::~FunctionException() throw()
{}

const char* FunctionException::what() const throw()
{
    return std::exception::what();
}


