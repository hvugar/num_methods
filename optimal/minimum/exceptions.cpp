#include "exceptions.h"

exceptions::exceptions()
{

}

exceptions::~exceptions()
{

}

FunctionException::~FunctionException() throw()
{}

const char* FunctionException::what() const throw()
{
    return std::exception::what();
}


