#include "diffequ.h"

OrdinaryDifferentialEquation::OrdinaryDifferentialEquation() {}

OrdinaryDifferentialEquation::OrdinaryDifferentialEquation(const OrdinaryDifferentialEquation &) {}

OrdinaryDifferentialEquation& OrdinaryDifferentialEquation::operator=(const OrdinaryDifferentialEquation &other)
{
    if (this == &other) { return *this; }
    return *this;
}

OrdinaryDifferentialEquation::~OrdinaryDifferentialEquation() {}

/**************************************************************************************************************/

ExceptionODE::ExceptionODE(unsigned int msgCode) NOEXCEPT : _msgCode(msgCode)  {}

ExceptionODE::~ExceptionODE() {}

const char* ExceptionODE::what() const NOEXCEPT
{
    if (_msgCode == 1) return "Error: Dimension of matrices do not matches!";
    if (_msgCode == 2) return "Error: Non-local condition matrix is not square matrix!";
    if (_msgCode == 3) return "Error: Non-local condition matrix dimension is not matches with count of differensial equations!";
    if (_msgCode == 4) return "Error: Rigth side of non-local conditions is not matches with count of differensial equations!";
    if (_msgCode == 5) return "Error: Equations count must be one.";
    return "Error: Unknown error!";
}

ExceptionPDE::ExceptionPDE(unsigned int msgCode) NOEXCEPT : _msgCode(msgCode)  {}

ExceptionPDE::~ExceptionPDE() {}

const char* ExceptionPDE::what() const NOEXCEPT
{
    if (_msgCode == 1) return "Error: Dimension of matrices do not matches!";
    if (_msgCode == 2) return "Error: Non-local condition matrix is not square matrix!";
    if (_msgCode == 3) return "Error: Non-local condition matrix dimension is not matches with count of differensial equations!";
    if (_msgCode == 4) return "Error: Rigth side of non-local conditions is not matches with count of differensial equations!";
    return "Error: Unknown error!";
}
