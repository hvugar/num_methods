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

#endif // EXCEPTIONS_H
