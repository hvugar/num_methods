#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include "global.h"
#include <exception>
#include <string>

class function_exception : public std::exception
{
public:
    MINIMUMSHARED_EXPORT explicit function_exception(char const* const _Message) NOEXCEPT;
    MINIMUMSHARED_EXPORT function_exception(char const* const _Message, int) NOEXCEPT;

    MINIMUMSHARED_EXPORT explicit function_exception(unsigned int msgCode = 0) NOEXCEPT;
    MINIMUMSHARED_EXPORT function_exception(const function_exception &) NOEXCEPT;
    MINIMUMSHARED_EXPORT virtual ~function_exception() NOEXCEPT;
    MINIMUMSHARED_EXPORT virtual const char* what() const NOEXCEPT;

    MINIMUMSHARED_EXPORT function_exception& operator= (const function_exception &) NOEXCEPT;

private:
    unsigned int msgCode;
};

class double_matrix_exception : public std::exception
{
public:
    double_matrix_exception(unsigned int msgCode = 0) NOEXCEPT;
    double_matrix_exception(const double_matrix_exception&) NOEXCEPT;
    virtual ~double_matrix_exception() NOEXCEPT;
    virtual const char* what() const NOEXCEPT;

    double_matrix_exception& operator= (const double_matrix_exception &other) NOEXCEPT;

private:
    unsigned int msgCode;
};

class double_vector_exception : public std::exception
{
public:
    double_vector_exception(unsigned int msgCode = 0) NOEXCEPT;
    double_vector_exception(const double_vector_exception&) NOEXCEPT;
    virtual ~double_vector_exception();
    virtual const char* what() const NOEXCEPT;

    double_vector_exception& operator= (const double_vector_exception &other) NOEXCEPT;

private:
    unsigned int msgCode;
};

class delta_grid_exception : public std::exception
{
public:
    delta_grid_exception(const std::string &msg = "") NOEXCEPT;
    delta_grid_exception(const delta_grid_exception&) NOEXCEPT;
    virtual ~delta_grid_exception() NOEXCEPT;
    virtual const char* what() const NOEXCEPT;

    delta_grid_exception& operator= (const delta_grid_exception &other) NOEXCEPT;

private:
    std::string message;
};

#endif // EXCEPTIONS_H
