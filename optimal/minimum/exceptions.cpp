#include "exceptions.h"

function_exception::function_exception(char const* const /*_Message*/) NOEXCEPT {}

function_exception::function_exception(char const* const /*_Message*/, int) NOEXCEPT {}

function_exception::function_exception(unsigned int msgCode) NOEXCEPT : msgCode(msgCode) {}

function_exception::function_exception(const function_exception & e) NOEXCEPT : msgCode(e.msgCode) {}

function_exception::~function_exception() NOEXCEPT {}

const char* function_exception::what() const NOEXCEPT
{
    return "";
}

function_exception& function_exception::operator =(const function_exception& other) NOEXCEPT
{
    if (this == &other) { return *this; }

    msgCode = other.msgCode;
    return *this;
}

/*********************************************/

double_matrix_exception::double_matrix_exception(unsigned int msgCode) NOEXCEPT : msgCode(msgCode) {}

double_matrix_exception::double_matrix_exception(const double_matrix_exception & e) NOEXCEPT : msgCode(e.msgCode) {}

double_matrix_exception::~double_matrix_exception() NOEXCEPT {}

const char* double_matrix_exception::what() const NOEXCEPT
{
    if (msgCode == 0) return "Error: Unknow error!";
    if (msgCode == 1) return "Error: Dimension of matrixs do not matches!";
    if (msgCode == 2) return "Error: Matrix is not square matrix!";
    if (msgCode == 3) return "Error: Matrix columns and row are not equals!";
    if (msgCode == 4) return "Error: Determinant of matrix is equal to zero.";
    if (msgCode == 5) return "Error: Row index out of range!";
    if (msgCode == 6) return "Error: Error: Column index out of range!";
    return "Error: Unknow error!";
}

double_matrix_exception& double_matrix_exception::operator =(const double_matrix_exception& other) NOEXCEPT
{
    if (this == &other) { return *this; }

    this->msgCode = other.msgCode;
    return *this;
}

/*********************************************/

double_vector_exception::double_vector_exception(unsigned int msgCode) NOEXCEPT : msgCode(msgCode) {}

double_vector_exception::double_vector_exception(const double_vector_exception & e) NOEXCEPT : msgCode(e.msgCode) {}

double_vector_exception::~double_vector_exception() NOEXCEPT {}

const char* double_vector_exception::what() const NOEXCEPT
{
    if (msgCode == 0) return "Error: Unknow error!";
    if (msgCode == 1) return "Error: Length of vertices do not matches!";
    if (msgCode == 2) return "Error: Index out of range!";
    return "Error: Unknow error!";
}

double_vector_exception& double_vector_exception::operator =(const double_vector_exception& other) NOEXCEPT
{
    if (this == &other) { return *this; }

    this->msgCode = other.msgCode;
    return *this;
}

/*********************************************/

delta_grid_exception::delta_grid_exception(const std::string &msg) NOEXCEPT : message(msg) {}

delta_grid_exception::delta_grid_exception(const delta_grid_exception&) NOEXCEPT {}

delta_grid_exception::~delta_grid_exception() NOEXCEPT {}

const char* delta_grid_exception::what() const noexcept
{
    //std::string msg = std::string("DeltaGridException: ") + message;
    return message.data();
}

delta_grid_exception& delta_grid_exception::operator =(const delta_grid_exception &other) NOEXCEPT
{
    if (this == &other) { return *this; }

    this->message = other.message;
    return *this;
}
