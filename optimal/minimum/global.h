#ifndef GLOBAL_H
#define GLOBAL_H

#ifdef _WIN32
#if defined(MINIMUM_LIBRARY)
#  define MINIMUMSHARED_EXPORT __declspec(dllexport)
#else
#  define MINIMUMSHARED_EXPORT __declspec(dllimport)
#endif
#else
#define MINIMUMSHARED_EXPORT
#endif

#define C_UNUSED(x) (void)x

#ifdef __GNUC__
#define UNUSED_PARAM __attribute__ ((unused))
#else
#define UNUSED_PARAM
#if !defined (M_PI)
#define M_PI 3.14159265358979323846
#endif
#endif

#if (defined(__cplusplus))

#if (__cplusplus == 1L)
#define NOEXCEPT
#endif

#if (__cplusplus == 199711L)
#define NOEXCEPT
#endif

#if (__cplusplus == 201103L)
#define NOEXCEPT noexcept
#endif

#if (__cplusplus == 201402L)
#define NOEXCEPT noexcept
#endif

#if (__cplusplus == 201500L)
#define NOEXCEPT
#endif

#if (__cplusplus == 201703L)
#define NOEXCEPT
#endif

#endif

#include <assert.h>

#define PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(ClassName) \
    public: \
    ClassName(); \
    ClassName(const ClassName &); \
    ClassName & operator =(const ClassName &); \
    virtual ~ClassName()

#define PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR_IMPL(ClassName) \
    ClassName::ClassName() {} \
    ClassName::ClassName(const InitialConditionODE &) {} \
    ClassName& ClassName::operator=(const ClassName &other) { \
        if (this == &other) { return *this; } \
        this->initialConditionType = other.initialConditionType; \
        this->value = other.value; \
        return *this; \
    } \
    InitialConditionODE::~InitialConditionODE() {}

//#define NDEBUG

#endif // GLOBAL
