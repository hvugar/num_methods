#ifndef DIFFENSIALEQUATION_H
#define DIFFENSIALEQUATION_H

#include "../global.h"
#include "../vector2d.h"
#include "../grid/grid.h"
#include "../linearequation.h"

/**
 * @brief Дифференциа́льное уравне́ние
 * Дифференциа́льное уравне́ние — уравнение, в которое входят производные функции,
 * и может входить сама функция, независимая переменная и параметры. Порядок входящих
 * в уравнение производных может быть различен (формально он ничем не ограничен).
 * Производные, функции, независимые переменные и параметры могут входить в уравнение
 * в различных комбинациях или могут отсутствовать вовсе, кроме хотя бы одной производной.
 * Не любое уравнение, содержащее производные неизвестной функции, является дифференциальным
 * уравнением.
 * @see OrdinaryDifferentialEquation
 * @see
 */
class MINIMUMSHARED_EXPORT DifferentialEquation
{
public:
    DifferentialEquation();
    virtual ~DifferentialEquation();
};

/**
 * @brief Обыкновенное дифференциальное уравнение
 * Обыкновенное дифференциальное уравне́ние (ОДУ) — это дифференциальное уравнение
 * для функции от одной переменной. (Этим оно отличается от уравнения в частных производных,
 * где неизвестная — это функция нескольких переменных.). Таким образом, ОДУ — это уравнения
 * вида F(x,y',y",...,y^(n)) = 0
 * где y (x) — неизвестная функция (возможно, вектор-функция, тогда F, как правило, тоже
 * вектор-функция со значениями в пространстве той же размерности; в этом случае говорят о
 * системе дифференциальных уравнений), зависящая от независимой
 * переменной x, штрих означает дифференцирование по x. Число n (порядок старшей производной,
 * входящей в данное уравнение) называется порядком дифференциального уравнения.
 * @see LinearODE
 * @see NonLinearODE
 */
class MINIMUMSHARED_EXPORT OrdinaryDifferentialEquation : public DifferentialEquation
{
public:
    OrdinaryDifferentialEquation();
    virtual ~OrdinaryDifferentialEquation();

    enum OdeSolverMethod
    {
        RK2,
        RK4,
        EULER,
        EULER_MOD
    };

    enum Direction
    {
        L2R, // Left to Right
        R2L  // Right to Left
    };

    const Dimension& dimension() const;
    void setDimension(const Dimension &dimension);

    virtual unsigned int count() const = 0;

private:
    Dimension _dimension;
};

/**
 * @brief Линейное дифференциальное уравнение с переменными коэффициентами
 * @see FirstOrderLinearODE
 * @see SecondOrderLinearODE
 */
class MINIMUMSHARED_EXPORT LinearODE : public OrdinaryDifferentialEquation
{};

/**
 * @brief The NonLinearODE class
 * @see FirstOrderNonLinearODE
 * @see SecondOrderNonLinearODE
 */
class MINIMUMSHARED_EXPORT NonLinearODE : public OrdinaryDifferentialEquation {};

class MINIMUMSHARED_EXPORT SystemDifferentialEquation
{
public:
    const Dimension& dimension() const;
    void setDimension(const Dimension &dimension);

protected:
    Dimension _dimension;
};

class ExceptionODE : std::exception
{
public:
    explicit ExceptionODE(unsigned int msgCode = 0) NOEXCEPT;
    virtual ~ExceptionODE();

    virtual const char* what() const NOEXCEPT;

private:
    unsigned int _msgCode;
};

class ExceptionPDE : std::exception
{
public:
    explicit ExceptionPDE(unsigned int msgCode = 0) NOEXCEPT;
    virtual ~ExceptionPDE();

    virtual const char* what() const NOEXCEPT;

private:
    unsigned int _msgCode;
};

#endif // DIFFENSIALEQUATION_H
