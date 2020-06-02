#ifndef DIFFENSIALEQUATION_H
#define DIFFENSIALEQUATION_H

#include "../global.h"
#include "../vector2d.h"
#include "../grid/grid.h"
#include "../grid/ibvp.h"
#include "../linearequation.h"

enum class ODESolverMethod
{
    EULER,
    EULER_MOD,
    RUNGE_KUTTA_2,
    RUNGE_KUTTA_4,
    RUNGE_KUTTA_6
};

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
class MINIMUMSHARED_EXPORT DifferentialEquation {};

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
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(OrdinaryDifferentialEquation);

    enum Direction
    {
        L2R, // Left to Right
        R2L  // Right to Left000
    };

    virtual auto dimension() const -> Dimension = 0;
    virtual auto count() const -> size_t = 0;
};

class MINIMUMSHARED_EXPORT CanonicalFormODE : public OrdinaryDifferentialEquation {};

/**
 * @brief Линейное дифференциальное уравнение с переменными коэффициентами
 * @see FirstOrderLinearODE
 * @see SecondOrderLinearODE
 */
class MINIMUMSHARED_EXPORT LinearODE : public CanonicalFormODE {};

class MINIMUMSHARED_EXPORT IFirstOrderLinearODE : public LinearODE
{
protected:
    /**
     * @brief A  A nxn dimensional matrix-function
     * @param node
     * @param row <= n
     * @param col <= n
     * @return
     */
    virtual auto A(const PointNodeODE &node, size_t row = 1,size_t col = 1) const -> double = 0;

    /**
     * @brief B n dimensional vector-function
     * @param node
     * @param row
     * @return
     */
    virtual auto B(const PointNodeODE &node, size_t row = 1) const -> double = 0;
};

class MINIMUMSHARED_EXPORT ISecondOrderLinearODE : public LinearODE
{
protected:
    /**
     * @brief A  A nxn dimensional matrix-function
     * @param node
     * @param row <= n
     * @param col <= n
     * @return
     */
    virtual auto A(const PointNodeODE &node, size_t row = 1, size_t col = 1) const -> double = 0;

    /**
     * @brief B B nxn dimensional matrix-function
     * @param node
     * @param row
     * @return
     */
    virtual auto B(const PointNodeODE &node, size_t row = 1, size_t col = 1) const -> double = 0;

    /**
     * @brief C B nxn dimensional vector-function
     * @param node
     * @param row
     * @return
     */
    virtual auto C(const PointNodeODE &node, size_t row = 1) const -> double = 0;
};


/**
 * @brief The NonLinearODE class
 * @see FirstOrderNonLinearODE
 * @see SecondOrderNonLinearODE
 */
class MINIMUMSHARED_EXPORT NonLinearODE : public CanonicalFormODE {};

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
