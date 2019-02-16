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
 * уравнением. Например,   f′(x) = f ( f ( x ) ) {\displaystyle \ f'(x)=f(f(x))} \ f'(x)=f(f(x)) не является дифференциальным уравнением
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
 *
 */
class MINIMUMSHARED_EXPORT LinearODE : public OrdinaryDifferentialEquation
{};

/**
 * @brief The NonLinearODE class
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

//class MINIMUMSHARED_EXPORT SystemDifferentialEquationODE : public SystemDifferentialEquation
//{};

//class MINIMUMSHARED_EXPORT SystemLinearODE : public SystemDifferentialEquationODE {};

//class MINIMUMSHARED_EXPORT SystemLinearODE1stOrder : public SystemLinearODE
//{
//protected:
//    virtual double A(double t, unsigned int k, unsigned int row = 0, unsigned int col = 0) const = 0;
//    virtual double B(double t, unsigned int k, unsigned int row = 0) const = 0;
//};

//class MINIMUMSHARED_EXPORT SystemNonLinearODE : public SystemDifferentialEquationODE
//{};

//class MINIMUMSHARED_EXPORT SystemNonLinearODE1stOrder : public SystemNonLinearODE
//{
//protected:
//    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const = 0;
//};

#endif // DIFFENSIALEQUATION_H
