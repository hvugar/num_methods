#ifndef GRADIENT_METHOD_H
#define GRADIENT_METHOD_H

#include "global.h"
#include "vector2d.h"
#include "function.h"
#include "exceptions.h"

class IPrinter;
class IProjection;
class IGradientPrinter;
class IVectorNormalizer;

/**
 * @brief The Abstract Gradient Method class
 */
class MINIMUMSHARED_EXPORT GradientMethod : protected R1Function
{
public:
    enum class MethodResult
    {
        NEXT_ITERATION = 0,
        FIRST_ITERATION = 1,
        BREAK_FIRST_ITERATION = 2,
        BREAK_GRADIENT_NORM_LESS = 3,
        BREAK_DISTANCE_LESS = 4
    };

public:
    GradientMethod();
    virtual ~GradientMethod();

    virtual void calculate(DoubleVector &x) = 0;

    virtual RnFunction* function() const;
    virtual void setFunction(RnFunction *function);

    virtual IGradient* gradient() const;
    virtual void setGradient(IGradient *gradient);

    /**
     * @brief optimalityTolerance
     * @return
     */
    auto optimalityTolerance() const -> double;
    /**
     * @brief setOptimalityTolerance
     * @param optimalityTolerance
     */
    auto setOptimalityTolerance(double optimalityTolerance) -> void;
    /**
     * @brief stepTolerance
     * @return
     */
    auto stepTolerance() const -> double;
    /**
     * @brief setStepTolerance
     * @param stepTolerance
     */
    auto setStepTolerance(double stepTolerance) -> void;
    /**
     * @brief functionTolerance
     * @return
     */
    auto functionTolerance() const -> double;
    /**
     * @brief setFunctionTolerance
     * @param functionTolerance
     */
    auto setFunctionTolerance(double functionTolerance) -> void;
    /**
     * @brief constraintTolerance
     * @return
     */
    auto constraintTolerance() const -> double;
    /**
     * @brief setConstraintTolerance
     * @param functionTolerance
     */
    auto setConstraintTolerance(double constraintTolerance) -> void;
    /**
     * @brief setTolerance
     * @param optimalityTolerance
     * @param stepTolerance
     * @param functionTolerance
     */
    auto setTolerance(double optimalityTolerance, double stepTolerance, double functionTolerance) -> void;
    /**
     * @brief setR1MinimizeEpsilon
     * @param step
     * @param epsilon1
     */
    void setR1MinimizeEpsilon(double step, double epsilon);
    /**
     * @brief count
     * @return
     */
    unsigned int count() const;
    /**
     * @brief Вывод информации о значениях параметров оптимизации на каждой итерации.
     * @param printer Интерфейс для вывода информации о значениях параметров оптимизации на каждой итерации.
     */
    void setPrinter(IPrinter *printer);    /**
     * @brief Метод проекции для значений параметров оптимизации.
     * @param projection Интерфейс для проекции значений параметров оптимизации.
     */
    void setProjection(IProjection *projection);
    /**
     * @brief setNormalize
     * @param normalize
     */
    void setNormalize(bool normalize);
    /**
     * @brief showExitMessage
     * @param shem
     */
    void showExitMessage(bool shem);
    /**
     * @brief setGradientNormalizer
     * @param normalizer
     */
    void setGradientNormalizer(IVectorNormalizer *normalizer);
    /**
     * @brief Максимальное количество повторений позволено.
     * @param maxIterations
     */
    auto setMaxIterationCount(unsigned int maxIterationCount) -> void;
    /**
     * @brief Максимальное количество повторений позволено.
     * @return
     */
    auto maxIterationCount() const -> unsigned int;
    /**
     * @brief Максимальное количество оценок функции позволено.
     * @param maxFunctionEvaluations
     */
    auto setMaxFunctionEvaluationCount(unsigned int maxFunctionEvaluationCount) -> void;
    /**
     * @brief Максимальное количество оценок функции позволено.
     * @return
     */
    auto maxFunctionEvaluationCount() const -> unsigned int;

    auto iterationNumber() const -> unsigned int;

    auto maxFunctionEvaluationNumber() const -> unsigned int;

protected:
    /**
     * @brief minimize
     * @param x
     * @param g
     * @return
     */
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) const = 0;
    /**
     * @brief          Метод прямого поиска. Установления границ интервала.
     * @param x        Произвольно выбранная начальная точка.
     * @param step     Величина шага.
     * @param a        Начальная точка отрезка.
     * @param b        Конечнная точка отрезка.
     * @param fxa      Величина функции в точке a.
     * @param fxb      Величина функции в точке b.
     * @param unimodal
     */
    void straightLineSearch(double x, double step, double &a, double &b, double &fxa, double &fxb, bool &unimodal) const;
    /**
     * @brief          Метод Свенна. Установления границ интервала.
     * @param x        Произвольно выбранная начальная точка.
     * @param step     Величина шага.
     * @param a        Начальная точка отрезка.
     * @param b        Конечнная точка отрезка.
     * @param fxa      Величина функции в точке a.
     * @param fxb      Величина функции в точке b.
     * @param unimodal
     */
    void swann(double x, double step, double &a, double &b, double &fx, double &fxb, bool &unimodal) const;
    /**
     * @brief          Метод золотого сечения.
     *                 Метод относится к последовательным стратегиям. Задается начальный интервал неопределенности и
     *                 требуемая точность. Алгоритм уменьшения интервала опирается на анализ значений функции в двух точках.
     *                 В качестве точек вычисления функции выбираются точки золотого сечения. Тогда с учетом свойств золотого
     *                 сечения на каждой итерации, кроме первой, требуется только одно новое вычисление функции. Условия
     *                 окончания процесса поиска стандартные: поиск заканчивается, когда длина текущего интервала
     *                 неопределенности оказывается меньше установленной величины.
     * @param x        Точка минимума.
     * @param a        Начальная точка отрезка.
     * @param b        Конечнная точка отрезка.
     * @param epsilon  Число эпсилон для останова метода.
     */
    void goldenSectionSearch(double &x, double &a, double &b, double epsilon) const;
    /**
     * @brief          Метод деления интервала пополам.
     *                 Метод относится к последовательным стратегиям и позволяет исключить из дальнейшего
     *                 рассмотрения на каждой итерации в точности половину текущего интервала неопределенности.
     *                 Задается начальный интервал неопределенности, а алгоритм уменьшения интервала, являясь, как и
     *                 в общем случае, "гарантирующим", основан на анализе величин функции в трех точках, равномерно
     *                 распределенных на текущем интервале (делящих его на четыре равные части). Условия окончания
     *                 процесса поиска стандартные: поиск заканчивается, когда длина текущего интервала неопределенности
     *                 оказывается меньше установленной величины.
     * @param x        Точка минимума.
     * @param a        Начальная точка отрезка.
     * @param b        Конечнная точка отрезка.
     * @param epsilon  Число эпсилон для останова метода.
     */
    void halphIntervalMethod(double &x, double &a, double &b, double epsilon) const;

protected:
    RnFunction *m_fn;
    IGradient *m_gr;
    IPrinter* m_printer;
    IProjection *m_projection;
    double m_optimalityTolerance;
    double m_functionTolerance;
    double m_stepTolerance;
    double m_constraintTolerance;
    double min_step;
    double min_epsilon;
    bool m_show_end_message;
    bool m_normalize;
    IVectorNormalizer *m_normalizer;
    unsigned int m_iterationNumber;
    unsigned int m_maxIterationCount;
    unsigned int m_functionEvaluationNumber;
    unsigned int m_maxFunctionEvaluationCount;
};

#endif // GRADIENT_METHOD_H
