#ifndef GRADIENT_BASED_METHOD_H
#define GRADIENT_BASED_METHOD_H

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
class MINIMUMSHARED_EXPORT GradientBasedMethod : protected R1Function
{
public:
    enum class MethodResult
    {
        NEXT_ITERATION = 0x00000000,
        FIRST_ITERATION = 0x00000001,
        BREAK_FIRST_ITERATION = 0x00000002,
        //BREAK_GRADIENT_NORM_LESS = 0x00000004,
        //BREAK_DISTANCE_LESS = 0x00000008,
        BREAK_OPTIMALITY_TOLERANCE = 0x00000010,
        BREAK_STEP_TOLERANCE = 0x00000020,
        BREAK_FUNCTION_TOLERANCE = 0x00000040,
        BREAK_CONSTRAAINT_TOLERANCE = 0x00000080,
        BREAK_ITERATION_NUMBER = 0x00000100,
        BREAK_FUNCTION_EVALUATION_NUMBER = 0x00000200
    };

public:
    GradientBasedMethod();
    virtual ~GradientBasedMethod();

    virtual void calculate(DoubleVector &x) = 0;

    virtual RnFunction* function() const;
    virtual void setFunction(RnFunction *function);

    virtual IGradient* gradient() const;
    virtual void setGradient(IGradient *gradient);

    /**
     * @brief optimalityTolerance
     * @return
     */
    double optimalityTolerance() const;
    /**
     * @brief setOptimalityTolerance
     * @param optimalityTolerance
     */
    void setOptimalityTolerance(double optimalityTolerance);
    /**
     * @brief stepTolerance
     * @return
     */
    double stepTolerance() const;
    /**
     * @brief setStepTolerance
     * @param stepTolerance
     */
    void setStepTolerance(double stepTolerance);
    /**
     * @brief functionTolerance
     * @return
     */
    double functionTolerance() const;
    /**
     * @brief setFunctionTolerance
     * @param functionTolerance
     */
    void setFunctionTolerance(double functionTolerance);
    /**
     * @brief constraintTolerance
     * @return
     */
    double constraintTolerance() const;
    /**
     * @brief setConstraintTolerance
     * @param functionTolerance
     */
    void setConstraintTolerance(double constraintTolerance);
    /**
     * @brief setTolerance
     * @param optimalityTolerance
     * @param stepTolerance
     * @param functionTolerance
     */
    void setTolerance(double optimalityTolerance, double stepTolerance, double functionTolerance);
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
    void setPrinter(IPrinter *printer);
    /**
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
    void setMaxIterationCount(unsigned int maxIterationCount);
    /**
     * @brief Максимальное количество повторений позволено.
     * @return
     */
    unsigned int maxIterationCount() const;
    /**
     * @brief Максимальное количество оценок функции позволено.
     * @param maxFunctionEvaluations
     */
    void setMaxFunctionEvaluationCount(unsigned int maxFunctionEvaluationCount);
    /**
     * @brief Максимальное количество оценок функции позволено.
     * @return
     */
    unsigned int maxFunctionEvaluationCount() const;

    unsigned int iterationNumber() const;

    unsigned int maxFunctionEvaluationNumber() const;

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

    bool checkForExit(double step_tolerance, double optimalityTolerance, double function_tolerance,
                      double f2, const DoubleVector &x, const DoubleVector &g, double alpha) const;

protected:
    RnFunction *m_fn;
    IGradient *m_gr;
    IPrinter* m_printer;
    IProjection *m_projection;
    double m_optimalityTolerance;
    double m_functionTolerance;
    double m_stepTolerance;
    double m_constraintTolerance;
public:
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

#endif // GRADIENT_BASED_METHOD_H
