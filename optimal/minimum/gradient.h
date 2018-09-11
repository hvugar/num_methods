#ifndef GRADIENT_H
#define GRADIENT_H

#include "global.h"

#include "r1minimize.h"
#include "exceptions.h"
#include "function.h"

class RnFunction;
class IGradient;
class IPrinter;
class IProjection;
struct IGradientPrinter;

/**
 * @brief The Abstract Gradient Method class
 */
class MINIMUMSHARED_EXPORT GradientMethod : protected R1Function
{
public:
    enum MethodResult
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
     * @brief Epsilon for gradient norm
     * @return
     */
    double epsilon1() const;
    void setEpsilon1(double epsilon);

    /**
     * @brief Epsilon for points distance
     * @return
     */
    double epsilon2() const;
    void setEpsilon2(double epsilon);

    /**
     * @brief Epsilon for function result
     * @return
     */
    double epsilon3() const;
    void setEpsilon3(double epsilon);

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
    int count() const;

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
     * @brief showEndMessage
     * @param showEndMessage
     */
    void showEndMessage(bool showEndMessage);

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

    RnFunction *m_fn;
    IGradient *m_gr;
    IPrinter* m_printer;
    IGradientPrinter *m_printer_gr;
    IProjection *m_projection;
    double m_epsilon1;
    double m_epsilon2;
    double m_epsilon3;
    double min_step;
    double min_epsilon;
    int m_iteration_count;
    bool m_normalize;
    bool m_show_end_message;
};

struct MINIMUMSHARED_EXPORT IGradientPrinter
{
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double fxResult, GradientMethod::MethodResult result, double alpha);
};

#endif // GRADIENT_H
