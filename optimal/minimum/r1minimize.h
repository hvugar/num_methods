#ifndef R1MINIMIZE_H
#define R1MINIMIZE_H

#include "global.h"
#include "function.h"

#ifdef __cplusplus
extern "C" {
#endif
void MINIMUMSHARED_EXPORT stranghLineSearch(double x, double step, double &a, double &b, R1Function *fn) noexcept;
double MINIMUMSHARED_EXPORT goldenSectionSearch(double &a, double &b, double &x, R1Function *f, double epsilon);
double MINIMUMSHARED_EXPORT goldenSectionSearch1(double &a, double &b, double &x, R1Function *f, double epsilon);
double MINIMUMSHARED_EXPORT goldenSectionSearch2(double &a, double &b, double &x, R1Function *f, double epsilon);
#ifdef __cplusplus
}
#endif

class MINIMUMSHARED_EXPORT R1FxMinimizer
{
public:

    struct MINIMUMSHARED_EXPORT Callback
    {
        friend class R1FxMinimizer;
        R1Function* function() const;

        virtual void straightLineSearchCallback(unsigned int iteration, double x, double a, double b, double fxa, double fxb, unsigned int fx_count) const;
        virtual void swannCallback(unsigned int iteration, double x, double a, double b, double fxa, double fxb, unsigned int fx_count) const;

        virtual void goldenSectionSearchCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;
        virtual void halphIntervalMethodCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;
        virtual void dichotomyCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;
        virtual void uniformLineSearchCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;
        virtual void fibonachiMethodCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;

    private:
        R1Function* mfunction;
    };

    enum UnimodalResultResult
    {
        UnimodalFunction = 0,
        NonUnimodalFunction = 1
    };

    void setFunction(R1Function *function);

    R1Function* function() const;

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

    /**
     * @brief          Метод дихотомии.
     *                 Метод относится к последовательным стратегиям. Задается начальный интервал неопределенности и
     *                 требуемая точность. Алгоритм опирается на анализ значений функции в двух точках. Для их нахождения
     *                 текущий интервал неопределенности делится пополам и в обе стороны от середины откладывается по — step/2,
     *                 где step - малое положительное число. Условия окончания процесса поиска стандартные: поиск заканчивается,
     *                 когда длина текущего интервала неопределенности оказывается меньше установленной величины.
     * @param x        Точка минимума.
     * @param a        Начальная точка отрезка.
     * @param b        Конечнная точка отрезка.
     * @param epsilon  Число эпсилон для останова метода.
     * @param step     Малое число.
     */
    void dichotomyMethod(double &x, double &a, double &b, double epsilon, double step) const;

    /**
     * @brief          Метод равномерного поиска.(метод перебора)
     *                 Метод относится к пассивным стратегиям. Задается начальный интервал
     *                 неопределенности [a, b] и количество вычислений функции n.
     *                 Вычисления производятся в n равноотстоящих друг от друга точках (при этом
     *                 интервал делится на n + 1 равных интервалов). Путем сравнения величин
     *                 f(xi), i = 1,...,n находится точка xк, в которой значение функции наименьшее.
     *                 Искомая точка минимума х* считается заключенной в интервале [хk-1, хk+1].
     * @param x        Точка минимума.
     * @param a        Начальная точка отрезка.
     * @param b        Конечнная точка отрезка.
     * @param n        Количество вычислений функции.
     */
    void uniformLineSearch(double &x, double &a, double &b, unsigned int n) const;

    /**
     * @brief          Метод Фибоначчи.
     *                 Метод относится к последовательным стратегиям. Задается начальный интервал непределенности и
     *                 количество N вычислений функции. Алгоритм уменьшения интервала опирается на анализ значений
     *                 функции в двух точках. Точки вычисления функции находятся с использованием последовательности
     *                 из n + 1 чисел Фибоначчи. Как в методе золотого сечения, на первой итерации требуются два
     *                 вычисления функции, а на каждой последующей - только по одному. Условия окончания процесса
     *                 поиска стандартные: поиск заканчивается, когда длина текущего интервала неопределенности
     *                 оказывается меньше установленной величины.
     * @param x        Точка минимума.
     * @param a        Начальная точка отрезка.
     * @param b        Конечнная точка отрезка.
     * @param epsilon  Число эпсилон для останова метода.
     * @param step     Малое число.
     */
    void fibonachiMethod(double &x, double &a, double &b, double epsilon, double step) const;

    void setCallback(R1FxMinimizer::Callback *callback);
    R1FxMinimizer::Callback* callback() const;

private:
    R1Function *mfunction = nullptr;
    R1FxMinimizer::Callback *mcallback = nullptr;
};

#endif // R1MINIMIZE_H
