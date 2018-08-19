#include "r1minimize.h"
#include <math.h>
#include <float.h>

#include <stdexcept>
#include <string>

/**
 * @brief         Метод прямого поиска. Установления границ интервала.
 * @param x       Произвольно выбранная начальная точка
 * @param step    Величина шага
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param fxa     Величина функции в точке a
 * @param fxb     Величина функции в точке b
 */
void R1FxMinimizer::straightLineSearch(double x, double step, double &a, double &b, double &fxa, double &fxb) const
{
    unsigned int iteration = 0;
    unsigned int fx_count = 0;

    if ( mfunction == NULL )
    {
        std::string msg = "in function \"stranghLineSearch\" function pointer is null.";
        throw std::runtime_error(msg);
    }

    if ( step <= 0.0 )
    {
        std::string msg = "in function \"stranghLineSearch\" step value is less than zero.";
        throw std::invalid_argument(msg);
    }

    double fstep = fabs(step);

    a = x - fstep;
    b = x + fstep;
    fxa = mfunction->fx(a);        fx_count++;
    fxb = mfunction->fx(b);        fx_count++;
    double fxx = mfunction->fx(x); fx_count++;

    if (mcallback) mcallback->straightLineSearchCallback(iteration, x, a, b, fxa, fxb, fx_count);

    // if fxa and fxb are both greater than fxx, then minimum point is inside a and b
    if (fxa >= fxx && fxb >= fxx) return;

    // if fxa and fxb are both lesser than fxx then there is not minimum. function is not unimodal
    if (fxa < fxx && fxx > fxb)
    {
        // if fxa equals fxb
        if (fabs(fxa - fxb) <= DBL_EPSILON)
        {
            //puts("y1 = y2");
            //a = x - fstep;
            //b = x + fstep;
            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
            //return;
        }
        else if (fxa < fxb)
        {
            //puts("y1 < y2");
            //a = x - fstep;
            b = x;
            fxb = mfunction->fx(b); fx_count++;

            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
            //return;
        }
        else if (fxa > fxb)
        {
            //puts("y1 > y2");
            a = x;
            //b = x + fstep;
            fxa = mfunction->fx(a); fx_count++;

            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
            //return;
        }
        //fprintf(stderr, "%f %f %f %f %f %f\n", fxa, fxx, fxb, x - fstep, x, x + fstep);
        fprintf(stderr, "%4d %8.4f %8.4f %8.4f %10.6f %10.6f %10.6f\n", -2, a, x, b, fxa, fxx, fxb);
        fputs("Function is not unimodal\n", stderr);
        return;
    }

    // from left to right -->>
    if ( fxa >= fxx && fxx >= fxb )
    {
        puts("---------------->>");
        //straightLineSearchCallback(iteration, x, a, b, fxa, fxb);
        while ( fxx > fxb)
        {
            x = x + fstep;
            a = x - fstep;
            b = x + fstep;

            fxa = fxx;
            fxx = fxb;
            fxb = mfunction->fx(b); fx_count++;

            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        return;
    }

    // from right to left <<--
    if ( fxa <= fxx && fxx <= fxb )
    {
        puts("<<----------------");
        //straightLineSearchCallback(iteration, x, a, b, fxa, fxb);
        while ( fxa < fxx )
        {
            x = x - fstep;
            a = x - fstep;
            b = x + fstep;

            fxb = fxx;
            fxx = fxa;
            fxa = mfunction->fx(a); fx_count++;

            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        return;
    }
}

/**
 * @brief          Метод Свенна. Установления границ интервала.
 * @param x        Произвольно выбранная начальная точка
 * @param step     Величина шага
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param fxa      Величина функции в точке a
 * @param fxb      Величина функции в точке b
 */
void R1FxMinimizer::swann(double x, double step, double &a, double &b, double &fxa, double &fxb) const
{
    unsigned int iteration = 0;
    unsigned int fx_count = 0;

    if ( mfunction == NULL )
    {
        std::string msg = "in function \"swann\" function pointer is null.";
        throw std::runtime_error(msg);
    }

    if ( step <= 0.0 )
    {
        std::string msg = "in function \"swann\" step value is less than zero.";
        throw std::invalid_argument(msg);
    }

    double fstep = fabs(step);

    a = x - step;
    b = x + step;
    fxa = mfunction->fx(a);        fx_count++;
    fxb = mfunction->fx(b);        fx_count++;
    double fxx = mfunction->fx(x); fx_count++;

    if (mcallback) mcallback->straightLineSearchCallback(iteration, x, a, b, fxa, fxb, fx_count);

    // if fxa and fxb are both greater than fxx, then minimum point is inside a and b
    // начальный интервал неопределенности найден
    if (fxa >= fxx && fxb >= fxx) return;

    // if fxa and fxb are both lesser than fxx then there is not minimum. function is not unimodal
    // функция не является
    if (fxa < fxx && fxx > fxb)
    {
        // if fxa equals fxb
        if (fabs(fxa - fxb) <= DBL_EPSILON)
        {
            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        else if (fxa < fxb)
        {
            b = x;
            fxb = mfunction->fx(b); fx_count++;

            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        else if (fxa > fxb)
        {
            a = x;
            fxa = mfunction->fx(a); fx_count++;

            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        //fprintf(stderr, "%f %f %f %f %f %f\n", fxa, fxx, fxb, x - fstep, x, x + fstep);
        fprintf(stderr, "%4d %8.4f %8.4f %8.4f %10.6f %10.6f %10.6f\n", -2, a, x, b, fxa, fxx, fxb);
        fputs("Function is not unimodal\n", stderr);
        return;
    }

    // from left to right -->>
    if ( fxa >= fxx && fxx >= fxb )
    {
        puts("---------------->>");
        unsigned int k = 1;
        while ( fxx > fxb )
        {
            a = x;
            x = b;
            b = x + pow(2.0, k) * fstep;
            k++;

            fxa = fxx;
            fxx = fxb;
            fxb = mfunction->fx(b); fx_count++;

            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        return;
    }

    // from right to left <<--
    if ( fxa <= fxx && fxx <= fxb )
    {
        puts("<<----------------");
        unsigned int k = 1;
        while ( fxa < fxx )
        {
            b = x;
            x = a;
            a = x - pow(2.0, k) * fstep;
            k++;

            fxb = fxx;
            fxx = fxa;
            fxa = mfunction->fx(a); fx_count++;

            if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        return;
    }
}

/**
 * @brief          Метод золотого сечения.
 *                 Метод относится к последовательным стратегиям. Задается начальный интервал неопределенности и
 *                 требуемая точность. Алгоритм уменьшения интервала опирается на анализ значений функции в двух точках.
 *                 В качестве точек вычисления функции выбираются точки золотого сечения. Тогда с учетом свойств золотого
 *                 сечения на каждой итерации, кроме первой, требуется только одно новое вычисление функции. Условия
 *                 окончания процесса поиска стандартные: поиск заканчивается, когда длина текущего интервала
 *                 неопределенности оказывается меньше установленной величины.
 * @param x        Точка минимума
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param epsilon  Число эпсилон для останова метода
 */
void R1FxMinimizer::goldenSectionSearch(double &x, double &a, double &b, double epsilon) const
{
    unsigned int iteration = 0;
    unsigned int fx__count = 0;

    if ( mfunction == NULL )
    {
        std::string msg = "in function \"goldenSectionSearch\" function pointer is null.";
        throw std::runtime_error(msg);
    }

    if ( epsilon <= DBL_EPSILON )
    {
        std::string msg = "in function \"goldenSectionSearch\" epsilon value is equal to or less than zero.";
        throw std::invalid_argument(msg);
    }

    double phi = 1.6180339887498948482045868343656;
    //double phi = 0.6180339887498948482045868343656;

    bool check_x1 = true;
    bool check_x2 = true;

    double x1 = 0.0;
    double x2 = 0.0;

    double y1 = 0.0;
    double y2 = 0.0;

    x = (a+b)/2.0;
    if (mcallback) mcallback->goldenSectionSearchCallback(iteration, x, a, b, NAN, NAN, NAN, fx__count);

    // Lazimi epsilon deqiqliyini alana qeder iterasiyalari davam edirik
    while ( fabs(b-a) > epsilon )
    {
        if (check_x1)
        {
            x1 = b - fabs(b-a)/phi;
            y1 = mfunction->fx(x1); fx__count++;
        }

        if (check_x2)
        {
            x2 = a + fabs(b-a)/phi;
            y2 = mfunction->fx(x2); fx__count++;
        }

        check_x1 = check_x2 = false;

        if (y1 >= y2)
        {
            a = x1;
            x1 = x2;         // Tapilmish x2 noqtesi ve bu noqtede funksiyanin qiymeti
            y1 = y2;         // sonraki iterasiyada x1 qiymeti kimi istifade olunacaq.
            check_x2 = true; // x2 novbeti iterasiyada axtarilacaq
        }
        else
        {
            b = x2;
            x2 = x1;         // Tapilmish x1 noqtesi ve bu noqtede funksiyanin qiymeti
            y2 = y1;         // sonraki iterasiyada x2 qiymeti kimi istifade olunacaq.
            check_x1 = true; // x1 novbeti iterasiyada axtarilacaq
        }

        iteration++;
        x = (b+a)/2.0;
        if (mcallback) mcallback->goldenSectionSearchCallback(iteration, x, a, b, NAN, NAN, NAN, fx__count);
    }

    //double c = (a+b)/2.0;
    //double fa = mfunction->fx(a);
    //double fb = mfunction->fx(b);
    //if (fa<fb)  c = a;
    //if (fa>fb)  c = b;
    //x = c;
    //return c;
}

/**
 * @brief          Метод деления интервала пополам.
 *                 Метод относится к последовательным стратегиям и позволяет исключить из дальнейшего
 *                 рассмотрения на каждой итерации в точности половину текущего интервала неопределенности.
 *                 Задается начальный интервал неопределенности, а алгоритм уменьшения интервала, являясь, как и
 *                 в общем случае, "гарантирующим", основан на анализе величин функции в трех точках, равномерно
 *                 распределенных на текущем интервале (делящих его на четыре равные части). Условия окончания
 *                 процесса поиска стандартные: поиск заканчивается, когда длина текущего интервала неопределенности
 *                 оказывается меньше установленной величины.
 * @param x        Точка минимума
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param epsilon  Число эпсилон для останова метода
 */
void R1FxMinimizer::halphIntervalMethod(double &x, double &a, double &b, double epsilon) const
{
    unsigned int iteration = 0;
    unsigned int fx__count = 0;

    if ( mfunction == NULL )
    {
        std::string msg = "in function \"halphIntervalMethod\" function pointer is null.";
        throw std::runtime_error(msg);
    }

    if ( epsilon <= DBL_EPSILON )
    {
        std::string msg = "in function \"halphIntervalMethod\" epsilon value is equal to or less than zero.";
        throw std::invalid_argument(msg);
    }

    double l = fabs(b - a);
    x = (a+b)/2.0;
    double f_xm = mfunction->fx(x); fx__count++;

    if (mcallback) mcallback->goldenSectionSearchCallback(iteration, x, a, b, NAN, NAN, NAN, fx__count);

    while ( l > epsilon )
    {
        double xm = (a+b)/2.0;
        double x1 = a + l/4.0;
        double x2 = b - l/4.0;

        double f_x1 = mfunction->fx(x1); fx__count++;
        double f_x2 = mfunction->fx(x2); fx__count++;

        if ( f_x1 < f_xm )
        {
            b = xm;
            xm = x1;
            f_xm = f_x1;
        }
        else if ( f_x2 < f_xm )
        {
            a = xm;
            xm = x2;
            f_xm = f_x2;
        }
        else
        {
            a = x1;
            b = x2;
        }
        l = fabs(b - a);

        iteration++;
        x = (b+a)/2.0;
        if (mcallback) mcallback->goldenSectionSearchCallback(iteration, x, a, b, NAN, NAN, NAN, fx__count);
    }

    //x = (b + a) /  2.0;
    //double fa = f->fx(a);
    //double fb = f->fx(b);

    //if (fa<fb) x = a;
    //if (fa>fb) x = b;
}

/**
 * @brief          Метод равномерного поиска.
 *                 Метод относится к пассивным стратегиям. Задается начальный интервал
 *                 неопределенности [a, b] и количество вычислений функции n.
 *                 Вычисления производятся в n равноотстоящих друг от друга точках (при этом
 *                 интервал делится на n + 1 равных интервалов). Путем сравнения величин
 *                 f(xi), i = 1,...,n находится точка xк, в которой значение функции наименьшее.
 *                 Искомая точка минимума х* считается заключенной в интервале [хk-1, хk+1]
 * @param x        Точка минимума
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param n        Количество вычислений функции
 */
void R1FxMinimizer::uniformLineSearch(double &x, double &a, double &b, unsigned int n) const
{
    //unsigned int iteration = 0;
    unsigned int fx__count = 0;

    if ( mfunction == NULL )
    {
        std::string msg = "in function \"uniformLineSearch\" function pointer is null.";
        throw std::runtime_error(msg);
    }

//    if ( epsilon <= DBL_EPSILON )
//    {
//        std::string msg = "in function \"uniformLineSearch\" epsilon value is equal to or less than zero.";
//        throw std::invalid_argument(msg);
//    }

//    if ( step <= 0.0 )
//    {
//        std::string msg = "in function \"uniformLineSearch\" step value is less than zero.";
//        throw std::invalid_argument(msg);
//    }

    double step = fabs(b - a) / n;
    DoubleVector x0(n+1);
    DoubleVector y0(n+1);

    for (unsigned int i=1; i<=n-1; i++)
    {
        x0[i] = a + i * step;
        y0[i] = mfunction->fx(x0[i]);
        fx__count++;
    }

    unsigned int k=1;
    double min = y0[k];
    for (unsigned int i=1; i<=n-1; i++)
    {
        if (y0[i] < min) { k = i; min = y0[i]; }
    }

    a = x0[k-1];
    b = x0[k+1];
    x = x0[k];

    x0.clear();
    y0.clear();
}

/**
 * @brief R1Minimize::setCallback
 * @param callback
 */
void R1FxMinimizer::setCallback(R1FxMinimizer::Callback *callback)
{
    this->mcallback = callback;
    this->mcallback->mfunction = mfunction;
}

/**
 * @brief R1Minimize::callback
 * @return
 */
R1FxMinimizer::Callback* R1FxMinimizer::callback() const
{
    return this->mcallback;
}

/**
 * @brief R1Minimize::function
 * @return
 */
R1Function* R1FxMinimizer::function() const
{
    return mfunction;
}

/**
 * @brief R1Minimize::setFunction
 * @param function
 */
void R1FxMinimizer::setFunction(R1Function *function)
{
    if (this->mcallback) mcallback->mfunction = function;
    this->mfunction = function;
}

/**
 * @brief R1Minimize::Callback::straightLineSearchCallback
 * @param iteration
 * @param x
 * @param a
 * @param b
 * @param fxa
 * @param fxb
 */
void R1FxMinimizer::Callback::straightLineSearchCallback(unsigned int iteration, double x, double a, double b, double fxa, double fxb, unsigned int fx_count) const
{
    C_UNUSED(iteration);
    C_UNUSED(a);
    C_UNUSED(b);
    C_UNUSED(fxa);
    C_UNUSED(fxb);
    C_UNUSED(fx_count);
}

/**
 * @brief R1Minimize::Callback::goldenSectionSearchCallback
 * @param iteration
 * @param x
 * @param a
 * @param b
 * @param fxx
 * @param fxa
 * @param fxb
 * @param fx_count
 */
void R1FxMinimizer::Callback::goldenSectionSearchCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const
{
    C_UNUSED(iteration);
    C_UNUSED(x);
    C_UNUSED(a);
    C_UNUSED(b);
    C_UNUSED(fxx);
    C_UNUSED(fxa);
    C_UNUSED(fxb);
    C_UNUSED(fx_count);
}

/**
 * @brief R1Minimize::Callback::function
 * @return
 */
R1Function* R1FxMinimizer::Callback::function() const
{
    return mfunction;
}

/**
 * @brief         Методы прямого поиска. Установления границ интервала.
 * @param x       Произвольно выбранная начальная точка
 * @param step    Величина шага
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param f       Целевая функция
 * @param fa      Величина функции в точке a
 * @param fb      Величина функции в точке b
 */
void stranghLineSearch(double x, double step, double &a, double &b, R1Function *f)
{
    if ( f == NULL )
    {
        std::string msg = "in function \"stranghLineSearch\" function pointer is null.";
        throw std::runtime_error(msg);
    }

    if ( step <= 0.0 )
    {
        std::string msg = "in function \"stranghLineSearch\" step value is less than zero.";
        throw std::invalid_argument(msg);
    }

    double fstep = fabs(step);

    double y0 = f->fx(x);
    double y1 = f->fx(x - fstep);
    double y2 = f->fx(x + fstep);

    // if y1 and y2 are both greater than y0 then minimum point is inside x1 and x2
    if (y1 >= y0 && y0 <= y2)
    {
        a = x - fstep;
        b = x + fstep;
        return;
    }

    // if y1 and y2 are both lesser than y0 then there is not minimum. function is not unimodal
    if (y1 <= y0 && y0 >= y2)
    {
        if (y1 <= y2)
        {
            a = x - fstep;
            b = x;
        }
        if (y1 >= y2)
        {
            a = x;
            b = x + fstep;
        }
        //a = x - fstep;
        //b = x + fstep;
        //        a = b = NAN;
        printf("%f %f %f %f %f %f\n", y1, y0, y2, x - fstep, x, x + fstep);
        fputs("Function is not unimodal\n", stderr);
        return;
    }

    if ( y1 >= y0 && y0 >= y2 )
    {
        while ( y2 < y0 )
        {
            x = x + fstep;
            y1 = y0;
            y0 = y2;
            y2 = f->fx(x + fstep);
        }
        a = x - fstep;
        b = x + fstep;
        return;
    }

    if ( y1 <= y0 && y0 <= y2 )
    {
        while ( y1 < y0 )
        {
            x = x - fstep;
            y2 = y0;
            y0 = y1;
            y1 = f->fx(x - fstep);
        }
        a = x - fstep;
        b = x + fstep;
        return;
    }

    a = x - fstep;
    b = x + fstep;
    printf("Test. stranghLineSearch %f %f %f %f %f %f\n", a, x, b, y1, y0, y2);
}

/**
 * @brief          Метод золотого сечения.
 *                 Метод относится к последовательным стратегиям. Задается начальный интервал неопределенности и
 *                 требуемая точность. Алгоритм уменьшения интервала опирается на анализ значений функции в двух точках.
 *                 В качестве точек вычисления функции выбираются точки золотого сечения. Тогда с учетом свойств золотого
 *                 сечения на каждой итерации, кроме первой, требуется только одно новое вычисление функции. Условия
 *                 окончания процесса поиска стандартные: поиск заканчивается, когда длина текущего интервала
 *                 неопределенности оказывается меньше установленной величины.
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param x
 * @param f        Целевая функция
 * @param epsilon  Число эпсилон для останова метода
 * @return
 */
double goldenSectionSearch(double &a, double &b, double &x, R1Function *f, double epsilon)
{
    if ( f == NULL )
    {
        std::string msg = "in function \"goldenSectionSearch\" function pointer is null.";
        throw std::runtime_error(msg);
    }

    if ( epsilon <= 0.0 )
    {
        std::string msg = "in function \"goldenSectionSearch\" epsilon value is equal to or less than zero.";
        throw std::invalid_argument(msg);
    }

    double phi = 1.6180339887498948482045868343656;

    double x1 = NAN;
    double x2 = NAN;

    double y1 = 0.0;
    double y2 = 0.0;

    // Lazimi epsilon deqiqliyini alana qeder iterasiyalari davam edirik
    while ( fabs(b-a) > epsilon )
    {
        if (isnan(x1))
        {
            x1 = b - fabs(b-a)/phi;
            y1 = f->fx(x1);
        }

        if (isnan(x2))
        {
            x2 = a + fabs(b-a)/phi;
            y2 = f->fx(x2);
        }

        if (y1 >= y2)
        {
            a = x1;
            x1 = x2;    // Tapilmish x2 noqtesi ve bu noqtede funksiyanin qiymeti
            y1 = y2;    // sonraki iterasiyada x1 qiymeti kimi istifade olunacaq.
            x2 = NAN;   // x2 novbeti iterasiyada axtarilacaq
        }
        else
        {
            b = x2;
            x2 = x1;    // Tapilmish x1 noqtesi ve bu noqtede funksiyanin qiymeti
            y2 = y1;    // sonraki iterasiyada x2 qiymeti kimi istifade olunacaq.
            x1 = NAN;   // x1 novbeti iterasiyada axtarilacaq
        }
    }

    double c = (a+b)/2.0;
    double fa = f->fx(a);
    double fb = f->fx(b);
    if (fa<fb)  c = a;
    if (fa>fb)  c = b;
    x = c;
    return c;
}

double goldenSectionSearch2(double &a, double &b, double &x, R1Function *f, double epsilon)
{
    double phi = 1.6180339887498948482045868343656;
    double x1 = b - fabs(b-a)/phi;
    double x2 = a + fabs(b-a)/phi;

    double y1 = 0.0;//f->fx(x1);
    double y2 = 0.0;//f->fx(x2);

    bool check1 = true;
    bool check2 = true;

    while (fabs(b-a) > epsilon)
    {
        if (check1)
        {
            x1 = b - fabs(b-a)/phi;
            y1 = f->fx(x1);
            check1 = false;
        }
        if (check2)
        {
            x2 = a + fabs(b-a)/phi;
            y2 = f->fx(x2);
            check2 = false;
        }

        if (y1 >= y2)
        {
            a = x1;
            x1 = x2;
            y1 = y2;
            //x2 = a + fabs(b-a)/phi;
            //y2 = f->fx(x2);
            check2 = true;
        }
        else
        {
            b = x2;
            x2 = x1;
            y2 = y1;
            //x1 = b - fabs(b-a)/phi;
            //y1 = f->fx(x1);
            check1 = true;
        }
    }

    double c = (a+b)/2.0;
    double fa = f->fx(a);
    double fb = f->fx(b);
    if (fa<fb) c = a;
    if (fa>fb) c = b;
    x = c;
    return c;
}

double goldenSectionSearch1(double &a, double &b, double &x, R1Function *f, double epsilon)
{
    double phi = 0.38196601125010515179541316563436;
    double x1 = a + phi * (b-a);
    double x2 = a + b  - x1;

    double fx1 = f->fx(x1);
    double fx2 = f->fx(x2);

    while (fabs(b-a) > epsilon)
    {
        if (fx1 <= fx2)
        {
            b = x2;
            x2 = x1;
            x1 = a + b - x1;
            fx2 = fx1;
            fx1 = f->fx(x1);
        }
        else
        {
            a = x1;
            x1 = x2;
            x2 = a + b - x2;
            fx1 = fx2;
            fx2 = f->fx(x2);
        }
    }

    double c = (a+b)/2.0;
    double fa = f->fx(a);
    double fb = f->fx(b);
    if (fa<fb) c = a;
    if (fa>fb) c = b;
    x = c;
    return c;
}

R1Minimize1::R1Minimize1() : m_f(0), m_x0(0.0), m_step(0.0), m_epsilon(0.0), m_a(0.0), m_b(0.0)
{}

R1Minimize1::~R1Minimize1()
{}

void R1Minimize1::setFunction(R1Function *f)
{
    m_f = f;
}

R1Function* R1Minimize1::function() const
{
    return m_f;
}

double R1Minimize1::x0() const
{
    return m_x0;
}

void R1Minimize1::setX0(double x0)
{
    m_x0 = x0;
}

double R1Minimize1::step() const
{
    return m_step;
}

void R1Minimize1::setStep(double step)
{
    m_step = step;
}

double R1Minimize1::epsilon() const
{
    return m_epsilon;
}

void R1Minimize1::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

double R1Minimize1::a() const
{
    return m_a;
}

double R1Minimize1::b() const
{
    return m_b;
}

double R1Minimize1::straightLineSearch()
{
    if ( m_f == NULL ) return NAN;
    if ( m_step == 0.0 ) return NAN;

    double y0 = m_f->fx(m_x0);
    double y1 = m_f->fx(m_x0 - m_step);
    double y2 = m_f->fx(m_x0 + m_step);

    // if y1 and y2 are both greater than y0 then minimum point is inside x1 and x2
    if (y1 >= y0 && y0 <= y2)
    {
        m_a = m_x0 - fabs(m_step);
        m_b = m_x0 + fabs(m_step);
        double c = (m_a + m_b) / 2.0;
        return c;
    }

    // if y1 and y2 are both lesser than y0 then there is not minimum. function is not unimodal
    if (y1 <= y0 && y0 >= y2)
    {
        fputs("Function is not unimodal", stderr);
        return NAN;
    }

    {
        if ( y1 >= y0 && y0 >= y2 )
        {
            while ( y2 < y0 )
            {
                m_x0 = m_x0 + fabs(m_step);
                y1 = y0;
                y0 = y2;
                y2 = m_f->fx(m_x0 + fabs(m_step));
            }
            m_a = m_x0 - fabs(m_step);
            m_b = m_x0 + fabs(m_step);
            double c = (m_a + m_b) / 2.0;
            return c;
        }

        if ( y1 <= y0 && y0 <= y2 )
        {
            while ( y1 < y0 )
            {
                m_x0 = m_x0 - fabs(m_step);
                y2 = y0;
                y0 = y1;
                y1 = m_f->fx(m_x0 - fabs(m_step));
            }
            m_a = m_x0 - fabs(m_step);
            m_b = m_x0 + fabs(m_step);
            double c = (m_a + m_b) / 2.0;
            return c;
        }
    }

    m_a = m_x0 - fabs(m_step);
    m_b = m_x0 + fabs(m_step);
    double c = (m_a + m_b) / 2.0;
    return c;
}

double R1Minimize1::halphIntervalMethod()
{
    if ( m_f == NULL ) return NAN;

    double L = m_b - m_a;

    while ( L > m_epsilon )
    {
        double xm = (m_a+m_b)/2.0;
        double x1 = m_a + L/4.0;
        double x2 = m_b - L/4.0;

        double f_xm = m_f->fx(xm);
        double f_x1 = m_f->fx(x1);
        double f_x2 = m_f->fx(x2);

        if (f_x1 < f_xm)
        {
            m_b = xm;
            xm = x1;
        }
        else
        {
            if ( f_x2 < f_xm )
            {
                m_a = xm;
                xm = x2;
            }
            else
            {
                m_a = x1;
                m_b = x2;
            }
        }
        L = m_b - m_a;
    }

    return L;
}

double R1Minimize1::goldenSectionSearch()
{
    double a = m_a;
    double b = m_b;
    //double sqrt_5 = 2.2360679774997896964091736687313
    //double phi = (sqrt(5) + 1.0) / 2.0;
    //double phi = (sqrt(5) - 1) / 2.0;
    double phi = 1.6180339887498948482045868343656;

    double x1 = NAN;
    double x2 = NAN;

    double y1 = 0.0;
    double y2 = 0.0;

    // Lazimi epsilon deqiqliyini alana qeder
    // iterasiyalari davam edirik
    while ( fabs(b-a) > m_epsilon )
    {
        if (isnan(x1))
        {
            x1 = b - fabs(b-a)/phi;
            y1 = m_f->fx(x1);
        }

        if (isnan(x2))
        {
            x2 = a + fabs(b-a)/phi;
            y2 = m_f->fx(x2);
        }

        if (y1 >= y2)
        {
            a = x1;
            x1 = x2;    // Tapilmish x2 noqtesi ve bu noqtede funkisiyanin qiymeti
            y1 = y2;    // sonraki iterasiyada x1 qiymeti kimi istifade olunacaq.
            x2 = NAN;   // x2 novbeti iterasiyada axtarilacaq
        }
        else
        {
            b = x2;
            x2 = x1;    // Tapilmish x1 noqtesi ve bu noqtede funkisiyanin qiymeti
            y2 = y1;    // sonraki iterasiyada x2 qiymeti kimi istifade olunacaq.
            x1 = NAN;   // x1 novbeti iterasiyada axtarilacaq
        }
    }

    double c = (a+b)/2.0;

    if (m_f->fx(a)<m_f->fx(b)) c = a;
    if (m_f->fx(a)>m_f->fx(b)) c = b;

    return c;
}

/**
 * @brief         Методы прямого поиска. Установления границ интервала.
 * @param x       Произвольно выбранная начальная точка
 * @param step    Величина шага
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param f       Целевая функция
 */
void R1Minimize1::StranghLineSearch(double x, double step, double &a, double &b, R1Function *f)
{
    if ( f == NULL ) return;
    if ( step == 0.0 ) return;

    double fstep = fabs(step);

    double y0 = f->fx(x);
    double y1 = f->fx(x - fstep);
    double y2 = f->fx(x + fstep);

    // if y1 and y2 are both greater than y0 then minimum point is inside x1 and x2
    if (y1 >= y0 && y0 <= y2)
    {
        a = x - fstep;
        b = x + fstep;
        return;
    }

    // if y1 and y2 are both lesser than y0 then there is not minimum. function is not unimodal
    if (y1 <= y0 && y0 >= y2)
    {
        a = b = NAN;
        fputs("Function is not unimodal\n", stderr);
        return;
    }

    if ( y1 >= y0 && y0 >= y2 )
    {
        while ( y2 < y0 )
        {
            x = x + fstep;
            y1 = y0;
            y0 = y2;
            y2 = f->fx(x + fstep);
        }
        a = x - fstep;
        b = x + fstep;
    }

    if ( y1 <= y0 && y0 <= y2 )
    {
        while ( y1 < y0 )
        {
            x = x - fstep;
            y2 = y0;
            y0 = y1;
            y1 = f->fx(x - fstep);
        }
        a = x - fstep;
        b = x + fstep;
    }

    a = x - fstep;
    b = x + fstep;
}

/**
 * @brief          Метод Свенна. Установления границ интервала.
 * @param x        Произвольно выбранная начальная точка
 * @param step     Величина шага
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param f        Целевая функция
 */
void R1Minimize1::Swann(double x, double step, double &a, double &b, R1Function *f)
{
    step = fabs(step);
    int k = 0;
    a = NAN;
    b = NAN;
    //double h = dx;

    if ( f == NULL ) return;
    if ( step  == 0.0  ) return;

    double y1 = f->fx(x - step);
    double y0 = f->fx(x);
    double y2 = f->fx(x + step);

    /* начальный интервал неопределенности найден */
    if ( y1 >= y0 && y0 <= y2 )
    {
        a = x - step;
        b = x + step;
        return;
    }

    /* функция не является  */
    if ( y1 <= y0 && y0 >= y2 )
    {
        fputs("Function is not unimodal\n", stderr);
        a = b = NAN;
        return;
    }

    if ( y1 >= y0 && y0 >= y2 )
    {
        step = fabs(step);
        a = x;
        x = x + step;
        k = 1;

        while ( f->fx(x+pow(2.0, k)*step) < f->fx(x) )
        {
            a = x;
            x = x + pow(2, k) * step;
            k++;
        }
        b = x+pow(2.0, k)*step;
    }

    if ( y1 <= y0 && y0 <= y2 )
    {
        step = -fabs(step);
        b = x;
        x = x - fabs(step);
        k = 1;

        while ( f->fx(x+pow(2.0, k)*step) < f->fx(x) )
        {
            b = x;
            x = x + pow(2, k) * step;
            k++;
        }
        a = x+pow(2.0, k)*step;
    }
}

/**
 * @brief          Метод золотого сечения.
 *                 Метод относится к последовательным стратегиям. Задается начальный интервал неопределенности и
 *                 требуемая точность. Алгоритм уменьшения интервала опирается на анализ значений функции в двух точках.
 *                 В качестве точек вычисления функции выбираются точки золотого сечения. Тогда с учетом свойств золотого
 *                 сечения на каждой итерации, кроме первой, требуется только одно новое вычисление функции. Условия
 *                 окончания процесса поиска стандартные: поиск заканчивается, когда длина текущего интервала
 *                 неопределенности оказывается меньше установленной величины.
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param x
 * @param f        Целевая функция
 * @param epsilon  Число эпсилон для останова метода
 * @return
 */
double R1Minimize1::GoldenSectionSearch(double &a, double &b, double &x, R1Function *f, double epsilon)
{
    double phi = 1.6180339887498948482045868343656;

    double x1 = NAN;
    double x2 = NAN;

    double y1 = 0.0;
    double y2 = 0.0;

    // Lazimi epsilon deqiqliyini alana qeder
    // iterasiyalari davam edirik
    while ( fabs(b-a) > epsilon )
    {
        if (isnan(x1))
        {
            x1 = b - fabs(b-a)/phi;
            y1 = f->fx(x1);
        }

        if (isnan(x2))
        {
            x2 = a + fabs(b-a)/phi;
            y2 = f->fx(x2);
        }

        if (y1 >= y2)
        {
            a = x1;
            x1 = x2;    // Tapilmish x2 noqtesi ve bu noqtede funkisiyanin qiymeti
            y1 = y2;    // sonraki iterasiyada x1 qiymeti kimi istifade olunacaq.
            x2 = NAN;   // x2 novbeti iterasiyada axtarilacaq
        }
        else
        {
            b = x2;
            x2 = x1;    // Tapilmish x1 noqtesi ve bu noqtede funkisiyanin qiymeti
            y2 = y1;    // sonraki iterasiyada x2 qiymeti kimi istifade olunacaq.
            x1 = NAN;   // x1 novbeti iterasiyada axtarilacaq
        }
    }

    double c = (a+b)/2.0;

    if (f->fx(a)<f->fx(b)) c = a;
    if (f->fx(a)>f->fx(b)) c = b;

    x = c;
    return c;
}

/**
 * @brief          Метод Фибоначчи.
 *                 Метод относится к последовательным стратегиям. Задается начальный интервал непределенности и
 *                 количество N вычислений функции. Алгоритм уменьшения интервала опирается на анализ значений
 *                 функции в двух точках. Точки вычисления функции находятся с использованием последовательности
 *                 из n + 1 чисел Фибоначчи. Как в методе золотого сечения, на первой итерации требуются два
 *                 вычисления функции, а на каждой последующей - только по одному. Условия окончания процесса
 *                 поиска стандартные: поиск заканчивается, когда длина текущего интервала неопределенности
 *                 оказывается меньше установленной величины.
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param x
 * @param f        Целевая функция
 * @param epsilon  Число эпсилон для останова метода
 * @return
 */
void R1Minimize1::FibonachiMethod(double &a, double &b, double &c, double step, double epsilon, R1Function *f)
{
    C_UNUSED(a);
    C_UNUSED(b);
    C_UNUSED(c);
    C_UNUSED(step);
    C_UNUSED(epsilon);
    C_UNUSED(f);

    //    double k = fabs(b-a)/epsilon;
    //    unsigned int FN = (unsigned int)ceil(k);
}

/**
 * @brief          Метод деления интервала пополам.
 *                 Метод относится к последовательным стратегиям и позволяет исключить из дальнейшего
 *                 рассмотрения на каждой итерации в точности половину текущего интервала неопределенности.
 *                 Задается начальный интервал неопределенности, а алгоритм уменьшения интервала, являясь, как и
 *                 в общем случае, "гарантирующим", основан на анализе величин функции в трех точках, равномерно
 *                 распределенных на текущем интервале (делящих его на четыре равные части). Условия окончания
 *                 процесса поиска стандартные: поиск заканчивается, когда длина текущего интервала неопределенности
 *                 оказывается меньше установленной величины.
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param x
 * @param f        Целевая функция
 * @param epsilon  Число эпсилон для останова метода
 * @return
 */
double R1Minimize1::HalphIntervalMethod(double &a, double &b, double &x, R1Function *f, double epsilon)
{
    if ( f == NULL ) return NAN;

    double L = b - a;

    while ( L > epsilon )
    {
        double xm = (a+b)/2.0;
        double x1 = a + L/4.0;
        double x2 = b - L/4.0;

        double f_xm = f->fx(xm);
        double f_x1 = f->fx(x1);
        double f_x2 = f->fx(x2);

        if (f_x1 < f_xm)
        {
            b = xm;
            xm = x1;
        }
        else
        {
            if ( f_x2 < f_xm )
            {
                a = xm;
                xm = x2;
            }
            else
            {
                a = x1;
                b = x2;
            }
        }
        L = b - a;
    }

    x = (b + a) /  2.0;
    double fa = f->fx(a);
    double fb = f->fx(b);

    if (fa<fb) x = a;
    if (fa>fb) x = b;

    return L;
}

/**
 * @brief          Метод равномерного поиска.
 *                 Метод относится к пассивным стратегиям. Задается начальный интервал
 *                 неопределенности [a, b] и количество вычислений функции n.
 *                 Вычисления производятся в n равноотстоящих друг от друга точках (при этом
 *                 интервал делится на n + 1 равных интервалов). Путем сравнения величин
 *                 f(xi), i = 1,...,n находится точка xк, в которой значение функции наименьшее.
 *                 Искомая точка минимума х* считается заключенной в интервале [хk-1, хk+1]
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param n        Количество вычислений функции
 * @param f        Целевая функция
 */
void R1Minimize1::UniformLineSearch(double &a, double &b, double &c, unsigned int n, R1Function *f)
{
    if ( f == NULL ) return;

    double step = (b - a) / (n+1);
    DoubleVector x(n);
    DoubleVector y(n);

    for (unsigned int i=1; i<=n; i++)
    {
        x[i-1] = a + i * step;
        y[i-1] = f->fx(x[i-1]);
    }

    unsigned int k=0;
    double min = y[k];
    for (unsigned int i=1; i<n; i++)
    {
        if (y[i] < min) { k = i; min = y[i]; }
    }

    a = x[k-1];
    b = x[k+1];
    c = x[k];

    x.clear();
    y.clear();
}

/**
 * @brief          Метод перебора.
 *                 Метод относится к пассивным стратегиям. Задается начальный интервал
 *                 неопределенности [a, b] и количество вычислений функции n.
 *                 Вычисления производятся в n равноотстоящих друг от друга точках (при этом
 *                 интервал делится на n + 1 равных интервалов). Путем сравнения величин
 *                 f(xi), i = 0,...,n+1 находится точка xк, в которой значение функции наименьшее.
 *                 Искомая точка минимума х* считается заключенной в интервале [хk-1, хk+1]
 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param n        Количество вычислений функции
 * @param f        Целевая функция
 */
void R1Minimize1::BruteForceLineSearch(double &a, double &b, double &c, unsigned int n, R1Function *f)
{
    if ( f == NULL ) return;

    double step = (b - a) / (n+1);
    DoubleVector x(n+2);
    DoubleVector y(n+2);

    for (unsigned int i=0; i<=n+1; i++)
    {
        x[i-1] = a + i * step;
        y[i-1] = f->fx(x[i-1]);
    }

    unsigned int k=0;
    double min = y[k];
    for (unsigned int i=1; i<=n+1; i++)
    {
        if (y[i] < min) { k = i; min = y[i]; }
    }

    if (k == 0)
    {
        a = x[0];
        b = x[1];
    }
    else if (k == n+1)
    {
        a = x[n];
        b = x[n+1];
    } else {
        a = x[k-1];
        b = x[k+1];
    }
    c = x[k];

    x.clear();
    y.clear();
}

/**
 * @brief          Метод дихотомии
 *                 Метод относится к последовательным стратегиям. Задается начальный интервал неопределенности и
 *                 требуемая точность. Алгоритм опирается на анализ значений функции в двух точках. Для их нахождения
 *                 текущий интервал неопределенности делится пополам и в обе стороны от середины откладывается по — step/2,
 *                 где step - малое положительное число. Условия окончания процесса поиска стандартные: поиск заканчивается,
 *                 когда длина текущего интервала неопределенности оказывается меньше установленной величины.

 * @param a        Начальная точка отрезка
 * @param b        Конечнная точка отрезка
 * @param c
 * @param f        Целевая функция
 * @param step     Малое число
 * @param epsilon  Число эпсилон для останова метода
 */
void R1Minimize1::DichotomyMethod(double &a, double &b, double &c, R1Function *f, double step, double epsilon)
{
    if ( f == NULL ) return;

    while ( fabs(b - a) > epsilon )
    {
        double y = (a + b - step) / 2.0;
        double z = (a + b + step) / 2.0;
        double f_y = f->fx(y);
        double f_z = f->fx(z);

        if (f_y <= f_z)
        {
            b = z;
        }
        else
        {
            a = y;
        }
    }
    c = (a + b)/2.0;
}
