#include "gradient.h"
#include "vectornormalizer.h"
#include "printer.h"
#include "function.h"
#include <math.h>
#include <float.h>

/******************************************** DefaultNormalizer ****************************************/

class DefaultNormalizer : public IVectorNormalizer
{
public:
    virtual ~DefaultNormalizer();
    inline virtual auto norm(const DoubleVector &v) const -> double { return EuclideanNorm(v); }
    inline virtual auto normalize(DoubleVector &v) const -> void { EuclideanNormalize(v); }
};

DefaultNormalizer::~DefaultNormalizer() {}

/******************************************** DefaultGradient *****************************************/

class DefaultGradient : public IGradient
{
public:
    virtual ~DefaultGradient();
    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;
};

DefaultGradient::~DefaultGradient() {}

auto DefaultGradient::gradient(const DoubleVector &, DoubleVector &) const -> void {}

/******************************************** DefaultPrinter *****************************************/

class DefaultPrinter : public IPrinter
{
public:
    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f,
                       double alpha, GradientMethod::MethodResult result) const -> void;
};

/*****************************************************************************************************/

GradientMethod::GradientMethod() : m_fn(NULL), m_gr(NULL), m_printer(NULL), m_projection(NULL),
    m_optimalityTolerance(0.1), m_functionTolerance(0.1), m_stepTolerance(0.1),
    min_step(0.1), min_epsilon(0.01),
    m_show_end_message(true), m_normalize(true), m_normalizer(NULL),
    m_iterationCount(0), m_maxIterations(UINT32_MAX)
{
    m_gr = new DefaultGradient();
    m_normalizer = new DefaultNormalizer();
}

GradientMethod::~GradientMethod() {}

/**
 * @brief Objective function
 * @param f
 */
void GradientMethod::setFunction(RnFunction *fn)
{
    m_fn = fn;
}

/**
 * @brief Objective function
 * @return
 */
RnFunction* GradientMethod::function() const
{
    return m_fn;
}

IGradient* GradientMethod::gradient() const
{
    return m_gr;
}

void GradientMethod::setGradient(IGradient *gr)
{
    m_gr = gr;
}

auto GradientMethod::optimalityTolerance() const -> double
{
    return m_optimalityTolerance;
}

auto GradientMethod::setOptimalityTolerance(double optimalityTolerance) -> void
{
    m_optimalityTolerance = optimalityTolerance;
}

auto GradientMethod::stepTolerance() const -> double
{
    return m_stepTolerance;
}

auto GradientMethod::setStepTolerance(double stepTolerance) -> void
{
    m_stepTolerance = stepTolerance;
}

auto GradientMethod::functionTolerance() const -> double
{
    return m_functionTolerance;
}

auto GradientMethod::setFunctionTolerance(double functionTolerance) -> void
{
    m_functionTolerance = functionTolerance;
}

auto GradientMethod::constraintTolerance() const -> double
{
    return m_constraintTolerance;
}

auto GradientMethod::setConstraintTolerance(double constraintTolerance) -> void
{
    m_constraintTolerance = constraintTolerance;
}

auto GradientMethod::setTolerance(double optimalityTolerance, double stepTolerance, double functionTolerance) -> void
{
    m_optimalityTolerance = optimalityTolerance;
    m_stepTolerance = stepTolerance;
    m_functionTolerance = functionTolerance;
}

void GradientMethod::setR1MinimizeEpsilon(double step, double epsilon)
{
    min_step = step;
    min_epsilon = epsilon;
}

unsigned int GradientMethod::count() const
{
    return m_iterationCount;
}

void GradientMethod::setPrinter(IPrinter *printer)
{
    this->m_printer = printer;
}

void GradientMethod::setProjection(IProjection *projection)
{
    this->m_projection = projection;
}

void GradientMethod::setNormalize(bool normalize)
{
    this->m_normalize = normalize;
}

void GradientMethod::showExitMessage(bool showEndMessage)
{
    m_show_end_message = showEndMessage;
}

void GradientMethod::setGradientNormalizer(IVectorNormalizer *normalizer)
{
    m_normalizer = normalizer;
}

auto GradientMethod::setMaxIterations(unsigned int maxIterations) -> void
{
    m_maxIterations = maxIterations;
}

auto GradientMethod::maxIterations() const -> unsigned int
{
    return m_maxIterations;
}

auto GradientMethod::setMaxFunctionEvaluations(unsigned int maxFunctionEvaluations) -> void
{
    m_maxFunctionEvaluations = maxFunctionEvaluations;
}

auto GradientMethod::maxFunctionEvaluations() const -> unsigned int
{
    return m_maxFunctionEvaluations;
}

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
void GradientMethod::straightLineSearch(double x, double step, double &a, double &b, double &fxa, double &fxb, bool &unimodal) const
{
    unsigned int iteration = 0;
    unsigned int fx_count = 0;

    if ( step <= 0.0 )
    {
        std::string msg = "in function \"stranghLineSearch\" step value is less than zero.";
        throw std::invalid_argument(msg);
    }

    double fstep = fabs(step);

    a = x - fstep;
    b = x + fstep;
    fxa = fx(a);        fx_count++;
    fxb = fx(b);        fx_count++;
    double fxx = fx(x); fx_count++;

    //if (mcallback) mcallback->straightLineSearchCallback(iteration, x, a, b, fxa, fxb, fx_count);

    // if fxa and fxb are both greater than fxx, then minimum point is inside a and b
    if (fxa >= fxx && fxb >= fxx)
    {
        unimodal = true;
        //if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        return;
    }

    if (fxa < fxx && fxx > fxb)
    {
        unimodal = false;
        //if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        //fprintf(stderr, "%4d %8.4f %8.4f %8.4f %10.6f %10.6f %10.6f\n", -2, a, x, b, fxa, fxx, fxb);
        //fputs("Function is not unimodal\n", stderr);
        return;
    }

    // from left to right -->>
    if ( fxa >= fxx && fxx >= fxb )
    {
        //puts("---------------->>");
        //straightLineSearchCallback(iteration, x, a, b, fxa, fxb);
        //fprintf(stderr, "%4d %8.4f %8.4f %8.4f %10.6f %10.6f %10.6f\n", iteration, a, x, b, fxa, fxx, fxb);
        while ( fxx > fxb)
        {
            x = x + fstep;
            a = x - fstep;
            b = x + fstep;

            fxa = fxx;
            fxx = fxb;
            fxb = fx(b); fx_count++;

            //fprintf(stderr, "%4d %8.4f %8.4f %8.4f %10.6f %10.6f %10.6f\n", ++iteration, a, x, b, fxa, fxx, fxb);
            //if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        unimodal = true;
        return;
    }

    // from right to left <<--
    if ( fxa <= fxx && fxx <= fxb )
    {
        //puts("<<----------------");
        //straightLineSearchCallback(iteration, x, a, b, fxa, fxb);
        //fprintf(stderr, "%4d %8.4f %8.4f %8.4f %10.6f %10.6f %10.6f\n", iteration, a, x, b, fxa, fxx, fxb);
        while ( fxa < fxx )
        {
            x = x - fstep;
            a = x - fstep;
            b = x + fstep;

            fxb = fxx;
            fxx = fxa;
            fxa = fx(a); fx_count++;

            //fprintf(stderr, "%4d %8.4f %8.4f %8.4f %10.6f %10.6f %10.6f\n", ++iteration, a, x, b, fxa, fxx, fxb);
            //if (mcallback) mcallback->straightLineSearchCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        unimodal = true;
        return;
    }
}

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
void GradientMethod::swann(double x, double step, double &a, double &b, double &fxa, double &fxb, bool &unimodal) const
{
    unsigned int iteration = 0;
    unsigned int fx_count = 0;

    if ( step <= 0.0 )
    {
        std::string msg = "in function \"swann\" step value is less than zero.";
        throw std::invalid_argument(msg);
    }

    double fstep = fabs(step);

    a = x - step;
    b = x + step;
    fxa = fx(a);        fx_count++;
    fxb = fx(b);        fx_count++;
    double fxx = fx(x); fx_count++;

    //if (mcallback) mcallback->swannCallback(iteration, x, a, b, fxa, fxb, fx_count);

    // if fxa and fxb are both greater than fxx, then minimum point is inside a and b
    // начальный интервал неопределенности найден
    if (fxa >= fxx && fxb >= fxx)
    {
        unimodal = true;
        //if (mcallback) mcallback->swannCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        return;
    }

    if (fxa < fxx && fxx > fxb)
    {
        unimodal = false;
        //if (mcallback) mcallback->swannCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        //fprintf(stderr, "%4d %8.4f %8.4f %8.4f %10.6f %10.6f %10.6f\n", -2, a, x, b, fxa, fxx, fxb);
        //fprintf(stderr, "e %4d a:%10.6f x:%10.6f b:%10.6f fxa:%10.6f fxx:%10.6f fxb:%10.6f fx_c:nt%4d\n", -2, a, x, b, fxa, fxx, fxb, fx_count);
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
            fxb = fx(b); fx_count++;

            //if (mcallback) mcallback->swannCallback(++iteration, x, a, b, fxa, fxb, fx_count);
        }
        unimodal = true;
        return;
    }

    // from right to left <<--
//    if ( fxa <= fxx && fxx <= fxb )
//    {
//        puts("<<----------------");
//        unsigned int k = 1;
//        while ( fxa < fxx )
//        {
//            b = x;
//            x = a;
//            a = x - pow(2.0, k) * fstep;
//            k++;

//            fxb = fxx;
//            fxx = fxa;
//            fxa = mfunction->fx(a); fx_count++;

//            if (mcallback) mcallback->swannCallback(++iteration, x, a, b, fxa, fxb, fx_count);
//        }
//        unimodal = true;
//        return;
//    }
}

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
void GradientMethod::goldenSectionSearch(double &x, double &a, double &b, double epsilon) const
{
    unsigned int iteration = 0;
    unsigned int fx__count = 0;

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
    //if (mcallback) mcallback->goldenSectionSearchCallback(iteration, NAN, a, b, NAN, NAN, NAN, fx__count);

    // Lazimi epsilon deqiqliyini alana qeder iterasiyalari davam edirik
    while ( fabs(b-a) > epsilon )
    {
        if (check_x1)
        {
            x1 = b - fabs(b-a)/phi;
            y1 = fx(x1); fx__count++;
        }

        if (check_x2)
        {
            x2 = a + fabs(b-a)/phi;
            y2 = fx(x2); fx__count++;
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

        x = (b+a)/2.0;
        iteration++;
        //if (mcallback) mcallback->goldenSectionSearchCallback(iteration, x, a, b, NAN, NAN, NAN, fx__count);
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
 * @param x        Точка минимума.
 * @param a        Начальная точка отрезка.
 * @param b        Конечнная точка отрезка.
 * @param epsilon  Число эпсилон для останова метода.
 */
void GradientMethod::halphIntervalMethod(double &x, double &a, double &b, double epsilon) const
{
    unsigned int iteration = 0;
    unsigned int fx__count = 0;

    if ( epsilon <= DBL_EPSILON )
    {
        std::string msg = "in function \"halphIntervalMethod\" epsilon value is equal to or less than zero.";
        throw std::invalid_argument(msg);
    }

    double l = fabs(b - a);
    x = (a+b)/2.0;
    double f_xm = fx(x); fx__count++;

    //if (mcallback) mcallback->halphIntervalMethodCallback(iteration, x, a, b, NAN, NAN, NAN, fx__count);

    while ( l > epsilon )
    {
        double xm = (a+b)/2.0;
        double x1 = a + l/4.0;
        double x2 = b - l/4.0;

        double f_x1 = fx(x1); fx__count++;
        double f_x2 = fx(x2); fx__count++;

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

        x = (b+a)/2.0;
        iteration++;
        //if (mcallback) mcallback->halphIntervalMethodCallback(iteration, x, a, b, NAN, NAN, NAN, fx__count);
    }

    //x = (b + a) /  2.0;
    //double fa = f->fx(a);
    //double fb = f->fx(b);

    //if (fa<fb) x = a;
    //if (fa>fb) x = b;
}
