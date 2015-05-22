#include "r1minimize.h"
#include <stdio.h>
#include <math.h>

R1Minimize::R1Minimize()
{
    m_a = m_b = m_x0 = 0.0;
    m_step = 0.1;
    m_epsilon = 0.01;
}

R1Minimize::~R1Minimize()
{
}

double R1Minimize::x0() const
{
    return m_x0;
}

void R1Minimize::setX0(double x0)
{
    m_x0 = x0;
}

double R1Minimize::step() const
{
    return m_step;
}
void R1Minimize::setStep(double step)
{
    m_step = step;
}

double R1Minimize::epsilon() const
{
    return m_epsilon;
}

void R1Minimize::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

double R1Minimize::a() const
{
    return m_a;
}

double R1Minimize::b() const
{
    return m_b;
}

double R1Minimize::straightLineSearch()
{
    if ( m_step == 0.0 )
    {
        return NAN;
    }

    double y0 = fx(m_x0);
    double y1 = fx(m_x0 - m_step);
    double y2 = fx(m_x0 + m_step);

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
                y2 = fx(m_x0 + fabs(m_step));
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
                y1 = fx(m_x0 - fabs(m_step));
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

double R1Minimize::halphIntervalMethod()
{
    double L = m_b - m_a;

    while ( L > m_epsilon )
    {
        double xm = (m_a+m_b)/2.0;
        double x1 = m_a + L/4.0;
        double x2 = m_b - L/4.0;

        double f_xm = fx(xm);
        double f_x1 = fx(x1);
        double f_x2 = fx(x2);

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

double R1Minimize::goldenSectionSearch()
{
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
    while ( fabs(m_b-m_a) > m_epsilon )
    {
        if (isnan(x1))
        {
            x1 = m_b - fabs(m_b-m_a)/phi;
            y1 = fx(x1);
        }

        if (isnan(x2))
        {
            x2 = m_a + fabs(m_b-m_a)/phi;
            y2 = fx(x2);
        }

        if (y1 >= y2)
        {
            m_a = x1;
            x1 = x2;    // Tapilmish x2 noqtesi ve bu noqtede funkisiyanin qiymeti
            y1 = y2;    // sonraki iterasiyada x1 qiymeti kimi istifade olunacaq.
            x2 = NAN;   // x2 novbeti iterasiyada axtarilacaq
        }
        else
        {
            m_b = x2;
            x2 = x1;    // Tapilmish x1 noqtesi ve bu noqtede funkisiyanin qiymeti
            y2 = y1;    // sonraki iterasiyada x2 qiymeti kimi istifade olunacaq.
            x1 = NAN;   // x1 novbeti iterasiyada axtarilacaq
        }
    }

    double c = (m_a+m_b)/2.0;

    if (fx(m_a)<fx(m_b)) c = m_a;
    if (fx(m_a)>fx(m_b)) c = m_b;

    return c;
}
