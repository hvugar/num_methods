#include "problem5ex1.h"

Problem5Ex1::Problem5Ex1()
{

}

double Problem5Ex1::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
#ifdef SAMPLE_1
    if (row == 1) { if (col == 1) { return +2.0; } if (col == 2) { return t; }      if (col == 3) { return -3.0; } }
    if (row == 2) { if (col == 1) { return +3.0; } if (col == 2) { return -4.0*t; } if (col == 3) { return -8.0; } }
    if (row == 3) { if (col == 1) { return +t; }   if (col == 2) { return +1.0; }   if (col == 3) { return -1.0; } }
#endif
}

double Problem5Ex1::B(double t, unsigned int k, unsigned int row) const
{
#ifdef SAMPLE_1
    if (row == 1) { return (2.0*t+1)   - (+2.0)*(t*t+t+2.0) + (t)       *(2.0*t-3.0) + (-3.0)*(t*t*t+t); }
    if (row == 2) { return (2.0)       - (+3.0)*(t*t+t+2.0) + (-4.0)*t*t*(2.0*t-3.0) + (-8.0)*(t*t*t+t); }
    if (row == 3) { return (3.0*t*t+1) - (+t)  *(t*t+t+2.0) + (+1.0)    *(2.0*t-3.0) + (-1.0)*(t*t*t+t); }
#endif
    return NAN;
}

double Problem5Ex1::B(double, unsigned int, unsigned int row, unsigned int col, unsigned int i) const
{
#ifdef SAMPLE_1
    if ( i == 0 )
    {
        if (row == 1) { if (col == 1) { return +2.0; } if (col == 2) { return +5.0; } if (col == 3) { return +3.0; } }
        if (row == 2) { if (col == 1) { return +4.0; } if (col == 2) { return +8.0; } if (col == 3) { return +1.0; } }
        if (row == 3) { if (col == 1) { return +1.0; } if (col == 2) { return +3.0; } if (col == 3) { return +4.0; } }
    }
    if ( i == 1 )
    {
        if (row == 1) { if (col == 1) { return +1.0; } if (col == 2) { return +3.0; } if (col == 3) { return +4.0; } }
        if (row == 2) { if (col == 1) { return +2.0; } if (col == 2) { return +3.0; } if (col == 3) { return +1.0; } }
        if (row == 3) { if (col == 1) { return +5.0; } if (col == 2) { return +2.0; } if (col == 3) { return +8.0; } }
    }
#endif
    return NAN;
}

double Problem5Ex1::X(double t, unsigned int k, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 1) { return t*t+t+2.0; }
    if (i == 2) { return 2.0*t-3.0; }
    if (i == 3) { return t*t*t+t; }
#endif
}

double Problem5Ex1::g(double t, unsigned int k, unsigned int row, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0)
    {
        if (row == 1) { return X(t,k,row)*X(0.3,k,row); }
        if (row == 2) { return X(t,k,row)*X(0.3,k,2)*X(0.3,k,row); }
        if (row == 3) { return X(t,k,row); }
    }
    if (i == 0)
    {
        if (row == 1) { return X(t,k,row)*X(t,k,row)*X(t,k,row); }
        if (row == 2) { return X(t,k,row); }
        if (row == 3) { return X(t,k,row)*X(t,k,row); }
    }
#endif
    return NAN;
}

double Problem5Ex1::C(double t, unsigned int k, unsigned int row) const
{
#ifdef SAMPLE_1
    return B(t,k,row,1,0)*g(0.3,k,1,0) + B(t,k,row,2,0)*g(0.3,k,2,0) + B(t,k,row,3,0)*g(0.3,k,3,0)+
           B(t,k,row,1,1)*g(0.6,k,1,1) + B(t,k,row,2,1)*g(0.6,k,2,1) + B(t,k,row,3,1)*g(0.6,k,3,1);
#endif
    return NAN;
}
