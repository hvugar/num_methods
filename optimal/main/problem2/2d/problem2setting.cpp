#include "problem2setting.h"

void P2Setting::toVector(DoubleVector &prms) const
{
    prms.clear();
    prms.resize(2*Lc*Lo + 2*Lo + 2*Lc);

    // k
    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            prms[i*Lo + j] = k[i][j];
        }
    }

    // z
    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            prms[i*Lo + j + Lc*Lo] = z[i][j];
        }
    }

    // xi
    for (unsigned int j=0; j<Lo; j++)
    {
        prms[2*j + 0 + 2*Lc*Lo] = xi[j].x;
        prms[2*j + 1 + 2*Lc*Lo] = xi[j].y;
    }

    // eta
    for (unsigned int i=0; i<Lc; i++)
    {
        prms[2*i + 0 + 2*Lo + 2*Lc*Lo] = eta[i].x;
        prms[2*i + 1 + 2*Lo + 2*Lc*Lo] = eta[i].y;
    }
}

void P2Setting::fromVector(const DoubleVector &prms)
{
    k.clear();
    k.resize(Lc, Lo);

    z.clear();
    z.resize(Lc, Lo);

    xi.clear();
    xi.resize(Lo);

    eta.clear();
    eta.resize(Lc);

    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            k[i][j] = prms[i*Lo + j];
        }
    }

    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            z[i][j] = prms[Lc*Lo + i*Lo + j];
        }
    }

    for (unsigned int j=0; j<Lo; j++)
    {
        xi[j].x = prms[2*Lc*Lo + 2*j];
        xi[j].y = prms[2*Lc*Lo + 2*j+1];
    }

    for (unsigned int i=0; i<Lc; i++)
    {
        eta[i].x = prms[2*Lc*Lo + 2*Lo + 2*i];
        eta[i].y = prms[2*Lc*Lo + 2*Lo + 2*i+1];
    }
}
