#include "islodenlcsv2.h"

void ISystemLinearODENonLocalContionsV2::calculateForward()
{
    unsigned int L = atimes.size();
    unsigned int N = mgrid.dimension().sizeN();

    unsigned int n0 = alphas.at(0).rows();

    for (unsigned int row = 0; row<n0; row++)
    {
        for (unsigned int s=0; s<L-1; s++)
        {

        }
    }
}

void ISystemLinearODENonLocalContionsV2::calculateBackward()
{}

void ISystemLinearODENonLocalContionsV2::addCondition(const DoubleMatrix &alpha, double time, unsigned int nmbr)
{
    alphas.push_back(alpha);
    atimes.push_back(time);
    anmbrs.push_back(nmbr);
}

void ISystemLinearODENonLocalContionsV2::addLoadPoint(const DoubleMatrix &betta, double time, unsigned int nmbr)
{
    bettas.push_back(betta);
    btimes.push_back(time);
    bnmbrs.push_back(nmbr);
}

void ISystemLinearODENonLocalContionsV2::setRightSize(const DoubleVector &gamma)
{
    this->gamma = gamma;
}

void ISystemLinearODENonLocalContionsV2::setGrid(const ODEGrid &grid)
{
    mgrid = grid;
}


