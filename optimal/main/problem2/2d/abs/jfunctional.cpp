#include "jfunctional.h"

JFunctional::JFunctional() : IFunctional() {}

double JFunctional::fx(const DoubleVector &prms) const
{
    JFunctional* functional = const_cast<JFunctional*>(this);

    double sum = 0.0;
    for (unsigned int i=0; i<fis.length(); i++)
    {
        functional->setIntTemperature(fis[i]);
        for (unsigned j=0; j<thetas.length(); j++)
        {
            functional->setEnvTemperature(thetas[j]);
            sum += IFunctional::fx(prms)*p_fis[i]*p_thetas[j];
        }
    }
    return sum;
}

void JFunctional::gradient(const DoubleVector &prms, DoubleVector &g)
{
    g.resize(prms.length(), 0.0);

    for (unsigned int i=0; i<fis.length(); i++)
    {
        setIntTemperature(fis[i]);
        for (unsigned j=0; j<thetas.length(); j++)
        {
            setEnvTemperature(thetas[j]);
            DoubleVector g0(prms.length(), 0.0);
            IFunctional::gradient(prms, g0);
            for (unsigned int k=0; k<prms.length(); k++)
            {
                g[k] += g0[k]*p_fis[i]*p_thetas[j];
            }
            g0.clear();
        }
    }
}

void JFunctional::setInitTemperatures(const DoubleVector &fis, const DoubleVector &p_fis)
{
    this->fis = fis;
    this->p_fis = p_fis;
}

void JFunctional::setEnvrTemperatures(const DoubleVector &thetas, const DoubleVector &p_thetas)
{
    this->thetas = thetas;
    this->p_thetas = p_thetas;
}
