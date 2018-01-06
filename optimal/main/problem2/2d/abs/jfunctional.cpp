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
    g.clear();
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

void JFunctional::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const
{
    IFunctional::print(i,x,g,f,result);

//    DoubleMatrix u;
//    std::vector<ExtendedSpaceNode2D> info;
//    forward->calculateMVD(u, info, true);

//    std::string f_name= "data_v" + std::to_string(i)+".txt";
//    Parameter prm(x, 2, 2);
//    FILE *file = fopen(f_name.data(), "w");
//    for (int j=0; j<mTimeDimension.sizeN(); j++)
//    {
//        double v1 = prm.k[0][0]*(info[0].value(j)-prm.z[0][0]) + prm.k[0][1]*(info[1].value(j)-prm.z[0][1]);
//        double v2 = prm.k[1][0]*(info[0].value(j)-prm.z[1][0]) + prm.k[1][1]*(info[1].value(j)-prm.z[1][1]);
//        fprintf(file, "%f %f\n", v1, v2);
//    }
//    fclose(file);
//    puts("---");
}
