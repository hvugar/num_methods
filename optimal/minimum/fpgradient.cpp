#include "fpgradient.h"

FastProximalGradient::FastProximalGradient() //: Gradient()
{
}

FastProximalGradient::~FastProximalGradient()
{
}

//void FastProximalGradient::setF(RnFunction *f)
//{
//    mf = f;
//}

//RnFunction* FastProximalGradient::f() const
//{
//    return mf;
//}

//void FastProximalGradient::gradient()
//{
//    double dx = 0.00001;
//    for (unsigned i=0; i<mx.size(); i++)
//    {
//        mx[i] = mx[i] - dx;
//        double f1 = mf->fx(mx);
//        mx[i] = mx[i] + 2*dx;
//        double f2 = mf->fx(mx);
//        mx[i] = mx[i] - dx;
//        mgrads[i] = (f2 - f1) / (2 * dx);
//    }
//}

//double FastProximalGradient::argmin(double alpha)
//{
//    return 0.0;
//}

//double FastProximalGradient::minimize()
//{
//    class Argmin : public R1Minimize
//    {
//    public:
//        std::vector<double>& x;
//        std::vector<double>& g;
//        Gradient* gm;
//        Argmin(std::vector<double>& x, std::vector<double>& g, Gradient* gm) : R1Minimize(), x(x), g(g), gm(gm) {}

//    protected:
//        double fx(double alpha) {
//            std::vector<double> x1;
//            for (unsigned int i=0; i < x.size(); i++) x1.push_back(x[i] - alpha * g[i]);
//            return gm->f()->fx(x);
//        }
//    };

//    Argmin r1(mx, mgrads, this);
//    r1.setX0(0.0);
//    r1.setStep(0.1);
//    r1.setEpsilon(0.000001);
//    r1.straightLineSearch();
//    double alpha = r1.goldenSectionSearch();
//    return alpha;
//}

//void FastProximalGradient::calculate()
//{
//    double grad_norm = 0.0;
//    double distance = 0.0;

//    do
//    {
//        gradient();

//        double alpha = minimize();

//        std::vector<double> x = mx;

//        grad_norm = 0.0;
//        for (unsigned int i=0; i<mx.size(); i++)
//        {
//            grad_norm = grad_norm + mx[i]*mx[i];
//        }
//        grad_norm = sqrt(grad_norm);

//        for (unsigned int i=0; i<mx.size(); i++)
//        {
//            mx[i] = mx[i] - alpha * mgrads[i];
//        }

//        distance = 0.0;
//        for (unsigned int i=0; i<mx.size(); i++)
//        {
//            distance = distance + (mx[i]-x[i])*(mx[i]-x[i]);
//        }
//        distance = sqrt(distance);


//    } while (grad_norm > epsilon() && distance > epsilon());
//}
