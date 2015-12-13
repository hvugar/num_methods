#include "gradient_sd.h"

SteepestDescentGradient::SteepestDescentGradient() : GradientMethod()
{
    setNormalize(true);
}

SteepestDescentGradient::~SteepestDescentGradient()
{}

void SteepestDescentGradient::calculate(DoubleVector &x)
{
    unsigned int n = x.size();

    iterationCount = 0;
    double distance = 0.0;
    double alpha = 0.0;

    DoubleVector g(n);

    do
    {
        /* calculating function gradient at current point */
        m_fn->gradient(x, g, grad_step);

        /* if gradinet norm at current point is less than epsilon then break. no minimize */
        double gradient_norm = g.L2Norm();
        if (gradient_norm < epsilon1())
        {
            puts("Optimisation ends, because L2 norm of gradient is less than epsilon...");
            break;
        }


        iterationCount++;

        /* Normalize vector */
        if (m_normalize) g.L2Normalize();

        /* R1 minimization in direct of antigradient */
        alpha = minimize(x, g);

        if (m_printer != NULL)
        {
            m_printer->print(iterationCount, x, g, alpha, function());
        }

        distance = 0.0;
        double f1 = m_fn->fx(x);
        for (unsigned int i=0; i < n; i++)
        {
            double cx = x[i];
            x[i] = x[i] - alpha * g[i];

            // calculating distance
            distance += (x[i]-cx)*(x[i]-cx);
        }
        distance = sqrt(distance);
        double f2 = m_fn->fx(x);

        /* calculating distance previous and new point */
        if (distance < epsilon2() && fabs(f2 - f1) < epsilon2())
        {
            puts("Optimisation ends, because distance beetween last and current point less than epsilon...");
            break;
        }

    } while (true);

    g.clear();
}

double SteepestDescentGradient::minimize(const DoubleVector &x, const DoubleVector &g)
{
    struct SteepestDescentR1Function : public R1Function
    {
        virtual double fx(double alpha)
        {
            DoubleVector cx(n);
            for (unsigned int i=0; i < n; i++)
            {
                cx[i] = x[i] - alpha * g[i];
            }
            return f->fx(cx);
        }

        SteepestDescentR1Function(const DoubleVector &x, const DoubleVector &g, RnFunction *f, unsigned int n) : x(x), g(g), f(f), n(n) {}
        const DoubleVector &x;
        const DoubleVector &g;
        RnFunction *f;
        unsigned int n;
    };

    SteepestDescentR1Function r1X(x, g, m_fn, x.size());

    double alpha0 = 0.0;
    R1Minimize r1;
    r1.setFunction(&r1X);
    r1.setX0(alpha0);
    r1.setStep(min_step);
    r1.setEpsilon(min_epsilon);
    r1.straightLineSearch();
    double alpha = r1.goldenSectionSearch();
    if (r1X.fx(alpha) > r1X.fx(alpha0)) alpha = alpha0;
    return alpha;

//    double alpha0 = 0.0;
//    double a,b,alpha;
//    R1Minimize::StranghLineSearch(alpha0, min_step, a, b, &r1X);
//    R1Minimize::GoldenSectionSearch(a, b, alpha, &r1X, min_epsilon);
//    if (r1X.fx(alpha) > r1X.fx(alpha0)) alpha = alpha0;
//    return alpha;
}
