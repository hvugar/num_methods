#include "gradient_cjt.h"

ConjugateGradient::ConjugateGradient() : GradientMethod()
{
    setNormalize(true);
}

ConjugateGradient::~ConjugateGradient()
{}

void ConjugateGradient::calculate(DoubleVector& x)
{
    unsigned int n1 = x.size();

    unsigned int n = 0;
    double gr0_mod = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    double distance = 0.0;
    double alpha = 0.0;

    DoubleVector g(n1);
    DoubleVector s(n1);
    do
    {
        // Gradient of objectiv function in current point
        //calculateGradient();
        m_fn->gradient(grad_step, x, g);

        double gradNorm = g.L2Norm();
        if (gradNorm < epsilon())
        {
            puts("Optimisation ends, because L2 norm of gradient is less than epsilon...");
            break;
        }

        iterationCount++;

        // Module of gradient
        gr0_mod = 0.0;
        for (unsigned int i=0; i<n1; i++) gr0_mod = gr0_mod + (g[i]*g[i]);
        //gr0_mod = sqrt(gr0_mod);

        if (n == 0)
        {
            gr1_mod = gr0_mod;
            // First direction is antigradient
            for (unsigned int i=0; i<n1; i++) s[i] = -g[i];
        }
        else
        {
            gr2_mod = gr0_mod;
            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;
            // Direction in next iteration (n != 0)
            for (unsigned int i=0; i<n1; i++) s[i] = -g[i] + s[i] * w;
        }

        /* Normalize vector */
        if (normalize) s.L2Normalize();

        alpha = minimize(x, s);

        if (printer != NULL) printer->print(iterationCount, x, s, alpha, function());

        distance = 0.0;
        for (unsigned int i=0; i<n1; i++)
        {
            double cx = x[i];
            x[i] = x[i] + alpha * s[i];

            distance += (x[i]-cx)*(x[i]-cx);
        }
        distance = sqrt(distance);

        if ( n == (x.size()) ) { n = 0; } else { n++; }

        /* calculating distance previous and new point */
    } while (distance > epsilon());
}

double ConjugateGradient::minimize(const DoubleVector &x, const DoubleVector &s)
{
    struct ConjugateR1Function : public R1Function
    {
        virtual double fx(double alpha)
        {
            DoubleVector cx(n);
            for (unsigned int i=0; i < n; i++)
            {
                cx[i] = x[i] + alpha * s[i];
            }
            return f->fx(cx);
        }

        ConjugateR1Function(const DoubleVector &x, const DoubleVector &s, RnFunction *f, unsigned int n) : x(x), s(s), f(f), n(n) {}
        const DoubleVector &x;
        const DoubleVector &s;
        RnFunction *f;
        unsigned int n;
    };

    ConjugateR1Function r1X(x, s, m_fn, x.size());

//    double alpha0 = 0.0;
//    R1Minimize r1;
//    r1.setFunction(&r1X);
//    r1.setX0(alpha0);
//    r1.setStep(min_step);
//    r1.setEpsilon(min_epsilon);
//    r1.straightLineSearch();
//    double alpha = r1.goldenSectionSearch();
//    if (r1X.fx(alpha) > r1X.fx(alpha0)) alpha = alpha0;
//    return alpha;

    double alpha0 = 0.0;
    double a,b,alpha;
    R1Minimize::StranghLineSearch(alpha0, min_step, a, b, &r1X);
    R1Minimize::GoldenSectionSearch(a, b, alpha, &r1X, min_epsilon);
    if (r1X.fx(alpha) > r1X.fx(alpha0)) alpha = alpha0;
    return alpha;
}
