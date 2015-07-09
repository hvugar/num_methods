#include "cjtgradient.h"

ConjugateGradient::ConjugateGradient() : Gradient()
{
}

ConjugateGradient::~ConjugateGradient()
{}

double ConjugateGradient::minimize()
{
    double alpha0 = 0.0;
    R1Minimize r1;
    r1.setF(this);
    r1.setX0(alpha0);
    r1.setStep(min_step);
    r1.setEpsilon(min_epsilon);
    r1.straightLineSearch();
    double alpha = r1.goldenSectionSearch();
    if ( fx(alpha) > fx(alpha0) ) alpha = alpha0;
    return alpha;
}

void ConjugateGradient::calculate()
{
    unsigned int n = 0;
    double gr0_mod = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    s = mx;

    do
    {
        // Gradient of objectiv function in current point
        calculateGradient();

        if (gradientNorm() < epsilon())
            break;

        printf("gradient: %.10f epsilon: %.10f\n", gradientNorm(), epsilon());

        k++;

        // Module of gradient
        gr0_mod = 0.0;
        for (unsigned int i=0; i<mg.size(); i++) gr0_mod = gr0_mod + mg[i]*mg[i];
        //gr0_mod = sqrt(gr0_mod);

        // First iteration
        if (n == 0)
        {
            gr1_mod = gr0_mod;
            // First direction is antigradient
            for (unsigned int i=0; i<mg.size(); i++) s[i] = -mg[i];
        }
        else
        {
            gr2_mod = gr0_mod;
            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;
            // Direction in next (k+1) iteration
            for (unsigned int i=0; i<mg.size(); i++) l  = -mg[i] + s[i] * w;
        }

        double sn = 0.0;
        for (unsigned int i=0; i<s.size(); i++) sn = sn + s[i]*s[i];
        sn = sqrt(sn);
        for (unsigned int i=0; i<s.size(); i++) s[i] = s[i] / sn;

        malpha = minimize();

        print();

        for (unsigned int i=0; i<mx.size(); i++)
        {
            mx[i] = mx[i] + malpha * s[i];
        }

        if ( n == mx.size() ) { n = 0; } else { n++; }

        printf("distance: %.10f epsilon: %.10f alpha: %.10f\n", distance(), epsilon(), malpha);

        /* calculating distance previous and new point */
    } while (distance() > epsilon());
}

double ConjugateGradient::distance() const
{
    double dist = 0.0;
    for (unsigned int i=0; i<mx.size(); i++)
    {
        dist = dist + (malpha * s[i]) * (malpha * s[i]);
    }
    dist = sqrt(dist);
    return dist;
}

void ConjugateGradient::print()
{
    if (k == 1)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|s1      \t|s2      \t|grad_norm  \t|alpha  \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+-------------\n");
    }

    double y = function()->fx(mx);
    double nr = gradientNorm();

    printf("%d\t", k);
    mx[0]>=0.0 ? printf("|+%.10f\t", fabs(mx[0])) : printf("|%.10f\t", mx[0]);
    mx[1]>=0.0 ? printf("|+%.10f\t", fabs(mx[1])) : printf("|%.10f\t", mx[1]);
    y>=0.0 ? printf("|%+10.6f\t", y) : printf("|%10.6f\t", y);
    s[0]>=0.0 ? printf("|%+10.6f\t", s[0]) : printf("|%10.6f\t", s[0]);
    s[1]>=0.0 ? printf("|%+10.6f\t", s[1]) : printf("|%10.6f\t", s[1]);
    nr>=0.0 ? printf("|%+10.6f\t", nr) : printf("|%10.6f\t", nr);
    malpha>=0.0 ? printf("|%+10.6f\t", malpha) : printf("|%10.6f\t", malpha);
    printf("\n");
}

double ConjugateGradient::fx(double alpha)
{
    std::vector<double> x1 = mx;
    for (unsigned int i=0; i < mx.size(); i++)
    {
        x1[i] = mx[i] + alpha * s[i];
    }
    return mfn->fx(x1);
}
