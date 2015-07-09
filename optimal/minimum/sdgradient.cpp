#include "sdgradient.h"

SteepestDescentGradient::SteepestDescentGradient() : Gradient()
{
}

SteepestDescentGradient::~SteepestDescentGradient()
{}

double SteepestDescentGradient::minimize()
{
    double alpha0 = 0.0;
    R1Minimize r1;
    r1.setF(this);
    r1.setX0(alpha0);
    r1.setStep(min_step);
    r1.setEpsilon(min_epsilon);
    r1.straightLineSearch();
    double alpha = r1.goldenSectionSearch();
    if ( this->fx(alpha) > this->fx(alpha0) ) alpha = alpha0;
    return alpha;
}

void SteepestDescentGradient::calculate()
{
    k = 0;
    do
    {
        /* calculating function gradient at current point */
        calculateGradient();

        /* if gradinet norm at current point is less than epsilon then break. no minimize */
        double gradient_norm = gradientNorm();
        if (gradient_norm < epsilon())
            break;

        k++;

        /* calculating unit vectors */
        for (unsigned int i=0; i<mg.size(); i++) mg[i] = mg[i] / gradient_norm;

        /* R1 minimization in direct of antigradient */
        malpha = minimize();

        print();

        for (unsigned int i=0; i<mx.size(); i++)
        {
            mx[i] = mx[i] - malpha * mg[i];
        }

        /* calculating distance previous and new point */
    } while (distance() > epsilon());
}

double SteepestDescentGradient::fx(double alpha)
{
    std::vector<double> x1 = mx;
    for (unsigned int i=0; i < mx.size(); i++)
    {
        x1[i] = mx[i] - alpha * mg[i];
    }
    return mfn->fx(x1);
}

void SteepestDescentGradient::print()
{
    if (k == 1)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|grad1      \t|grad2      \t|grad_norm  \t|alpha  \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+-------------\n");
    }

    double y = function()->fx(mx);
    double nr = gradientNorm();

    printf("%d\t", k);
    mx[0]>=0.0 ? printf("|+%.10f\t", fabs(mx[0])) : printf("|%.10f\t", mx[0]);
    mx[1]>=0.0 ? printf("|+%.10f\t", fabs(mx[1])) : printf("|%.10f\t", mx[1]);
    y>=0.0 ? printf("|%+10.6f\t", y) : printf("|%10.6f\t", y);
    mg[0]>=0.0 ? printf("|%+10.6f\t", mg[0]) : printf("|%10.6f\t", mg[0]);
    mg[1]>=0.0 ? printf("|%+10.6f\t", mg[1]) : printf("|%10.6f\t", mg[1]);
    nr>=0.0 ? printf("|%+10.6f\t", nr) : printf("|%10.6f\t", nr);
    malpha>=0.0 ? printf("|%+10.6f\t", malpha) : printf("|%10.6f\t", malpha);
    printf("\n");
}
