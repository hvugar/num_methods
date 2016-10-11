#include "cgradient.h"

void conjucate_gradient(struct c_gradient_t *gr)
{
    unsigned int n = gr->n;

    unsigned int k = 0;
    double gr0_mod = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    double distance = 0.0;
    double alpha = 0.0;
    double f1 = 0.0;
    double f2 = 0.0;
    int firstIteration = 1;


//    DoubleVector g(n);
//    DoubleVector s(n);

//    mx = &x;
//    ms = &s;

    double minimize(double *x UNUSED_PARAM, unsigned int n UNUSED_PARAM)
    {
//        double alpha0 = 0.0;
//        double a,b,alpha;

//        stranghLineSearch(alpha0, min_step, a, b, this);
//        goldenSectionSearch(a, b, alpha, this, min_epsilon);

//        if (this->fx(alpha) > this->fx(alpha0)) alpha = alpha0;

//        return alpha;
        return 1.0;
    }

    gr->count = 0;
    double *g = (double*)malloc(sizeof(double)*n);
    double *s = (double*)malloc(sizeof(double)*n);
    do
    {
        // Gradient of objectiv function in current point
        gr->grad(gr->x, g, n);

        double gradient_norm = c_vector_L2Norm(g, n );
        if (gradient_norm < gr->grad_norm_esp)
        {
            if (gr->showMsg) puts("Optimisation ends, because L2 norm of gradient is less than epsilon...");
            break;
        }

        gr->count++;

        // Module of gradient
        gr0_mod = 0.0;
        for (unsigned int i=0; i<n; i++) gr0_mod = gr0_mod + (g[i]*g[i]);
        //gr0_mod = sqrt(gr0_mod);

        if (k == 0)
        {
            gr1_mod = gr0_mod;
            // First direction is antigradient
            for (unsigned int i=0; i<n; i++) s[i] = -g[i];
        }
        else
        {
            gr2_mod = gr0_mod;
            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;
            // Direction in next iteration (k != 0)
            for (unsigned int i=0; i<n; i++) s[i] = -g[i] + s[i] * w;
        }

        /* Normalize vector */
        if (gr->normalize) c_vector_L2Normalize(s, n);

        /* R1 minimization in direct of antigradient */
//        alpha = minimize(x, s);

        f1 = f2;
        if (firstIteration)
        {
            f1 = gr->func(gr->x, n);
            firstIteration = 0;
        }

        distance = 0.0;
        for (unsigned int i=0; i<n; i++)
        {
            double cx = gr->x[i];
            gr->x[i] = gr->x[i] + alpha * s[i];

//            if (m_projection != NULL) m_projection->project(x, i);

            distance += (gr->x[i]-cx)*(gr->x[i]-cx);
        }
        distance = sqrt(distance);
        f2 = gr->func(gr->x, n);

        if ( k == n ) { k = 0; } else { k++; }

//        if (m_printer != NULL) m_printer->print(iterationCount, x, g, alpha, m_fn);

        /* calculating distance previous and new point */
        if (distance < gr->distance_eps && fabs(f2 - f1) < gr->func_diff_eps)
        {
            if (gr->showMsg) puts("Optimisation ends, because distance between last and current point less than epsilon...");
            break;
        }
    } while (1);
}

void steepest_descent_gradient(struct c_gradient_t *gr)
{
    C_UNUSED(gr);
}

double c_vector_L2Norm(double *x, unsigned int n)
{
    C_UNUSED(x);
    C_UNUSED(n);
    return 0.0;
}

void c_vector_L2Normalize(double *x, unsigned int n)
{
    C_UNUSED(x);
    C_UNUSED(n);

}
