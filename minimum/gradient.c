#include "gradient.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double L2Norm(double *vctr, int n)
{
    double sum = 0.0;
    int i;
    for (i=0; i<n; i++)
    {
        sum += vctr[i] * vctr[i];
    }
    return sqrt(sum);
}

void Gradient(RnFunction f, double dx, double *x, int n, double *gr)
{
    int i = 0;
    for (i=0; i<n; i++)
    {
        double c = x[i];
        x[i] = x[i] - dx;
        double f1 = f(x, n);
        x[i] = c;
        x[i] = x[i] + dx;
        double f2 = f(x, n);
        x[i] = x[i] - dx;
        x[i] = c;
        gr[i] = (f2 - f1) / (2 * dx);
    }
}

void SteepestDescentMethod(RnFunction f, GFunction gradient, double *x, int n, double line_step, double gold_eps, double grad_eps, double epsilon1, double epsilon2)
{
    int i = 0;

    int iterationCount = 0;
    double distance = 0.0;
    double alpha = 0.0;

    double *gr = (double*) malloc(sizeof(double) * n);
    double *xc = (double*) malloc(sizeof(double) * n);

    do
    {
        /* calculating function gradient at current point */
        gradient(f, grad_eps, x, n, gr);

        /* if gradinet norm at current point is less than epsilon then break. no minimize */
        double gradient_norm = L2Norm(gr, n);
        if (gradient_norm < epsilon1)
        {
            puts("Optimisation ends, because L2 norm of gradient is less than epsilon...");
            break;
        }

        iterationCount++;

        /* Normalize vector */
        for (i=0; i<n; i++) gr[i] = gr[i] / gradient_norm;

        /* R1 minimization in direct of antigradient */
        double argmin(double alpha)
        {
            int i;
            for (i=0; i<n; i++)
                xc[i] = x[i] - alpha * gr[i];
            return f(xc, n);
        }

        double alpha0 = 0.0;
        double a = 0.0;
        double b = 0.0;
        straight_line_search_metod(argmin, alpha0, line_step, &a, &b);
        alpha = golden_section_search_min(argmin, a, b, gold_eps);
        if ( argmin(alpha) > argmin(alpha0) ) alpha = alpha0;
        /* R1 minimization in direct of antigradient */

        /* printing... */

        distance = 0.0;
        double f1 = f(x, n);
        for (i=0; i<n; i++)
        {
            double cx = x[i];
            x[i] = x[i] - alpha * gr[i];

            // calculating distance
            distance += (x[i]-cx)*(x[i]-cx);
        }
        distance = sqrt(distance);
        double f2 = f(x, n);

        if (f2 > f1)
        {
            puts("VValue of the function in the previous point less than in the current point...");
            break;
        }

        /* calculating distance previous and new point */
        if (distance < epsilon2 && fabs(f2 - f1) < epsilon2)
        {
            puts("Optimisation ends, because distance beetween last and current point less than epsilon...");
            break;
        }
    }
    while ( 1 );

    printf("%.6f %.6f %.6f\n", x[0], x[1], f(x, n));

    free(xc);
    free(gr);
}

void ConjugateGradientMethod(RnFunction f, GFunction gradient, double *x, int n, double line_step, double gold_eps, double grad_eps, double epsilon1, double epsilon2)
{
    int i = 0;
    int k = 0;
    int iterationCount = 0;

    // Gradient of x
    double* gr = (double*) malloc(sizeof(double) * n);
    // Direction
    double *s  = (double*) malloc(sizeof(double) * n);
    // copy of x
    double *xc = (double*) malloc(sizeof(double) * n);

    double gr0_mod = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;

    do
    {
        printf("%.6f %.6f %.6f\n", x[0], x[1], f(x, n));
        // Gradient of objectiv function in current point
        gradient(f, grad_eps, x, n, gr);

        double gradient_norm = L2Norm(gr, n);
        if (gradient_norm < epsilon1)
        {
            puts("Optimisation ends, because L2 norm of gradient is less than epsilon...");
            break;
        }

        iterationCount++;

        // Module of gradient
        gr0_mod = 0.0;
        for (i=0; i<n; i++) gr0_mod = gr0_mod + gr[i]*gr[i];
        //gr0_mod = sqrt(gr0_mod);

        // First iteration
        if (k == 0)
        {
            gr1_mod = gr0_mod;
            // First direction is antigradient
            for (i=0; i<n; i++) s[i] = -gr[i];
        }
        else
        {
            gr2_mod = gr0_mod;
            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;
            // Direction in next iteration (k != 0)
            for (i=0; i<n; i++) s[i] = -gr[i] + s[i] * w;
        }

        /* Normalize vector */
        //        double sn = 0.0;
        //        for (i=0; i<n; i++) sn = sn + s[i]*s[i];
        //        for (i=0; i<n; i++) s[i] = s[i] / sn;

        /* R1 minimization in direct of antigradient */
        double argmin(double alpha)
        {
            int i;
            for (i=0; i<n; i++)
                xc[i] = x[i] + alpha * s[i];
            double result = f(xc, n);
            return result;
        }

        double alpha0 = 0.0;
        double a = 0.0;
        double b = 0.0;
        straight_line_search_metod(argmin, alpha0, line_step, &a, &b);
        double alpha = golden_section_search_min(argmin, a, b, gold_eps);
        if ( argmin(alpha) > argmin(alpha0) ) alpha = alpha0;
        /* R1 minimization in direct of antigradient */

        /* printing... */

        double distance = 0.0;
        double f1 = f(x, n);
        for (i=0; i<n; i++)
        {
            double cx = x[i];
            x[i] = x[i] + alpha * s[i];

            // calculating distance
            distance += (x[i]-cx)*(x[i]-cx);
        }
        distance = sqrt(distance);
        double f2 = f(x, n);

        if ( k == n ) { k = 0; } else { k++; }

        if (f2 > f1)
        {
            puts("VValue of the function in the previous point less than in the current point...");
            break;
        }

        /* calculating distance previous and new point */
        if (distance < epsilon2 && fabs(f2 - f1) < epsilon2)
        {
            puts("Optimisation ends, because distance beetween last and current point less than epsilon...");
            break;
        }

    } while ( 1 );

    free(xc);
    free(s);
    free(gr);
}
