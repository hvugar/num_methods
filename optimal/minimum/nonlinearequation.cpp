#include "nonlinearequation.h"

void NonLinearEquation::calculateSimpleIdetartion(const DoubleVector &x0, DoubleVector &x, double epsilon)
{
    unsigned int n = x0.size();
    x = x0;
    DoubleVector x1 = x0;
    DoubleVector dx(n);
    while (true)
    {
        for (unsigned int i=0; i<n; i++)
        {
            x1[i] = fx(x, i);
            dx[i] = x1[i] - x[i];
        }
        x = x1;
        if (dx.L2Norm() <= epsilon) break;
    }
}

void NonLinearEquation::calculateNewtonMethod(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double epsilon)
{
    unsigned int n = x0.size();

    rx = x0;
    DoubleMatrix W(n, n);
    DoubleMatrix WI(n,n);
    DoubleVector e(x0.size());

    do {

        // Forming Jacobian matrix
        DoubleVector x2;
        DoubleVector x1;
        for (unsigned int row=0; row<n; row++)
        {
            for (unsigned int col=0; col<n; col++)
            {
                x2 = x1 = rx;
                x2[col] += diffEspilon;
                x1[col] -= diffEspilon;
                W.at(row, col) = (fx(x2, row) - fx(x1, row))/(2.0*diffEspilon);
            }
        }

        WI = W;
        WI.inverse();

        DoubleVector x3 = rx;

         for (unsigned int row=0; row<n; row++)
        {
            e[row] = 0.0;
            for (unsigned int col=0; col<n; col++)
            {
                e[row] += WI[row][col]*fx(x3,col);
            }
            rx[row] -= e[row];
        }

    } while (e.LInfNorm() > epsilon);
}

void NonLinearEquation::calculateNewtonMethodMod(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double epsilon)
{
    unsigned int n = x0.size();
    rx = x0;

    DoubleMatrix W(n, n);
    DoubleVector e(n);
    do {

        // Forming Jacobian matrix
        for (unsigned int row=0; row<n; row++)
        {
            for (unsigned int col=0; col<n; col++)
            {
                DoubleVector x2 = rx;
                DoubleVector x1 = rx;
                x2[col] += diffEspilon;
                x1[col] -= diffEspilon;

                W.at(row, col) = (fx(x2, row) - fx(x1, row))/(2.0*diffEspilon);
            }
        }
        // Regulization
        DoubleMatrix V = W + 0.0*DoubleMatrix::IdentityMatrix(n);

        V.inverse();

        double alpha = minimize(W, V, rx, n);

        DoubleVector x = rx;

        for (unsigned int row=0; row<n; row++)
        {
            e[row] = 0.0;
            for (unsigned int col=0; col<n; col++)
            {
                e[row] += V[row][col]*fx(rx,col);
            }
            x[row] -= alpha*e[row];
        }

        rx = x;

    } while (e.LInfNorm() > epsilon);
}

double NonLinearEquation::minimize(const DoubleMatrix &W, const DoubleMatrix &V, const DoubleVector& rx, unsigned int n)
{
    class FindMinimum : public R1Function
    {
    public:
        FindMinimum(const NonLinearEquation &p, const DoubleMatrix& W, const DoubleMatrix& V, const DoubleVector& rx) : p(p), W(W), V(V), rx(rx) {}

        virtual double fx(double alpha) const
        {
            double SUM = 0.0;
            for (unsigned int row=0; row<n; row++)
            {
                double norm = 0.0;
                for (unsigned int col=0; col<n; col++)
                {
                    norm += W[row][col] * W[row][col];
                }

                DoubleVector x = rx;
                for (unsigned int col=0; col<n; col++)
                {
                    x[row] -= alpha*V[row][col] * p.fx(rx, col);
                }

                SUM += p.fx(x, row) / norm;
            }
            return SUM;
        }

        const NonLinearEquation &p;
        const DoubleMatrix &W;
        const DoubleMatrix &V;
        const DoubleVector &rx;
        unsigned int n;
    };

    FindMinimum fm(*this, W,V,rx);
    fm.n = n;

    double alpha0 = 0.0;
    double min_step = 0.01;
    double min_epsilon = 0.0001;
    double a,b,alpha;
    stranghLineSearch(alpha0, min_step, a, b, &fm);
    goldenSectionSearch(a, b, alpha, &fm, min_epsilon);
    if (fm.fx(alpha) > fm.fx(alpha0)) alpha = alpha0;

    return alpha;
}
