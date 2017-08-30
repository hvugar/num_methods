#include "nonlinearfunctionsystem.h"
#include <math.h>

void NonLinearFunction::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    NonLinearFunction nlfs;

    //    DoubleVector x0;
    //    x0 << 3.5 << 2.2;
    //    DoubleVector xs;
    //    nlfs.calculateSimpleIdetartion(x0, xs, 0.001);
    //    IPrinter::print(xs, xs.size());

    DoubleVector x1;
    x1 << 3.4 << 2.2;
    IPrinter::print(x1, x1.size(), 10, 4);
    DoubleVector rx;
    nlfs.calculateNewtonMethodMod(x1, rx, 0.0001, 0.0001);
    IPrinter::print(rx, rx.size(), 10, 4);

}

double NonLinearFunction::fx(const DoubleVector &x, unsigned int num) const
{
    double x1 = x[0];
    double x2 = x[1];
    //    if (num == 0) return sqrt(x1 + 3.0*log10(x1));
    //    if (num == 1) return sqrt((x1*(x2+5.0)-1)/2.0);

    if (num == 0) return x1+3.0*log10(x1)-x2*x2;
    if (num == 1) return 2.0*x1*x1 - x1*x2 - 5.0*x1 + 1.0;

    return NAN;
}

void INonLinearFunction::calculateSimpleIdetartion(const DoubleVector &x0, DoubleVector &x, double epsilon)
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

void INonLinearFunction::calculateNewtonMethod(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double epsilon)
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

void INonLinearFunction::calculateNewtonMethodMod(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double epsilon)
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

double INonLinearFunction::minimize(const DoubleMatrix &W, const DoubleMatrix &V, const DoubleVector& rx, unsigned int n)
{
    class FindMinimum : public R1Function
    {
    public:
        FindMinimum(const INonLinearFunction &p, const DoubleMatrix& W, const DoubleMatrix& V, const DoubleVector& rx) : p(p), W(W), V(V), rx(rx) {}

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

        const INonLinearFunction &p;
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
