#include "nonlinearequation.h"

NonLinearEquation::~NonLinearEquation()
{}

void NonLinearEquation::calculateSimpleIdetartion(const DoubleVector &x0, DoubleVector &x, double epsilon)
{
    unsigned int n = x0.length();
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
    unsigned int n = x0.length();

    rx = x0;
    DoubleMatrix W(n, n);
    DoubleMatrix WI(n,n);
    DoubleVector e(x0.length());

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
//      }   while (e.L2Norm() > epsilon);
}

void NonLinearEquation::calculateNewtonMethodMod(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double epsilon)
{
    unsigned int n = x0.length();
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
                //printf("%f %f %f\n", fx(x2, row), fx(x1, row),(fx(x2, row) - fx(x1, row))/(2.0*diffEspilon));
                W[row][col] = (fx(x2, row) - fx(x1, row))/(2.0*diffEspilon);
            }
        }
//IPrinter::printSeperatorLine();
//        IPrinter::print(W);
//        IPrinter::printSeperatorLine();

        // Regulization
        DoubleMatrix V = W + 0.00*DoubleMatrix::IdentityMatrix(n);

        V.inverse();

        double alpha = 1.0;minimize(W, V, rx, n);

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

void NonLinearEquation::calculateNewtonMethodMod2(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double epsilon)
{
    unsigned int n = x0.length();
    rx = x0;
    DoubleMatrix W(2, 2);
    DoubleVector e(2);
    unsigned int iter = 0;
    do {
        printf("%4d %20.10f %20.10f\n", iter++, rx[0], rx[1]);

        // Forming Jacobian matrix
        DoubleVector x2;
        DoubleVector x1;

        x1 = x2 = rx; x2[0] += diffEspilon; x1[0] -= diffEspilon;
        W[0][0] = (fx(x2, 0) - fx(x1, 0))/(2.0*diffEspilon);

        x1 = x2 = rx; x2[1] += diffEspilon; x1[1] -= diffEspilon;
        W[0][1] = (fx(x2, 0) - fx(x1, 0))/(2.0*diffEspilon);

        x1 = x2 = rx; x2[0] += diffEspilon; x1[0] -= diffEspilon;
        W[1][0] = (fx(x2, 1) - fx(x1, 1))/(2.0*diffEspilon);

        x1 = x2 = rx; x2[1] += diffEspilon; x1[1] -= diffEspilon;
        W[1][1] = (fx(x2, 1) - fx(x1, 1))/(2.0*diffEspilon);

        // Regulization
        DoubleMatrix V = W;// + 0.00001*DoubleMatrix::IdentityMatrix(n);

        V.inverse();

        double alpha = minimize(W, V, rx, n);// printf("alpha %.10f\n", alpha);

        //DoubleVector x = rx;

        e[0] = V[0][0]*fx(rx,0)+V[0][1]*fx(rx,1);
        e[1] = V[1][0]*fx(rx,0)+V[1][1]*fx(rx,1);

        rx[0] -= alpha*e[0];
        rx[1] -= alpha*e[1];

    } while (e.LInfNorm() > epsilon);
    printf("%4d %20.10f %20.10f\n", iter, rx[0], rx[1]);
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

//            DoubleVector x;
//            double s;

//            x = rx;
//            s = V(0,0)*p.fx(rx,0)+V(0,1)*p.fx(rx,1);
//            x[0] = rx[0] - alpha*s;
//            SUM += (p.fx(x, 0)*p.fx(x, 0)) / (W[0][0]*W[0][0] + W[0][1]*W[0][1]);

//            x = rx;
//            s = V(1,0)*p.fx(rx,0)+V(1,1)*p.fx(rx,1);
//            x[1] = rx[1] - alpha*(V[1][0]*p.fx(rx,0) + V[1][1]*p.fx(rx,1));
//            SUM += (p.fx(x, 1)*p.fx(x, 1)) / (W[1][0]*W[1][0] + W[1][1]*W[1][1]);

            for (unsigned int row=0; row<n; row++)
            {
                DoubleVector x;
                x = rx;
                double s = 0.0;
                for (unsigned int col=0; col<n; col++)
                {
                    s += alpha*V[row][col] * p.fx(rx, col);
                }
                x[row] = rx[row] - alpha*s;

                double norm = 0.0;
                for (unsigned int col=0; col<n; col++)
                {
                    norm += W[row][col] * W[row][col];
                }

                SUM += (p.fx(x, row)*p.fx(x, row)) / norm;
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

    double alpha0 = 1.0;
    double min_step = 0.1;
    double min_epsilon = 0.02;
    double a,b,alpha;
    stranghLineSearch(alpha0, min_step, a, b, &fm);
    goldenSectionSearch(a, b, alpha, &fm, min_epsilon);
    if (fm.fx(alpha) > fm.fx(alpha0)) alpha = alpha0;

    return alpha;
}
