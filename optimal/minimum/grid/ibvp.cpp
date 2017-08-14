#include "ibvp.h"

const UniformGrid& InitialValueProblem::timeGrid() const
{
    return timegrid;
}

void InitialValueProblem::setTimeGrid(const UniformGrid &grid)
{
    timegrid = grid;
}

unsigned int InitialBoundaryValueProblemPDE::dimSize()
{
    return mspaceDimension.size();
}

void InitialBoundaryValueProblemPDE::setTimeDimension(const Dimension &dimension)
{
    mtimeDimension = dimension;
}

const Dimension &InitialBoundaryValueProblemPDE::timeDimension() const
{
    return mtimeDimension;
}

void InitialBoundaryValueProblemPDE::addSpaceDimension(const Dimension &dimension)
{
    mspaceDimension.push_back(dimension);
}

const Dimension &InitialBoundaryValueProblemPDE::spaceDimension(Dimension::SpaceDimension dim) const
{
    return mspaceDimension.at(dim);
}

void qovma_first_row(unsigned int K, unsigned int N,
             DoubleVector &betta, double eta,
             const DoubleMatrix &alpha,
             const DoubleMatrix &phi,
             const DoubleVector &psi,
             DoubleVector &nx)
{
    for (unsigned int k=0; k<=N-K; k++)
    {
        for (unsigned int i=1; i<=K; i++)
        {
            betta[k+i] = betta[k+i] + alpha[k][i]*betta[k];
        }
        eta = eta - alpha[k][0]*betta[k];
    }

    DoubleMatrix M(K+1,K+1);
    DoubleVector A(K+1);
    DoubleVector x(K+1);

    M[0][0] = 0.0;
    for (unsigned int j=1; j<=K; j++)
    {
        M[0][j] = betta[N-(K-j)];
    }
    A[0] = eta;

    for (unsigned int i=1; i<=K; i++)
    {
        for (unsigned int j=0; j<=K; j++)
            M[i][j] = phi[i-1][j];
        A[i] = psi[i-1];
    }

    GaussianElimination(M, A, x);

    nx.resize(N+1);
    for (unsigned int i=0; i<=K; i++)
    {
        nx[N-i] = x[K-i];
    }

    for (unsigned int i=1; i<=K; i++)
    {
        betta[N-(K-i)] -= alpha[N-K][i]*betta[N-K];
    }
    eta += alpha[N-K][0]*betta[N-K];

    for (unsigned int k=N-K; k>0; k--)
    {
        for (unsigned int j=0; j<K; j++)
        {
            betta[k+j] -= alpha[k-1][j+1]*betta[k-1];
        }
        eta += alpha[k-1][0]*betta[k-1];

        nx[k-1] = eta;
        for (unsigned int i=k; i<=N; i++)
        {
            nx[k-1] -= betta[i]*nx[i];
        }
        nx[k-1] /= betta[k-1];
    }

    x.clear();
    A.clear();
    M.clear();
}


