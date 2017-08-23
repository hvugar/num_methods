#ifndef INITIALBOUNDARYVALUEPROBLEM_H
#define INITIALBOUNDARYVALUEPROBLEM_H

#include "bvp.h"
#include "grid.h"
#include "uniformgrid.h"

void MINIMUMSHARED_EXPORT qovma_first_row(unsigned int K, unsigned int N,
             DoubleVector &betta, double eta,
             const DoubleMatrix &alpha,
             const DoubleMatrix &phi, const DoubleVector &psi, DoubleVector &nx);

enum SweepMethodDirection
{
    ForwardSweep = 1,
    BackwardSweep = 2,
    Centered = 3
};

class MINIMUMSHARED_EXPORT InitialValueProblem
{
public:
    enum Method
    {
        RK2,
        RK4,
        EULER,
        EULER_MOD
    };

    enum Direction
    {
        L2R, // Left to Right
        R2L  // Right to Left
    };

    const UniformODEGrid& timeGrid() const;
    void setTimeGrid(const UniformODEGrid &grid);

private:
    UniformODEGrid timegrid;
};

/**
 * @brief The InitialBoundaryValueProblemPDE class
 * @see ParabolicIBVP
 */
class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemPDE : protected BoundaryValueProblemPDE, protected InitialValueProblem
{
public:
    void setTimeDimension(const Dimension &dimension);
    const Dimension& timeDimension() const;
    void addSpaceDimension(const Dimension &dimension);
    const Dimension& spaceDimension(Dimension::SpaceDimension dim) const;
    unsigned int dimSize();

protected:



    Dimension mtimeDimension;
    std::vector<Dimension> mspaceDimension;
};

#endif // INITIALBOUNDARYVALUEPROBLEM_H
