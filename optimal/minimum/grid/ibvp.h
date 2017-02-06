#ifndef INITIALBOUNDARYVALUEPROBLEM_H
#define INITIALBOUNDARYVALUEPROBLEM_H

#include "bvp.h"
#include "grid.h"

enum SweepMethodDirection
{
    ForwardSweep = 1,
    BackwardSweep = 2,
    Centered = 3
};

class MINIMUMSHARED_EXPORT InitialValueProblem {};

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
