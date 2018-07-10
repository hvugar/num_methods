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
class InitialBoundaryValueProblemPDE : public BoundaryValueProblemPDE, public InitialValueProblem
{
public:
    MINIMUMSHARED_EXPORT void setTimeDimension(const Dimension &dimension);
    MINIMUMSHARED_EXPORT const Dimension& timeDimension() const;
    MINIMUMSHARED_EXPORT void addSpaceDimension(const Dimension &dimension);
    MINIMUMSHARED_EXPORT const Dimension& spaceDimension(Dimension::SpaceDimension dim) const;
    MINIMUMSHARED_EXPORT unsigned int dimSize();

protected:
    Dimension mtimeDimension;
    std::vector<Dimension> mspaceDimension;
};

#endif // INITIALBOUNDARYVALUEPROBLEM_H
