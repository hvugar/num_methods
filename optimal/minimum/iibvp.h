#ifndef IPDEINITIALBOUNDARYVALUEPROBLEM_H
#define IPDEINITIALBOUNDARYVALUEPROBLEM_H

#include "global.h"
#include "vector2d.h"
#include "matrix2d.h"
#include "matrix3d.h"
#include "cmethods.h"
#include "printer.h"
#include <stdio.h>
#include <math.h>

class IPDEInitialBoundaryValueProblem
{
public:
    enum Boundary
    {
        Left = 0,
        Right = 1
    };

    enum BoundaryCondition
    {
        Dirichlet = 1,
        Neumann = 2,
        Mixed = 3
    };
};

#endif // IPDEINITIALBOUNDARYVALUEPROBLEM_H
