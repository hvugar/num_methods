#ifndef BOUNDARYVALUEPROBLEM_H
#define BOUNDARYVALUEPROBLEM_H

#include "../vector2d.h"
#include "../matrix2d.h"
#include "../matrix2d.h"
#include "../cmethods.h"
#include "../printer.h"
#include <math.h>

struct MINIMUMSHARED_EXPORT BoundaryValueProblem
{
    enum BoundaryCondition
    {
        Dirichlet = 1,
        Neumann = 2,
        Mixed = 3
    };

    enum BoundaryType
    {
        Left = 0,
        Right = 1
    };
};


struct MINIMUMSHARED_EXPORT BoundaryValueProblemODE : protected BoundaryValueProblem
{
   virtual double boundary(BoundaryType bound) const = 0;
};

struct MINIMUMSHARED_EXPORT BoundaryValueProblemPDE : protected BoundaryValueProblem
{
    virtual double boundary(unsigned int m, BoundaryType bound) const = 0;
};

struct MINIMUMSHARED_EXPORT Grid
{
    Grid();

    double ht;
    double hx1;
    double hx2;
    double hx3;
    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    unsigned int N3;
};

#endif // BOUNDARYVALUEPROBLEM_H
