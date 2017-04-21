#ifndef ART_PROBLEM1L2_H
#define ART_PROBLEM1L2_H

#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <cmethods.h>
#include <gradient/igradient.h>

#include "../iproblem1.h"

class ArtProblem1L2 : public IProjection, public IProblem1
{
public:
    ArtProblem1L2();
    virtual ~ArtProblem1L2() {}

    void initialize();
    void startOptimize();

    void optimize(DoubleVector &x0) const;

    virtual void project(DoubleVector &x, int index);

    void table1Generate();
    void table2Generate();
    void table3Generate();

    void image1Generate();
    void image2Generate();
    void image3Generate();

    void imager2L();
    void imager3L();

    static void Main(int argc, char* argv[]);
};

#endif // ART_PROBLEM1L2_H
