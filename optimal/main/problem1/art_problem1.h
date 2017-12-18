#ifndef ART_PROBLEM1_H
#define ART_PROBLEM1_H

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

#include "iproblem1.h"

class ArtProblem1 : public IProjection, public IProblem1
{
public:
    ArtProblem1();
    virtual ~ArtProblem1() {}

    void initialize();
    void startOptimize();

    void optimize(DoubleVector &x0) const;

    virtual void project(DoubleVector &x, unsigned int index);

    void table1Generate();
    void table2Generate();
    void table3Generate();

    void image1Generate();
    void image2Generate();
    void image3Generate();

    void imager2L();
    void imager3L();

    void image1L();

    static void Main(int argc, char* argv[]);
};

#endif // ART_PROBLEM1_H
