#include "problem0h_exporter.h"
#include "problem0h_solver.h"

h0p::ProblemSolver fw1;

void init_problem()
{
//    fw1.external_source = {SpacePoint(0.25, 0.36), 0.10, 0.05, 0.05, 0.01};

//    fw1.source_number = 2;

//    fw1.waveSpeed = 1.0;
//    fw1.waveDissapation = 0.0;

//    fw1.eps1 = 1.0;
//    fw1.eps2 = 1.0;

//    unsigned int Nt = 4000;
//    fw1.setDimension(Dimension(0.01, 0, static_cast<int>(Nt)), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
//    fw1.optimalParameters[0].distribute(SpacePoint(0.38, 0.68));
//    fw1.optimalParameters[1].distribute(SpacePoint(0.78, 0.28));

//    DoubleVector x;
//    fw1.parameterToVector(x);
}

void get_vector_size(double size)
{
    size = 10;
}

void init_strt_vector(double* x0)
{
    DoubleVector x;
    fw1.parameterToVector(x);

    unsigned int size = x.length();
    x0 = new double[size];
    for (unsigned int i=0; i<size; i++) x0[i] = x[i];


}

double call_fx(double *x)
{
    return fw1.fx(DoubleVector(x, 16));
}

void call_gr(double *x, double *g, unsigned int size)
{
    //DoubleVector px(size); for (unsigned int i=0; i<size; i++) px[i] = x[i];
    DoubleVector gr(size); for (unsigned int i=0; i<size; i++) gr[i] = 0.0;

    fw1.gradient(DoubleVector(x, size), gr);
    for (unsigned int i=0; i<size; i++) g[i] = gr[i];

    //px.clear();
    gr.clear();
}

void setPenaltyR(double r)
{
//    prob.r = r;
}
