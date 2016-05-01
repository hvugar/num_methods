#include "heatexample1.h"

void HeatExample1::main(int argc, char **argv)
{
    HeatExample1 h1;
    DoubleMatrix u;
    h1.calculateN(u, h1.hx, h1.ht, h1.N, h1.M);

    for (unsigned int i=0; i<=h1.M; i++)
    {
        char buffer[20];
        int n = 0;
        if (i<10) n = sprintf(buffer, "data/0000000%d.txt", i);
        if (i<100 && i>=10) n = sprintf(buffer, "data/000000%d.txt", i);
        if (i<1000 && i>=100) n = sprintf(buffer, "data/00000%d.txt", i);
        if (i<10000 && i>=1000) n = sprintf(buffer, "data/0000%d.txt", i);
        buffer[n] = '\0';
        FILE *file = fopen(buffer, "w");
        IPrinter::printVector(u[i], NULL, h1.N, 0, 0, file);
        fclose(file);
    }

//    IPrinter::printVector(u[h1.M]);

//    DoubleVector v(h1.N+1);
//    for (unsigned int i=0; i<=h1.N; i++)
//    {
//        double x = i*h1.hx;
//        v[i] = 0.5*x*x + 1.0;
//    }

    //IPrinter::printVector(v);
}

HeatExample1::HeatExample1() : IParabolicEquation()
{
    a = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    x0 = 0.0;
    x1 = 1.0;
    ht = 0.001;
    hx = 0.001;
    N  = 1000;
    M  = 1000;
}

double HeatExample1::initial(unsigned int i) const
{
    double x = i*hx;
    return 0.0;
}

double HeatExample1::boundary(Boundary type, unsigned int j) const
{
    double t = j*ht;
    if (type == Left)  return 0.0;
    if (type == Right) return 5.0;
    return 0.0;
}

double HeatExample1::f(unsigned int i, unsigned int j) const
{
    double t = j*ht;
    return 0.0;//2.0*t - 2.0*a*a;
}
