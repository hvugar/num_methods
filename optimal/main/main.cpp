#include "headers.h"

#include "control/example1.h"
#include "control/example2.h"

struct A : public CauchyProblem
{
    double f(double x, const DoubleVector &y) const
    {
        double y1 = y[0];
        double y2 = y[1];
        double y3 = y[2];
        return x*y1 + y2 + 2.0*y3 + (2.0*x - 2.0*x*x*x - 2.0*sin(x));
    }
};

struct B : public CauchyProblem
{
    double f(double x, const DoubleVector &y) const
    {
        double y1 = y[0];
        double y2 = y[1];
        double y3 = y[2];
        return y1 + 2.0*x*x;
    }
};

struct C : public CauchyProblem
{
    double f(double x, const DoubleVector &y) const
    {
        double y1 = y[0];
        double y2 = y[1];
        double y3 = y[2];
        return y1 + 2.0*y2 + y3 + (cos(x) - x*x - 2.0*x*x*x - sin(x));
    }
};

class P : public IParabolicEquation
{
public:
    double initial(unsigned int i) const
    {
        double x = i*hx;
        return x*x*x;
    }
    double boundary(Boundary type, unsigned int j) const
    {
        double t = j*ht;
        if (type == Left) return t*t;
        if (type == Right) return t*t+1.0;
        return 0.0;
    }
    double f(unsigned int i, unsigned int j) const
    {
        double x = i*hx;
        double t = j*ht;
        return 2.0*t - 6.0*x*a*a;
    }
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double a;
};

class P1 : public IParabolicEquation
{
public:
    double initial(unsigned int i) const
    {
        double x = i*hx;
        return x*x*x*x;
    }
    double boundary(Boundary type, unsigned int j) const
    {
        double t = j*ht;
        if (type == Left) return 0.0;
        if (type == Right) return 4.0;
        return 0.0;
    }
    double f(unsigned int i, unsigned int j) const
    {
        double x = i*hx;
        double t = j*ht;
        return 4.0*t*t*t - 12.0*a*a*x*x;
    }
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double a;
};

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    Matrix2D *m = new Matrix2D(100, 200);
    printf("%d %d\n", m->rows(), m->cols());
    printf("%f\n", (*m)(0,0));
    (*m)(0,0) = 10;
    printf("%f\n", (*m)(0,0));
    m->at(10, 20) = 12.0f;
    printf("%f\n", m->at(10,20));

    //    {
    //        DoubleMatrix u;
    //        P p;
    //        p.hx = 0.01;
    //        p.N = 100;

    //        p.ht = 0.00005;
    //        p.M = 20000;

    //        p.a = 1.0;
    //        p.calculateL(u, p.hx, p.ht, p.N, p.M, p.a);
    //        IPrinter::printMatrix(u);
    ////        IPrinter::printVector(u[10]);
    //    }

    //    {
    //        P1 p;
    //        p.hx = 0.001;
    //        p.N = 1000;

    //        p.ht = 0.001;
    //        p.M = 1000;

    //        p.a = 1.0;
    //        DoubleMatrix u1;
    //        p.calculateN1(u1, p.hx, p.ht, p.N, p.M, p.a);
    //        IPrinter::printMatrix(u1);
    //        puts("---");
    //        DoubleMatrix u2;
    //        p.calculateN(u2, p.hx, p.ht, p.N, p.M, p.a);
    //        IPrinter::printMatrix(u2);
    //    }

    //    {
    //        std::vector<CauchyProblem*> cps;
    //        cps.resize(3);
    //        cps[0] = new A;
    //        cps[0]->x0 = 0.0;
    //        cps[0]->y0 = 0.0;
    //        cps[1] = new B;
    //        cps[1]->x0 = 0.0;
    //        cps[1]->y0 = 0.0;
    //        cps[2] = new C;
    //        cps[2]->x0 = 0.0;
    //        cps[2]->y0 = 0.0;
    //        DoubleMatrix m;
    //        CauchyProblem::rungeKutta(cps, 0.0, 0.001, 1000, m);
    //        IPrinter::printVector(m[0]);
    //        IPrinter::printVector(m[1]);
    //        IPrinter::printVector(m[2]);
    //    }

    //    Example2::main(argc, argv);

    //    ControlFunction4::main(argc, argv);
    //    SampleMain();
    //    SampleLoaderBorder::main();
    //    HyperbolicControl2D21::main(argc, argv);
    //    Rosenbrock::main(argc, argv);
    return 0;
}
