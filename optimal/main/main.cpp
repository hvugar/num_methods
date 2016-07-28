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

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    //    DoubleMatrix u;
    //    P p;
    //    p.hx = 0.001;
    //    p.ht = 0.001;
    //    p.N = 1000;
    //    p.M = 1000;
    //    p.a = 1.0;
    //    p.calculateL(u, p.hx, p.ht, p.N, p.M, p.a);
    //    IPrinter::printMatrix(u);

    //    Vector a(1, 1.0);
    //    printf("%f %f %f\n", a[0], a[1], a[2]);
    //    double &r = a[0];
    //    r = 2.0;
    //    printf("%f %f %f\n", a[0], a[1], a[2]);
    //    printf("size: %d\n", a.size());
    //    int in;
    //    scanf("%d", &in);
    //    a.resize(100*1024*1024, 10.0);
    //    scanf("%d", &in);
    //    a.resize(1, 10.0);
    //    printf("size: %d\n", a.size());
    //    printf("%f %f %f %f %f\n", a[0], a[1], a[2], a[3], a[4]);
    //    scanf("%d", &in);

    Matrix m(3, 4, 2.0);
    printf("%f %f %f %f\n", m.row(0)[0], m.row(0)[1], m.row(0)[2], m.row(0)[3]);
    printf("%f %f %f %f\n", m.row(1)[0], m.row(1)[1], m.row(1)[2], m.row(1)[3]);
    printf("%f %f %f %f\n", m.row(2)[0], m.row(2)[1], m.row(2)[2], m.row(2)[3]);


    //    std::vector<CauchyProblem*> cps;
    //    cps.resize(3);
    //    cps[0] = new A;
    //    cps[0]->x0 = 0.0;
    //    cps[0]->y0 = 0.0;
    //    cps[1] = new B;
    //    cps[1]->x0 = 0.0;
    //    cps[1]->y0 = 0.0;
    //    cps[2] = new C;
    //    cps[2]->x0 = 0.0;
    //    cps[2]->y0 = 0.0;
    //    DoubleMatrix m;
    //    CauchyProblem::rungeKutta(cps, 0.0, 0.001, 1000, m);
    //    IPrinter::printVector(m[0]);
    //    IPrinter::printVector(m[1]);
    //    IPrinter::printVector(m[2]);

    //    Example2::main(argc, argv);

    //    ControlFunction4::main(argc, argv);
    //    SampleMain();
    //    SampleLoaderBorder::main();
    //    HyperbolicControl2D21::main(argc, argv);
    //    Rosenbrock::main(argc, argv);
    return 0;
}
