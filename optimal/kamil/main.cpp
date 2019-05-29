#include <iostream>

using namespace std;

#include <grid/bvp.h>
#include <ode/lode1o.h>

class FirstOrderLinearODEY0 : public FirstOrderLinearODE
{
public:
    virtual ~FirstOrderLinearODEY0() {}

protected:
    virtual double A(const PointNodeODE &node, unsigned int row, unsigned int col) const
    {
        double t = node.x;
        if (row == 1 && col == 1) return +1.0; if (row == 1 && col == 2) return -exp(-t*t);
        if (row == 2 && col == 1) return +t*t; if (row == 2 && col == 2) return +2.0*t;
        return NAN;
    }

    virtual double B(const PointNodeODE &node, unsigned int row) const
    {
        double t = node.x;
        if (row == 1) return exp(t);
        if (row == 2) return -t*t*(t*exp(t)+1.0);
        return NAN;
    }

protected:
    virtual auto initial(InitialCondition, unsigned int row) const -> double
    {
        if (row == 1) return +1.0-exp(-1.0);
        if (row == 2) return +1.0-exp(+2.0);
        return NAN;
    }

    virtual auto boundary(const PointNodeODE &, BoundaryConditionODE &, unsigned int) const -> double { return NAN; }
protected:
    virtual unsigned int count() const { return 2; }
};

class FirstOrderLinearODEY1 : public FirstOrderLinearODE
{
public:
    virtual ~FirstOrderLinearODEY1() {}

protected:
    virtual double A(const PointNodeODE &node, unsigned int row, unsigned int col) const
    {
        double t = node.x;
        if (row == 1 && col == 1) return +1.0; if (row == 1 && col == 2) return -exp(-t*t);
        if (row == 2 && col == 1) return +t*t; if (row == 2 && col == 2) return +2.0*t;
        return NAN;
    }

    virtual double B(const PointNodeODE &node, unsigned int row) const
    {
        double t = node.x;
        if (row == 1) return +0.0;
        if (row == 2) return +0.0;
        return NAN;
    }

protected:
    virtual auto initial(InitialCondition, unsigned int row) const -> double
    {
        if (row == 1) return +0.0;
        if (row == 2) return -1.0;
        return NAN;
    }

    virtual auto boundary(const PointNodeODE &, BoundaryConditionODE &, unsigned int) const -> double { return NAN; }
protected:
    virtual unsigned int count() const { return 2; }
};

int main()
{
    FirstOrderLinearODEY0 f1;
    f1.setDimension(Dimension(0.0001, -10000, 10000));
    std::vector<DoubleVector> y0;
    f1.solveInitialValueProblem(y0);

    for (unsigned int i=0; i<=y0.size(); i++) { if (i%2000==0) printf("%12.6f ", y0[i][0]); } puts("");
    for (unsigned int i=0; i<=y0.size(); i++) { if (i%2000==0) printf("%12.6f ", y0[i][1]); } puts("");

    cout << "-------" << endl;

    FirstOrderLinearODEY1 f2;
    f2.setDimension(Dimension(0.0001, -10000, 10000));
    std::vector<DoubleVector> y1;
    f2.solveInitialValueProblem(y1);

    for (unsigned int i=0; i<=y1.size(); i++) { if (i%2000==0) printf("%12.6f ", y1[i][0]); } puts("");
    for (unsigned int i=0; i<=y1.size(); i++) { if (i%2000==0) printf("%12.6f ", y1[i][1]); } puts("");

    cout << "-------" << endl;

    double y02_m1 = y0[0][1];
    double y02_00 = y0[10000][1];
    double y02_p1 = y0[20000][1];
    double y12_00 = y1[10000][1];
    double y12_p1 = y1[20000][1];

    double a = 1.0 + y12_00*y12_00 - 2.0*y12_p1*y12_p1;
    double b = -2.0*y02_00*y12_00
               -2.0*y12_00*y12_00*y02_m1
               +4.0*y02_p1*y12_p1
               +4.0*y12_p1*y12_p1*y02_m1;
    double c = y02_00*y02_00
             + 2.0*y02_00*y12_00*y02_m1
             + y12_00*y12_00*y02_m1*y02_m1
             - 2.0*y02_p1*y02_p1
             - 4.0*y02_p1*y12_p1*y02_m1
             - 2.0*y12_p1*y12_p1*y02_m1*y02_m1
             - y02_m1;

    double D = b*b - 4.0*a*c;
    double z1_1 = (-b - sqrt(D))/(2.0*a);
    double z1_2 = (-b + sqrt(D))/(2.0*a);
    double g = y02_m1 - z1_1;
    double z2 = y02_00 + y12_00*g;
    double z3 = y02_p1 + y12_p1*g;
    printf("%f %f %f %f %f\n", z1_1, z1_2, D, z2, z3);

    std::vector<DoubleVector> x(y1.size());
    std::vector<DoubleVector> x0(y1.size());
    for (unsigned int i=0; i<y1.size(); i++)
    {
        x[i].resize(2);
        x0[i].resize(2);
        x[i] = y0[i] + y1[i]*g;

        double t = (static_cast<int>(i)-10000)*0.0001;
        x0[i][0] = t*exp(t) + 1.0;
        x0[i][1] = exp(t*t);
    }

    for (unsigned int i=0; i<=y1.size(); i++) { if (i%2000==0) printf("%12.6f ", x[i][0]); } puts("");
    for (unsigned int i=0; i<=y1.size(); i++) { if (i%2000==0) printf("%12.6f ", x0[i][0]); } puts("");
    for (unsigned int i=0; i<=y1.size(); i++) { if (i%2000==0) printf("%12.6f ", x[i][1]); } puts("");
    for (unsigned int i=0; i<=y1.size(); i++) { if (i%2000==0) printf("%12.6f ", x0[i][1]); } puts("");

    cout << "Hello World!" << endl;
    return 0;
}
