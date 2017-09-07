#ifndef PROBLEM4EX2_H
#define PROBLEM4EX2_H

#include <function.h>
#include <nonlinearequation.h>
#include <grid/grid.h>
#include <ode/lode1o.h>

class NonLinearEquationEx2;
class Problem4Ex2Zetta0;
class Problem4Ex2Zettai;

#define SAMPLE_2

class Problem4Ex2
{
public:
    static void Main(int agrc, char *argv[]);

    Problem4Ex2();
    virtual ~Problem4Ex2();

    void initialize();
    double X(double t, unsigned int num) const;
    //virtual double g(unsigned int num, unsigned int row) const;
    virtual double g(const DoubleVector &x, unsigned int num, unsigned int row) const;

    void printResult1(const DoubleVector &x);

    UniformODEGrid grid;

    std::vector<LinearODE1stOrder::Condition> cs;
    DoubleVector betta;

protected:
    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int k, unsigned int row) const;
    virtual double C(double t, unsigned int k, unsigned int num, unsigned int row, unsigned int col) const;

private:
    std::vector<DoubleVector> zm0;
    std::vector<std::vector<DoubleVector>> zm1;
    std::vector<std::vector<DoubleVector>> zm2;

    friend class NonLinearEquationEx2;
    friend class Problem4Ex2Zetta0;
    friend class Problem4Ex2Zettai;
};

class NonLinearEquationEx2 : public NonLinearEquation
{
public:
    NonLinearEquationEx2(const Problem4Ex2 &p);
//protected:
    virtual double fx(const DoubleVector &x, unsigned int num = 0) const;
private:
    const Problem4Ex2 &p;
};

class Problem4Ex2Zetta0 : public LinearODE1stOrder
{
public:
    Problem4Ex2Zetta0(const Problem4Ex2 &p);
    virtual unsigned int equationsNumber() const;
protected:
    virtual double A(double x, unsigned int i, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(double x, unsigned int i, unsigned int row = 0) const;
private:
    const Problem4Ex2 &p;
};

class Problem4Ex2Zettai : public LinearODE1stOrder
{
public:
    Problem4Ex2Zettai(const Problem4Ex2 &p4, unsigned int i);
    virtual unsigned int equationsNumber() const;

    void calculateM(const std::vector<LinearODE1stOrder::Condition> &cs, const DoubleMatrix &betta, std::vector<std::vector<DoubleVector>> &zmi);

protected:
    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int k, unsigned int row) const;
    virtual double C(double t, unsigned int k, unsigned int row, unsigned int col) const;
private:
    const Problem4Ex2 &p;
    unsigned int i;
    unsigned int cur_col;
};


#endif // PROBLEM4EX2_H
