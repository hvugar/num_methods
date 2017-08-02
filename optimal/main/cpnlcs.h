#ifndef CAUCHYPROBLEMNONLOCALCONTIONS_H
#define CAUCHYPROBLEMNONLOCALCONTIONS_H

#include <vector2d.h>
#include <matrix2d.h>
#include <vector>
#include <grid/cauchyp.h>

#define SAMPLE_1

using namespace std;

class CauchyProblemNonLocalContions : public SystemLinearODE1stOrder
{
public:
    static void Main(int agrc, char *argv[]);

    struct NonSeparatedCondition
    {
        double time;
        unsigned int nmbr;
        DoubleMatrix alpha;
        unsigned int rows;
        unsigned int cols;
    };

    struct SeparatedCondition
    {
        double time;
        unsigned int nmbr;
        DoubleMatrix alpha;
        unsigned int rows;
        unsigned int cols;
    };

    CauchyProblemNonLocalContions(const ODEGrid &grid);

    virtual void initialize();

    /**
     * @brief nlcs Неразделенные условия.
     */
    std::vector<NonSeparatedCondition> nscs;
    /**
     * @brief nlcs Разделенные условия.
     */
    std::vector<SeparatedCondition> scs;
    DoubleVector betta;

    unsigned int n = 2;
    unsigned int L = 3;

    double x1(unsigned int k) const;
    double x2(unsigned int k) const;

    void calculate1();

    void calculateInterval(unsigned int s, unsigned int r);
    void calculateCondition(unsigned int r);
    \
public:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const;
};

class CauchyProblemM1stOrderA : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderA(CauchyProblemNonLocalContions &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}

protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
    double S0(double t, const DoubleVector &x, unsigned int k) const;
private:
    CauchyProblemNonLocalContions &p;
};

class CauchyProblemM1stOrderB : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderB(CauchyProblemNonLocalContions &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}
protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
private:
    CauchyProblemNonLocalContions &p;
};

#endif // CAUCHYPROBLEMNONLOCALCONTIONS_H
